from anndata import AnnData
import numpy as np
from typing import Tuple, Dict
import pandas as pd
from flask import request
import bisect
from utils.colors import get_color_map
import plotly.express as px
from plotly.graph_objs import Figure

class AlignmentData:
    """
        A class to manage the alignment and transformation of 3D scatter data from an AnnData object.
        
        Attributes:
            _adata (AnnData): The AnnData object containing the dataset.
            _xfield (str): The field used for the x-axis in alignment.
            _yfield (str): The field used for the y-axis in alignment.
            _zfield (str): The field used for the z-axis in alignment.
            _x_aligned_field (str): The field used for the x-aligned coordinates.
            _y_aligned_field (str): The field used for the y-aligned coordinates.
            _orifig (Figure): The original figure for visualization.
            _alifig (Figure): The aligned figure for visualization.
            _rangescale (tuple): Range and scale for the transformation.
            _slicelist (list): List of slice values for the z-field.
            _coord_trans_mtx (dict): Dictionary storing transformation matrix for each user.
            _actref_slice (dict): Dictionary storing active and reference slices for each user.
            _alistatus (dict): Dictionary to track the status of the alignment process. Includes status for each step, the creator, and the thread.
    """
    def __init__(self):
        """
            Initializes the AlignmentData object with default values.
        """
        self._adata = None
        self._xfield = None
        self._yfield = None
        self._x_aligned_field = None
        self._y_aligned_field = None
        self._zfield = None
        self._orifig = None
        self._alifig = None
        self._rangescale = None #(oriRange, oriScale)
        self._slicelist = None
        self._coord_trans_mtx = {} # {'usrname': [dx, dy, deg] }
        self._actref_slice = {} # {'usrname': [actSlice, refSlice, reverse:bool] }
        self._alistatus = {
            'creator':None,
            'thread':None,
            'steps':{
                1:{'running':False, 'complete':False, 'failed':False, 'experror':False},
                2:{'running':False, 'complete':False, 'failed':False},
                3:{'running':False, 'complete':False, 'failed':False, 'percent':0},
                4:{'running':False, 'complete':False, 'failed':False},
                5:{'running':False, 'complete':False, 'failed':False},
            }
        }

    def set_alistatus_thread(self, thread):
        """
        Sets the thread associated with the alignment process.

        Args:
            thread: The thread identifier or object.
        """
        self._alistatus['thread'] = thread

    def set_alistatus_creator(self, creator:str):
        """
        Sets the creator of the alignment process.

        Args:
            creator (str): The name or identifier of the creator.
        """
        self._alistatus['creator'] = creator
    
    def set_alistatus_steps(self, step:int, running:bool=None, complete:bool=None, failed:bool=None, percent:float=None, experror:bool=None):
        """
        Updates the status of a specific step in the alignment process.

        Args:
            step (int): The step number to update.
            running (bool, optional): Indicates if the step is running.
            complete (bool, optional): Indicates if the step is complete.
            failed (bool, optional): Indicates if the step has failed.
            percent (float, optional): The progress percentage (for step 3).
            experror (bool, optional): Indicates if there is an gene expression value error.
        """
        s = self.get_alistatus_steps()[step]
        if running is not None:
            s['running'] = running
        if complete is not None:
            s['complete'] = complete
        if failed is not None:
            s['failed'] = failed
        if percent is not None:
            s['percent'] = percent
        if experror is not None:
            s['experror'] = experror

    def reset_alistatus(self):
        """
        Resets the alignment status to its default state.
        Clears all progress and status fields for the alignment process.
        """
        self._alistatus = {
            'creator':None,
            'thread':None,
            'steps':{
                1:{'running':False, 'complete':False, 'failed':False, 'experror':False},
                2:{'running':False, 'complete':False, 'failed':False},
                3:{'running':False, 'complete':False, 'failed':False, 'percent':0},
                4:{'running':False, 'complete':False, 'failed':False},
                5:{'running':False, 'complete':False, 'failed':False},
            }
        }

    def get_alistatus_thread(self):
        """
        Returns the thread identifier associated with the alignment process.

        Returns:
            The thread associated with the alignment process.
        """
        return self._alistatus['thread']
    
    def get_alistatus_creator(self):
        """
        Returns the creator of the alignment process.

        Returns:
            The creator identifier or name.
        """
        return self._alistatus['creator']
    
    def get_alistatus_steps(self)->Dict:
        """
        Returns the status of all steps in the alignment process.

        Returns:
            dict: The status dictionary for each step (running, complete, failed, etc.).
        """
        return self._alistatus['steps']
    
    def test(self, fpath:str, xfield:str, yfield:str, zfield:str)->None:
        """
            A test function to initialize the class with a given AnnData file and fields for x, y, and z axes.

            Args:
                fpath (str): The path to the AnnData file.
                xfield (str): The field name for the x-axis.
                yfield (str): The field name for the y-axis.
                zfield (str): The field name for the z-axis.
        """
        import scanpy as sc
        adata = sc.read_h5ad(fpath)
        self.set_adata(adata)
        self.set_xfield(xfield)
        self.set_yfield(yfield)
        self.set_zfield(zfield)
        self.plot_3d_scatter(origin_fig=True)
        self.set_x_aligned_field('x_aligned')
        self.set_y_aligned_field('y_aligned')
        self.plot_3d_scatter()

    def reset_props(self)->None:
        """
            Resets all internal properties of the object to their default values.
        """
        self._adata = None
        self._xfield = None
        self._yfield = None
        self._zfield = None
        self._x_aligned_field = None
        self._y_aligned_field = None
        self._orifig = None
        self._alifig = None
        self._rangescale = None
        self._slicelist = None
        self._coord_trans_mtx = {}
        self._cur_actref_slice = {}

    def delete_user_data(self, usrname:str)->None:
        """
            Deletes the transformation data and slice information for a given user.

            Args:
                usrname (str): The username for which data is to be deleted.
        """
        if usrname in self._coord_trans_mtx:
            del self._coord_trans_mtx[usrname]
        if usrname in self._actref_slice:
            del self._actref_slice[usrname]

    def get_movestep_size(self, activeSlice:float, referenceSlice:float)->float:
        """
            Calculates the average move step size between two slices based on their aligned coordinates.

            Args:
                activeSlice (float): The active slice value.
                referenceSlice (float): The reference slice value.

            Returns:
                float: The calculated move step size.
        """
        adata = self._adata
        z_column = self._zfield
        activeObs = adata[adata.obs[z_column]==activeSlice].obs
        referenceObs = adata[adata.obs[z_column]==referenceSlice].obs
        act_x_min = np.min(activeObs['x_aligned'])
        act_x_max = np.max(activeObs['x_aligned'])
        ref_x_min = np.min(referenceObs['x_aligned'])
        ref_x_max = np.max(referenceObs['x_aligned'])
        act_y_min = np.min(activeObs['y_aligned'])
        act_y_max = np.max(activeObs['y_aligned'])
        ref_y_min = np.min(referenceObs['y_aligned'])
        ref_y_max = np.max(referenceObs['y_aligned'])
        step = np.mean([
            abs(act_x_min-ref_x_min), 
            abs(act_x_max-ref_x_max), 
            abs(act_y_min-ref_y_min),
            abs(act_y_max-ref_y_max)
        ])
        return step

    def get_graph_view(self, actSlice:float, refSlice:float)->dict:
        """
            Determines the 3D camera view for the plot based on the active and reference slices.

            Args:
                actSlice (float): The active slice value.
                refSlice (float): The reference slice value.

            Returns:
                dict: A dictionary containing the camera view for 'bottom', 'top', and 'center'.
        """
        eye = {
            'bottom': {'x': -3.1294811103329714e-20, 'y': -1.3257190114690003e-16, 'z': -2.1650635094610964},
            'top': {'x': 1.2149787566408422e-05, 'y': -0.006897060070669372, 'z': 2.1650525237080886},
            'center': {'x': 1.25, 'y': 1.25, 'z': 1.25}
        }
        if not actSlice or not refSlice:
            return eye['center']
        if float(actSlice)>=float(refSlice):
            return eye['top']
        return eye['bottom']

    def set_actref_slice(self, usrname:str, act_slice=None, ref_slice=None, reverse=None)->None:
        """
            Sets the active and reference slices for a given user.

            Args:
                usrname (str): The username.
                act_slice (float, optional): The active slice value.
                ref_slice (float, optional): The reference slice value.
                reverse (bool, optional): Whether to reverse the transformation.
        """
        if usrname not in self._actref_slice:
            self._actref_slice[usrname] = [None, None, None]
        if act_slice:
            self._actref_slice[usrname][0] = float(act_slice)
        if ref_slice:
            self._actref_slice[usrname][1] = float(ref_slice)
        if reverse!=None:
            self._actref_slice[usrname][2] = reverse

    def get_actref_slice(self, usrname:str)->list:
        """
            Retrieves the active and reference slices for a given user.

            Args:
                usrname (str): The username.

            Returns:
                list: The list containing active slice, reference slice, and reverse flag.
        """
        if usrname not in self._actref_slice:
            return None
        return self._actref_slice[usrname]

    def set_coord_trans_mtx(self, usrname:str, dx=None, dy=None, reg=None)->None:
        """
            Sets the coordinate transformation matrix for a given user.

            Args:
                usrname (str): The username.
                dx (float, optional): The translation along the x-axis.
                dy (float, optional): The translation along the y-axis.
                reg (float, optional): The rotation angle in degrees.
        """
        if not dx and not dy and not reg:
            self._coord_trans_mtx[usrname] = [0, 0, 0]
            return
        if usrname not in self._coord_trans_mtx:
            self._coord_trans_mtx[usrname] = [0, 0, 0]
        reverse = self.get_actref_slice(usrname)[2]
        if dx:
            self._coord_trans_mtx[usrname][0]+=dx
        if dy:
            self._coord_trans_mtx[usrname][1]+=-dy if reverse else dy
        if reg:
            self._coord_trans_mtx[usrname][2]+=-reg if reverse else reg

    def get_coord_trans_mtx(self, usrname:str)->list:
        """
            Retrieves the coordinate transformation matrix for a user.

            Args:
                usrname (str): Username.

            Returns:
                list: A list containing the translation and rotation matrix for the user, or None if not set.
        """
        if usrname in self._coord_trans_mtx:
            return self._coord_trans_mtx[usrname]
        return None

    def get_slice_index(self, slice:float):
        """
            Retrieves the index of the given slice from the slice list.

            Args:
                slice (float): Slice value.

            Returns:
                int: Index of the slice in the slice list, or None if not found.
        """
        slices = self.get_slice_list()
        if slices is None:
            return None
        idx = bisect.bisect_left(slices, slice)
        return idx

    def get_slice_list(self)->list:
        """
            Retrieves the list of unique slices from the data.

            Returns:
                list: A sorted list of slice values, or None if no slices are available.
        """
        return self._slicelist
    
    def set_slice_list(self)->None:
        """
            Sets the slice list based on the unique values of the z-field in the AnnData object.
        """
        ad = self._adata
        zfield = self._zfield
        if ad is None or zfield is None:
            return
        slices = list(ad.obs[zfield].unique())
        slices.sort()
        self._slicelist = slices

    def get_zfield_min(self)->float:
        """
            Retrieves the minimum value of the z-field in the AnnData object.

            Returns:
                float: Minimum value of the z-field, or None if the AnnData or z-field is not set.
        """
        if self._adata is None or self._zfield is None:
            return None
        col = self._adata.obs[self._zfield]
        return col.min()
    
    def get_zfield_max(self)->float:
        """
            Retrieves the maximum value of the z-field in the AnnData object.

            Returns:
                float: Maximum value of the z-field, or None if the AnnData or z-field is not set.
        """
        if self._adata is None or self._zfield is None:
            return None
        col = self._adata.obs[self._zfield]
        return col.max()

    def set_range_scale(self, ori_range:float, ori_scale:float)->None:
        """
            Sets the original range and scale for adjusting plot scaling.

            Args:
                ori_range (float): Original range for scaling.
                ori_scale (float): Original scale factor.
        """
        self._rangescale = (ori_range, ori_scale)
    
    def get_range_scale(self)->tuple:
        """
            Retrieves the original range and scale for adjusting plot scaling.

            Returns:
                tuple: A tuple containing the original range and scale values.
        """
        return self._rangescale
    
    def set_orifig(self, fig:Figure)->None:
        """
            Sets the original figure for plotting.

            Args:
                fig (Figure): The plotly figure to be set as the original figure.
        """
        self._orifig = fig
    
    def get_orifig(self)->Figure:
        """
            Retrieves the original figure for plotting.

            Returns:
                Figure: The original plotly figure.
        """
        return self._orifig
    
    def set_alifig(self, fig:Figure)->None:
        """
            Sets the aligned figure for plotting.

            Args:
                fig (Figure): The plotly figure to be set as the aligned figure.
        """
        self._alifig = fig
    
    def get_alifig(self)->Figure:
        """
            Retrieves the aligned figure for plotting.

            Returns:
                Figure: The aligned plotly figure.
        """
        return self._alifig
    
    def has_orifig(self)->bool:
        """
            Checks if the original figure is set.

            Returns:
                bool: True if the original figure is set, otherwise False.
        """
        return self._orifig!=None
    
    def has_alifig(self)->bool:
        """
            Checks if the aligned figure is set.

            Returns:
                bool: True if the aligned figure is set, otherwise False.
        """
        return self._alifig!=None
    
    def set_adata(self, adata:AnnData)->None:
        """
            Sets the AnnData object containing the data.

            Args:
                adata (AnnData): The AnnData object to be set.
        """
        self._adata = adata
    
    def get_adata(self)->AnnData:
        """
            Retrieves the AnnData object containing the data.

            Returns:
                AnnData: The AnnData object.
        """
        return self._adata
    
    def set_x_aligned_field(self, xfield:str)->None:
        """
            Sets the aligned x-field for plotting.

            Args:
                xfield (str): The name of the aligned x-field in the AnnData object.
        """
        self._x_aligned_field = xfield
    
    def get_x_aligned_field(self)->str:
        """
            Retrieves the aligned x-field for plotting.

            Returns:
                str: The name of the aligned x-field.
        """
        return self._x_aligned_field
    
    def set_y_aligned_field(self, yfield:str)->None:
        """
            Sets the aligned y-field for plotting.

            Args:
                yfield (str): The name of the aligned y-field in the AnnData object.
        """
        self._y_aligned_field = yfield
    
    def get_y_aligned_field(self)->str:
        """
            Retrieves the aligned y-field for plotting.

            Returns:
                str: The name of the aligned y-field.
        """
        return self._y_aligned_field
    def set_xfield(self, xfield:str)->None:
        """
            Sets the x-field for plotting and ensures its values are numeric.

            Args:
                xfield (str): The name of the x-field in the AnnData object.
        """
        self._adata.obs[xfield] = pd.to_numeric(self._adata.obs[xfield], errors='raise')
        self._xfield = xfield
    
    def get_xfield(self)->str:
        """
            Retrieves the x-field for plotting.

            Returns:
                str: The name of the x-field.
        """
        return self._xfield
    
    def set_yfield(self, yfield:str)->None:
        """
            Sets the y-field for plotting and ensures its values are numeric.

            Args:
                yfield (str): The name of the y-field in the AnnData object.
        """
        self._adata.obs[yfield] = pd.to_numeric(self._adata.obs[yfield], errors='raise')
        self._yfield = yfield
    
    def get_yfield(self)->str:
        """
            Retrieves the y-field for plotting.

            Returns:
                str: The name of the y-field.
        """
        return self._yfield
    
    def set_zfield(self, zfield:str)->None:
        """
            Sets the z-field for plotting and ensures its values are numeric, then updates the slice list.

            Args:
                zfield (str): The name of the z-field in the AnnData object.
        """
        self._adata.obs[zfield] = pd.to_numeric(self._adata.obs[zfield], errors='raise')
        self._zfield = zfield
        self.set_slice_list()
    
    def get_zfield(self)->str:
        """
            Retrieves the z-field for plotting.

            Returns:
                str: The name of the z-field.
        """
        return self._zfield

    def plot_3d_scatter(self, origin_fig:bool=False, marker_size:int=5, boarder_width:int=1, boarder_color:str='#0d0015', grid_color:str = '#5F9EA0')->Figure:
        """
            Plots a 3D scatter plot based on the data in the AnnData object.

            Args:
                origin_fig (bool, optional): If True, plots the original figure; otherwise, plots the aligned figure.
                marker_size (int, optional): Size of the plot markers. Default is 5.
                boarder_width (int, optional): Border width of the markers. Default is 1.
                boarder_color (str, optional): Color of the marker borders. Default is '#0d0015'.
                grid_color (str, optional): Color of the gridlines. Default is '#5F9EA0'.

            Returns:
                Figure: The generated plotly figure.
        """
        df = self.get_adata().obs.copy()
        if origin_fig:
            x = self._xfield
            y = self._yfield
        else:
            x = self._x_aligned_field
            y = self._y_aligned_field
        z = self._zfield
        df['layers'] = df[z].astype(str)
        df['index'] = df.index
        df = df.sort_values(by=z)
        cmap = get_color_map(set(df['layers']))
        fig = px.scatter_3d(
            df, x=x, y=y, z=z, 
            color='layers', 
            custom_data='index',
            color_discrete_map=cmap,
        )

        x_min = np.min(df[x])
        x_max = np.max(df[x])
        y_min = np.min(df[y])
        y_max = np.max(df[y])

        min_val = min(x_min, y_min)
        max_val = max(x_max, y_max)
        dist = max_val-min_val
        max_val += dist/2
        min_val -= dist/2
        dist = max_val-min_val

        if origin_fig or self.get_range_scale()==None:
            step = dist/20
            x_scale = dist/(x_max-x_min)
            y_scale = dist/(y_max-y_min)
            scale = (x_scale+y_scale)/2
            if origin_fig:
                self.set_range_scale(dist, scale)
        else:
            ori_range, ori_scale = self.get_range_scale()
            multy_factor = dist/ori_range
            step = dist/20/multy_factor
            scale = multy_factor*ori_scale

        fig.update_traces(
            marker=dict(
                size=marker_size,
                line=dict(color=boarder_color, width=boarder_width)
            )
        )
        fig.update_layout(
            autosize=True,
            legend=dict(
                itemsizing='constant',
                traceorder = 'reversed',
                title=None
            ),
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            scene=dict(
                aspectratio=dict(x=scale, y=scale, z=1),
                xaxis=dict(visible=False, showgrid=True, showbackground=False, showticklabels=False, gridcolor=grid_color, tickfont=dict(color=grid_color), title='', dtick=step, range=[min_val, max_val]), 
                yaxis=dict(visible=False, showgrid=True, showbackground=False, showticklabels=False, gridcolor=grid_color, tickfont=dict(color=grid_color), title='', dtick=step, range=[min_val, max_val]), 
                zaxis=dict(visible=False, showgrid=False, showbackground=False, showticklabels=False, gridcolor=grid_color),
                camera=dict(projection=dict(type='orthographic'))
            ),
        )
        if origin_fig:
            self.set_orifig(fig)
        else:
            self.set_alifig(fig)
        return fig
    
    def get_obs_fields(self)->list:
        """
            Retrieves the list of column names (fields) from the 'obs' attribute of the AnnData object.

            Returns:
                list: A list of column names in the 'obs' dataframe, or an empty list if the AnnData object is None.
        """
        ad = self.get_adata()
        if ad is None:
            return []
        return ad.obs.columns.tolist()
    
    def get_gene_list(self)->list:
        """
            Retrieves a sorted list of gene names from the 'var_names' attribute of the AnnData object.

            Returns:
                list: A sorted list of gene names, or an empty list if the AnnData object is None.
        """
        ad = self.get_adata()
        if ad is None:
            return []
        genes = list(ad.var_names)
        genes.sort()
        return genes
    
alidata = AlignmentData()