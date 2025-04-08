import scanpy as sc
from dataManager.alignment_d import alidata
from controller.notice import set_head_notice
from controller.auth import verify_modify_permission
import pandas as pd
import threading
import feffery_antd_components as fac
from typing import Tuple, List
import numpy as np
from dash import set_props, html, Patch
from utils.colors import get_color_map, get_scale_colors
from api.alignment import paste1, paste2, ensure_numeric_fields, validate_exp_data
import os

# alidata.test('/data1/zhouyb/public/data/stMouseEmbryo/embryo_1-2-E7.75_min400_Ann_HC0.5.h5ad', 'x', 'y', 'z', 'x', 'y')

graphStyle = {'display':'block', 'height':'90vh', 'width':'100%'}
grapHidden = {'display':'none', 'height':'90vh', 'width':'100%'}

def update_alignment_progress(progress_only:bool=False):
    steps = alidata.get_alistatus_steps()
    success = False
    for i in range(1, 6):
        id = 'step'+str(i)
        if steps[i]['failed']:
            set_props(id, {'icon':'fc-high-priority'})
            set_props('alignment-interval', {'disabled':True})
            set_props('alignment-button-alignSlices', {'disabled':False})
            if i==1 and steps[1]['experror']:
               set_head_notice('Gene expression values cannot be negative. Please check and resubmit your data !', 'warning') 
            alistatus = alidata.get_alistatus()
            e = alistatus['exception']
            if e:
                set_props('alignment-button-showBug', dict(style={'backgroundColor':'#bb5548'}))
            return
        elif steps[i]['complete']:
            set_props(id, {'icon':'fc-ok'})
            if i==3:
                set_props('alignment-alipercent', {'percent':100})
            if i==5:
                success = True
        elif steps[i]['running']:
            set_props(id, {'icon':'antd-loading'})
            if i==3:
                set_props('alignment-alipercent', {'percent':steps[3]['percent']})
            return
    if success:
        alidata.set_alistatus_thread(thread=None)
        set_props('alignment-interval', {'disabled':True})
        set_props('alignment-button-alignSlices', {'disabled':False})
        set_props('alignment-paragraph-fieldnotice', {'style':{'marginTop':'-20px'}})
        if not progress_only:
            alidata.set_x_aligned_field('x_aligned')
            alidata.set_y_aligned_field('y_aligned')
            fig = alidata.plot_3d_scatter()
            set_props('alignment-graph-aligned', {'figure':fig, 'style':graphStyle})
            obs_fields = alidata.get_obs_fields()
            genes = alidata.get_gene_list()
            set_props('graphSetting-selecter-colorField', {'options':obs_fields, 'value':None})
            set_props('graphSetting-selecter-colorGene', {'options':genes})
            slices = alidata.get_slice_list()
            if slices is not None:
                set_props('graphSetting-selecter-activeSlice', {'options':slices})
                set_props('graphSetting-selecter-referenceSlice', {'options':slices})
                set_props('graphSetting-table-applySlices', {'data':[{'slices':slice} for slice in slices]})

def align_slices_with_paste(model:str, use_gpu:bool, creator:str)->None:
    """
    Align slices using a given model and apply the alignment to the slices.

    Args:
        model (str): The alignment model to use ('paste1' or 'paste2').
        use_gpu (bool): Whether to use GPU for alignment computation.
    
    """
    permission = verify_modify_permission()
    if not permission:
        return
    status = check_alignment_data()
    if not status:
        return
    has_orifig = check_orifig()
    if not has_orifig:
        return
    ad = alidata.get_adata()
    x = alidata.get_xfield()
    y = alidata.get_yfield()
    z = alidata.get_zfield()
    alidata.reset_alistatus()
    alidata.set_alistatus_creator(creator)
    alistatus = alidata.get_alistatus()
    if model == 'paste1':
        thread = threading.Thread(target=paste1, args=(ad, x, y, z), kwargs={'use_gpu': use_gpu, 'alistatus': alistatus})
    else:
        thread = threading.Thread(target=paste2, args=(ad, x, y, z), kwargs={'alistatus': alistatus})
    alidata.set_alistatus_thread(thread)
    set_props('alignment-interval', {'disabled':False})
    reset_alistatus_progress(creator)
    thread.start()
def reset_alistatus_progress(creator:str)->None:
    set_props('alignment-button-alignSlices', {'disabled':True})
    set_props('alignment-timeline-creator', {'children':creator})
    set_props('alignment-button-showBug', dict(style={'backgroundColor':'#bb5548', 'display':'none'}))
    set_props('alignment-alipercent', {'percent':0})
    for i in range(1, 6):
        id = 'step'+str(i)
        set_props(id, {'icon':'md-schedule'})
    set_props('alignment-paragraph-fieldnotice', {'style':{'marginTop':'-20px', 'display':'none'}})
    set_props('alignment-popupcard-alignTask', {'visible':True})

def get_transformed_coord(usrname:str, slices:set, x='x_aligned', y='y_aligned', storeNewCoord=False)->dict:
        """
        Get transformed coordinates and optionally store them.

        Args:
            usrname (str): Username used to retrieve a specific transformation matrix.
            slices (set): Set of slices to be processed.
            x (str): X-axis field name, default is 'x_aligned'.
            y (str): Y-axis field name, default is 'y_aligned'.
            storeNewCoord (bool): Whether to store the new coordinates, default is False.

        Returns:
            dict: A dictionary mapping each slice to its transformed coordinates.
                If storing new coordinates, returns the updated coordinates.
        """
        transMtx = alidata.get_coord_trans_mtx(usrname)
        adata = alidata.get_adata()
        zColumn = alidata.get_zfield()
        if transMtx==None or not adata or zColumn==None:
            return None
        aliFig = alidata.get_alifig()
        newPoints = {}
        for slice in slices:
            idx = alidata.get_slice_index(slice)
            obsIndex = aliFig['data'][idx]['customdata'].flatten()
            points_ori = np.array(adata.obs.loc[obsIndex, [x, y]])
            points_new = transform_points(points_ori,transMtx)
            newPoints[slice] = points_new
            if storeNewCoord:
                adata.obs.loc[obsIndex, [x, y]] = points_new
                aliFig['data'][idx]['x'] = points_new[:, 0]
                aliFig['data'][idx]['y'] = points_new[:, 1]

        if storeNewCoord:
            adata.obsm['X_spatial_registered'] = np.array(adata.obs[[x, y, zColumn]], dtype=np.float64)

        return newPoints

def transform_points(points:np.ndarray, coordTransMtx:list)->np.ndarray:
        """
        Transform coordinates based on a given transformation matrix.

        Args:
            points (np.ndarray): Original coordinates, shape (n, 2).
            coordTransMtx (list): Transformation matrix containing translation (dx, dy) and rotation angle (angle_deg).

        Returns:
            np.ndarray: Transformed coordinates.
        """
        dx, dy, angle_deg = coordTransMtx
        center = np.mean(points, axis=0)    
        points_centered = points - center
        angle_rad = np.radians(angle_deg)
        rotation_matrix = np.array([
            [np.cos(angle_rad), -np.sin(angle_rad)],
            [np.sin(angle_rad),  np.cos(angle_rad)]
        ])
        rotated_points = points_centered @ rotation_matrix
        transformed_points = rotated_points + center + np.array([dx, dy])
        return transformed_points

def get_field_legends(field:str, sliceList:list, startIndex:float, colorType='tab20', fetchSize=30)-> Tuple[List, bool]:
    """
    Retrieve legends for a specific field (e.g., labels and colors).

    Args:
        field (str): Field name used to generate legends.
        sliceList (list): List of slices.
        startIndex (float): Starting index for retrieving legends.
        colorType (str): Color scheme type, default is 'tab20'.
        fetchSize (int): Number of legend items to fetch per request, default is 30.

    Returns:
        Tuple[List, bool]: A list of legend items and a flag indicating whether more legends are available.
    """
    if not alidata.has_alifig():
        return [], False
    aliFig = alidata.get_alifig()
    ad = alidata.get_adata()
    legend_items = set()
    labels = set(ad.obs[field].unique())
    colorMap = get_color_map(labels, type=colorType)
    for slice in sliceList:
        idx = alidata.get_slice_index(slice)
        obsIndex = aliFig['data'][idx]['customdata'].flatten()
        fieldValues = ad.obs.loc[obsIndex, field]
        legend_items = legend_items.union(fieldValues)
    legend_items = list(legend_items)
    legend_items.sort()
    legendColor = {}
    length = len(legend_items)
    for i in range(startIndex, startIndex+fetchSize):
        if i<length:
            label = legend_items[i]
            legendColor[label] = colorMap[label]
        else:
            break
    hasLegend = True
    if startIndex+fetchSize>=length:
        hasLegend = False
    legends = [get_legend_item(legendColor[key], key) for key in legendColor]
    return legends, hasLegend

def get_legend_item(color:str, label:str)->fac.AntdTooltip:
    """
    Create a single legend item.

    Args:
        color (str): Color associated with the legend item.
        label (str): Label for the legend item.

    Returns:
        AntdTooltip: A component containing the color and label for the legend.
    """
    return fac.AntdTooltip(
        html.Div(
            children=[
                html.Div(
                    style={
                        'width': '10px',
                        'height': '10px',
                        'borderRadius': '50%',
                        'backgroundColor': color,
                        'marginRight': '10px',
                    }
                ),
                html.Div(
                    label, 
                    style={
                        'marginBottom': '4px',
                        'whiteSpace': 'nowrap',
                        'overflow': 'hidden', 
                        'textOverflow': 'ellipsis', 
                        'maxWidth': '150px', 
                    }
                ),
            ],
            style={
                'display': 'flex',
                'alignItems': 'center',
            }
        ), 
        title=fac.AntdText(label, style={'color':'white'}),
        color=color
    )

def reset_slicecolor(patch:Patch, sliceList:list)->None:
    """
    Reset the color of slices to their original state.

    Args:
        patch (dict): Dictionary holding the data for slices and their markers.
        sliceList (list): List of slices to reset the colors for.

    """
    aliFig = alidata.get_alifig()
    for slice in sliceList:
        if slice:
            idx = alidata.get_slice_index(slice)
            patch['data'][idx]['marker']['color'] = aliFig['data'][idx]['marker']['color']

def set_color_bygene(patch:Patch, gene:str, sliceList:list, colorType='C1')->None:
    """
    Set the color of the points in slices based on gene expression.

    Args:
        patch (dict): Dictionary holding the data for slices and their markers.
        gene (str): The gene used to determine the color based on its expression level.
        sliceList (list): List of slices to apply the color changes.
        colorType (str): Type of color scale to use, default is 'C1'.

    """
    aliFig = alidata.get_alifig()
    ad = alidata.get_adata()

    if colorType=='C1':
        colorType='red'
    elif colorType=='C2':
        colorType='green'
    elif colorType=='C3':
        colorType='blue'

    minValue = float(np.min(ad[:, gene].X))
    maxValue = float(np.max(ad[:, gene].X))

    for slice in sliceList:
        idx = alidata.get_slice_index(slice)
        colorList = None
        obsIndex = aliFig['data'][idx]['customdata'].flatten()
        geneExp = ad[obsIndex, gene].X
        colorList = get_scale_colors(geneExp, minValue, maxValue, colorType=colorType)
        if colorList==None:
            patch['data'][idx]['marker']['color'] = aliFig['data'][idx]['marker']['color']
        else:
            patch['data'][idx]['marker']['color'] = colorList

def set_color_byfield(patch:Patch, field:str, sliceList:list, colorType='tab20')->None:
    """
    Set the color of the points in slices based on a specific field.

    Args:
        patch (dict): Dictionary holding the data for slices and their markers.
        field (str): The field based on which colors will be set.
        sliceList (list): List of slices to apply the color changes.
        colorType (str): The color scale to use, default is 'tab20'.
    
    """
    aliFig = alidata.get_alifig()
    ad = alidata.get_adata()
    colorMap = None
    legend_items = set()
    labels = set(ad.obs[field].unique())
    colorMap = get_color_map(labels, type=colorType)
    for slice in sliceList:
        idx = alidata.get_slice_index(slice)
        if colorMap==None:
            patch['data'][idx]['marker']['color'] = aliFig['data'][idx]['marker']['color']
        else:
            obsIndex = aliFig['data'][idx]['customdata'].flatten()
            fieldValues = ad.obs.loc[obsIndex, field]
            legend_items = legend_items.union(fieldValues)
            patch['data'][idx]['marker']['color'] = [colorMap[val] for val in fieldValues]

def get_below_slice(act_slice:float)->float:
    """
    Get the slice that is below the given slice in the list.

    Args:
        act_slice (float): The current active slice.

    Returns:
        float: The slice below the given slice, or None if no such slice exists.
    
    """
    idx = alidata.get_slice_index(act_slice)
    below_slice = alidata.get_slice_list()[idx-1] if idx is not None and idx!=0 else None
    return below_slice

# def align_slices_with_paste(model:str, use_gpu:bool)->None:
#     """
#     Align slices using a given model and apply the alignment to the slices.

#     Args:
#         model (str): The alignment model to use ('paste1' or 'paste2').
#         use_gpu (bool): Whether to use GPU for alignment computation.
    
#     """
#     permission = verify_modify_permission()
#     if not permission:
#         return
#     status = check_alignment_data()
#     if not status:
#         return
#     has_orifig = check_orifig()
#     if not has_orifig:
#         return
#     ad = alidata.get_adata()
#     x = alidata.get_xfield()
#     y = alidata.get_yfield()
#     z = alidata.get_zfield()
#     if model == 'paste1':
#         success = paste1(ad, x, y, z, use_gpu=use_gpu)
#     else:
#         success = paste2(ad, x, y, z)
#     if not success:
#         set_head_notice('Gene expression values cannot be negative. Please check and resubmit your data !', 'warning')
#         return
#     alidata.set_x_aligned_field('x_aligned')
#     alidata.set_y_aligned_field('y_aligned')
#     fig = alidata.plot_3d_scatter()
#     set_props('alignment-graph-aligned', {'figure':fig, 'style':graphStyle})
#     obs_fields = alidata.get_obs_fields()
#     genes = alidata.get_gene_list()
#     set_props('graphSetting-selecter-colorField', {'options':obs_fields, 'value':None})
#     set_props('graphSetting-selecter-colorGene', {'options':genes})
#     slices = alidata.get_slice_list()
#     if slices is not None:
#         set_props('graphSetting-selecter-activeSlice', {'options':slices})
#         set_props('graphSetting-selecter-referenceSlice', {'options':slices})
#         set_props('graphSetting-table-applySlices', {'data':[{'slices':slice} for slice in slices]})

def plot_origin_fig(x:str, y:str, z:str):
    """
    Plot the original 3D figure using the provided x, y, z fields.

    Args:
        x (str): The field used for the x-axis.
        y (str): The field used for the y-axis.
        z (str): The field used for the z-axis.
    """
    permission = verify_modify_permission()
    if not permission:
        return
    status = check_alignment_data()
    if not status:
        return
    if x is None or y is None or z is None:
        set_head_notice('Please select the x, y, z field !', 'warning')
        return
    try:
        alidata.set_xfield(x)
    except Exception as e:
        set_head_notice('The x field is not numeric !', 'warning')
        return
    try:
        alidata.set_yfield(y)
    except Exception as e:
        set_head_notice('The y field is not numeric !', 'warning')
        return
    try:
        alidata.set_zfield(z)
        zfield_min = alidata.get_zfield_min()
        zfield_max = alidata.get_zfield_max()
        if zfield_min is not None and zfield_max is not None:
            props = dict()
            props['min'] = zfield_min
            props['max'] = zfield_max
            props['value'] = [zfield_min, zfield_max]
            props['tooltipPrefix'] = f'{alidata.get_zfield()}: '
            set_props('alignmentGraphSetting-slider-z', props)

    except Exception as e:
        set_head_notice('The z field is not numeric !', 'warning')
        return
    fig = alidata.plot_3d_scatter(origin_fig=True)
    set_props('alignment-graph-origion', {'figure':fig, 'style':graphStyle})

def check_orifig()->bool:
    """
    Check if the original figure exists.

    Returns:
        bool: True if the original figure exists, False otherwise.
    """
    if alidata.get_orifig() is None:
        set_head_notice('The original figure is empty, Please click the Plot Figure button first !', 'warning')
        return False
    return True

def check_alignment_data()->bool:
    """
    Check if the alignment data is available.

    Returns:
        bool: True if alignment data exists, False otherwise.
    """
    status = alidata.get_adata() is not None
    if not status:
        set_head_notice('There was no data exist, Please importing previsously !', 'warning')
    return status
    
def export_alignment_file(path:str)->bool:
    """
    Export the alignment data to a specified file.

    Args:
        path (str): The file path where the data should be exported.

    Returns:
        bool: True if export is successful, False otherwise.
    """
    access_path = path
    if not os.path.exists(path):
        access_path = os.path.dirname(path)
    if not os.access(access_path, os.W_OK):
        set_head_notice('You have no permission to write data to this folder !', 'warning')
        return False
    status = False
    try:
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        alidata.get_adata().write_h5ad(path, compression='gzip')
    except Exception as e:
        set_head_notice('export failed!', 'error')
        raise e
    else:
        status = True
        set_head_notice('export successfully!', 'success')
    return status

def read_alignment_file(path:str)->bool:
    """
    Read alignment data from a specified file.

    Args:
        path (str): The file path from which the alignment data should be loaded.

    Returns:
        bool: True if the data is successfully read, False otherwise.
    """
    alidata.reset_props()
    status = False
    try:
        adata = sc.read_h5ad(path)
        alidata.set_adata(adata)
        obs_fields = alidata.get_obs_fields()
    except Exception as e:
        set_head_notice('import failed!', 'error')
        raise e
    else:
        status = True
        set_props('alignment-select-x', {'options':obs_fields, 'value':None})
        set_props('alignment-select-y', {'options':obs_fields, 'value':None})
        set_props('alignment-select-z', {'options':obs_fields, 'value':None})
        set_props('graphSetting-selecter-colorField', {'options':[], 'value':None})
        set_props('graphSetting-selecter-colorGene', {'options':[], 'value':None})
        set_props('graphSetting-selecter-activeSlice', {'options':[], 'value':None})
        set_props('graphSetting-selecter-referenceSlice', {'options':[], 'value':None})
        set_props('graphSetting-table-applySlices', {'data':[], 'selectedRowKeys':[]})
        set_props('alignment-graph-origion', {'style':grapHidden})
        set_props('alignment-graph-aligned', {'style':grapHidden})
        props = dict()
        props['min'] = 0
        props['max'] = 0
        props['value'] =[]
        props['tooltipPrefix'] = ''
        set_props('alignmentGraphSetting-slider-z', props)
        set_head_notice('import successfully!', 'success')
    return status

def restore_initial_data()->None:
    """
    Restore the initial data and settings to their default values.
    """
    obs_fields = alidata.get_obs_fields()
    if obs_fields:
        set_props('alignment-select-x', {'options':obs_fields, 'value':alidata.get_xfield()})
        set_props('alignment-select-y', {'options':obs_fields, 'value':alidata.get_yfield()})
        set_props('alignment-select-z', {'options':obs_fields, 'value':alidata.get_zfield()})
    orifig = alidata.get_orifig()
    alifig = alidata.get_alifig()
    if orifig:
        set_props('alignment-graph-origion', {'figure':orifig, 'style':graphStyle})
    if alifig:
        set_props('graphSetting-selecter-colorField', {'options':obs_fields})
        set_props('alignment-graph-aligned', {'figure':alifig, 'style':graphStyle})
        genes = alidata.get_gene_list()
        set_props('graphSetting-selecter-colorGene', {'options':genes})
        slices = alidata.get_slice_list()
        if slices is not None:
            set_props('graphSetting-selecter-activeSlice', {'options':slices})
            set_props('graphSetting-selecter-referenceSlice', {'options':slices})
            set_props('graphSetting-table-applySlices', {'data':[{'slices':slice} for slice in slices]})

    zfield_min = alidata.get_zfield_min()
    zfield_max = alidata.get_zfield_max()
    if zfield_min is not None and zfield_max is not None:
        props = dict()
        props['min'] = zfield_min
        props['max'] = zfield_max
        props['value'] = [zfield_min, zfield_max]
        props['tooltipPrefix'] = f'{alidata.get_zfield()}: '
        set_props('alignmentGraphSetting-slider-z', props)

    thread = alidata.get_alistatus_thread()
    if thread is not None:
        creator = alidata.get_alistatus_creator()
        set_props('alignment-popupcard-alignTask', {'visible':True})
        set_props('alignment-interval', {'disabled':False})
        set_props('alignment-button-alignSlices', {'disabled':True})
        set_props('alignment-timeline-creator', {'children':creator})

