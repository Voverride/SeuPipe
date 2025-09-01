from dataManager.workspace import *
from controller.auth import get_host
from dataManager.users import search_user
import scanpy as sc
import re
import plotly.express as px
from dash import Patch
import shutil
import pickle
import plotly.graph_objects as go
from utils.commonfuc import *

class MaskViewerData:
    def __init__(self):
        self._userdata = {} # {usrname: Anndata}
    def get_contrast_mask_contour_figure(self, showMask=True, showContour=False):
        """
        获取两种算法分割结果对比图
        """
        ad = self.get_user_adata()
        stain = ad.layers['stain']
        cellpose_mask = ad.uns['cellpose_mask']
        cellpose_contour = ad.uns['cellpose_contour']
        watershed_mask = ad.uns['watershed_mask']
        watershed_contour = ad.uns['watershed_contour']
        cellpose_fig = self.get_mask_contour_figure(stain, cellpose_mask, cellpose_contour, title='Cellpose Mask', showcontour=showContour, showmask=showMask)
        watershed_fig = self.get_mask_contour_figure(stain, watershed_mask, watershed_contour, title='Watershed Mask', showcontour=showContour, showmask=showMask)
        return cellpose_fig, watershed_fig
    def get_mask_contour_figure(self, stain, mask, contour, title=None, showmask=True, showcontour=False, mask_opacity=0.7):
        """
        获取stain， mask，contour， 叠加交互图
        """
        if stain.ndim == 2:
            stain_rgb = np.stack([stain]*3, axis=-1)
        else:
            stain_rgb = stain[..., :3]
        
        fig = get_imgfig_withplotly(stain_rgb, title=title)

        fig.add_trace(go.Image(
            source=f"data:image/png;base64,{array_to_base64(mask)}",
            opacity=mask_opacity,
            visible=showmask,
        ))

        fig.add_trace(go.Image(
            source=f"data:image/png;base64,{array_to_base64(contour)}",
            visible=showcontour,
        ))
        return fig
    
    def get_user_adata(self):
        usrname = self.get_request_usrname()
        if usrname in self._userdata:
            return self._userdata[usrname]
        else:
            return None

    def get_registration_figure(self):
        """
        获取配准前与配准后图像
        """
        ad = self.get_user_adata()
        before_fig = get_imgfig_withplotly(ad.uns['before'], title='Before Aligned')
        aligned_fig = get_imgfig_withplotly(ad.uns['aligned'], title='After Aligned')
        return before_fig, aligned_fig
    def read_slice_adata(self, taskname, slice):
        """
        读取切片h5ad
        """
        seg_path = get_segmentation_workspace()
        task_path = os.path.join(seg_path, taskname)
        slice_path = os.path.join(task_path, 'slices', str(slice)+'.h5ad')
        adata_path = sc.read_h5ad(slice_path)
        usrname = self.get_request_usrname()
        self._userdata[usrname] = adata_path

    def get_request_usrname(self):
        """
        获取请求的用户名
        """
        host = get_host()
        user = search_user(usrhost=host)
        if len(user)==0:
            return None
        username = user[0]['usrname']
        return username
    def get_task_slices(self, taskname):
        """
        获取任务切片列表
        """
        seg_path = get_segmentation_workspace()
        task_path = os.path.join(seg_path, taskname)
        if os.path.exists(task_path):
            slice_path = os.path.join(task_path, 'slices')
            if os.path.exists(slice_path):
                slices = [s.split(".")[0] for s in os.listdir(slice_path)]
                slices = sorted(slices, key=lambda x: int(re.search(r'\d+', x).group()))
                return slices
        return []

    def get_exist_tasks(self):
        """
        获取已创建的任务列表
        """
        seg_path = get_segmentation_workspace()
        tasks = os.listdir(seg_path)
        tasks = [task for task in tasks if os.path.isdir(os.path.join(seg_path, task))]
        tasks.sort()
        return tasks

maskData = MaskViewerData()