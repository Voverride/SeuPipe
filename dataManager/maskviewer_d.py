from dataManager.workspace import *
from controller.auth import get_host
from dataManager.users import search_user
import scanpy as sc
import re
import plotly.express as px
from dash import Patch
import shutil
from dash import no_update
import pickle
from dataManager.segmentation_d import segData
import plotly.graph_objects as go
from utils.commonfuc import *

class MaskViewerData:
    def __init__(self):
        self.export_task = {} # key:username, value:{taskname:'', type:''}

    def set_export_task(self, taskname, type):
        """
        设置用户导出任务
        """
        username = self.get_username()
        self.export_task[username] = {
            'taskname': taskname,
            'type': type
        }

    def get_export_task(self):
        """
        获取当前用户导出任务
        """
        username = self.get_username()
        if username in self.export_task:
            return self.export_task[username]
        return {}

    def get_username(self):
        """
        获取当前用户名
        """
        userInfo = search_user(usrhost = get_host())[0]
        return userInfo['usrname']
    
    def get_graph(self, task, slice, graphName, showMask=True, showContour=False):
        """
        获取分割图像
        """
        if task and slice and graphName:
            slice_index = slice.split("_")[1]
            stain_path = segData.get_seg_stain_figure_path(task, slice_index)
            figure_folder = segData.get_seg_figure_folder(task, slice_index)
            mask_path = os.path.join(figure_folder, f'{graphName}_mask.pkl')
            contour_path = os.path.join(figure_folder, f'{graphName}_contour.pkl')
            if os.path.exists(stain_path) and os.path.exists(contour_path) and os.path.exists(mask_path):
                mask = read_pkl(mask_path)
                contour = read_pkl(contour_path)
                stain = read_pkl(stain_path)
                return self.get_mask_contour_figure(stain, mask, contour, title=graphName, showmask=showMask, showcontour=showContour)
            return no_update

        return no_update
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
    
    def get_figure_types(self, taskname, slice):
        """
        获取任务切片的图像类型
        """
        slice_index = slice.split("_")[1]
        figure_folder = segData.get_seg_figure_folder(taskname, slice_index)
        figure_types = set()
        for fname in os.listdir(figure_folder):
            if '_' in fname:
                figure_type = fname.split("_")[0]
                figure_types.add(figure_type)
        figure_types = list(figure_types)
        figure_types.sort()
        return figure_types

    def get_registration_figure(self, taskname, slice):
        """
        获取配准前与配准后图像
        """
        slice_index = slice.split("_")[1]
        before_path = segData.get_registration_before_path(taskname, slice_index)
        aligned_path = segData.get_registration_aligned_path(taskname, slice_index)
        before = read_pkl(before_path)
        aligned = read_pkl(aligned_path)
        before_fig = get_imgfig_withplotly(before, title='Before Aligned')
        aligned_fig = get_imgfig_withplotly(aligned, title='After Aligned')
        return before_fig, aligned_fig

    def get_task_slices(self, taskname):
        """
        获取任务切片列表
        """
        slice_folder = segData.get_seg_slice_folder(taskname)
        if os.path.exists(slice_folder):
            slices = os.listdir(slice_folder)
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