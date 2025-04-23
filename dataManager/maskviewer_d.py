from dataManager.workspace import *
from controller.auth import get_host
from dataManager.users import search_user
import scanpy as sc
import re
import plotly.express as px
from dash import Patch
import shutil
import pickle

class MaskViewerData:
    def __init__(self):
        self._userdata = {} # {usrname: Anndata}

    def get_registration_figure(regmtx):
        """
        获取图像与基因表达配准图像
        """
        fig = px.imshow(regmtx)
        fig.update_layout(
            autosize=True,
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
        )
        return fig
    def read_slice_adata(self, taskname, z):
        """
        读取切片h5ad
        """
        seg_path = get_segmentation_workspace()
        task_path = os.path.join(seg_path, taskname)
        slice_path = os.path.join(task_path, 'slices', str(z)+'.h5ad')
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