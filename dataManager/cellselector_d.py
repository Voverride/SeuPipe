from dataManager.workspace import *
from controller.auth import get_request_usrname
from dataManager.segmentation_d import segData
import shutil
from utils.observer import observer
from websocket.websocket import ws
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
import cachetools
import numpy as np
from utils.typing import StepIcon
from utils.commonfuc import write_json, read_json, read_pkl

class CellSelectorData:
    def __init__(self):
        self._projects = cachetools.TTLCache(maxsize=10, ttl=1800)
        self._positions = cachetools.TTLCache(maxsize=10, ttl=1800)
        self._iswritting = set()
        self._exportProject = dict()
        self._clustering = dict()

    def is_clustering(self, project):
        """
        判断项目是否有聚类任务正在运行
        """
        if project in self._clustering:
            return self._clustering[project].get('running', False)
        return False
    def add_clustering(self, project, metadata):
        """
        添加项目聚类任务
        """
        observered_metadata = observer.observe(metadata, ws.notifyUpdateClusteringStatus, project=project)
        self._clustering[project] = observered_metadata
        return observered_metadata

    def remove_clustering(self, project):
        """
        删除项目聚类任务
        """
        if project in self._clustering:
            observer.disobserve(self._clustering[project], ws.notifyUpdateClusteringStatus)
            del self._clustering[project]

    def set_export_project(self, project):
        """
        设置用户导出任务
        """
        username = get_request_usrname()
        self._exportProject[username] = project
    
    def get_export_project(self):
        """
        获取用户导出任务
        """
        username = get_request_usrname()
        return self._exportProject.get(username, None)

    def set_slice_writting(self, project, slice):
        """
        设置切片正在写入
        """
        self._iswritting.add(f'{project}_{slice}')
    
    def is_slice_writting(self, project, slice):
        """
        判断切片是否正在写入
        """
        return f'{project}_{slice}' in self._iswritting
    
    def remove_slice_writting(self, project, slice):
        """
        删除切片正在写入
        """
        self._iswritting.remove(f'{project}_{slice}')

    def get_position_path(self, project, slice):
        """
        获取细胞位置路径
        """
        slice_folder = self.get_slice_folder(project, slice)
        position_path = os.path.join(slice_folder, 'position.csv')
        return position_path
    def generate_cell_position(self, project, slice, mask):
        """
        生成细胞空间位置，并转换为左下角为原点的笛卡尔坐标系（用于绘图）。
        """
        if issparse(mask):
            rows, cols = mask.nonzero()
            cell_ids = mask.data
        else:
            rows, cols = np.where(mask != 0)
            cell_ids = mask[rows, cols]

        df_coords = pd.DataFrame({
            'cell_id': cell_ids,
            'row': rows,
            'col': cols
        })

        df_coords = df_coords[df_coords['cell_id'] != 0]

        df_coords['x'] = df_coords['col']
        df_coords['y'] = df_coords['row']

        centroids = df_coords.groupby('cell_id').agg(
            x=('x', 'mean'),
            y=('y', 'mean')
        ).reset_index()

        index_labels = [f'z{slice}_{int(cell_id)}' for cell_id in centroids['cell_id']]

        df = pd.DataFrame({
            'x': centroids['x'].values,
            'y': centroids['y'].values
        }, index=index_labels)
        df['selected'] = True
        df['cluster'] = 0
        self.add_cell_position(project, slice, df)
        position_path = self.get_position_path(project, slice)
        df.to_csv(position_path, sep='\t')

        return df
    
    def add_cell_position(self, project, slice, position):
        """
        添加细胞位置
        """
        key = f'{project}_{slice}'            
        self._positions[key] = position

    def get_cell_position(self, project, slice):
        """
        获取细胞位置
        """
        key = f'{project}_{slice}'
        if key in self._positions:
            position = self._positions[key]
            del self._positions[key]
        else:
            position_path = self.get_position_path(project, slice)
            if not os.path.exists(position_path):
                return None
            df = pd.read_csv(position_path, sep='\t', index_col=0)
            position = df
        self._positions[key] = position
        return position
    def get_seg_mask(self, project_name, slice):
        """
        获取分割结果
        """
        project_info = self.get_project_info(project_name)
        result_type = project_info['result_field']
        adata_path = project_info['slices'][str(slice)]['adata_path']
        try:
            adata = sc.read_h5ad(adata_path)
            return adata.layers[result_type]
        except:
            return None
    
    def get_stain_fig(self, project_name, slice):
        """
        获取原图
        """
        stain_path = segData.get_stain_path(project_name, slice)
        fig = read_pkl(stain_path)
        return fig
    
    def get_gem_path(self, project_name, slice):
        """
        获取gem路径
        """
        project_info = self.get_project_info(project_name)
        return project_info['slices'][str(slice)]['gem_path']
    
    def get_project_info(self, project_name: str):
        """
        获取细胞筛选项目信息
        """
        project_info = dict(self._clustering.get(project_name, {}))
        if project_info: 
            return project_info
        elif project_name in self._projects:
            project_info = self._projects[project_name]
            del self._projects[project_name]
        else:
            project_folder = self.get_project_folder(project_name)
            metadata_path = os.path.join(project_folder, 'metadata.json')
            project_info = read_json(metadata_path)
        self._projects[project_name] = project_info
        return project_info
    
    def get_project_slices(self, project_name: str):
        """
        获取细胞筛选项目切片列表
        """
        project_info = self.get_project_info(project_name)
        slices = project_info['slices'].keys()
        slices = [int(slice) for slice in slices]
        slices.sort()
        return slices
    def get_project_folder(self, project_name: str):
        """
        获取细胞筛选项目路径
        """
        sel_path = get_cellselector_workspace()
        project_path = os.path.join(sel_path, project_name)
        return project_path
    def delete_project(self, project: str):
        """
        删除细胞筛选项目
        """
        project_folder = self.get_project_folder(project)
        if os.path.exists(project_folder):
            shutil.rmtree(project_folder)

    def get_slice_folder(self, project: str, slice_index):
        """
        获取细胞筛选切片目录
        """
        project_folder = self.get_project_folder(project)
        slice_folder = os.path.join(project_folder, f'z_{slice_index}')
        check_path(slice_folder)
        return slice_folder

    def create_project(self, seg_project: str, result: str, resolution: float, iteration: int):
        """
        创建细胞筛选项目
        """
        project_name = f'{seg_project}_{result}'
        project_folder = self.get_project_folder(project_name)
        check_path(project_folder)
        slices_metadata = {}
        segProjectInfo = segData.get_project_info(seg_project)
        for slice in segProjectInfo['slices']:
            slice_index = slice['z']
            slice_folder = os.path.join(project_folder, f'z_{slice_index}')
            check_path(slice_folder)
            stain_path = slice['img']
            gem_path = slice['gem']
            adata_path = segData.get_result_path(seg_project, slice_index)
            slices_metadata[slice_index] = {
                'stain_path': stain_path,
                'gem_path': gem_path,
                'adata_path': adata_path
            }
        steps = {
            'step1': StepIcon.SCHEDULE,
            'step2': StepIcon.SCHEDULE,
            'step3': StepIcon.SCHEDULE,
            'step4': StepIcon.SCHEDULE,
            'step5': StepIcon.SCHEDULE,
            'percent': 0,
        }
        metadata = {
            'creator': get_request_usrname(),
            'project_name': seg_project,
            'result_field': f'{result}_mask',
            'resolution': resolution,
            'iteration': iteration,
            'exception': None,
            'steps': steps,
            'slices': slices_metadata,
        }
        write_json(os.path.join(project_folder, 'metadata.json'), metadata)         
        return project_name, metadata

    def get_exist_projects(self):
        """
        获取细胞筛选项目列表
        """
        sel_path = get_cellselector_workspace()
        projects = os.listdir(sel_path)
        projects = [p for p in projects if os.path.isdir(os.path.join(sel_path, p))]
        projects.sort()
        return projects
        
    

selData = CellSelectorData()