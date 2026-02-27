from dataManager.workspace import *
from utils.commonfuc import *
from utils.typing import StepStatus
from websocket.websocket import ws
from utils.observer import observer
from controller.auth import get_request_usrname
import scanpy as sc
import cachetools
import shutil
import bisect

class AlignmentData:
    def __init__(self):
        self._alignment = {}
        self._projects = cachetools.TTLCache(maxsize=20, ttl=1800)
        self._export_project = cachetools.TTLCache(maxsize=200, ttl=1800)

    def get_syncslice_index(self, project, slice):
        """
        获取同步切片索引
        """
        znums = self.get_znums(project)
        idx = bisect.bisect_left(znums, slice)
        return len(znums) - idx -1
    
    def get_znums(self, projectname):
        """
        获取切片列表
        """
        metadata = self.get_project_info(projectname)
        znums = metadata.get('znums', [])
        return znums

    def create_alignment_task(self, projectname):
        """
        创建对齐任务
        """
        metadata = self.get_project_info(projectname)
        metadata['steps'] = {
            'preprocess': StepStatus.WAIT,
            'alignment': StepStatus.WAIT,
            'postprocess': StepStatus.WAIT,
            'percent': 0,
        }
        metadata['exception'] = None
        observed_metadata = observer.observe(metadata, ws.notifyUpdateAlignmentStatus, project=projectname)
        self._alignment[projectname] = observed_metadata
        return observed_metadata

    def remove_alignment_task(self, projectname):
        """
        移除对齐任务
        """
        if projectname in self._alignment:
            observer.disobserve(self._alignment[projectname], ws.notifyUpdateAlignmentStatus)
            del self._alignment[projectname]
    
    def is_running(self, projectname):
        """
        检查是否存在对齐任务
        """
        if projectname not in self._alignment:
            return False
        metadata = self.get_project_info(projectname)
        return metadata.get('running', False)
    def get_exist_projects(self):
        """
        获取存在的对齐项目列表
        """
        ali_path = get_alignment_workspace()
        projects = os.listdir(ali_path)
        projects = [p for p in projects if os.path.isdir(os.path.join(ali_path, p))]
        projects.sort()
        return projects

    def get_project_info(self, project_name: str):
        """
        获取对齐项目信息
        """
        project_info = dict(self._alignment.get(project_name, {}))
        if project_info: 
            return project_info
        elif project_name in self._projects:
            project_info = self._projects[project_name]
            del self._projects[project_name]
        else:
            metadata_path = self.get_metadata_path(project_name)
            project_info = read_json(metadata_path)
        self._projects[project_name] = project_info
        return project_info

    def get_metadata_path(self, project_name: str):
        """
        获取对齐项目元数据路径
        """
        project_folder = self.get_project_folder(project_name)
        metadata_path = os.path.join(project_folder, 'metadata.json')
        return metadata_path

    def get_project_folder(self, project_name: str):
        """
        获取对齐项目文件夹路径
        """
        ali_path = get_alignment_workspace()
        project_path = os.path.join(ali_path, project_name)
        check_path(project_path)
        return project_path
    
    def delete_project(self, project: str):
        """
        删除对齐选项目
        """
        project_folder = self.get_project_folder(project)
        if project in self._projects:
            del self._projects[project]
        if os.path.exists(project_folder):
            shutil.rmtree(project_folder)

    def get_initial_folder(self, project):
        """
        获取初始数据文件夹
        """
        project_folder = self.get_project_folder(project)
        initial_folder = os.path.join(project_folder, 'initial')
        check_path(initial_folder)
        return initial_folder
    
    def get_initialdata_path(self, project):
        """
        获取原始数据路径
        """
        initial_folder = self.get_initial_folder(project)
        initial_path = os.path.join(initial_folder, 'data.h5ad')
        return initial_path
    
    def get_coordinate_folder(self, project):
        """
        获取对齐项目坐标文件夹路径
        """
        project_folder = self.get_project_folder(project)
        coordinate_folder = os.path.join(project_folder, 'coordinate')
        check_path(coordinate_folder)
        return coordinate_folder
    
    def get_coordinate_path(self, project):
        """
        获取对齐项目坐标路径
        """
        coordinate_folder = self.get_coordinate_folder(project)
        coordinate_path = os.path.join(coordinate_folder, 'coordinate.csv')
        return coordinate_path
    
    def update_coordinate(self, project, coordinate):
        """
        更新对齐项目坐标
        """
        coordinate_path = self.get_coordinate_path(project)
        coordinate.to_csv(coordinate_path)
    
    def get_figure_folder(self, project):
        """
        获取对齐项目可视化文件夹路径
        """
        project_folder = self.get_project_folder(project)
        figure_folder = os.path.join(project_folder, 'figure')
        check_path(figure_folder)
        return figure_folder
    
    def get_resultfig_path(self, project):
        """
        获取对齐项目可视化结果路径
        """
        figure_folder = self.get_figure_folder(project)
        resultfig_path = os.path.join(figure_folder, 'resultfig.pkl')
        return resultfig_path
    
    def update_resultfig(self, project, coordinate, x, y, z):
        """
        更新对齐项目可视化结果
        """
        resultfig_path = self.get_resultfig_path(project)
        resultfig = self.plot_3d_scatter(project, coordinate, x, y, z)
        write_pkl(resultfig, resultfig_path)
    
    def get_resultfig(self, project):
        """
        获取对齐项目可视化结果
        """
        resultfig_path = self.get_resultfig_path(project)
        if os.path.exists(resultfig_path):
            try:
                resultfig = read_pkl(resultfig_path)
            except:
                resultfig = None
        else:
            resultfig = None
        return resultfig
    
    def update_initialfig(self, project, coordinate, x, y, z):
        """
        更新初始数据可视化结果
        """
        initialfig_path = self.get_initialfig_path(project)
        initialfig = self.plot_3d_scatter(project, coordinate, x, y, z, initial=True)
        write_pkl(initialfig, initialfig_path)
        return initialfig
    
    def get_initialfig_path(self, project):
        """
        获取初始数据可视化结果路径
        """
        figure_folder = self.get_figure_folder(project)
        initialfig_path = os.path.join(figure_folder, 'initialfig.pkl')
        return initialfig_path
    
    def get_initialfig(self, project):
        """
        获取初始数据可视化结果
        """
        initialfig_path = self.get_initialfig_path(project)
        if os.path.exists(initialfig_path):
            try:
                initialfig = read_pkl(initialfig_path)
            except:
                initialfig = None
        else:
            initialfig = None
        return initialfig
    
    def update_project_info(self, project_name: str, metadata: dict):
        """
        更新对齐项目信息
        """
        metadata_path = self.get_metadata_path(project_name)
        write_json(metadata_path, metadata)

    def create_project(self, importdata, x, y, z, projectname):
        """
        创建对齐项目
        """
        initialdata_path = self.get_initialdata_path(projectname)
        shutil.copy(importdata, initialdata_path)
        data = sc.read_h5ad(initialdata_path, backed='r')
        coordinate = data.obs.copy()
        data.file.close()
        self.update_coordinate(projectname, coordinate)
        zmin = coordinate[z].min()
        zmax = coordinate[z].max()
        znums = list(set(coordinate[z].unique()))
        znums.sort()
        znum = len(znums)
        x_min = np.min(coordinate[x])
        x_max = np.max(coordinate[x])
        y_min = np.min(coordinate[y])
        y_max = np.max(coordinate[y])
        min_val = min(x_min, y_min)
        max_val = max(x_max, y_max)
        dist = max_val-min_val
        max_val += dist/2
        min_val -= dist/2
        iniRange = max_val-min_val
        x_scale = iniRange/(x_max-x_min)
        y_scale = iniRange/(y_max-y_min)
        iniScale = (x_scale+y_scale)/2
        metadata = {
            'creator': get_request_usrname(),
            'projectName': projectname,
            'date': get_current_date(),
            'exception': None,
            'running': False,
            'initialData': importdata,
            'x': x,
            'y': y,
            'z': z,
            'zmin': zmin,
            'zmax': zmax,
            'znum': znum,
            'znums': znums,
            'iniRange': iniRange,
            'iniScale': iniScale,
            'steps': {
                'preprocess': StepStatus.WAIT,
                'alignment': StepStatus.WAIT,
                'postprocess': StepStatus.WAIT,
                'percent': 0,
            }
        }
        metadata_path = self.get_metadata_path(projectname)
        write_json(metadata_path, metadata)
        self.update_resultfig(projectname, coordinate, x, y, z)

    def get_graph_view(self, actSlice:float, refSlice:float)->dict:
        """
        获取图像视角
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

    def plot_3d_scatter(self, project, coordinate, x, y, z, initial=False, marker_size:int=5, boarder_width:int=1, boarder_color:str='#0d0015', grid_color:str = '#5F9EA0'):
        """
        绘制3D散点图
        """
        df = coordinate.copy()
        df = df.sort_values(by=z)
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
        if initial:
            step = dist/20
            x_scale = dist/(x_max-x_min)
            y_scale = dist/(y_max-y_min)
            scale = (x_scale+y_scale)/2
        else:
            metadata = self.get_project_info(project)
            iniRange = metadata['iniRange']
            iniScale = metadata['iniScale']
            multy_factor = dist/iniRange
            step = dist/20/multy_factor
            scale = multy_factor*iniScale

        df[z] = df[z].astype(str)
        df['index'] = df.index
        cmap = get_color_map(set(df[z]), type='COLORS_60')
        fig = px.scatter_3d(
            df, x=x, y=y, z=z, 
            color=z, 
            custom_data='index',
            color_discrete_map=cmap,
        )
        fig.update_traces(
            marker=dict(
                size=marker_size,
                line=dict(color=boarder_color, width=boarder_width)
            )
        )
        fig.update_layout(
            autosize=True,
            uirevision=True,
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
        return fig

alidata = AlignmentData()