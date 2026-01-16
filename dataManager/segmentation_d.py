from dataManager.workspace import *
from controller.auth import get_request_usrname
from utils.commonfuc import read_json, write_json
from websocket.websocket import ws
from utils.observer import observer
import shutil

class SegmentationData:
    def __init__(self):
        self._project_metadata = {}
        self._running_projects = {}

    def get_result_path(self, project, z):
        """
        获取分割结果路径
        """
        slice_folder = self.get_slice_folder(project, z)
        result_path = os.path.join(slice_folder, f'z_{z}.h5ad')
        return result_path
    
    def get_slice_folder(self, project, z):
        """
        获取分割结果切片文件夹
        """
        seg_workspace = get_segmentation_workspace()
        slice_folder = os.path.join(seg_workspace, project, f'z_{z}')
        check_path(slice_folder)
        return slice_folder
    
    def add_running_project(self, projectName: str, project_info: dict):
        """
        添加正在运行的项目
        """
        project_info['progress'] = 0
        project_info['exception'] = None
        slices = project_info['slices']
        default_status = {
            "status": "warning",
            "text": "waiting"
        }
        for slice in slices:
            slice['segmentation'] = default_status.copy()
            slice['postprocess'] = default_status.copy()  
        observed_info = observer.observe(project_info, ws.notifyUpdateSegmentationStatus, project=projectName)
        self._running_projects[projectName] = observed_info
        return observed_info

    def remove_running_project(self, projectName: str):
        """
        删除正在运行的项目
        """
        if projectName in self._running_projects:
            observer.disobserve(self._running_projects[projectName], ws.notifyUpdateSegmentationStatus)
            del self._running_projects[projectName]

    def has_running_project(self, projectName: str):
        """
        检查项目是否正在运行
        """
        if projectName in self._running_projects:
            project_info = self._running_projects[projectName]
            if project_info.get('running', False):
                return True
        return False
    def get_project_info(self, project_name: str):
        """
        获取项目信息
        """
        if project_name in self._running_projects:
            project_info = self._running_projects[project_name]
            return dict(project_info)
        project_folder = self.get_project_folder(project_name)
        metadata_path = os.path.join(project_folder, 'metadata.json')
        if not os.path.exists(metadata_path):
            return None
        project_info = read_json(metadata_path)
        return project_info
    
    def get_project_folder(self, project_name: str):
        """
        获取项目文件夹路径
        """
        return os.path.join(get_segmentation_workspace(), project_name)

    def delete_project(self, projectName):
        """
        从磁盘删除项目
        """
        project_folder = self.get_project_folder(projectName)
        if os.path.exists(project_folder):
            shutil.rmtree(project_folder)
        result_fig_folder = os.path.join(get_maskviewer_workspace(), projectName)
        if os.path.exists(result_fig_folder):
            shutil.rmtree(result_fig_folder)

    def write_metadata_to_disk(self, project_name: str, metadata: dict):
        """
        将项目元数据写入磁盘
        """
        project_folder = self.get_project_folder(project_name)
        check_path(project_folder)
        output_path = os.path.join(project_folder, 'metadata.json')
        write_json(output_path, metadata)

    def save_project_metadata(self, metadata:dict):
        """
        保存项目元数据
        """
        username = get_request_usrname()
        project_metadata = self._project_metadata[username]
        metadata.update(project_metadata)
        project_name = metadata['name']
        self.write_metadata_to_disk(project_name, metadata)
        del self._project_metadata[username]

    def get_exist_projects(self):
        """
        获取已创建的项目列表
        """
        seg_path = get_segmentation_workspace()
        projects = os.listdir(seg_path)
        projects.sort()
        return [project for project in projects if os.path.isdir(os.path.join(seg_path, project))]

    def reset_project_metadata(self):
        """
        重置用户创建的项目元数据
        """
        username = get_request_usrname()
        if username in self._project_metadata:
            self._project_metadata[username]['slices'] = None

    def set_project_metadata(self, slices:list):
        """
        设置用户创建的项目元数据
        """
        username = get_request_usrname()
        if username not in self._project_metadata:
            self._project_metadata[username] = {}
        self._project_metadata[username]['slices'] = slices
    
    def get_stain_path(self, project: str, slice: str):
        """
        获取 stain 图像路径
        """
        slice_folder = self.get_slice_folder(project, slice)
        stain_path = os.path.join(slice_folder, 'stain.pkl')
        return stain_path

segData = SegmentationData()
