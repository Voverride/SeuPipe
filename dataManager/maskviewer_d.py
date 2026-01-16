from dataManager.workspace import *
from controller.auth import get_request_usrname
from dash import no_update
from utils.commonfuc import *

class MaskViewerData:
    def __init__(self):
        self.export_project = {}

    def get_project_folder(self, project):
        """
        获取项目文件夹
        """
        workspace = get_maskviewer_workspace()
        project_folder = os.path.join(workspace, project)
        check_path(project_folder)
        return project_folder
    
    def get_slice_folder(self, project, slice):
        """
        获取切片文件夹
        """
        project_folder = self.get_project_folder(project)
        slice_folder = os.path.join(project_folder, f'z_{slice}')
        check_path(slice_folder)
        return slice_folder
    
    def get_overlay_path(self, project, slice):
        """
        获取叠加图像路径
        """
        slice_folder = self.get_slice_folder(project, slice)
        overlay_path = os.path.join(slice_folder, 'overlay.pkl')
        return overlay_path
    
    def get_expansion_path(self, project, slice):
        """
        获取expansion图像路径
        """
        slice_folder = self.get_slice_folder(project, slice)
        cellpose_path = os.path.join(slice_folder, 'expansion.pkl')
        return cellpose_path
    
    def get_cellpose_path(self, project, slice):
        """
        获取cellpose图像路径
        """
        slice_folder = self.get_slice_folder(project, slice)
        cellpose_path = os.path.join(slice_folder, 'cellpose.pkl')
        return cellpose_path
    
    def get_watershed_path(self, project, slice):
        """
        获取watershed图像路径
        """
        slice_folder = self.get_slice_folder(project, slice)
        watershed_path = os.path.join(slice_folder, 'watershed.pkl')
        return watershed_path
    
    def set_export_project(self, project, type):
        """
        设置用户导出任务
        """
        username = get_request_usrname()
        self.export_project[username] = {
            'project': project,
            'type': type
        }

    def get_export_project(self):
        """
        获取当前用户导出项目
        """
        username = get_request_usrname()
        if username in self.export_project:
            return self.export_project[username]
        return {}
    
    def get_graph(self, project, slice, figureType, showMask=True, showContour=False):
        """
        获取分割图像
        """
        if project and slice and figureType:
            slice_folder = self.get_slice_folder(project, slice)
            figure_path = os.path.join(slice_folder, f'{figureType}.pkl')
            fig = read_pkl(figure_path)
            if figureType!='overlay':
                fig['data'][1]['visible'] = showMask
                fig['data'][2]['visible'] = showContour
            return fig
        return no_update
    
    def get_figure_types(self, project, slice):
        """
        获取切片的图像类型
        """
        slice_folder = self.get_slice_folder(project, slice)
        files = os.listdir(slice_folder)
        figure_types = set()
        for fname in files:
            if '.pkl' in fname:
                figure_type = fname.split(".pkl")[0]
                figure_types.add(figure_type)
        remove_overlay = False
        if 'overlay' in figure_types:
            figure_types.remove('overlay')
            remove_overlay = True
        figure_types = list(figure_types)
        figure_types.sort()
        if remove_overlay:
            figure_types.insert(0, 'overlay')
        return figure_types

    def get_project_slices(self, project):
        """
        获取项目切片列表
        """
        project_folder = self.get_project_folder(project)
        if os.path.exists(project_folder):
            files = os.listdir(project_folder)
            slices = []
            for file in files:
                if file.startswith('z_') and os.path.isdir(os.path.join(project_folder, file)):
                    slice = int(file.split("_")[1])
                    slices.append(slice)
            slices.sort()
            return slices
        return []

    def get_exist_projects(self):
        """
        获取已创建的任务列表
        """
        maskviewer_folder = get_maskviewer_workspace()
        projects = os.listdir(maskviewer_folder)
        projects = [project for project in projects if os.path.isdir(os.path.join(maskviewer_folder, project))]
        projects.sort()
        return projects

maskData = MaskViewerData()