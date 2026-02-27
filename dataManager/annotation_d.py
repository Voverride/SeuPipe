from dataManager.workspace import *
from utils.commonfuc import *
from utils.typing import StepStatus
from controller.auth import get_request_usrname
from websocket.websocket import ws
from utils.observer import observer
import cachetools
import shutil

class AnnotationData:
    def __init__(self):
        self._annotation = {}
        self._projects = cachetools.TTLCache(maxsize=20, ttl=1800)
        self._export_project = cachetools.TTLCache(maxsize=200, ttl=1800)

    def set_export_project(self, projectname):
        """
        设置用户导出项目
        """
        username = get_request_usrname()
        self._export_project[username] = projectname

    def get_export_project(self):
        """
        获取用户导出项目
        """
        username = get_request_usrname()
        return self._export_project.get(username, None)

    def create_annotation_task(self, projectname):
        """
        创建细胞注释任务
        """
        metadata = self.get_project_info(projectname)
        metadata['steps'] = {
            'preprocess': StepStatus.WAIT,
            'training': StepStatus.WAIT,
            'postprocess': StepStatus.WAIT,
            'percent': 0,
        }
        metadata['exception'] = None
        observed_metadata = observer.observe(metadata, ws.notifyUpdateAnnotationStatus, project=projectname)
        self._annotation[projectname] = observed_metadata
        return observed_metadata
    
    def remove_annotation_task(self, projectname):
        """
        移除细胞注释任务
        """
        if projectname in self._annotation:
            observer.disobserve(self._annotation[projectname], ws.notifyUpdateAnnotationStatus)
            del self._annotation[projectname]
    
    def is_running(self, projectname):
        """
        检查是否存在细胞注释任务
        """
        metadata = self.get_project_info(projectname)
        return metadata.get('running', False)

    def get_project_info(self, project_name: str):
        """
        获取细胞注释项目信息
        """
        project_info = dict(self._annotation.get(project_name, {}))
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
    
    def update_project_info(self, project_name: str, metadata: dict):
        """
        更新细胞注释项目信息
        """
        metadata_path = self.get_metadata_path(project_name)
        write_json(metadata_path, metadata)
    
    def get_metadata_path(self, project_name: str):
        """
        获取细胞注释项目元数据路径
        """
        project_folder = self.get_project_folder(project_name)
        metadata_path = os.path.join(project_folder, 'metadata.json')
        return metadata_path
    
    def get_result_folder(self, project_name: str):
        """
        获取细胞注释项目结果文件夹路径
        """
        project_folder = self.get_project_folder(project_name)
        result_folder = os.path.join(project_folder, 'result')
        check_path(result_folder)
        return result_folder
    
    def set_parameter(self, project_name: str, parameter: dict):
        """
        设置注释结果参数
        """
        parameter_path = self.get_parameter_path(project_name)
        write_json(parameter_path, parameter)
    
    def get_parameter_path(self, project_name: str):
        """
        获取注释结果参数路径
        """
        result_folder = self.get_result_folder(project_name)
        parameter_path = os.path.join(result_folder, 'parameter.json')
        return parameter_path
    
    def get_parameter(self, project_name: str):
        """
        获取注释结果参数
        """
        parameter_path = self.get_parameter_path(project_name)
        if not os.path.exists(parameter_path):
            return None
        parameter = read_json(parameter_path)
        return parameter
    
    def get_classes(self, project_name: str):
        """
        获取注释结果类别数
        """
        parameter = self.get_parameter(project_name)
        if parameter is None:
            return None
        classes = parameter.get('classes', None)
        return classes
    
    def get_z_range(self, project_name: str):
        """
        获取注释结果z轴范围
        """
        parameter = self.get_parameter(project_name)
        if parameter is None:
            return None, None
        z_min = parameter.get('z_min', None)
        z_max = parameter.get('z_max', None)
        return z_min, z_max
    
    def get_diffgene_fig_path(self, project_name: str):
        """
        获取细胞注释差异基因表达图路径
        """
        result_folder = self.get_result_folder(project_name)
        diffgene_fig_path = os.path.join(result_folder, 'diffgene_fig.pkl')
        return diffgene_fig_path
    
    def get_diffgene_fig(self, project_name: str):
        """
        获取细胞注释差异基因表达图
        """
        diffgene_fig_path = self.get_diffgene_fig_path(project_name)
        if not os.path.exists(diffgene_fig_path):
            return None
        diffgene_fig = read_pkl(diffgene_fig_path)
        return diffgene_fig
    
    def set_diffgene_fig(self, project_name: str, diffgene_fig):
        """
        设置细胞注释差异基因表达图
        """
        diffgene_fig_path = self.get_diffgene_fig_path(project_name)
        write_pkl(diffgene_fig, diffgene_fig_path)

    def get_annotation_fig_path(self, project_name: str):
        """
        获取细胞注释图路径
        """
        result_folder = self.get_result_folder(project_name)
        annotation_fig_path = os.path.join(result_folder, 'annotation_fig.pkl')
        return annotation_fig_path
    
    def get_annotation_fig(self, project_name: str):
        """
        获取细胞注释图
        """
        annotation_fig_path = self.get_annotation_fig_path(project_name)
        if not os.path.exists(annotation_fig_path):
            return None
        annotation_fig = read_pkl(annotation_fig_path)
        return annotation_fig
    
    def set_annotation_fig(self, project_name: str, annotation_fig):
        """
        设置细胞注释图
        """
        annotation_fig_path = self.get_annotation_fig_path(project_name)
        write_pkl(annotation_fig, annotation_fig_path)
    
    def get_initial_folder(self, project_name: str):
        """
        获取细胞注释项目原始数据文件夹路径
        """
        project_folder = self.get_project_folder(project_name)
        initial_folder = os.path.join(project_folder, 'initial')
        check_path(initial_folder)
        return initial_folder
    
    def get_refdata_path(self, project_name: str):
        """
        获取细胞注释项目参考数据路径
        """
        initial_folder = self.get_initial_folder(project_name)
        refdata_path = os.path.join(initial_folder, 'ref.h5ad')
        return refdata_path
    
    def get_querydata_path(self, project_name: str):
        """
        获取细胞注释项目查询数据路径
        """
        initial_folder = self.get_initial_folder(project_name)
        querydata_path = os.path.join(initial_folder, 'query.h5ad')
        return querydata_path

    def create_project(self, refdata, querydata, label, x, y, z, rm_mt, rm_ribo, rm_hb, use_hvg, n_layers, n_hiddens, n_latent, epochs, batch_size, dropout, projectname):
        """
        创建细胞注释项目
        """
        project_folder = self.get_project_folder(projectname)
        check_path(project_folder)
        refdata_path = self.get_refdata_path(projectname)
        querydata_path = self.get_querydata_path(projectname)
        shutil.copy(refdata, refdata_path)
        shutil.copy(querydata, querydata_path)
        metadata = {
            'creator': get_request_usrname(),
            'project_name': projectname,
            'date': get_current_date(),
            'exception': None,
            'refdata': refdata,
            'querydata': querydata,
            'label': label,
            'x': x,
            'y': y,
            'z': z,
            'rm_mt': rm_mt,
            'rm_ribo': rm_ribo,
            'rm_hb': rm_hb,
            'use_hvg': use_hvg,
            'n_layers': n_layers,
            'n_hiddens': n_hiddens,
            'n_latent': n_latent,
            'epochs': epochs,
            'batch_size': batch_size,
            'dropout': dropout,
            'steps': {
                'preprocess': StepStatus.WAIT,
                'training': StepStatus.WAIT,
                'postprocess': StepStatus.WAIT,
                'percent': 0,
            }
        }
        metadata_path = self.get_metadata_path(projectname)
        write_json(metadata_path, metadata)
    
    def get_project_folder(self, project_name: str):
        """
        获取注释选项目路径
        """
        ann_path = get_annotation_workspace()
        project_path = os.path.join(ann_path, project_name)
        return project_path
    
    def delete_project(self, project: str):
        """
        删除注释选项目
        """
        if project in self._projects:
            del self._projects[project]
        project_folder = self.get_project_folder(project)
        if os.path.exists(project_folder):
            shutil.rmtree(project_folder)

    def get_exist_projects(self):
        """
        获取注释选项目列表
        """
        ann_path = get_annotation_workspace()
        projects = os.listdir(ann_path)
        projects = [p for p in projects if os.path.isdir(os.path.join(ann_path, p))]
        projects.sort()
        return projects

annData = AnnotationData()