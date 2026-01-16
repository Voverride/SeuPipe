import re
import os
import threading
from controller.notice import set_head_notice, set_aside_notice
from dash import set_props
from api.segmentation import *
from controller.auth import get_request_usrname
from dataManager.segmentation_d import segData
from utils.commonfuc import get_current_date

def raise_runtime_bug(project):
    """
    抛出运行时的异常
    """
    project_info = segData.get_project_info(project)
    exception = project_info['exception']
    set_aside_notice('Runtime Error', exception, 'error')

def start_segmentation_project(projectname):
    """
    启动分割项目
    """
    if projectname in segData._running_projects:
        set_head_notice(f'{projectname} is running, please wait for complete !', type='warning')
        return
    project_info = segData.get_project_info(projectname)
    if project_info is None or len(project_info['slices'])==0:
        set_head_notice(f'There was no slice to execute in {projectname} !', type='warning')
    else:
        observed_info = segData.add_running_project(projectname, project_info)
        thread = threading.Thread(
            target=run_segmentation_project,
            args=(observed_info,)
        )
        thread.start()
        set_head_notice(f'{projectname} has been started !', type='success')

def restore_init_data(projectName):
    """
    恢复初始数据
    """
    set_props('seg-select-project', dict(value=projectName))    

def update_table_with_project(projectName):
    """
    更新项目数据表格
    """
    project_info = segData.get_project_info(projectName)
    if project_info is None:
        set_props('seg-select-project', dict(value=None))
        set_head_notice(f'{projectName} related data may been removed from your disk !', type='warning')
        set_props('seg-table-metadata', dict(data=[]))
        set_props('seg-table-slices', dict(data=[]))
        set_props('seg-bug-panel', dict(style={'display':'none'}))
    else:
        slices = project_info['slices']
        project_info.pop('slices')
        
        project_info['gpu'] = str(project_info['gpu'])
        if project_info['diameter']==0:
            project_info['diameter'] = 'auto'
        exception = project_info['exception']
        if exception is None:
            set_props('seg-bug-panel', dict(style={'display':'none'}))
        else:
            set_props('seg-bug-panel', dict(style={'display':'flex'}))
        set_props('seg-table-metadata', dict(data=[project_info]))
        set_props('seg-table-slices', dict(data=slices))
        is_running = segData.has_running_project(projectName)
        if is_running:
            set_props('seg-start-project', dict(disabled=True))
            set_props('seg-delete-project', dict(disabled=True))
        else:
            set_props('seg-start-project', dict(disabled=False))
            set_props('seg-delete-project', dict(disabled=False))

def delete_project_from_disk(projectName):
    """
    删除项目
    """
    segData.delete_project(projectName)
    set_head_notice(f'{projectName} has been removed from your disk !', type='success')
    set_props('seg-table-metadata', dict(data=[]))
    set_props('seg-table-slices', dict(data=[]))
    set_props('seg-select-project', dict(value=None))

def process_submited_project(projectName, fileStatus, modelType, diameter, batchsize, useGPU):
    """
    处理提交的项目，并将项目持久化到本地
    """
    if projectName is None or projectName.strip()=='':
        set_head_notice('Project Name cannot be empty', type='error')
        return False
    projectList = segData.get_exist_projects()
    if projectName in projectList:
        set_head_notice(f'Project Name {projectName} already exists', type='warning')
        return False
    if fileStatus!='success':
        set_head_notice(f'Please upload your file first for project {projectName}', type='warning')
        return False
    username = get_request_usrname()
    current_date = get_current_date()

    segData.save_project_metadata({
        'name': projectName,
        'creator': username,
        'date': current_date,
        'model': modelType,
        'diameter': diameter,
        'batchsize': batchsize,
        'gpu': useGPU,
        'progress': 0,
        'exception': None
    })

    set_head_notice(f'Project {projectName} created successfully!', type='success')
    set_props('seg-modal-newproject', dict(visible=False))
    set_props('seg-select-project', dict(value=projectName))

def read_project_metadata_file(fpath:str):
    """
    读取项目元数据文件
    """
    success = True
    try:
        with open(fpath, 'r') as f:
            parse_project_metadata(f.readlines())
    except Exception as e:
        success = False

    filename = os.path.basename(fpath)
    if success:
        set_props('seg-project-filename', dict(type='success', children=filename))
        set_head_notice(filename+' import successfully!', type='success')
    else:
        set_props('seg-project-filename', dict(type='secondary', children='No file'))
        set_head_notice(filename+' import failed, please check your file format', type='error')
    return success
def parse_project_metadata(lines:list):
    """
    解析项目元数据列表
    """
    try:
        slices = []
        for line in lines:
            z, img, gem, *extra = re.split(r'[,]+', line.strip())
            slices.append({
                'z': int(z),
                'img': img.strip(),
                'gem': gem.strip(),
                'segmentation': {'status': 'warning', 'text': 'waiting'},
                'postprocess': {'status': 'warning', 'text': 'waiting'},
            })
        segData.set_project_metadata(slices)
    except Exception as e:
        segData.reset_project_metadata()
        raise e