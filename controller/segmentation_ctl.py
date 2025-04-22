from dataManager.segmentation_d import segData
from dash import set_props
from controller.auth import *
from api.segmentation import *
from controller.notice import *
import threading
import re


def raise_runtime_bug(taskName):
    """
    抛出运行时的异常
    """
    taskinfo = segData.read_taskinfo(taskName)
    exception = taskinfo['exception']
    set_aside_notice('Task Error', exception, 'error')

def start_segmentation_task(taskName):
    """
    启动分割任务
    """
    task_info = segData.read_taskinfo(taskName)
    if task_info is None or len(task_info['data'])==0:
        set_head_notice('There was no task to execute !', type='warning')
    else:
        segData.reset_taskinfo(task_info)
        segData.save_taskinfo(task_info)
        segData.create_running_task(task_info)
        thread = threading.Thread(
            target=run_segmentation_task,
            args=(task_info,)
        )
        task_info['thread'] = thread
        thread.start()
        set_head_notice(taskName+' has been started !', type='success')
        set_props('seg-table-metadata', dict(data=[task_info['metadata']]))
        set_props('seg-table-tasklist', dict(data=task_info['data']))


def delete_task_from_disk(taskName):
    """
    从磁盘删除任务
    """
    segData.delete_task(taskName)
    set_head_notice(taskName+' has been removed from your disk !', type='success')
    set_props('seg-table-metadata', dict(data=[]))
    set_props('seg-table-tasklist', dict(data=[]))
    set_props('seg-select-taskname', dict(value=None))
def restore_initial_data(lastTaskName):
    """
    恢复网页初始数据
    """
    if lastTaskName:
        set_props('seg-select-taskname', dict(value=lastTaskName))

def update_table_with_tasklist(taskName):
    """
    根据下拉列表选项更新任务列表
    """
    set_props('seg-store-taskname', dict(data=taskName))
    if segData.has_userVirtualTask(taskName):
        task_info = segData.get_userVirtualTask(taskName)
    else:
        task_info = segData.read_taskinfo(taskName)
    if task_info is None:
        # set_props('seg-select-taskname', dict(value=None))
        set_head_notice(taskName+' related data may have been removed from your disk !', type='warning')
        set_props('seg-table-metadata', dict(data=[]))
        set_props('seg-table-tasklist', dict(data=[]))
        set_props('seg-bug-panel', dict(style={'display':'none'}))
    else:
        metadata = task_info['metadata']
        metadata['GPU'] = str(metadata['GPU'])
        if metadata['diameter']==0:
            metadata['diameter'] = 'auto'
        exception = task_info['exception']
        if exception is None:
            set_props('seg-bug-panel', dict(style={'display':'none'}))
        else:
            set_props('seg-bug-panel', dict(style={'display':'flex'}))
            
        set_props('seg-table-metadata', dict(data=[metadata]))
        set_props('seg-table-tasklist', dict(data=task_info['data']))

def process_submited_tasklist(taskName, fileStatus, modelType, diameter, batchsize, useGPU, username):
    """
    处理提交的任务列表，并将任务持久化到本地
    """
    if taskName.strip()=='' or taskName is None:
        set_head_notice('Task Name cannot be empty', type='error')
        return False
    taskList = segData.get_exist_tasks()
    if taskName in taskList:
        set_head_notice('Task Name already exists', type='warning')
        return False
    if fileStatus!='success':
        set_head_notice('Please upload your file first', type='warning')
        return False
    
    segData.set_temptask_metadata(taskName,{
        'creator': username,
        'model': modelType,
        'diameter': diameter,
        'batchsize': batchsize,
        'GPU': useGPU,
        'progress': 0
    })

    segData.save_temptask()
    set_head_notice('Task '+taskName+' created successfully!', type='success')
    set_props('seg-modal-newtask', dict(visible=False))
    set_props('seg-select-taskname', dict(options=list(taskList), value=taskName))

def parse_tasklist(lines:list):
    """
    解析任务列表
    """
    # data=[
    #     {
    #         'z': 1,
    #         'image': 'image2',
    #         'gem': 'gem1',
    #         'registration': {'status': 'success', 'text': 'success'},
    #         'segmentation': {'status': 'success', 'text': 'success'},
    #     }
    # ],
    try:
        data = []
        for line in lines:
            z, image, gem, *extra = re.split(r'[,\s]+', line.strip())
            data.append({
                'z': int(z),
                'image': image.strip(),
                'gem': gem.strip(),
                'registration': {'status': 'warning', 'text': 'waiting'},
                'segmentation': {'status': 'warning', 'text': 'waiting'},
            })
        segData.set_temptask_data(data)
    except Exception as e:
        segData.reset_temptask_data()
        raise e
def read_tasklist_file(fpath:str):
    """
    读取任务列表文件
    """
    success = True
    try:
        with open(fpath, 'r') as f:
            parse_tasklist(f.readlines())
    except Exception as e:
        success = False

    filename = os.path.basename(fpath)
    if success:
        set_props('seg-tasklist-filename', dict(type='success', children=filename))
        set_head_notice(filename+' import successfully!', type='success')
    else:
        set_props('seg-tasklist-filename', dict(type='secondary', children='No file'))
        set_head_notice(filename+' import failed, please check your file format', type='error')
    return success