from controller.auth import *
from dataManager.expansion_d import expData
from dataManager.segmentation_d import segData
from websocket.websocket import ws
from controller.notice import *
import threading

from api.expansion import *

def start_expansion_with_taskname(taskName):
    """
    启动胞域扩增任务
    """
    segInfo = segData.get_project_info(taskName)
    progress = segInfo.get('progress', 0)
    if progress < 1:
        set_head_notice('Please wait for the segmentation task to complete !', type='warning')
    else:
        expData.create_running_task(taskName, ws.notifyUpdateExpansionStatus)
        thread = threading.Thread(
            target=run_scs, 
            args=(taskName,)
        )
        thread.start()

def raise_runtime_bug(taskName):
    """
    抛出运行时的异常
    """
    taskinfo = expData.read_taskinfo(taskName)
    exception = taskinfo['exception']
    set_aside_notice('Task Error', exception, 'error')

def update_table_metadata_with_taskname(taskName):
    """
    更新任务列表元数据
    """
    metadata = expData.read_taskinfo(taskName)
    slices = metadata.pop('slices', [])
    if metadata['exception'] is not None:
        set_props('exp-bug-panel', dict(style={'display':'flex'}))
    else:
        set_props('exp-bug-panel', dict(style={'display':'none'}))
    set_props('exp-table-metadata', dict(data=[metadata]))
    set_props('exp-table-tasklist', dict(data=slices))

def delete_exptask_from_disk(taskName):
    """
    从磁盘删除胞域扩增任务
    """
    expData.delete_task(taskName)
    set_head_notice(taskName+' has been removed from your disk !', type='success')
    set_props('exp-table-metadata', dict(data=[]))
    set_props('exp-table-tasklist', dict(data=[]))
    set_props('expand-select-taskname', dict(value=None))
def get_expandable_task_options()->list:
    """
    获取可胞域扩增任务列表
    """
    segTasks = segData.get_exist_projects()
    expTasks = set(expData.get_exist_tasks())
    return [task for task in segTasks if task not in expTasks]

def submit_expansion_task(userName, taskName, mode, patchsize, binsize, epochs, diameter, neighbor):
    """
    提交胞域扩增任务
    """
    expData.save_expansion_task(userName, taskName, mode, patchsize, binsize, epochs, diameter, neighbor)
    set_head_notice(f'Expansion task {taskName} submitted successfully!', type='success')