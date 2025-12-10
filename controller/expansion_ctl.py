from controller.auth import *
from dataManager.expansion_d import expData
from dataManager.segmentation_d import segData

def update_table_metadata_with_taskname(taskName):
    """
    更新任务列表元数据
    """
    metadata = expData.read_taskinfo(taskName)
    slices = metadata.pop('slices', [])
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
    segTasks = segData.get_exist_tasks()
    expTasks = set(expData.get_exist_tasks())
    return [task for task in segTasks if task not in expTasks]

def submit_expansion_task(userName, taskName, mode, patchsize, binsize, epochs, diameter, neighbor):
    """
    提交胞域扩增任务
    """
    expData.save_expansion_task(userName, taskName, mode, patchsize, binsize, epochs, diameter, neighbor)
    set_head_notice(f'Expansion task {taskName} submitted successfully!', type='success')