from controller.expansion_ctl import *
from dash import callback, Input, Output, State, no_update, Patch
from dash.exceptions import PreventUpdate
from pages.components.fileSelecter import fileSelecter
from controller.auth import *


@callback(
    Input('expand-select-taskname', 'value'),
)
def update_table_metadata(taskName):
    """
    更新任务列表元数据
    """
    if taskName:
        update_table_metadata_with_taskname(taskName)

@callback(
    Input('exp-delete-task-confirm', 'confirmCounts'),
    State('expand-select-taskname', 'value'),
    running=[
        (Output('exp-delete-task', 'loading'), True, False),
        (Output('exp-start-task', 'disabled'), True, False),
        (Output('expand-select-taskname', 'disabled'), True, False)
    ],
    prevent_initial_call=True
)
def delete_task(nc, taskName):
    """
    删除任务
    """
    if nc and taskName and verify_modify_permission():
        # 如果任务正在运行禁止删除
        delete_exptask_from_disk(taskName)
@callback(
    Output('expand-select-taskname', 'options', allow_duplicate=True),
    Input('expand-select-taskname-tooltip', 'open'),
    prevent_initial_call=True
)
def update_expTask_options(open):
    """
    更新胞域扩增任务列表
    """
    if open:
        expTasks = expData.get_exist_tasks()
        return expTasks
    return no_update

@callback(
    Output('expandtask-select-taskname', 'options'),
    Input('expandtask-select-taskname-tooltip', 'open'),
    prevent_initial_call=True
)
def update_tasklist_options(open):
    """
    获取可扩增分割任务列表
    """
    if open:
        return get_expandable_task_options()
    return no_update

@callback(
    Output('expand-modal-newtask', 'visible', allow_duplicate=True),
    Output('expand-select-taskname', 'options'),
    Output('expand-select-taskname', 'value'),
    Input('expand-button-submitTask', 'nClicks'),
    State('main-title-username', 'children'),
    State('expandtask-select-taskname', 'value'),
    State('expand-select-mode', 'value'),
    State('expand-input-patchsize', 'value'),
    State('expand-input-binsize', 'value'),
    State('expand-input-epochs', 'value'),
    State('expand-input-diameter', 'value'),
    State('expand-input-neighbor', 'value'),
    prevent_initial_call=True
)
def submit_task(nClicks, userName, taskName, mode, patchsize, binsize, epochs, diameter, neighbor):
    """
    提交胞域扩增任务
    """
    if nClicks:
        if not taskName:
            set_head_notice('Please select a task', type='warning')
        else:
            submit_expansion_task(userName, taskName, mode, patchsize, binsize, epochs, diameter, neighbor)
            expTasks = expData.get_exist_tasks()
            return False, expTasks, taskName
    return no_update, no_update, no_update

@callback(
    Output('expand-input-patchsize', 'disabled'),
    Output('expand-input-patchsize', 'value'),
    Input('expand-select-mode', 'value')
)
def update_patchsize_disabled(mode):
    """
    根据分割模式选择是否禁用patchSize输入框
    """
    if mode.startswith('p'):
        return False, 300
    else:
        return True, 0
    
@callback(
    Output('expand-modal-newtask', 'visible'),
    Output('expandtask-select-taskname', 'value', allow_duplicate=True),
    Input('expansion-button-newtask', 'nClicks'),
    prevent_initial_call=True
)
def show_newtask_modal(nClicks):
    """
    显示新建任务对话框
    """
    if nClicks and verify_modify_permission():
        return True, None
    return no_update, no_update