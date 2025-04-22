from controller.segmentation_ctl import *
from dash import callback, Input, Output, State, no_update, Patch
from pages.components.fileSelecter import fileSelecter
from dash.exceptions import PreventUpdate


@callback(
    Input('seg-button-showBug', 'nClicks'),
    State('seg-select-taskname', 'value'),
    prevent_initial_call=True
)
def show_bug(nc, taskName):
    if nc:
        raise_runtime_bug(taskName)
@callback(
    Output('seg-table-metadata', 'data', allow_duplicate='True'),
    Output('seg-table-tasklist', 'data', allow_duplicate='True'),
    Output('seg-bug-panel', 'style'),
    Input('segmentation-interval', 'n_intervals'),
    State('seg-select-taskname', 'value'),
    prevent_initial_call=True
)
def update_segment_status(n_intervals, taskname):
    """
    同步后台进程至前端
    """
    if n_intervals:
        try:
            metadata = no_update
            data = no_update
            style = no_update
            update = False
            pmetadata, pdata, exception = segData.get_diff(taskname)
            if pmetadata is not None:
                metadata = pmetadata
                update = True
            if pdata is not None:
                data = pdata
                update = True
            if exception is not None:
                update = True
                style = {'display':'flex'}
            if update:
                return metadata, data, style
        except Exception as e:
            pass
    raise PreventUpdate

@callback(
    Output('segmentation-interval', 'disabled'),
    Input('segmentation-event-loop', 'n_intervals'),
    State('seg-select-taskname', 'value'),
    prevent_initial_call=True
)
def event_loop(n_intervals, taskname):
    """
    监听分割任务状态
    """
    if n_intervals:
        task_running = segData.has_running_task(taskname)
        task_virtual = segData.has_userVirtualTask(taskname)
        if task_running or task_virtual:
            return False
        return True
    raise PreventUpdate

@callback(
    Input('seg-start-task', 'nClicks'),
    State('seg-select-taskname', 'value'),
    prevent_initial_call=True
)
def start_segment_task(nc, taskname):
    """
    启动分割任务
    """
    if nc and verify_modify_permission():
        if taskname:
            task_running = segData.has_running_task(taskname)
            task_virtual = segData.has_userVirtualTask(taskname)
            if task_running or task_virtual:
                set_head_notice(taskname+' is running, please wait for it to complete !', type='warning')
                return
            start_segmentation_task(taskname)
        else:
            set_head_notice('Please select a task to start !', type='warning')

@callback(
    Input('seg-delete-task-confirm', 'confirmCounts'),
    State('seg-select-taskname', 'value'),
    running=[
        (Output('seg-delete-task', 'loading'), True, False),
        (Output('seg-start-task', 'disabled'), True, False),
        (Output('seg-select-taskname', 'disabled'), True, False)
    ],
    prevent_initial_call=True
)
def delete_task(nc, taskName):
    """
    删除任务
    """
    if nc and taskName and verify_modify_permission():
        task_running = segData.has_running_task(taskName)
        task_virtual = segData.has_userVirtualTask(taskName)
        if task_running or task_virtual:
            set_head_notice(taskName+' is running, please wait for it to complete !', type='warning')
            return
        delete_task_from_disk(taskName)
@callback(
    Input('init-restore-segmentation', 'n_intervals'),
    State('seg-store-taskname', 'data')
)
def restore_segmentation(n_intervals, lastTaskName):
    """
    初始回调函数，恢复网页状态
    """
    if n_intervals:
        restore_initial_data(lastTaskName)
@callback(
    Input('seg-select-taskname', 'value'),
)
def update_table_tasklist(taskName):
    """
    根据下拉列表选项更新任务表格
    """
    if taskName:
        update_table_with_tasklist(taskName)

@callback(
    Output('seg-select-taskname', 'options'),
    Input('seg-select-taskname-tooltip', 'open'),
    prevent_initial_call=True
)
def update_tasklist_options(open):
    if open:
        tasks = list(segData.get_exist_tasks())
        tasks.sort()
        return tasks
    return no_update

@callback(
    Input('seg-button-submitTaskList', 'nClicks'),
    State('seg-input-taskname', 'value'),
    State('seg-tasklist-filename', 'type'),
    State('seg-select-modelType', 'value'),
    State('seg-input-diameter', 'value'),
    State('seg-input-batchsize', 'value'),
    State('seg-checkbox-useGPU', 'checked'),
    State('main-title-username', 'children'),
    prevent_initial_call=True
)
def submit_tasklist(nc, taskName, fileStatus, modelType, diameter, batchsize, useGPU, username):
    """
    检查提交的任务列表，并将任务持久化到本地
    """
    if nc:
        process_submited_tasklist(taskName, fileStatus, modelType, diameter, batchsize, useGPU, username)

@callback(
    Output('seg-tasklist-filename', 'children', allow_duplicate='True'),
    Output('seg-tasklist-filename', 'type', allow_duplicate='True'),
    Input('seg-dragger-upload', 'lastUploadTaskRecord'),
    prevent_initial_call=True
)
def upload_status(lastRecord):
    """
    监听文件上传状态
    """
    if lastRecord['taskStatus']=='success':
        return lastRecord['fileName'], 'success'
    else:
        set_head_notice(lastRecord['fileName']+' upload failed, please check file format!', type='error')
        return 'No file', 'secondary'
@callback(
    Input('seg-button-importTaskList', 'nClicks'),
    prevent_initial_call=True  
)
def open_import_box(nc):
    """
    打开导入任务列表文件窗口
    """
    if nc:
        fileSelecter.open_import_box()

@callback(
    Output('seg-modal-newtask', 'visible'),
    Output('seg-tasklist-filename', 'children'),
    Output('seg-tasklist-filename', 'type'),
    Input('segmentation-button-newtask', 'nClicks'),
    prevent_initial_call=True
)
def show_newtask_modal(nClicks):
    """
    显示创建任务对话框
    """
    if nClicks and verify_modify_permission():
        return True, 'No file', 'secondary'
    return no_update, no_update, no_update
