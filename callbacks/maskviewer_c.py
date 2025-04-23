from controller.maskviewer_ctl import *
from dash import callback, Input, Output, State, no_update, Patch
from dash.exceptions import PreventUpdate

@callback(
    Input('maskviewer-select-taskname', 'value')
)
def change_task(taskname):
    """
    切换项目更新列表
    """
    if taskname:
        update_taskname(taskname)

@callback(
    Input('maskviewer-select-graph', 'value'),
)
def graph_change(graph):
    """
    切换不同图
    """
    if graph:
        update_graph(graph)
@callback(
    Input('init-restore-maskviewer', 'n_intervals'),
    State('maskviewer-store-taskname', 'data')
)
def restore_segmentation(n_intervals, lastTaskName):
    """
    初始回调函数，恢复网页状态
    """
    if n_intervals:
        restore_initial_data(lastTaskName)

@callback(
    Output('maskviewer-select-slice', 'options'),
    Input('maskviewer-select-slice-tooltip', 'open'),
    State('maskviewer-select-taskname', 'value'),
    prevent_initial_call=True
)
def update_tasklist_options(open, taskname):
    """
    每次点击选择切片下拉列表时更新切片状态
    """
    if open and taskname:
        tasks = maskData.get_task_slices(taskname)
        return tasks
    return no_update
@callback(
    Output('maskviewer-select-taskname', 'options'),
    Input('maskviewer-select-taskname-tooltip', 'open'),
    prevent_initial_call=True
)
def update_tasklist_options(open):
    """
    每次点击下拉列表时更新任务状态
    """
    if open:
        tasks = maskData.get_exist_tasks()
        return tasks
    return no_update
