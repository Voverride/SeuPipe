from controller.maskviewer_ctl import *
from dash import callback, Input, Output, State, no_update, Patch
from dash.exceptions import PreventUpdate
import feffery_antd_components as fac

@callback(
    Output('mv-graph-left', 'figure', allow_duplicate=True),
    Output('mv-graph-right', 'figure', allow_duplicate=True),
    Input('mv-checkbox-mask', 'checked'),
    Input('mv-checkbox-contour', 'checked'),
    prevent_initial_call=True
)
def switch_mask_contour(checked_mask, checked_contour):
    """
    开关mask与contour
    """
    patch = Patch()
    patch['data'][1]['visible'] = checked_mask
    patch['data'][2]['visible'] = checked_contour
    return patch, patch
@callback(
    Output("mv-graph-left", "figure", allow_duplicate=True),
    Input("mv-graph-right", "relayoutData"),
    State("mv-graph-right", "figure"),
    State('maskviewer-select-graph', 'value'),
    prevent_initial_call=True,
)
def update_relayout(rightLayout, rightFig, graph):
    """
    同步左图状态
    """
    if graph=='registration' or rightFig is None or rightLayout is None:
        raise PreventUpdate
    patch = Patch()
    patch['layout']['xaxis'] = rightFig['layout']['xaxis']
    patch['layout']['yaxis'] = rightFig['layout']['yaxis']
    return patch

@callback(
    Output("mv-graph-right", "figure", allow_duplicate=True),
    Input("mv-graph-left", "relayoutData"),
    State("mv-graph-left", "figure"),
    State('maskviewer-select-graph', 'value'),
    prevent_initial_call=True,
)
def update_relayout(leftLayout, leftFig, graph):
    """
    同步右图状态
    """
    if graph=='registration' or leftFig is None or leftLayout is None:
        raise PreventUpdate
    patch = Patch()
    patch['layout']['xaxis'] = leftFig['layout']['xaxis']
    patch['layout']['yaxis'] = leftFig['layout']['yaxis']
    return patch

@callback(
    Input('maskviewer-select-graph', 'value'),
    State('maskviewer-select-slice', 'value'),
    State('mv-checkbox-mask', 'checked'),
    State('mv-checkbox-contour', 'checked'),
    running=[
        (Output('maskviewer-select-graph', 'disabled'), True, False),
        (Output('mv-select-graph-label', 'children'), 'Loading...', 'Select Graph'),
    ]
)
def graph_change(graph, slice, showMask, showContour):
    """
    切换不同图
    """
    if graph:
        update_graph_with_type(slice, graph, showMask=showMask, showContour=showContour)
@callback(
    Input('maskviewer-select-slice', 'value'),
    State('maskviewer-select-graph', 'value'),
    State('maskviewer-select-taskname', 'value'),
    State('mv-checkbox-mask', 'checked'),
    State('mv-checkbox-contour', 'checked'),
    running=[
        (Output('maskviewer-select-slice', 'disabled'), True, None),
        (Output('mv-select-slice-label', 'children'), 'Loading...', 'Select Slice'),
    ]
)
def slice_change(slice, graph, taskname, showMask, showContour):
    """
    切换图类别和切片时，更新图像
    """
    if taskname:
        update_graph_with_slice(taskname, slice, graph, showMask, showContour)

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
