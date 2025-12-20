from controller.maskviewer_ctl import *
from dash import callback, Input, Output, State, no_update, Patch
from dash.exceptions import PreventUpdate
import feffery_antd_components as fac
from pages.components.fileSelecter import fileSelecter


@callback(
    Input('maskviewer-button-exportData', 'nClicks'),
    State('maskviewer-select-taskname', 'value'),
    State('maskviewer-select-export-type', 'value')
)
def export_data(nc, taskname, type):
    """
    显示导出数据窗口
    """
    if nc:
        if not taskname:
            set_head_notice('Please select task!', 'warning')
        elif not type:
            set_head_notice('Please select export type!', 'warning')
        else:
            maskData.set_export_task(taskname, type)
            fileSelecter.open_export_box()

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
    Output('mv-graph-right', 'figure', allow_duplicate=True),
    Input('maskviewer-right-graph', 'value'),
    State('maskviewer-select-graph', 'value'),
    State('maskviewer-select-taskname', 'value'),
    State('maskviewer-select-slice', 'value'),
    State('mv-checkbox-mask', 'checked'),
    State('mv-checkbox-contour', 'checked'),
    running=[
        (Output('maskviewer-right-graph', 'disabled'), True, False),
        (Output('mv-right-graph-label', 'children'), 'Loading...', 'Right Graph'),
    ],
    prevent_initial_call=True,
)
def right_graph_change(graphName, graphType, taskname, slice, showMask, showContour):
    """
    切换右图
    """
    if graphType != 'registration':
        if taskname and slice:
            graph = maskData.get_graph(taskname, slice, graphName, showMask=showMask, showContour=showContour)
            return graph
    return no_update

@callback(
    Output('mv-graph-left', 'figure', allow_duplicate=True),
    Input('maskviewer-left-graph', 'value'),
    State('maskviewer-select-graph', 'value'),
    State('maskviewer-select-taskname', 'value'),
    State('maskviewer-select-slice', 'value'),
    State('mv-checkbox-mask', 'checked'),
    State('mv-checkbox-contour', 'checked'),
    running=[
        (Output('maskviewer-left-graph', 'disabled'), True, False),
        (Output('mv-left-graph-label', 'children'), 'Loading...', 'Left Graph'),
    ],
    prevent_initial_call=True,
)
def left_graph_change(graphName, graphType, taskname, slice, showMask, showContour):
    """
    切换左图
    """
    if graphType != 'registration':
        if taskname and slice:
            graph = maskData.get_graph(taskname, slice, graphName, showMask=showMask, showContour=showContour)
            return graph
    return no_update

@callback(
    Input('maskviewer-select-graph', 'value'),
    State('maskviewer-select-taskname', 'value'),
    State('maskviewer-select-slice', 'value'),
    State('mv-checkbox-mask', 'checked'),
    State('mv-checkbox-contour', 'checked'),
    State('maskviewer-left-graph', 'value'),
    State('maskviewer-right-graph', 'value'),
    running=[
        (Output('maskviewer-select-graph', 'disabled'), True, False),
        (Output('mv-select-graph-label', 'children'), 'Loading...', 'Select Graph'),
    ]
)
def graph_change(graph, taskname, slice, showMask, showContour, leftGraph, rightGraph):
    """
    切换不同图
    """
    if graph and slice and taskname:
        update_graph_with_type(taskname, slice, graph, showMask, showContour, leftGraph, rightGraph)
@callback(
    Input('maskviewer-select-slice', 'value'),
    State('maskviewer-select-graph', 'value'),
    State('maskviewer-select-taskname', 'value'),
    State('mv-checkbox-mask', 'checked'),
    State('mv-checkbox-contour', 'checked'),
    State('maskviewer-left-graph', 'value'),
    State('maskviewer-right-graph', 'value'),
    running=[
        (Output('maskviewer-select-slice', 'disabled'), True, None),
        (Output('mv-select-slice-label', 'children'), 'Loading...', 'Select Slice'),
    ]
)
def slice_change(slice, graph, taskname, showMask, showContour, leftGraph, rightGraph):
    """
    切换图类别和切片时，更新图像
    """
    if taskname and slice:
        update_graph_with_slice(taskname, slice, graph, showMask, showContour, leftGraph, rightGraph)

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

@callback(
    Output('maskviewer-select-export-type', 'options'),
    Input('maskviewer-select-export-type-tooltip', 'open'),
    State('maskviewer-select-taskname', 'value'),
    prevent_initial_call=True
)
def update_export_type_options(open, taskname):
    """
    更新导出结果选项
    """
    if open and taskname:
        types = set()
        for slice in maskData.get_task_slices(taskname):
            curTypes = maskData.get_figure_types(taskname, slice)
            curTypes = set(curTypes)
            if len(types) > 0:
                types = types.intersection(curTypes)
            else:
                types = types.union(curTypes)
        types = list(types)
        types.sort()
        return types
    return no_update

@callback(
    Output('maskviewer-left-graph', 'options'),
    Input('maskviewer-left-graph-tooltip', 'open'),
    State('maskviewer-select-taskname', 'value'),
    State('maskviewer-select-slice', 'value'),
    prevent_initial_call=True
)
def update_leftGraph_options(open, taskname, slice):
    """
    更新左图选项
    """
    if open and taskname and slice:
        types = maskData.get_figure_types(taskname, slice)
        return types
    return no_update

@callback(
    Output('maskviewer-right-graph', 'options'),
    Input('maskviewer-right-graph-tooltip', 'open'),
    State('maskviewer-select-taskname', 'value'),
    State('maskviewer-select-slice', 'value'),
    prevent_initial_call=True
)
def update_rightGraph_options(open, taskname, slice):
    """
    更新右图选项
    """
    if open and taskname and slice:
        types = maskData.get_figure_types(taskname, slice)
        return types
    return no_update
