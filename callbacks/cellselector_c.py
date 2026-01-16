from controller.cellselector_ctl import *
from controller.auth import *
from websocket.websocket import ms
from pages.components.fileSelecter import fileSelecter
from dash.exceptions import PreventUpdate
from dash import callback, Input, Output, State, no_update, callback_context, Patch

@callback(
    Input('sel-cluster-showBug', 'nClicks'),
    State('sel-select-project', 'value'),
    prevent_initial_call=True
)
def show_bug(nc, project):
    """
    显示聚类异常信息
    """
    if nc and project:
        show_bug_info(project)


@callback(
    Input('sel-button-exportData', 'nClicks'),
    State('sel-select-project', 'value'),
)
def export_data(nc, project):
    """
    显示导出数据窗口
    """
    if nc:
        if not project:
            set_head_notice('Please select a project first!', 'warning')
        else:
            selData.set_export_project(project)
            fileSelecter.open_export_box()

@callback(
    Input('refresh-sel-status', 'children'),
    State('sel-select-project', 'value'),
    State('sel-select-slice', 'value'),
    prevent_initial_call=True
)
def refresh_cell_selector_status(_, project, slice):
    """
    刷新细筛选状态
    """
    if project and not slice:
        set_props('sel-select-project', dict(value=project))
    if project and slice:
        set_props('sel-select-slice', dict(value=slice))

@callback(
     Input('sel-button-retain', 'nClicks'),
     Input('sel-button-remove', 'nClicks'),
     Input('sel-button-reset', 'nClicks'),
     State('sel-select-project', 'value'),
     State('sel-select-slice', 'value'),
     State('sel-cell-graph', 'selectedData')
)
def retain_cells(retain, remove, reset, project, slice, selectedData):
    """
    筛选细胞
    """
    triggered = callback_context.triggered
    if not triggered:
        return
    prop_id = triggered[0]['prop_id']
    if 'sel-button-retain' in prop_id:
        action = 'retain'
    elif 'sel-button-remove' in prop_id:
        action = 'remove'
    elif 'sel-button-reset' in prop_id:
        action = 'reset'
    else:
        action = None
    if action is None or not verify_modify_permission():
        return
    if action=='reset' and project and slice:
        set_props('sel-button-retain', dict(disabled=True))
        set_props('sel-button-remove', dict(disabled=True))
        select_cell(project, slice, None, operation=action)
        return
    if selectedData is None:
        return
    points = selectedData.get('points', None)
    if points is None or len(points) == 0:
        return
    pointIndex = [point['customdata'] for point in points]
    set_props('sel-button-retain', dict(disabled=True))
    set_props('sel-button-remove', dict(disabled=True))
    select_cell(project, slice, pointIndex, operation=action)

@callback(
    Output('sel-button-retain', 'disabled'),
    Output('sel-button-remove', 'disabled'),
    Input('sel-cell-graph', 'selectedData')
)
def update_operation_button(selectedData):
    """
    当选中细胞时，更新操作按钮禁用状态
    """
    if selectedData is None:
        return True, True
    points = selectedData.get('points', None)
    if points is None or len(points) == 0:
        return True, True
    return False, False

@callback(
    Output('sel-cell-graph', 'figure', allow_duplicate=True),
    Input('sel-checkbox-image', 'checked'),
    Input('sel-checkbox-spot', 'checked'),
    Input('sel-select-spotSize', 'value'),
    Input('sel-colorPicker-spotColor', 'value'),
    State('sel-select-slice', 'value'),
    State('sel-store-clusters', 'data'),
    prevent_initial_call=True
)
def update_show_image_scatter(show_image, show_scatter, spotSize, spotColor, slice, clusters):
    """
    更新是否显示图像和散点图
    """
    if not slice:
        raise PreventUpdate
    fig = Patch()
    fig = update_image_style(fig, show_image, show_scatter, clusters, scatter_color=spotColor, scatter_size=spotSize)
    return fig

@callback(
    Output('sel-cell-graph', 'figure'),
    Input('sel-select-slice', 'value'),
    State('sel-select-project', 'value'),
    State('sel-checkbox-image', 'checked'),
    State('sel-checkbox-spot', 'checked'),
    State('sel-select-spotSize', 'value'),
    State('sel-colorPicker-spotColor', 'value'),
    running=[
        (Output('sel-select-slice', 'disabled'), True, False)
    ]
)
def update_slice_graph(slice, project, show_image, show_scatter, spotSize, spotColor):
    """
    更新切片图像
    """
    if slice is None:
        return
    ms.clinetSend(ms._updateParams_c2s, params=dict(
        project=project,
        slice=slice
    ))
    fig, clusters = get_slice_graph(project, slice)
    set_props('sel-store-clusters', dict(data=clusters))
    if clusters == 1:
        set_props('sel-colorPicker-spotColor', dict(disabled=False))
    else:
        set_props('sel-colorPicker-spotColor', dict(disabled=True))
    fig = update_image_style(fig, show_image, show_scatter, clusters, scatter_color=spotColor, scatter_size=spotSize)
    if fig is None:
        return no_update
    return fig

@callback(
    Output('sel-select-slice', 'value'),
    Output('sel-select-slice', 'options'),
    Input('sel-select-project', 'value')
)
def update_slice_options(project):
    """
    更新切片选项
    """
    ms.clinetSend(ms._updateParams_c2s, params=dict(
        project=project,
        slice=None
    ))
    if project is None:
        return no_update, no_update
    show_clustering = update_cluster_status(project)
    if show_clustering:
        return None, []
    slices = selData.get_project_slices(project)
    return None, slices

@callback(
    Output('sel-select-project', 'options'),
    Input('sel-select-project-tooltip', 'open')
)
def update_project_options(open):
    """
    更新项目选项
    """
    if open:
        projects = selData.get_exist_projects()
        return projects
    else:
        raise PreventUpdate

@callback(
    Output('sel-cell-graph', 'figure', allow_duplicate=True),
    Output('sel-select-project', 'value', allow_duplicate=True),
    Output('sel-select-slice', 'value', allow_duplicate=True),
    Output('sel-select-slice', 'options', allow_duplicate=True),
    Input('sel-delete-project-confirm', 'confirmCounts'),
    State('sel-select-project', 'value'),
    running=[
        (Output('sel-delete-project', 'loading'), True, False),
    ],
    prevent_initial_call=True
)
def delete_project(nc, project_name):
    """
    删除项目
    """
    if nc and project_name and verify_modify_permission():
        if selData.is_clustering(project_name):
            set_head_notice(f'{project_name} is clustering, please wait...!', type='warning')
            return no_update, no_update, no_update, no_update
        selData.delete_project(project_name)
        set_head_notice('Delete successfully!', type='success')
        return None, None, None, []
    return no_update, no_update, no_update, no_update

@callback(
    Output('sel-select-project', 'value'),
    Output('sel-modal-newproject', 'visible', allow_duplicate=True),
    Input('sel-button-submitProject-modal', 'nClicks'),
    State('sel-select-project-modal', 'value'),
    State('sel-select-result-modal', 'value'),
    State('sel-input-resolution-modal', 'value'),
    State('sel-input-iteration-modal', 'value'),
    running=[
        (Output('sel-button-submitProject-modal', 'loading'), True, False),
    ],
    prevent_initial_call=True
)
def submit_project(nc, project, result, resolution, iteration):
    """
    提交项目
    """
    if not nc:
        return no_update, no_update
    if not project:
        set_head_notice('Please select a project first!', type='warning')
        return no_update, True
    if not result:
        set_head_notice('Please select a result first!', type='warning')
        return no_update, True
    project_name = create_project(project, result, resolution, iteration)
    set_head_notice('Created successfully!', type='success')
    return project_name, False

@callback(
    Output('sel-input-resolution-modal', 'value'),
    Output('sel-input-resolution-modal', 'disabled'),
    Output('sel-input-iteration-modal', 'value'),
    Output('sel-input-iteration-modal', 'disabled'),
    Input('sel-switch-clustering-modal', 'checked')
)
def update_resolution_input_disabled(checked):
    """
    更新分辨率输入框是否禁用
    """
    if checked:
        return 0.8, False, 0, False
    else:
        return None, True, None, True

@callback(
    Output('sel-select-result-modal', 'options'),
    Input('sel-select-project-modal', 'value')
)
def update_project_result_options(project):
    """
    更新分割结果选项
    """
    if project is None:
        return no_update
    result_options = get_result_options(project)
    return result_options

@callback(
    Output('sel-modal-newproject', 'visible'),
    Output('sel-select-project-modal', 'options'),
    Output('sel-select-project-modal', 'value'),
    Output('sel-select-result-modal', 'value'),
    Output('sel-switch-clustering-modal', 'checked'),
    Input('sel-button-newProject', 'nClicks')
)
def open_newproject_modal(nClicks):
    """
    打开新建项目对话框
    """
    if nClicks and verify_modify_permission():
        projects = get_seg_projects()
        return True, projects, None, None, False
    else:
        return no_update, no_update, no_update, no_update, False  