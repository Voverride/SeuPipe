from controller.annotation_ctl import *
from dash import callback, Input, Output, State, no_update, Patch
from websocket.message import ms
from pages.components.fileSelecter import fileSelecter
from controller.auth import *

@callback(
    Input('ann-button-exportData', 'nClicks'),
    State('ann-select-project', 'value'),
)
def export_data(nc, project_name):
    """
    导出数据
    """
    if nc:
        if not project_name:
            set_head_notice('Please select a project!', type='warning')
            return
        parameter_path = annData.get_parameter_path(project_name)
        if not os.path.exists(parameter_path):
            set_head_notice(f'Project {project_name} has not been annotated yet.', type='warning')
            return
        annData.set_export_project(project_name)
        fileSelecter.open_export_box()

@callback(
    Output('ann-graph-result', 'figure', allow_duplicate=True),
    Input('ann-button-slicer', 'nClicks'),
    State('ann-select-project', 'value'),
    State('ann-slider-slicer', 'value'),
    State('ann-select-spotSize', 'value'),
    State('ann-select-borderWidth', 'value'),
    State('ann-colorPicker-boarderColor', 'value'),
    prevent_initial_call=True
)
def slicer(nc, project_name, slicer_range, spot_size, border_width, border_color):
    """
    筛选切片结果
    """
    if nc:
        z_min = slicer_range[0]
        z_max = slicer_range[1]
        fig = get_slicer_result(project_name, z_min, z_max, spot_size, border_width, border_color)
        return fig
        
@callback(
    Input('ann-select-spotSize', 'value'),
    Input('ann-select-borderWidth', 'value'),
    Input('ann-colorPicker-boarderColor', 'value'),
    State('ann-select-project', 'value'),
)
def update_spot_property(spot_size, border_width, border_color, project_name):
    """
    更新散点属性
    """
    if not project_name:
        return
    classes = annData.get_classes(project_name)
    if classes is not None:
        patch = Patch()
        for i in range(classes):
            patch['data'][i]['marker']['size'] = spot_size
            patch['data'][i]['marker']['line']['width'] = border_width
            patch['data'][i]['marker']['line']['color'] = border_color
        set_props('ann-graph-result', dict(figure=patch))

@callback(
    Input('ann-button-showBug', 'nClicks'),
    State('ann-select-project', 'value'),
    prevent_initial_call=True
)
def show_bug(nc, project):
    """
    显示注释异常信息
    """
    if nc and project:
        show_bug_info(project)

@callback(
    Output('ann-select-project', 'value', allow_duplicate=True),
    Input('ann-delete-project-confirm', 'confirmCounts'),
    State('ann-select-project', 'value'),
    running=[
        (Output('ann-delete-project', 'loading'), True, False),
    ],
    prevent_initial_call=True
)
def delete_project(nc, project_name):
    """
    删除项目
    """
    if nc and project_name and verify_modify_permission():
        if annData.is_running(project_name):
            set_head_notice(f'{project_name} is running, please wait...!', type='warning')
            return no_update
        annData.delete_project(project_name)
        set_head_notice('Delete successfully!', type='success')
        return None
    return no_update

@callback(
    Input('refresh-ann-status', 'children'),
    State('ann-select-project', 'value'),
    prevent_initial_call=True
)
def refresh_annotation_status(_, project):
    """
    刷新细胞注释状态
    """
    if project:
        set_props('ann-select-project', dict(value=project))

@callback(
    Input('ann-start-project', 'nClicks'),
    State('ann-select-project', 'value'),
)
def start_project(nc, projectname):
    """
    启动注释项目
    """
    if nc and verify_modify_permission():
        if not projectname:
            set_head_notice('Please Select Project Name', type='warning')
            return
        start_ann_project(projectname)

@callback(
    Input('ann-select-project', 'value'),
    State('ann-select-spotSize', 'value'),
    State('ann-select-borderWidth', 'value'),
    State('ann-colorPicker-boarderColor', 'value')
)
def update_project_tabledata(projectname, spot_size, border_width, border_color):
    """
    更新项目元数据
    """
    ms.clinetSend(ms._updateParams_c2s, params=dict(project=projectname))
    update_project_metadata(projectname, spot_size, border_width, border_color)

@callback(
    Input('ann-button-submit', 'nClicks'),
    State('ann-store-refdata', 'data'),
    State('ann-store-querydata', 'data'),
    State('ann-select-label', 'value'),
    State('ann-select-x', 'value'),
    State('ann-select-y', 'value'),
    State('ann-select-z', 'value'),
    State('ann-remove-mt', 'checked'),
    State('ann-remove-ribo', 'checked'),
    State('ann-remove-hb', 'checked'),
    State('ann-use-hvg', 'checked'),
    State('ann-nlayers', 'value'),
    State('ann-nhiddens', 'value'),
    State('ann-nlatent', 'value'),
    State('ann-epochs', 'value'),
    State('ann-batchsize', 'value'),
    State('ann-dropout', 'value'),
    State('ann-projectname', 'value'),
    running=[
        (Output('ann-button-submit', 'loading'), True, False)
    ],
    prevent_initial_call=True
)
def submit_project(nc, refdata, querydata, label, x, y, z, rm_mt, rm_ribo, rm_hb, use_hvg, n_layers, n_hiddens, n_latent, epochs, batch_size, dropout, projectname):
    """
    提交创建项目
    """
    if nc and verify_modify_permission():
        if not refdata:
            set_head_notice('Please Upload Reference Data', type='warning')
            return
        if not querydata:
            set_head_notice('Please Upload Query Data', type='warning')
            return
        if not label:
            set_head_notice('Please Select Label', type='warning')
            return
        if not x:
            set_head_notice('Please Select X', type='warning')
            return
        if not y:
            set_head_notice('Please Select Y', type='warning')
            return
        if not z:
            set_head_notice('Please Select Z', type='warning')
            return
        if not projectname:
            set_head_notice('Please Input Project Name', type='warning')
            return
        create_project(refdata, querydata, label, x, y, z, rm_mt, rm_ribo, rm_hb, use_hvg, n_layers, n_hiddens, n_latent, epochs, batch_size, dropout, projectname)

@callback(
    Output('ann-select-project', 'options'),
    Input('ann-select-project-tooltip', 'open'),
    prevent_initial_call=True
)
def update_project_list(open):
    """
    更新项目列表
    """
    if open:
        project_list = annData.get_exist_projects()
        return project_list
    return no_update

@callback(
    Input('ann-button-refdata', 'nClicks'),
    prevent_initial_call=True
)
def upload_refdata(nc):
    """
    打开上传参考数据窗口
    """
    if nc and verify_modify_permission():
        fileSelecter.set_annotask('ref')
        fileSelecter.open_import_box()
    
@callback(
    Input('ann-button-querydata', 'nClicks'),
    prevent_initial_call=True
)
def upload_querydata(nc):
    """
    打开上传查询数据窗口
    """
    if nc and verify_modify_permission():
        fileSelecter.set_annotask('query')
        fileSelecter.open_import_box()


@callback(
    Output('ann-projectname', 'status'),
    Input('ann-projectname', 'value'),
    prevent_initial_call=True
)
def update_projectname_status(value):
    """
    更新项目名称状态
    """
    if value is None or value == '':
        return 'error'
    return None

@callback(
    Input('ann-button-newproject', 'nClicks'),
    prevent_initial_call=True
)
def open_create_project_modal(nc):
    """
    显示新建项目面板
    """
    if nc and verify_modify_permission():
        open_new_project_modal()