from controller.annotation_ctl import *
from dash import callback, Input, Output, State, no_update, Patch
from dash.exceptions import PreventUpdate
from pages.components.fileSelecter import fileSelecter
from controller.auth import *

@callback (
    Output('annotask-button-showBug', 'id'),
    Input('annotask-button-showBug', 'nClicks'),
    prevent_initial_call=True
)
def show_bug(nc):
    """
    显示注释进程运行bug
    """
    anntstatus = annData.get_anntstatus()
    e = anntstatus['exception']
    if nc and e:
        raise e
    return no_update
@callback(
    Input('annotation-button-exportData', 'nClicks')
)
def export_data(nc):
    """
    导出数据
    """
    if nc is not None:
        status = check_result_data()
        if status:
            fileSelecter.open_export_box()
@callback(
    Output('annotation-graph-result', 'figure', allow_duplicate=True),
    Input('annotation-button-slicer', 'nClicks'),
    State('annotation-slider-slicer', 'value'),
    prevent_initial_call=True
)
def slice_graph(nc, values):
    """
    裁剪切片
    """
    if nc is None:
        return no_update
    result = annData.get_resultfig()
    if result is None:
        return no_update
    min_value, max_value = values
    fig = plot_3d_scatter(annData, min_z=min_value, max_z=max_value)
    return fig

@callback(
    Output('annotation-graph-result', 'figure', allow_duplicate=True),
    # Output('annotation-graph-refumap', 'figure', allow_duplicate=True),
    Output('annotation-graph-queryumap', 'figure', allow_duplicate=True),
    Input('annotation-select-spotSize', 'value'),
    Input('annotation-select-borderWidth', 'value'),
    Input('annotation-colorPicker-boarderColor', 'value'),
    prevent_initial_call=True
)
def update_spot_style(spotSize, borderWidth, borderColor):
    """
    调整散点大小和边框
    """
    # 'marker': {'color': '#6495ED', 'line': {'color': '#0d0015', 'width': 1}, 'size': 5, 'symbol': 'circle'}
    resultfig = annData.get_resultfig()
    # refumap = annData.get_refumap()
    queryumap = annData.get_queryumap()
    result_patch = no_update
    # ref_patch = no_update
    query_patch = no_update
    def update_patch(fig, patch, use_border=False):
        for i in range(len(fig['data'])):
            patch['data'][i]['marker']['size'] = spotSize
            if use_border:
                patch['data'][i]['marker']['line']['width'] = borderWidth
                patch['data'][i]['marker']['line']['color'] = borderColor

    if resultfig is not None:
        result_patch = Patch()
        update_patch(resultfig, result_patch, use_border=True)

    # if refumap is not None:
    #     ref_patch = Patch()
    #     update_patch(refumap, ref_patch)

    if queryumap is not None:
        query_patch = Patch()
        update_patch(queryumap, query_patch)
    
    return result_patch, query_patch


@callback(
    Output('annotation-interval', 'disabled', allow_duplicate=True),
    Input('annotation-event-loop', 'n_intervals'),
    State('annotation-interval', 'disabled'),
    prevent_initial_call=True
)
def annotation_event_loop(_, disabled):
    """
    当有用户开始训练任务时，其他用户实时同步状态
    """
    if disabled:
        anntstatus = annData.get_anntstatus()
        thread = anntstatus['thread']
        creator = anntstatus['creator']
        if thread is None or not thread.is_alive():
            return True
        reset_annstatus_progress(creator)
        return False
    return no_update

@callback(
    Output('annotation-interval', 'disabled', allow_duplicate=True),
    Input('annotation-interval', 'n_intervals'),
    prevent_initial_call=True
)
def update_annotation_status(_):
    """
    轮回查询注释状态
    """
    anntstatus = annData.get_anntstatus()
    thread = anntstatus['thread']
    update_annotation_progress()
    if thread is None or not thread.is_alive():
        return True
    return no_update

@callback(
    Output('annotation-interval', 'disabled'),
    Input('anntask-button-train', 'nClicks'),
    State('annotask-remove-mt', 'checked'),
    State('annotask-remove-ribo', 'checked'),
    State('annotask-remove-hb', 'checked'),
    State('annotask-use-hvg', 'checked'),
    State('annotask-nlayers', 'value'),
    State('annotask-nhiddens', 'value'),
    State('annotask-nlatent', 'value'),
    State('annotask-epochs', 'value'),
    State('annotask-batchsize', 'value'),
    State('annotask-dropout', 'value'),
    State('main-title-username', 'children'),
    prevent_initial_call=True
)
def start_training(nc, rm_mt, rm_ribo, rm_hb, use_hvg, n_layers, n_hiddens, n_latent, epochs, batch_size, dropout, usrname):
    """
    开始训练
    """
    if nc:
        access = verify_modify_permission()
        if access:
            run_scvi(rm_mt, rm_ribo, rm_hb, use_hvg, n_layers, n_hiddens, n_latent, epochs, batch_size, dropout, usrname)
            return False
    return no_update

# tosica代码
# @callback(
#     Output('annotation-interval', 'disabled'),
#     Input('anntask-button-train', 'nClicks'),
#     State('annotask-remove-mt', 'checked'),
#     State('annotask-remove-ribo', 'checked'),
#     State('annotask-remove-hb', 'checked'),
#     State('annotask-use-hvg', 'checked'),
#     State('annotask-epoch', 'value'),
#     State('annotask-depth', 'value'),
#     State('annotask-batchsize', 'value'),
#     State('annotask-lr', 'value'),
#     State('annotask-gmt', 'value'),
#     State('main-title-username', 'children'),
#     prevent_initial_call=True
# )
# def start_training(nc, rm_mt, rm_ribo, rm_hb, use_hvg, epoch, depth, batchsize, lr, gmt, usrname):
#     """
#     开始训练
#     """
#     if nc:
#         access = verify_modify_permission()
#         if access:
#             pass
#             # run_tosica(rm_mt, rm_ribo, rm_hb, use_hvg, epoch, depth, batchsize, lr, gmt, usrname)
#             return False
#     return no_update

@callback(
    Input('init-restore-annotation', 'n_intervals'),
)
def update_init_component(_):
    """
        恢复初始数据
    """
    restore_initial_data()
@callback(
    Output('annotation-slider-slicer', 'marks'),
    Input('annotation-slider-slicer', 'value'),
)
def update_slider_markers(value):
    """
        更新slider显示数值
    """
    if value:
        return {val:val for val in value}
    raise PreventUpdate
@callback(
    Output('annotask-select-z', 'value'),
    Input('annotask-select-z', 'value'),
    prevent_initial_call=True
)
def set_zfield(value):
    """
    设置z坐标并检查数据类型重置slicer滑动条数据
    """
    zfield = annData.get_zfield()
    if zfield==value:
        raise PreventUpdate
    access = verify_modify_permission()
    if value is not None and access:
        status = check_queryfield_type(value)
        if status:
            annData.set_zfield(value)
            zfield_min = annData.get_zfield_min()
            zfield_max = annData.get_zfield_max()
            if zfield_min is not None and zfield_max is not None:
                props = dict()
                props['min'] = zfield_min
                props['max'] = zfield_max
                props['value'] = [zfield_min, zfield_max]
                props['tooltipPrefix'] = f'{annData.get_zfield()}: '
                set_props('annotation-slider-slicer', props)
            return no_update
    return zfield

@callback(
    Output('annotask-select-y', 'value'),
    Input('annotask-select-y', 'value'),
    prevent_initial_call=True
)
def set_yfield(value):
    """
    设置y坐标并检查数据类型
    """
    yfield = annData.get_yfield()
    if yfield==value:
        raise PreventUpdate
    access = verify_modify_permission()
    if value is not None and access:
        status = check_queryfield_type(value)
        if status:
            annData.set_yfield(value)
            return no_update
    return yfield

@callback(
    Output('annotask-select-x', 'value'),
    Input('annotask-select-x', 'value'),
    prevent_initial_call=True
)
def set_xfield(value):
    """
    设置x坐标并检查数据类型
    """
    xfield = annData.get_xfield()
    if xfield==value:
        raise PreventUpdate
    access = verify_modify_permission()
    if value is not None and access:
        status = check_queryfield_type(value)
        if status:
            annData.set_xfield(value)
            return no_update
    return xfield

@callback(
    Output('annotask-select-label', 'value'),
    Input('annotask-select-label', 'value'),
    prevent_initial_call=True
)
def set_labelfield(value):
    """
    设置label坐标
    """
    labelfield = annData.get_labelfield()
    if labelfield==value:
        raise PreventUpdate
    access = verify_modify_permission()
    if value is not None and access:
        annData.set_labelfield(value)
        return no_update
    return labelfield

@callback(
    Input('annotask-button-querydata', 'nClicks'),
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
    Input('annotask-button-refdata', 'nClicks'),
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
    Output('annotask-drawer', 'visible'),
    Input('annotation-button-newtask', 'nClicks'),
    prevent_initial_call=True
)
def open_annotask_drawer(nc):
    """
    弹出新建注释任务抽屉
    """
    if nc:
        return True
    raise PreventUpdate