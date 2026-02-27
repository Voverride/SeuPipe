from pages.components.fileSelecter import fileSelecter
from controller.alignment_ctl import *
from websocket.websocket import ms
from dash import dcc
from dash.exceptions import PreventUpdate
from dash import Input, Output, State, no_update, callback, Patch, ctx
from controller.auth import verify_modify_permission

@callback(
    Output('ma-icon-autoPickAbove', 'icon'),
    Input('ma-button-autoPickAbove', 'nClicks'),
    State('ma-icon-autoPickAbove', 'icon')
)
def update_autoPickAbove_icon(click, icon):
    """
    更新自动选择上方切片图标
    """
    if click:
        if icon=='antd-check':
            return 'antd-close'
        return 'antd-check'
    return no_update

@callback(
    Output('ma-icon-autoPickBelow', 'icon'),
    Input('ma-button-autoPickBelow', 'nClicks'),
    State('ma-icon-autoPickBelow', 'icon'),
)
def update_autoPickBelow_icon(click, icon):
    """
    更新自动选择下方切片图标
    """
    if click:
        if icon=='antd-check':
            return 'antd-close'
        return 'antd-check'
    return no_update


@callback(
    Input('ma-button-pickAbove', 'nClicks'),
    Input('ma-button-pickBelow', 'nClicks'),
    State('ma-selecter-activeSlice', 'value'),
    State('ma-table-SyncSlices', 'selectedRowKeys'),
    State('ali-select-project', 'value'),
    State('ma-slider-znum', 'data'),
)
def update_selected_SyncSlices(pickAbove, pickBelow, activeSlice, selectedRowKeys, project, znum):
    """
    选中activeSlice上方或下方切片
    """
    if not pickAbove and not pickBelow:
        return
    if not activeSlice:
        set_head_notice('Please ensure that the active slice was selected !', type='warning')
        return
    tid = ctx.triggered_id
    idx = alidata.get_syncslice_index(project, float(activeSlice))
    if selectedRowKeys==None:
        selectedRowKeys = []
    if tid=='ma-button-pickAbove':
        for i in range(0, idx+1):
            rowkey = str(i)
            if rowkey not in selectedRowKeys:
                selectedRowKeys.append(rowkey)
    if tid=='ma-button-pickBelow':
        for i in range(idx, znum):
            rowkey = str(i)
            if rowkey not in selectedRowKeys:
                selectedRowKeys.append(rowkey)
    set_props('ma-table-SyncSlices', dict(selectedRowKeys=selectedRowKeys))

@callback(
    Output('ali-graph-left', 'figure', allow_duplicate=True),
    Output('ali-graph-right', 'figure', allow_duplicate=True),
    Input("ali-graph-left", "relayoutData"),
    Input("ali-graph-right", "relayoutData"),
    State('ali-icon-add-contrast', 'icon'),
    State('ali-store-figureScale', 'data'),
    prevent_initial_call=True
)
def update_relayout(aliLayout, oriLayout, icon, figureScale):
    """
    同步两图状态
    """
    if icon=='antd-add':
        raise PreventUpdate

    oriScale = figureScale.get('initScale', 1)
    aliScale = figureScale.get('resultScale', 1)

    tid = ctx.triggered_id
    patchOri = Patch()
    patchAli = Patch()
    if tid == 'ali-graph-right':
        mark = False
        if 'scene.camera' in oriLayout:
            patchAli['layout']['scene']['camera']['eye'] = oriLayout['scene.camera']['eye']
            patchOri['layout']['scene']['camera']['eye'] = oriLayout['scene.camera']['eye']
            mark = True
        if 'scene.aspectratio' in oriLayout:
            factor = aliScale/oriScale
            patchAli['layout']['scene']['aspectmode'] = 'manual'
            patchOri['layout']['scene']['aspectmode'] = 'manual'
            aspectratio = oriLayout['scene.aspectratio'].copy()
            aspectratio['x']*=factor
            aspectratio['y']*=factor
            patchAli['layout']['scene']['aspectratio'] = aspectratio
            patchOri['layout']['scene']['aspectratio'] = oriLayout['scene.aspectratio']
            mark = True
        if mark:
            return patchAli, patchOri
    if tid == 'ali-graph-left':
        mark = False
        if 'scene.camera' in aliLayout:
            patchOri['layout']['scene']['camera']['eye'] = aliLayout['scene.camera']['eye']
            patchAli['layout']['scene']['camera']['eye'] = aliLayout['scene.camera']['eye']
            # z_dire = aliLayout['scene.camera']['eye']['z']
            # alidata.set_actref_slice(usrname, reverse=False if z_dire>=0 else True)
            mark = True
        if 'scene.aspectratio' in aliLayout:
            factor = oriScale/aliScale
            patchOri['layout']['scene']['aspectmode'] = 'manual'
            patchAli['layout']['scene']['aspectmode'] = 'manual'
            aspectratio = aliLayout['scene.aspectratio'].copy()
            aspectratio['x']*=factor
            aspectratio['y']*=factor
            patchOri['layout']['scene']['aspectratio'] = aspectratio
            patchAli['layout']['scene']['aspectratio'] = aliLayout['scene.aspectratio']
            mark = True
        if mark:
            return patchAli, patchOri
    raise PreventUpdate

@callback(
    Output('ali-graph-left', 'figure', allow_duplicate=True),
    Output('ma-button-showGrid', 'style'),
    Input('ma-button-showGrid', 'nClicks'),
    State('ma-button-showGrid', 'style'),
    prevent_initial_call=True
)
def switch_grid(nClicks, style):
    """
    开关图像网格
    """
    if not nClicks:
        raise PreventUpdate
    patch = no_update
    status = ['#949495', '#ca8269'] # 0 -> 关， 1-> 开
    patch = Patch()
    if style['backgroundColor']==status[0]:
        style['backgroundColor'] = status[1]
        patch['layout']['scene']['xaxis']['visible'] = True
        patch['layout']['scene']['yaxis']['visible'] = True
        patch['layout']['scene']['zaxis']['visible'] = True
    else:
        style['backgroundColor'] = status[0]
        patch['layout']['scene']['xaxis']['visible'] = False
        patch['layout']['scene']['yaxis']['visible'] = False
        patch['layout']['scene']['zaxis']['visible'] = False
    return patch, style

@callback(
    Output('ali-graph-left', 'figure', allow_duplicate=True),
    Input('ma-button-reFocus', 'nClicks'),
    State('ma-selecter-activeSlice', 'value'),
    State('ma-selecter-referenceSlice', 'value'),
    prevent_initial_call=True
)
def reset_view(nClicks, actSlice, refSlice):
    """
    重置图像视角
    """
    patch = no_update
    if nClicks:
        patch = Patch()
        patch['layout']['scene']['camera']['eye'] = alidata.get_graph_view(actSlice, refSlice)
    return patch

@callback(
    Input('ma-selecter-referenceSlice', 'value'),
    Input('ma-ColorPicker-refSlice', 'value'),
    State('ali-icon-add-contrast', 'icon'),
    State('ali-select-project', 'value'),
    State('ma-selecter-activeSlice', 'value'),
)
def reference_slice_change(slice, color, icon, project, actSlice):
    """
    更新参考切片显示配置
    """
    if slice:
        update_reference_slice(project, slice, actSlice, color, icon)

@callback(
    Input('ma-selecter-activeSlice', 'value'),
    Input('ma-ColorPicker-actSlice', 'value'),
    State('ali-icon-add-contrast', 'icon'),
    State('ali-select-project', 'value'),
    State('ma-selecter-referenceSlice', 'value'),
    State('ma-icon-autoPickAbove', 'icon'),
    State('ma-icon-autoPickBelow', 'icon'),
    State('ma-slider-znum', 'data'),
)
def active_slice_change(slice, color, icon, project, refSlice, autoPickAbove, autoPickBelow, znum):
    """
    更新当前切片显示配置
    """
    tid = ctx.triggered_id
    if tid=='ma-selecter-activeSlice':
        if not slice:
            selectedRowKeys = []
        else:
            idx = alidata.get_syncslice_index(project, float(slice))
            selectedRowKeys = [str(idx)]
            if autoPickAbove=='antd-check':
                for i in range(0, idx):
                    rowkey = str(i)
                    selectedRowKeys.append(rowkey)
            if autoPickBelow=='antd-check':
                for i in range(idx+1, znum):
                    rowkey = str(i)
                    selectedRowKeys.append(rowkey)
        set_props('ma-table-SyncSlices', dict(selectedRowKeys=selectedRowKeys))

    if slice:
        update_active_slice(project, slice, refSlice, color, icon)

@callback(
    Input('ma-button-slicer', 'nClicks'),
    State('ma-slider-znums', 'data'),
    State('ma-slider-z', 'value'),
    State('ali-icon-add-contrast', 'icon'),
)
def update_slicer(nc, znums, zvalue, icon):
    """
    更新切片配置
    """
    if nc:
        zmin = min(zvalue)
        zmax = max(zvalue)
        slicer_layers(zmin, zmax, znums, icon)
        

@callback(
    Input('ma-selecter-spotSize', 'value'),
    Input('ma-selecter-borderWidth', 'value'),
    Input('ma-ColorPicker-boarderColor', 'value'),
    State('ma-slider-znum', 'data')
)
def update_spot(spotSize, borderWidth, borderColor, znum):
    """
    更新点属性
    """
    update_spot_properties(spotSize, borderWidth, borderColor, znum)

@callback(
    Output('ali-splitter-figure', 'items'),
    Output('ali-icon-add-contrast', 'icon'),
    Output('ali-content-tabs', 'activeKey', allow_duplicate=True),
    Input('ali-button-add-contrast', 'nClicks'),
    State('ali-select-project', 'value'),
    State('ali-icon-add-contrast', 'icon'),
    State('ma-selecter-spotSize', 'value'),
    State('ma-selecter-borderWidth', 'value'),
    State('ma-ColorPicker-boarderColor', 'value'),
    State('ma-slider-z', 'value'),
    running=[
        (Output('ali-button-add-contrast', 'loading'), True, False),
    ],
    prevent_initial_call=True
)
def add_contrast(nc, project, icon, spotSize, borderWidth, borderColor, zvalue):
    """
    添加对比图像
    """
    if nc:
        if not project:
            set_head_notice('Please Select One Project', type='warning')
            return no_update, no_update, no_update
        patch = Patch()
        activeKey = 'Figure'
        if icon == 'antd-plus':
            icon = 'antd-minus'
            figure = get_contrast_figure(project)
            scale = figure['layout']['scene']['aspectratio']['x']
            scalePatch = Patch()
            scalePatch['initScale'] = scale
            set_props('ali-store-figureScale', scalePatch)
            update_figure_spot(figure, spotSize, borderWidth, borderColor)
            zmin = min(zvalue)
            zmax = max(zvalue)
            for trace in figure['data']:
                trace['visible'] = float(trace['z'][0]) >= zmin and float(trace['z'][0]) <= zmax
            patch.append({
                'children': fac.AntdCenter(
                    dcc.Graph(
                        figure=figure,
                        config={'displaylogo':False}, 
                        id='ali-graph-right', 
                        style={'height': '100%', 'width': '100%'},
                    ),
                    style={'height': '100%', 'width': '100%'}
                ),
                'collapsible': True,
            })
        else:
            icon = 'antd-plus'
            del patch[-1]
        return patch, icon, activeKey
    return no_update, no_update, no_update

@callback(
    Input('ali-button-showBug', 'nClicks'),
    State('ali-select-project', 'value'),
    prevent_initial_call=True
)
def show_bug(nc, project):
    """
    显示对齐异常信息
    """
    if nc and project:
        show_bug_info(project)

@callback(
    Input('refresh-ali-status', 'children'),
    State('ali-select-project', 'value'),
    State('ma-selecter-spotSize', 'value'),
    State('ma-selecter-borderWidth', 'value'),
    State('ma-ColorPicker-boarderColor', 'value'),
    prevent_initial_call=True
)
def refresh_alignment_status(_, project, spotSize, borderWidth, borderColor):
    """
    刷新对齐状态
    """
    if project:
        update_project_metadata(project, spotSize, borderWidth, borderColor)

@callback(
    Input('ali-start-project', 'nClicks'),
    State('ali-select-project', 'value'),
    State('ali-select-model', 'value'),
    State('ali-radio-device', 'value'),
)
def start_project(nc, projectname, model, device):
    """
    启动对齐项目
    """
    if nc and verify_modify_permission():
        if not projectname:
            set_head_notice('Please Select One Project', type='warning')
            return
        start_alignment_project(projectname, model, device)

@callback(
    Output('ali-select-project', 'value', allow_duplicate=True),
    Input('ali-delete-project-confirm', 'confirmCounts'),
    State('ali-select-project', 'value'),
    running=[
        (Output('ali-delete-project', 'loading'), True, False),
    ],
    prevent_initial_call=True
)
def delete_project(nc, project_name):
    """
    删除项目
    """
    if nc and project_name and verify_modify_permission():
        if alidata.is_running(project_name):
            set_head_notice(f'{project_name} is running, please wait...!', type='warning')
            return no_update
        alidata.delete_project(project_name)
        set_head_notice('Delete successfully!', type='success')
        return None
    return no_update

@callback(
    Input('ali-select-project', 'value'),
    State('ali-icon-add-contrast', 'icon'),
    State('ma-selecter-spotSize', 'value'),
    State('ma-selecter-borderWidth', 'value'),
    State('ma-ColorPicker-boarderColor', 'value'),
    running=[
        (Output('ali-select-project', 'disabled'), True, False)
    ],
)
def update_project_tabledata(project, icon, spotSize, borderWidth, borderColor):
    """
    更新项目元数据
    """
    ms.clinetSend(ms._updateParams_c2s, params=dict(project=project))
    if icon == 'antd-minus':
        if project:
            figure = get_contrast_figure(project)
            scale = figure['layout']['scene']['aspectratio']['x']
            scalePatch = Patch()
            scalePatch['initScale'] = scale
            set_props('ali-store-figureScale', scalePatch)
            update_figure_spot(figure, spotSize, borderWidth, borderColor)
        else:
            figure = None
        set_props('ali-graph-right', dict(figure=figure))
    update_project_metadata(project, spotSize, borderWidth, borderColor, userTriggered=True)

@callback(
    Input('ali-button-importdata', 'nClicks'),
    prevent_initial_call=True
)
def open_import_data(nc):
    """
    打开数据导入面板
    """
    if nc and verify_modify_permission():
        fileSelecter.open_import_box()


@callback(
    Input('ali-button-newproject', 'nClicks'),
    prevent_initial_call=True
)
def open_new_project(nc):
    """
    打开新建项目面板
    """
    if nc and verify_modify_permission():
        open_new_project_modal()

@callback(
    Output('ali-projectname', 'status'),
    Input('ali-projectname', 'value'),
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
    Input('ali-button-submit', 'nClicks'),
    State('ali-store-importdata', 'data'),
    State('ali-select-x', 'value'),
    State('ali-select-y', 'value'),
    State('ali-select-z', 'value'),
    State('ali-projectname', 'value'),
    running=[(Output('ali-button-submit', 'loading'), True, False)],
    prevent_initial_call=True
)
def submit_project(nc, importdata, x, y, z, projectname):
    """
    提交新建项目
    """
    if nc and verify_modify_permission():
        if not projectname:
            set_head_notice('Please Input Project Name', type='warning')
            return
        if not importdata:
            set_head_notice('Please Import Data First', type='warning')
            return
        if not x or not y or not z:
            set_head_notice('Please Select X, Y, Z Fields', type='warning')
            return
        create_project(importdata, x, y, z, projectname)


@callback(
    Output('ali-select-project', 'options'),
    Input('ali-select-project-tooltip', 'open'),
    prevent_initial_call=True
)
def update_project_options(open):
    """
    更新项目下拉列表选项
    """
    if open:
        return update_project_list()
    return no_update

@callback(
    Output('ma-slider-z', 'marks'),
    Input('ma-slider-z', 'value'),
)
def update_slider_markers(value):
    """
    更新slider显示数值
    """
    if value:
        return {val:val for val in value}
    raise PreventUpdate

@callback(
    Output('gs-drawer-manual', 'visible'),
    Output('ali-content-tabs', 'activeKey'),
    Input('ali-button-manualAdjust', 'nClicks'),
    State('ali-select-project', 'value'),
)
def open_manual_adjust_drawer(nc, project):
    """
    打开手动调整抽屉
    """
    if nc:
        if not project:
            set_head_notice('Please Select One Project', type='warning')
            return no_update, no_update
        if alidata.is_running(project):
            set_head_notice(f'{project} is running, please wait...!', type='warning')
            return no_update, no_update
        return True, 'Figure'
    return no_update, no_update