from pages.components.fileSelecter import fileSelecter
from controller.alignment_ctl import *
from websocket.websocket import ms
from dash.exceptions import PreventUpdate
from dash import Input, Output, State, no_update, callback, Patch, ctx
from controller.auth import verify_modify_permission
import plotly.express as px
import pandas as pd

hiddenDf = pd.DataFrame({
    'x': [5],
    'y': [3],
    'z': [8]
})

hiddenGraph = px.scatter_3d(hiddenDf, x='x', y='y', z='z')
hiddenGraph.update_layout(
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)',
)

@callback(
    Input('ma-selecter-colorMode', 'value'),
    State('ma-selecter-activeSlice', 'value'),
    State('ma-selecter-referenceSlice', 'value'),
    State('ma-ColorPicker-actSlice', 'value'),
    State('ma-ColorPicker-refSlice', 'value'),
    State('ma-selecter-colorField', 'value'),
    State('ma-selecter-colorGene', 'value'),
    State('ma-selecter-colorRange', 'value'),
    State('ali-select-project', 'value'),
    State('ali-icon-add-contrast', 'icon'),
    running=[
        (Output('ma-selecter-colorMode', 'disabled'), True, False)
    ],
    prevent_initial_call=True
)
def colormode_change(colorMode, actSlice, refSlice, actColor, refColor, colorField, colorGene, colorRange, project, icon):
    """
    颜色模式改变
    """
    patch = Patch()
    if colorMode=='Field' and colorField:
        set_color_byfield(patch, project, colorField)
        set_props('ali-graph-left', dict(figure=patch))
        if icon=='antd-minus':
            set_props('ali-graph-right', dict(figure=patch))
        return
    if colorMode=='Gene' and colorGene:
        set_color_bygene(patch, project, colorGene, colorRange)
        set_props('ali-graph-left', dict(figure=patch))
        if icon=='antd-minus':
            set_props('ali-graph-right', dict(figure=patch))
        return
    elif colorMode=='Custom':
        if refSlice and refColor:
            i = alidata.get_slice_index(project, refSlice)
            patch['data'][i]['marker']['color'] = refColor
        if actSlice and actColor:
            i = alidata.get_slice_index(project, actSlice)
            patch['data'][i]['marker']['color'] = actColor
        set_props('ali-graph-left', dict(figure=patch))
        if icon=='antd-minus':
            set_props('ali-graph-right', dict(figure=patch))

@callback(
    Input('ma-selecter-colorGene', 'value'),
    Input('ma-selecter-colorRange', 'value'),
    State('ali-select-project', 'value'),
    State('ali-icon-add-contrast', 'icon'),
    running=[
        (Output('ma-selecter-colorGene', 'disabled'), True, False),
        (Output('ma-selecter-colorRange', 'disabled'), True, False),
    ]
)
def update_colorGene(colorGene, colorRange, project, icon):
    """
    更新基因颜色
    """
    if not colorGene or not colorRange:
        return
    patch = Patch()
    set_color_bygene(patch, project, colorGene, colorRange)
    set_props('ali-graph-left', dict(figure=patch))
    if icon=='antd-minus':
        set_props('ali-graph-right', dict(figure=patch))

@callback(
    Output('ma-div-legend', 'children', allow_duplicate=True),
    Output('ma-div-loadMore', 'style', allow_duplicate=True),
    Input('ma-button-loadMore', 'nClicks'),
    State('ma-div-loadMore', 'style'),
    State('ma-selecter-colorMode', 'value'),
    State('ma-selecter-colorField', 'value'),
    State('ali-select-project', 'value'),
    prevent_initial_call=True
)
def load_more_legend(nClicks, loadMore, colorMode, colorField, project):
    """
    加载更多图例
    """
    if nClicks and colorMode=='Field' and colorField:
        legends, hasLegend = get_field_legends(project, colorField, loadMore['legendNum'])
        patch_style = Patch()
        patch_legend = Patch()
        displayLoadMore = 'none'
        if hasLegend:
            displayLoadMore = 'block'
        patch_style['display'] = displayLoadMore
        patch_style['legendNum'] += len(legends)
        patch_legend+=legends
        return patch_legend, patch_style
    raise PreventUpdate

@callback(
    Output('ma-div-legend', 'children', allow_duplicate=True),
    Output('ma-div-loadMore', 'style', allow_duplicate=True),
    Input('ma-selecter-colorMode', 'value'),
    Input('ma-selecter-colorField', 'value'),
    State('ali-select-project', 'value'),
    prevent_initial_call=True
)
def update_expand_fieldLegend_panel(colorMode, colorField, project):
    """
    刷新field域对应图例内容
    """
    if colorMode=='Field' and colorField:
        legends, hasLegend = get_field_legends(project, colorField, 0)
        patch = Patch()
        displayLoadMore = 'none'
        if hasLegend:
            displayLoadMore = 'block'
        patch['display'] = displayLoadMore
        patch['legendNum'] = len(legends)
        return legends, patch
    raise PreventUpdate

@callback(
    Output('ma-card-legend', 'visible', allow_duplicate=True),
    Output('ma-card-legend', 'title', allow_duplicate=True),
    Input('ma-selecter-colorMode', 'value'),
    Input('ma-selecter-colorField', 'value'),
    Input('gs-drawer-manual', 'visible'),
    prevent_initial_call=True
)
def switch_expand_fieldLegend_panel(colorMode, colorField, drawerVisable):
    """
    开关field域对应图例
    """
    if colorMode=='Field' and drawerVisable and colorField:
        return True, colorField
    return False, no_update

@callback(
    Input('ma-selecter-colorField', 'value'),
    State('ali-select-project', 'value'),
    State('ali-icon-add-contrast', 'icon'),
    running=[
        (Output('ma-selecter-colorField', 'disabled'), True, False),
    ],
)
def update_sliceColor_byField(field, project, icon):
    """
    基于数据域显示切片颜色
    """
    if not field:
        return
    patch = Patch()
    set_color_byfield(patch, project, field)
    set_props('ali-graph-left', dict(figure=patch))
    if icon=='antd-minus':
        set_props('ali-graph-right', dict(figure=patch))

@callback(
    Output('ma-ColorPicker-actSlice', 'disabled'),
    Output('ma-ColorPicker-refSlice', 'disabled'),
    Output('ma-selecter-colorField', 'disabled'),
    Output('ma-selecter-colorGene', 'disabled'),
    Output('ma-selecter-colorRange', 'disabled'),
    Input('ali-button-manualAdjust', 'nClicks'),
    Input('ma-selecter-colorMode', 'value'),
)
def disabled_color_mode(_, value):
    """
    基于colorMode禁用颜色组件
    """
    if value=='Field':
        return True, True, False, True, True
    elif value=='Gene':
        return True, True, True, False, False
    return False, False, True, True, True

@callback(
    Input('ali-button-exportData', 'nClicks'),
    State('ali-select-project', 'value'),
)
def export_data(nc, project_name):
    """
    导出数据
    """
    if nc:
        if not project_name:
            set_head_notice('Please select a project!', type='warning')
            return
        alidata.set_export_project(project_name)
        fileSelecter.open_export_box()

@callback(
    Input('ma-button-up', 'nClicks'),
    Input('ma-button-down', 'nClicks'),
    Input('ma-button-left', 'nClicks'),
    Input('ma-button-right', 'nClicks'),
    Input('ma-button-clockwise', 'nClicks'),
    Input('ma-button-unclockwise', 'nClicks'), 
    State('ma-selecter-activeSlice', 'value'),
    State('ma-selecter-referenceSlice', 'value'),
    State('ma-inputNum-stepSize', 'value'),
    State('ma-inputNum-rotationAngle', 'value'),
    State('ali-select-project', 'value'),
    State('ma-table-SyncSlices', 'selectedRowKeys'),
)
def move_slice(up, down, left, right, clockwise, unclockwise, actSlice, refSlice, stepSize, rotationAngle, project, selectedRowKeys):
    """
    根据方向按钮调整切片位置
    """
    if alidata.is_running(project):
        set_head_notice(f'Project {project} is running, please wait !', type='warning')
        return
    modifyName = alidata.is_project_modify(project)
    if modifyName:
        # set_head_notice(f'{modifyName} is in operation, please wait !', type='warning')
        return
    try:
        alidata.set_project_modify(project)
        buttonPressed = up or down or left or right or clockwise or unclockwise
        if not buttonPressed:
            return
        permission = verify_modify_permission()
        if not permission:
            set_head_notice(f'permission denied !', type='warning')
            return
        elif not actSlice:
            set_head_notice(f'The active slice cannot be empty !', type='warning')
            return
        elif not refSlice:
            set_head_notice(f'The reference slice cannot be empty !', type='warning')
            return
        elif not stepSize:
            set_head_notice(f'Step size cannot be empty !', type='warning')
            return
        elif not rotationAngle:
            set_head_notice(f'Rotation angle cannot be empty !', type='warning')
            return
        tid = ctx.triggered_id
        if tid=='ma-button-up':
            transMtx = alidata.get_coord_trans_mtx(dy=stepSize)
        elif tid=='ma-button-down':
            transMtx = alidata.get_coord_trans_mtx(dy=-stepSize)
        elif tid=='ma-button-left':
            transMtx = alidata.get_coord_trans_mtx(dx=-stepSize)
        elif tid=='ma-button-right':
            transMtx = alidata.get_coord_trans_mtx(dx=stepSize)
        elif tid=='ma-button-clockwise':
            transMtx = alidata.get_coord_trans_mtx(reg=rotationAngle)
        elif tid=='ma-button-unclockwise':
            transMtx = alidata.get_coord_trans_mtx(reg=-rotationAngle)
        if transMtx is None:
            return
        syncSlices = get_sync_slices(project, selectedRowKeys)
        actSlice = float(actSlice)
        syncSlices.add(actSlice)
        get_transformed_coord(project, transMtx, syncSlices)
    except Exception as e:
        raise e
    finally:
        alidata.unset_project_modify(project)

@callback(
    Input('key-pressed-events', 'data'),
    State('ma-selecter-activeSlice', 'value'),
    State('ma-selecter-referenceSlice', 'value'),
    State('ma-inputNum-stepSize', 'value'),
    State('ma-inputNum-rotationAngle', 'value'),
    State('ali-select-project', 'value'),
    State('ma-table-SyncSlices', 'selectedRowKeys'),
    State('gs-drawer-manual', 'visible'),
)
def detect_pressed_key(keyboard, actSlice, refSlice, stepSize, rotationAngle, project, selectedRowKeys, visible):
    """
    绑定键盘事件调整切片
    37 : left
    38 : up
    39 : right
    40 : down
    1037 : ctrl+left
    1038 : ctrl+up
    1039 : ctrl+right
    1040 : ctrl+down
    """
    if not visible:
        return
    if alidata.is_running(project):
        set_head_notice(f'Project {project} is running, please wait !', type='warning')
        return
    modifyName = alidata.is_project_modify(project)
    if modifyName:
        set_head_notice(f'{modifyName} is in operation, please wait !', type='warning')
        return
    try:
        alidata.set_project_modify(project)
        keyCode = int(keyboard) if keyboard.isdigit() else 0
        detectedKeys = {37, 38, 39, 40, 1037, 1038, 1039, 1040}
        if not actSlice or not refSlice or not stepSize or not rotationAngle or keyCode not in detectedKeys:
            raise PreventUpdate
        permission = verify_modify_permission()
        if not permission:
            return
        if keyCode==38:
            transMtx = alidata.get_coord_trans_mtx(dy=stepSize)
        elif keyCode==40:
            transMtx = alidata.get_coord_trans_mtx(dy=-stepSize)
        elif keyCode==37:
            transMtx = alidata.get_coord_trans_mtx(dx=-stepSize)
        elif keyCode==39:
            transMtx = alidata.get_coord_trans_mtx(dx=stepSize)
        elif keyCode==1037 or keyCode==1038:
            transMtx = alidata.get_coord_trans_mtx(reg=-rotationAngle)
        elif keyCode==1039 or keyCode==1040:
            transMtx = alidata.get_coord_trans_mtx(reg=rotationAngle)
        if transMtx is None:
            return
        syncSlices = get_sync_slices(project, selectedRowKeys)
        actSlice = float(actSlice)
        syncSlices.add(actSlice)
        get_transformed_coord(project, transMtx, syncSlices)
    except Exception as e:
        raise e
    finally:
        alidata.unset_project_modify(project)

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
    Input('ali-store-leftLayout', 'data'),
    State('ali-icon-add-contrast', 'icon'),
    State('ali-store-figureScale', 'data'),
)
def update_right_relayout(aliLayout, icon, figureScale):
    """
    同步右图状态
    """
    if aliLayout is None:
        return
    
    patchAli = Patch()
    update_relayoutfig(patchAli, aliLayout)
    set_props('ali-graph-left', dict(figure=patchAli))

    if icon=='antd-plus':
        return

    oriScale = figureScale.get('initScale', 1)
    aliScale = figureScale.get('resultScale', 1)

    patch = Patch()
    
    if 'scene.camera' in aliLayout:
        patch['layout']['scene']['camera']['eye'] = aliLayout['scene.camera']['eye']
    if 'scene.aspectratio' in aliLayout:
        factor = oriScale/aliScale
        patch['layout']['scene']['aspectmode'] = 'manual'
        aspectratio = aliLayout['scene.aspectratio'].copy()
        aspectratio['x']*=factor
        aspectratio['y']*=factor
        patch['layout']['scene']['aspectratio'] = aspectratio
    set_props('ali-graph-right', dict(figure=patch))

@callback(
    Input('ali-store-leftLayout', 'data'),
    State('ali-icon-add-contrast', 'icon'),
    State('ali-store-figureScale', 'data'),
)
def update_right_relayout(aliLayout, icon, figureScale):
    """
    同步右图状态
    """
    if aliLayout is None:
        return
    
    patchAli = Patch()
    update_relayoutfig(patchAli, aliLayout)
    set_props('ali-graph-left', dict(figure=patchAli))

    if icon=='antd-plus':
        return

    oriScale = figureScale.get('initScale', 1)
    aliScale = figureScale.get('resultScale', 1)

    patch = Patch()
    
    if 'scene.camera' in aliLayout:
        patch['layout']['scene']['camera']['eye'] = aliLayout['scene.camera']['eye']
    if 'scene.aspectratio' in aliLayout:
        factor = oriScale/aliScale
        patch['layout']['scene']['aspectmode'] = 'manual'
        aspectratio = aliLayout['scene.aspectratio'].copy()
        aspectratio['x']*=factor
        aspectratio['y']*=factor
        patch['layout']['scene']['aspectratio'] = aspectratio

    set_props('ali-graph-right', dict(figure=patch))

@callback(
    Input('ali-store-rightLayout', 'data'),
    State('ali-store-figureScale', 'data'),
)
def update_left_relayout(oriLayout, figureScale):
    """
    同步左图状态
    """
    if oriLayout is None:
        return

    oriScale = figureScale.get('initScale', 1)
    aliScale = figureScale.get('resultScale', 1)

    patch = Patch()
    patchOri = Patch()
    update_relayoutfig(patchOri, oriLayout) 

    if 'scene.camera' in oriLayout:
        patch['layout']['scene']['camera']['eye'] = oriLayout['scene.camera']['eye']
    if 'scene.aspectratio' in oriLayout:
        factor = aliScale/oriScale
        patch['layout']['scene']['aspectmode'] = 'manual'
        aspectratio = oriLayout['scene.aspectratio'].copy()
        aspectratio['x']*=factor
        aspectratio['y']*=factor
        patch['layout']['scene']['aspectratio'] = aspectratio
    set_props('ali-graph-left', dict(figure=patch))  
    set_props('ali-graph-right', dict(figure=patchOri))

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
    State('ma-selecter-colorMode', 'value'),
    State('ma-selecter-colorField', 'value'),
    State('ma-selecter-colorGene', 'value')
)
def reference_slice_change(slice, color, icon, project, actSlice, colorMode, colorField, colorGene):
    """
    更新参考切片显示配置
    """
    if slice:
        update_reference_slice(project, slice, actSlice, color, icon, colorMode, colorField, colorGene)

@callback(
    Input('ma-selecter-activeSlice', 'value'),
    Input('ma-ColorPicker-actSlice', 'value'),
    State('ali-icon-add-contrast', 'icon'),
    State('ali-select-project', 'value'),
    State('ma-selecter-referenceSlice', 'value'),
    State('ma-icon-autoPickAbove', 'icon'),
    State('ma-icon-autoPickBelow', 'icon'),
    State('ma-slider-znum', 'data'),
    State('ma-selecter-colorMode', 'value'),
    State('ma-selecter-colorField', 'value'),
    State('ma-selecter-colorGene', 'value')
)
def active_slice_change(slice, color, icon, project, refSlice, autoPickAbove, autoPickBelow, znum, colorMode, colorField, colorGene):
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
        update_active_slice(project, slice, refSlice, color, icon, colorMode, colorField, colorGene)

@callback(
    Input('ma-button-slicer', 'nClicks'),
    State('ma-slider-znums', 'data'),
    State('ma-slider-z', 'value'),
    State('ali-icon-add-contrast', 'icon')
)
def update_slicer(nc, znums, zvalue, icon):
    """
    裁剪切片
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
    Output('ma-ColorPicker-boarderColor', 'value', allow_duplicate=True),
    Input('ali-icon-add-contrast', 'icon'),
    State('ma-ColorPicker-boarderColor', 'value'),
    prevent_initial_call=True
)
def refresh_right_graph(icon, borderColor):
    """
    刷新右侧图像
    """
    if icon == 'antd-minus':
        if borderColor=='#0D0015':
            return '#0D0016'
        return '#0D0015'
    return no_update

@callback(
    Output('ali-icon-add-contrast', 'icon'),
    Output('ali-content-tabs', 'activeKey', allow_duplicate=True),
    Input('ali-button-add-contrast', 'nClicks'),
    State('ali-select-project', 'value'),
    State('ali-icon-add-contrast', 'icon'),
    State('ma-selecter-spotSize', 'value'),
    State('ma-selecter-borderWidth', 'value'),
    State('ma-ColorPicker-boarderColor', 'value'),
    State('ma-slider-z', 'value'),
    State('ma-selecter-colorMode', 'value'),
    State('ma-selecter-colorField', 'value'),
    State('ma-selecter-colorGene', 'value'),
    State('ma-selecter-colorRange', 'value'),
    running=[
        (Output('ali-button-add-contrast', 'loading'), True, False),
    ],
    prevent_initial_call=True
)
def add_contrast(nc, project, icon, spotSize, borderWidth, borderColor, zvalue, colorMode, colorField, colorGene, colorRange):
    """
    添加对比图像
    """
    if nc:
        if not project:
            set_head_notice('Please Select One Project', type='warning')
            return no_update, no_update
        patch = Patch()
        activeKey = 'Figure'
        if icon == 'antd-plus':
            icon = 'antd-minus'
            figure = get_contrast_figure(project)
            scale = figure['layout']['scene']['aspectratio']['x']
            scalePatch = Patch()
            scalePatch['initScale'] = scale
            set_props('ali-store-figureScale', dict(data=scalePatch))
            update_figure_spot(figure, spotSize, borderWidth, borderColor)
            zmin = min(zvalue)
            zmax = max(zvalue)
            for trace in figure['data']:
                if float(trace['z'][0]) >= zmin and float(trace['z'][0]) <= zmax:
                    trace['visible'] = True
                else:
                    trace['visible'] = 'legendonly'
            if colorMode=='Field' and colorField:
                set_color_byfield(figure, project, colorField)
            if colorMode=='Gene' and colorGene:
                set_color_bygene(figure, project, colorGene, colorRange)
            set_props('ali-graph-right', dict(figure=figure))
            patch[0]['size'] = '50%'
            patch[1]['size'] = '50%'
        else:
            icon = 'antd-plus'
            patch[0]['size'] = '100%'
            patch[1]['size'] = '0%'
            set_props('ali-graph-right', dict(figure=hiddenGraph))
        set_props('ali-splitter-figure', dict(items=patch))
        return icon, activeKey
    return no_update, no_update

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
            set_props('ali-store-figureScale', dict(data=scalePatch))
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
    State('ali-select-project', 'value')
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

@callback(
    Input('ali-content-tabs', 'activeKey'),
)
def close_live_userspace(activeKey):
    """
    关闭实时用户区域
    """
    patch = Patch()
    if activeKey == 'ProjectInfo':
        patch['opacity'] = 0
        set_props('ali-live-userspace', dict(style=patch))