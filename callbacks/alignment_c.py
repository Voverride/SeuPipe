from pages.components.fileSelecter import fileSelecter
from controller.alignment_ctl import *
from dash.exceptions import PreventUpdate
from dash import Input, Output, State, no_update, callback, Patch, ctx

@callback(
    Input('alignment-interval', 'n_intervals')
)
def update_align_status(n):
    """
        轮回查询对齐状态
    """
    update_alignment_progress()

@callback(
    # Output('main-loading-area', 'children', allow_duplicate=True),
    Input('alignment-button-alignSlices', 'nClicks'),
    State('alignment-select-model', 'value'),
    State('alignment-radio-device', 'value'),
    State('main-title-username', 'children'),
    prevent_initial_call=True,
)
def align_slices(nc, model, device, usrname):
    """
        对齐切片
    """
    if nc is not None:
        use_gpu = True
        if device == 'CPU':
            use_gpu = False
        align_slices_with_paste(model, use_gpu, creator=usrname)
    # return no_update

@callback(
  Output("alignment-graph-origion", "figure", allow_duplicate=True),
  Output("alignment-graph-aligned", "figure", allow_duplicate=True),
  Input("alignment-graph-origion", "restyleData"),
  Input("alignment-graph-aligned", "restyleData"),
  prevent_initial_call=True
)
def update_legend(rd1, rd2):
    fig1 = alidata.get_orifig()
    fig2 = alidata.get_alifig()
    if not fig1 or not fig2:
        raise PreventUpdate
    tid = ctx.triggered_id
    if tid == 'alignment-graph-origion':
        patch = Patch()
        legendState = rd1[0]['visible']
        legendIndex = rd1[1]
        for idx, status in zip(legendIndex, legendState):
            patch['data'][idx]['visible'] = status
        return no_update, patch
    else:
        patch = Patch()
        legendState = rd2[0]['visible']
        legendIndex = rd2[1]
        for idx, status in zip(legendIndex, legendState):
            patch['data'][idx]['visible'] = status
        return patch, no_update

@callback(
    Output("alignment-graph-origion", "figure", allow_duplicate=True),
    Output("alignment-graph-aligned", "figure", allow_duplicate=True),
    Input("alignment-graph-origion", "relayoutData"),
    Input("alignment-graph-aligned", "relayoutData"),
    State('main-title-username', 'children'),
    prevent_initial_call=True,
)
def update_relayout(oriLayout, aliLayout, usrname):
    """
        同步两图状态
    """
    fig1 = alidata.get_orifig()
    fig2 = alidata.get_alifig()
    if not fig1 or not fig2:
        raise PreventUpdate
    
    oriScale = fig1['layout']['scene']['aspectratio']['x']
    aliScale = fig2['layout']['scene']['aspectratio']['x']

    tid = ctx.triggered_id
    patchOri = Patch()
    patchAli = Patch()
    if tid == 'alignment-graph-origion':
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
            return patchOri, patchAli
    if tid == 'alignment-graph-aligned':
        mark = False
        if 'scene.camera' in aliLayout:
            patchOri['layout']['scene']['camera']['eye'] = aliLayout['scene.camera']['eye']
            patchAli['layout']['scene']['camera']['eye'] = aliLayout['scene.camera']['eye']
            z_dire = aliLayout['scene.camera']['eye']['z']
            alidata.set_actref_slice(usrname, reverse=False if z_dire>=0 else True)
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
            return patchOri, patchAli
    raise PreventUpdate

@callback(
    Output('main-notice-area', 'children', allow_duplicate=True),
    Output('alignment-graph-aligned', 'figure', allow_duplicate=True),
    Input('graphSetting-button-apply', 'nClicks'),
    State('graphSetting-table-applySlices', 'selectedRowKeys'),
    State('main-title-username', 'children'),
    prevent_initial_call=True
)
def apply_to_other_slices(nClicks, selectedRowKeys, usrname):
    """
        将当前变换应用到其他切片
    """
    if not nClicks:
        raise PreventUpdate
    permission = verify_modify_permission()
    if not permission:
        return no_update, no_update
    notice = no_update
    patch = no_update

    mtx = alidata.get_coord_trans_mtx(usrname)
    if not alidata.has_alifig():
        notice = fac.AntdMessage(content=f'There was no aligned graph exists !', type='warning')
    elif mtx==None or mtx==[0, 0, 0]:
        notice = fac.AntdMessage(content=f'No modifications have been made to the current active slice !', type='warning')
    elif selectedRowKeys==None or len(selectedRowKeys)==0:
        notice = fac.AntdMessage(content=f'You have not select any slice to apply !', type='warning')
    else:
        patch = Patch()
        slices = alidata.get_slice_list()
        select_slices = [slices[int(i)] for i in selectedRowKeys]
        coordTransed = get_transformed_coord(usrname, select_slices, storeNewCoord=True)

        alidata.set_coord_trans_mtx(usrname)

        for key in selectedRowKeys:
            idx = int(key)
            slice = slices[idx]
            coord = coordTransed[slice]
            patch['data'][idx]['x'] = coord[:, 0]
            patch['data'][idx]['y'] = coord[:, 1]
        notice = fac.AntdMessage(content=f'Operation completed successfully !', type='success')
    return notice, patch

@callback(
    Output('main-notice-area', 'children', allow_duplicate=True),
    Output('graphSetting-table-applySlices', 'selectedRowKeys'),
    Input('graphSetting-button-pickAbove', 'nClicks'),
    Input('graphSetting-button-pickBelow', 'nClicks'),
    State('graphSetting-selecter-activeSlice', 'value'),
    State('graphSetting-table-applySlices', 'selectedRowKeys'),
    prevent_initial_call=True
)
def update_selected_applySlices(pickAbove, pickBelow, activeSlice, selectedRowKeys):
    """
        选中activeSlice上方或下方切片
    """
    if not pickAbove and not pickBelow:
        raise PreventUpdate
    if not activeSlice:
        return fac.AntdMessage(content=f'Please ensure that the active slice was selected !', type='warning'), no_update
    tid = ctx.triggered_id
    slices_length = len(alidata.get_slice_list())
    idx = alidata.get_slice_index(float(activeSlice))
    if selectedRowKeys==None:
        selectedRowKeys = []
    if tid=='graphSetting-button-pickAbove':
        for i in range(idx+1, slices_length):
            rowkey = str(i)
            if rowkey not in selectedRowKeys:
                selectedRowKeys.append(rowkey)
    if tid=='graphSetting-button-pickBelow':
        for i in range(0, idx):
            rowkey = str(i)
            if rowkey not in selectedRowKeys:
                selectedRowKeys.append(rowkey)
    return no_update, selectedRowKeys

@callback(
    Output('main-notice-area', 'children', allow_duplicate=True),
    Output('alignment-graph-aligned', 'figure', allow_duplicate=True),
    Output('graphSetting-dialog-transMtxAlert', 'visible', allow_duplicate=True),
    Input('graphSetting-modify-discard', 'nClicks'),
    Input('graphSetting-modify-apply', 'nClicks'),
    State('graphSetting-table-applySlices', 'selectedRowKeys'),
    State('graphSetting-selecter-activeSlice', 'value'),
    State('graphSetting-selecter-referenceSlice', 'value'),
    State('main-title-username', 'children'),
    prevent_initial_call=True
)
def apply_discard_modify(discard, apply, selectedRowKeys, actSlice, refSlice, usrname):
    tid = ctx.triggered_id
    notice = no_update
    if discard or apply:
        ad = alidata.get_adata()
        patch = Patch()
        if tid=='graphSetting-modify-apply':
            notice = fac.AntdMessage(content=f'The changes have been successfully applied !', type='success')
            slices = alidata.get_slice_list()
            select_slices = [slices[int(i)] for i in selectedRowKeys]
            coordTransed = get_transformed_coord(usrname, select_slices, storeNewCoord=True)
            for key in selectedRowKeys:
                idx = int(key)
                slice = slices[idx]
                coord = coordTransed[slice]
                patch['data'][idx]['x'] = coord[:, 0]
                patch['data'][idx]['y'] = coord[:, 1]

        alidata.set_coord_trans_mtx(usrname)

        if tid=='graphSetting-modify-discard':
            notice = fac.AntdMessage(content=f'The data has been restored !', type='success')
            idx = alidata.get_slice_index(actSlice)
            obsIndex = alidata.get_alifig()['data'][idx]['customdata'].flatten()
            df = ad.obs.loc[obsIndex, ['x_aligned', 'y_aligned']]
            patch['data'][idx]['x'] = df['x_aligned']
            patch['data'][idx]['y'] = df['y_aligned']
            
        return notice, patch, False
        
    raise PreventUpdate

@callback(
    Output('main-notice-area', 'children', allow_duplicate=True),
    Output("alignment-graph-aligned", "figure", allow_duplicate=True),
    Output('graphSetting-button-showGrid', 'style'),
    Input('graphSetting-button-showGrid', 'nClicks'),
    State('graphSetting-button-showGrid', 'style'),
    prevent_initial_call=True
)
def switch_grid(nClicks, style):
    """
        开关图像网格
    """
    if not nClicks:
        raise PreventUpdate
    notice = no_update
    patch = no_update
    status = ['#949495', '#ca8269'] # 0 -> 关， 1-> 开
    if not alidata.has_alifig():
        notice = fac.AntdMessage(content=f'There was no aligned graph exists !', type='warning')
        style = no_update
    else:
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
    return notice, patch, style

@callback(
    Output('main-notice-area', 'children', allow_duplicate=True),
    Output("alignment-graph-aligned", "figure", allow_duplicate=True),
    Input('graphSetting-button-reFocus', 'nClicks'),
    State('graphSetting-selecter-activeSlice', 'value'),
    State('graphSetting-selecter-referenceSlice', 'value'),
    prevent_initial_call=True
)
def reset_view(nClicks, actSlice, refSlice):
    """
        重置图像视角
    """
    if not nClicks:
        raise PreventUpdate
    notice = no_update
    patch = no_update
    if not alidata.has_alifig():
        notice = fac.AntdMessage(content=f'There was no aligned graph exists !', type='warning')
    else:
        patch = Patch()
        patch['layout']['scene']['camera']['eye'] = alidata.get_graph_view(actSlice, refSlice)
        # patch['layout']['scene']['aspectratio'] = alidata.get_alifig()['layout']['scene']['aspectratio']
    return notice, patch

@callback(
    Output('alignment-graph-aligned', 'figure', allow_duplicate=True),
    Input('key-pressed-events', 'data'),
    State('graphSetting-selecter-activeSlice', 'value'),
    State('graphSetting-selecter-referenceSlice', 'value'),
    State('graphSetting-inputNum-stepSize', 'value'),
    State('graphSetting-inputNum-rotationAngle', 'value'),
    State('main-title-username', 'children'),
    prevent_initial_call=True
)
def detect_pressed_key(keyboard, actSlice, refSlice, stepSize, rotationAngle, usrname):
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
    keyCode = int(keyboard) if keyboard.isdigit() else 0
    detectedKeys = {37, 38, 39, 40, 1037, 1038, 1039, 1040}
    if not actSlice or not refSlice or not stepSize or not rotationAngle or keyCode not in detectedKeys:
        raise PreventUpdate
    permission = verify_modify_permission()
    if not permission:
        return no_update
    patch = Patch()
    if keyCode==38:
        alidata.set_coord_trans_mtx(usrname=usrname, dy=stepSize)
    elif keyCode==40:
        alidata.set_coord_trans_mtx(usrname=usrname, dy=-stepSize)
    elif keyCode==37:
        alidata.set_coord_trans_mtx(usrname=usrname, dx=-stepSize)
    elif keyCode==39:
        alidata.set_coord_trans_mtx(usrname=usrname, dx=stepSize)
    elif keyCode==1037 or keyCode==1038:
        alidata.set_coord_trans_mtx(usrname=usrname, reg=-rotationAngle)
    elif keyCode==1039 or keyCode==1040:
        alidata.set_coord_trans_mtx(usrname=usrname, reg=rotationAngle)
    actSlice = float(actSlice)
    coordTransed = get_transformed_coord(usrname, {actSlice})
    slice_index = alidata.get_slice_index(actSlice)
    coord = coordTransed[actSlice]
    patch['data'][slice_index]['x'] = coord[:,0]
    patch['data'][slice_index]['y'] = coord[:,1]
    return patch

@callback(
    Output('main-notice-area', 'children', allow_duplicate=True),
    Output('alignment-graph-aligned', 'figure', allow_duplicate=True),
    Input('graphSetting-button-up', 'nClicks'),
    Input('graphSetting-button-down', 'nClicks'),
    Input('graphSetting-button-left', 'nClicks'),
    Input('graphSetting-button-right', 'nClicks'),
    Input('graphSetting-button-clockwise', 'nClicks'),
    Input('graphSetting-button-unclockwise', 'nClicks'), 
    State('graphSetting-selecter-activeSlice', 'value'),
    State('graphSetting-selecter-referenceSlice', 'value'),
    State('graphSetting-inputNum-stepSize', 'value'),
    State('graphSetting-inputNum-rotationAngle', 'value'),
    State('main-title-username', 'children'),
    prevent_initial_call=True
)
def move_slice(up, down, left, right, clockwise, unclockwise, actSlice, refSlice, stepSize, rotationAngle, usrname):
    """
        根据方向按钮调整切片位置
    """
    buttonPressed = up or down or left or right or clockwise or unclockwise
    if not buttonPressed:
        raise PreventUpdate
    permission = verify_modify_permission()
    if not permission:
        return no_update, no_update
    notice = None
    if not alidata.has_alifig():
        notice = fac.AntdMessage(content=f'There was no aligned graph to modify !', type='warning')
    elif not actSlice:
        notice = fac.AntdMessage(content=f'The active slice cannot be empty !', type='warning')
    elif not refSlice:
        notice = fac.AntdMessage(content=f'The reference slice cannot be empty !', type='warning')
    elif not stepSize:
        notice = fac.AntdMessage(content=f'Step size cannot be empty !', type='warning')
    elif not rotationAngle:
        notice = fac.AntdMessage(content=f'Rotation angle cannot be empty !', type='warning')
    
    if notice!=None:
        return notice, no_update
    
    notice = no_update

    patch  = Patch()

    tid = ctx.triggered_id
    if tid=='graphSetting-button-up':
        alidata.set_coord_trans_mtx(usrname, dy=stepSize)
    elif tid=='graphSetting-button-down':
        alidata.set_coord_trans_mtx(usrname, dy=-stepSize)
    elif tid=='graphSetting-button-left':
        alidata.set_coord_trans_mtx(usrname, dx=-stepSize)
    elif tid=='graphSetting-button-right':
        alidata.set_coord_trans_mtx(usrname, dx=stepSize)
    elif tid=='graphSetting-button-clockwise':
        alidata.set_coord_trans_mtx(usrname, reg=rotationAngle)
    elif tid=='graphSetting-button-unclockwise':
        alidata.set_coord_trans_mtx(usrname, reg=-rotationAngle)
    actSlice = float(actSlice)
    coordTransed = get_transformed_coord(usrname, {actSlice})
    slice_index = alidata.get_slice_index(actSlice)
    coord = coordTransed[actSlice]
    patch['data'][slice_index]['x'] = coord[:,0]
    patch['data'][slice_index]['y'] = coord[:,1]
    
    return notice, patch

@callback(
    Output('alignment-graph-aligned', 'figure', allow_duplicate=True),
    Output('graphSetting-dialog-transMtxAlert', 'visible'),
    Output('GraphSetting-drawer-setting', 'visible', allow_duplicate=True),
    Output('graphSetting-selecter-activeSlice', 'value', allow_duplicate=True),
    Output('graphSetting-table-applySlices', 'selectedRowKeys', allow_duplicate=True),
    Input('GraphSetting-drawer-setting', 'visible'),
    Input('graphSetting-selecter-activeSlice', 'value'),
    Input('graphSetting-selecter-referenceSlice', 'value'),
    State('main-title-username', 'children'),
    prevent_initial_call = True
)
def check_transMtx_and_direction(visible, actSlice, refSlice, usrname):
    """
        关闭抽屉或切换active slice时检查是否有未保存操作，
        基于活动和参考切片设置当前方向状态
    """
    tid = ctx.triggered_id
    dialogVisible = no_update
    graphSettingVisible = no_update
    selectedRowKeys = no_update
    actOri = actSlice
    patch = no_update
    transMtx = alidata.get_coord_trans_mtx(usrname)
    hasTransMtx = True if transMtx and transMtx!=[0, 0, 0] else False
    if tid=='GraphSetting-drawer-setting':
        if not visible:
            if hasTransMtx:
                dialogVisible = True
                graphSettingVisible = True
            else:
                if alidata.has_alifig():
                    patch = Patch()
                    reset_slicecolor(patch, alidata.get_slice_list())
    else:
        if tid=='graphSetting-selecter-activeSlice':
            if hasTransMtx:
                dialogVisible = True
                actSlice = alidata.get_actref_slice(usrname)[0]
            elif actSlice:
                idx = alidata.get_slice_index(float(actSlice))
                selectedRowKeys = [str(idx)] 
                act = float(actSlice)
                alidata.set_actref_slice(usrname=usrname, act_slice=act)
                if refSlice:      
                    ref = float(refSlice)
                    alidata.set_actref_slice(usrname=usrname, ref_slice=ref, reverse=False if act>=ref else True)    
        else:
            if actSlice and refSlice:
                act = float(actSlice)
                ref = float(refSlice)
                alidata.set_actref_slice(usrname, act, ref, False if act>=ref else True)

    if actOri==actSlice:
        actSlice = no_update 

    return patch, dialogVisible, graphSettingVisible, actSlice, selectedRowKeys

@callback(
    Output('graphSetting-div-legend', 'children', allow_duplicate=True),
    Output('graphSetting-div-loadMore', 'style', allow_duplicate=True),
    Input('graphSetting-button-loadMore', 'nClicks'),
    State('graphSetting-div-loadMore', 'style'),
    State('graphSetting-selecter-colorMode', 'value'),
    State('graphSetting-selecter-colorField', 'value'),
    State('graphSetting-selecter-activeSlice', 'value'),
    State('graphSetting-selecter-referenceSlice', 'value'),
    prevent_initial_call=True
)
def load_more_legend(nClicks, loadMore, colorMode, colorField, actSlice, refSlice):
    """
        加载更多图例
    """
    if nClicks and colorMode=='Field' and colorField:
        sliceList = [actSlice, refSlice]
        if not actSlice or not refSlice:
            sliceList = alidata.get_slice_list()
        legends, hasLegend = get_field_legends(colorField, sliceList, loadMore['legendNum'])
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
    Output('graphSetting-div-legend', 'children', allow_duplicate=True),
    Output('graphSetting-div-loadMore', 'style', allow_duplicate=True),
    Input('graphSetting-selecter-colorMode', 'value'),
    Input('graphSetting-selecter-colorField', 'value'),
    Input('graphSetting-selecter-activeSlice', 'value'),
    Input('graphSetting-selecter-referenceSlice', 'value'),
    prevent_initial_call=True
)
def update_expand_fieldLegend_panel(colorMode, colorField, actSlice, refSlice):
    """
        刷新field域对应图例内容
    """
    if colorMode=='Field' and colorField:
        sliceList = [actSlice, refSlice]
        if not actSlice or not refSlice:
            sliceList = alidata.get_slice_list()
        legends, hasLegend = get_field_legends(colorField, sliceList, 0)
        patch = Patch()
        displayLoadMore = 'none'
        if hasLegend:
            displayLoadMore = 'block'
        patch['display'] = displayLoadMore
        patch['legendNum'] = len(legends)
        return legends, patch
    raise PreventUpdate

@callback(
    Output('graphSetting-card-legend', 'visible', allow_duplicate=True),
    Output('graphSetting-card-legend', 'title', allow_duplicate=True),
    Input('graphSetting-selecter-colorMode', 'value'),
    Input('graphSetting-selecter-colorField', 'value'),
    Input('GraphSetting-drawer-setting', 'visible'),
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
    Output("alignment-graph-aligned", "figure", allow_duplicate=True),
    Input('graphSetting-selecter-colorMode', 'value'),
    Input('alignment-button-graphSetting', 'nClicks'),
    State('graphSetting-selecter-colorField', 'value'),
    State('graphSetting-selecter-colorGene', 'value'),
    State('graphSetting-selecter-colorRange', 'value'),
    State('graphSetting-ColorPicker-actSlice', 'value'),
    State('graphSetting-ColorPicker-refSlice', 'value'),
    State('graphSetting-selecter-activeSlice', 'value'),
    State('graphSetting-selecter-referenceSlice', 'value'),
    prevent_initial_call=True
)
def modify_sliceColor_by_colorMode(value, nClicks, field, gene, colorType, actColor, refColor, actSlice, refSlice):
    """
        基于colorMode修改切片颜色
    """
    if not value or not nClicks:
        raise PreventUpdate
    
    if not alidata.has_alifig():
        raise PreventUpdate
    
    patch = Patch()
    sliceList = [actSlice, refSlice]
    if not actSlice or not refSlice:
        sliceList = alidata.get_slice_list()

    if value=='Field':
        if field is not None:
            set_color_byfield(patch, field, sliceList)
        return patch
    
    elif value=='Gene':
        if gene is not None:
            set_color_bygene(patch, gene, sliceList, colorType)
        return patch
    
    if actSlice and refSlice:
        patch['data'][alidata.get_slice_index(actSlice)]['marker']['color'] = actColor
        patch['data'][alidata.get_slice_index(refSlice)]['marker']['color'] = refColor

    return patch

@callback(
    Output("alignment-graph-aligned", "figure", allow_duplicate=True),
    Input('graphSetting-selecter-colorRange', 'value'),
    State('graphSetting-selecter-colorGene', 'value'),
    State('graphSetting-selecter-activeSlice', 'value'),
    State('graphSetting-selecter-referenceSlice', 'value'),
    prevent_initial_call=True
)
def update_sliceColor_byColorType(colorType, gene, actSlice, refSlice):
    """
        基于colortype显示基因表达
    """
    if not gene or not alidata.has_alifig():
        raise PreventUpdate
    sliceList = [actSlice, refSlice]
    if not actSlice or not refSlice:
        sliceList = alidata.get_slice_list()
    patch = Patch()
    set_color_bygene(patch, gene, sliceList, colorType)
    return patch

@callback(
    Output("alignment-graph-aligned", "figure", allow_duplicate=True),
    Input('graphSetting-selecter-colorGene', 'value'),
    State('graphSetting-selecter-activeSlice', 'value'),
    State('graphSetting-selecter-referenceSlice', 'value'),
    State('graphSetting-selecter-colorRange', 'value'),
    prevent_initial_call=True
)
def update_sliceColor_byGene(gene, actSlice, refSlice, colorType):
    """
        基于基因表达显示切片颜色
    """
    if not gene:
        raise PreventUpdate
    patch = Patch()
    sliceList = [actSlice, refSlice]
    if not actSlice or not refSlice:
        sliceList = alidata.get_slice_list()
    set_color_bygene(patch, gene, sliceList, colorType)
    return patch

@callback(
    Output("alignment-graph-aligned", "figure", allow_duplicate=True),
    Input('graphSetting-selecter-colorField', 'value'),
    State('graphSetting-selecter-activeSlice', 'value'),
    State('graphSetting-selecter-referenceSlice', 'value'),
    prevent_initial_call=True
)
def update_sliceColor_byField(field, actSlice, refSlice):
    """
        基于数据域显示切片颜色
    """
    if not field:
        raise PreventUpdate
    sliceList = [actSlice, refSlice]
    if not actSlice or not refSlice:
        sliceList = alidata.get_slice_list()
    patch = Patch()
    set_color_byfield(patch, field, sliceList)
    return patch

@callback(
    Output("alignment-graph-aligned", "figure", allow_duplicate=True),
    Input('graphSetting-ColorPicker-actSlice', 'value'),
    Input('graphSetting-ColorPicker-refSlice', 'value'),
    State('graphSetting-selecter-activeSlice', 'value'),
    State('graphSetting-selecter-referenceSlice', 'value'),
    prevent_initial_call = True
)
def show_adjusting_slice(actColor, refColor, actSlice, refSlice):
    """
        更改活动和参考切片颜色
    """
    tid = ctx.triggered_id
    if not alidata.has_alifig() or not actSlice or not refSlice:
        raise PreventUpdate
    patch = Patch()
    if tid=='graphSetting-ColorPicker-actSlice':
        idx = alidata.get_slice_index(actSlice)
        patch['data'][idx]['marker']['color'] = actColor
    else:
        idx = alidata.get_slice_index(refSlice)
        patch['data'][idx]['marker']['color'] = refColor

    return patch

@callback(
    Output('graphSetting-ColorPicker-actSlice', 'disabled'),
    Output('graphSetting-ColorPicker-refSlice', 'disabled'),
    Output('graphSetting-selecter-colorField', 'disabled'),
    Output('graphSetting-selecter-colorGene', 'disabled'),
    Output('graphSetting-selecter-colorRange', 'disabled'),
    Input('graphSetting-selecter-colorMode', 'value'),
    prevent_initial_call=True
)
def disabled_color_mode(value):
    """
        基于colorMode禁用颜色组件
    """
    if not value:
        raise PreventUpdate
    
    if value=='Field':
        return True, True, False, True, True
    elif value=='Gene':
        return True, True, True, False, False
    return False, False, True, True, True

@callback(
    Output("alignment-graph-aligned", "figure", allow_duplicate=True),
    Output('graphSetting-inputNum-stepSize', 'value'),
    Input('graphSetting-selecter-activeSlice', 'value'),
    Input('graphSetting-selecter-referenceSlice', 'value'),
    State('main-title-username', 'children'),
    State('graphSetting-selecter-colorMode', 'value'),
    State('graphSetting-selecter-colorField', 'value'),
    State('graphSetting-ColorPicker-actSlice', 'value'),
    State('graphSetting-ColorPicker-refSlice', 'value'),
    State('graphSetting-selecter-colorGene', 'value'),
    State('graphSetting-selecter-colorRange', 'value'),
    prevent_initial_call = True
)
def show_adjusting_slice(activeSlice, referenceSlice, usrname, colorMode, field, actColor, refColor, gene, colorType):
    """
        只显示需要调整切片，并转换到直角坐标系视角，同时计算出默认的移动步长
    """
    transMtx = alidata.get_coord_trans_mtx(usrname)
    hasTransMtx = True if transMtx and transMtx!=[0, 0, 0] else False
    if ctx.triggered_id=='graphSetting-selecter-activeSlice' and hasTransMtx:
        raise PreventUpdate
    if activeSlice and referenceSlice and alidata.has_alifig():
        patch = Patch()
        patch['layout']['scene']['camera']['eye'] = alidata.get_graph_view(activeSlice, referenceSlice)
        slices = alidata.get_slice_list()
        for i, slice in enumerate(slices):
            if slice==activeSlice or slice==referenceSlice:
                patch['data'][i]['visible'] = True
            else:
                patch['data'][i]['visible'] = 'legendonly'

        sliceList = [activeSlice, referenceSlice]
        if colorMode=='Custom':
            for slice in sliceList:
                i = alidata.get_slice_index(slice)
                patch['data'][i]['marker']['color'] = actColor if slice==activeSlice else refColor
        elif colorMode=='Field':
            if field is not None:
                set_color_byfield(patch, field, sliceList)
        else:
            if gene is not None:
                set_color_bygene(patch, gene, sliceList, colorType)
        step = alidata.get_movestep_size(activeSlice, referenceSlice)
        return patch, step
    raise PreventUpdate

@callback(
    Output('graphSetting-selecter-referenceSlice', 'value'),
    Input('graphSetting-selecter-activeSlice', 'value'),
    prevent_initial_call=True
)
def update_refSlice_by_actSlice(act_slice):
    """
        调整actSlice时默认绑定refSlice为下面的切片
    """
    if not act_slice:
        raise PreventUpdate
    below_slice = get_below_slice(float(act_slice))
    return below_slice

@callback(
    Output('alignment-graph-origion', 'figure', allow_duplicate=True),
    Output('alignment-graph-aligned', 'figure', allow_duplicate=True),
    Input('alignmentGraphSetting-button-slicer', 'nClicks'),
    State('alignmentGraphSetting-slider-z', 'value'),
    prevent_initial_call=True
)
def slicer_data(nClicks, values):
    """
        根据滑动条筛选切片
    """
    def slicer(values):
        values.sort()
        min, max = values
        fig = Patch()
        slices = alidata.get_slice_list()
        if slices is None:
            return no_update
        size = len(slices)
        for i in range(size):
            z = slices[i]
            if z>=min and z<=max:
                fig['data'][i]['visible'] = True
            else:
                fig['data'][i]['visible'] = 'legendonly'
        return fig

    figOri = no_update
    figAli = no_update

    if alidata.has_orifig() or alidata.has_orifig():
        fig = slicer(values)
        if alidata.has_alifig():
            figAli = fig
        if alidata.has_orifig():
            figOri = fig

    return figOri, figAli

@callback(
    Output('alignmentGraphSetting-slider-z', 'marks'),
    Input('alignmentGraphSetting-slider-z', 'value'),
)
def update_slider_markers(value):
    """
        更新slider显示数值
    """
    if value:
        return {val:val for val in value}
    raise PreventUpdate

@callback(
    Output("alignment-graph-origion", "figure", allow_duplicate=True),
    Output("alignment-graph-aligned", "figure", allow_duplicate=True),
    Input('graphSetting-selecter-borderWidth', 'value'),
    Input('graphSetting-ColorPicker-boarderColor', 'value'),
    prevent_initial_call=True
)
def update_boarder_width(width, color):
    """
        调整散点边框大小和颜色
    """
    patch1 = no_update
    patch2 = no_update
    orifig = alidata.get_orifig()
    alifig = alidata.get_alifig()
    if orifig:
        patch1 = Patch()
        for i in range(len(orifig['data'])):
            patch1['data'][i]['marker']['line']['color'] = color
            patch1['data'][i]['marker']['line']['width'] = width
    if alifig:
        patch2 = Patch()
        for i in range(len(alifig['data'])):
            patch2['data'][i]['marker']['line']['color'] = color
            patch2['data'][i]['marker']['line']['width'] = width
    return patch1, patch2

@callback(
    Output("alignment-graph-origion", "figure", allow_duplicate=True),
    Output("alignment-graph-aligned", "figure", allow_duplicate=True),
    Input('graphSetting-selecter-spotSize', 'value'),
    prevent_initial_call=True
)
def update_spot_size(size):
    """
        调整图像散点大小
    """
    patch1 = no_update
    patch2 = no_update
    orifig = alidata.get_orifig()
    alifig = alidata.get_alifig()
    if orifig:
        patch1 = Patch()
        for i in range(len(orifig['data'])):
            patch1['data'][i]['marker']['size'] = size
    if alifig:
        patch2 = Patch()
        for i in range(len(alifig['data'])):
            patch2['data'][i]['marker']['size'] = size
    return patch1, patch2

@callback(
    Output('GraphSetting-drawer-setting', 'visible'),
    Input('alignment-button-graphSetting', 'nClicks')
)
def open_graph_setting_box(nc):
    if nc:
        return True
    raise PreventUpdate

@callback(
    Output('main-loading-area', 'children', allow_duplicate=True),
    Input('alignment-button-plotFig', 'nClicks'),
    State('alignment-select-x', 'value'),
    State('alignment-select-y', 'value'),
    State('alignment-select-z', 'value'),
    prevent_initial_call=True,
)
def alignment_plotFig(nClicks, x, y, z):
    """
        绘制原始3d图像
    """
    if nClicks:
        plot_origin_fig(x, y, z)
    return no_update

@callback(
    # Output('main-loading-area', 'children', allow_duplicate=True),
    Input('main-title-header', 'children'),
    # prevent_initial_call=True
)
def update_init_component(head):
    """
        恢复初始数据
    """
    if head=='Alignment':
        restore_initial_data()
    # return no_update


@callback(
    Input('alignment-button-importData', 'nClicks'),
)
def alignment_importData(nc):
    """
        显示导入数据窗口
    """
    if nc is not None:
        fileSelecter.open_import_box()

@callback(
    Input('alignment-button-exportData', 'nClicks'),
)
def alignment_exportData(nc):
    """
        显示导出数据窗口
    """
    if nc is not None:
        status = check_alignment_data()
        if status:
            fileSelecter.open_export_box()