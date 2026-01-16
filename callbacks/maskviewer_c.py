from controller.maskviewer_ctl import *
from dash import callback, Input, Output, State, no_update, Patch, set_props
from dash.exceptions import PreventUpdate
from pages.components.fileSelecter import fileSelecter

@callback(
    Input('mv-button-exportData', 'nClicks'),
    State('mv-select-export-project', 'value'),
    State('mv-select-export-type', 'value')
)
def export_data(nc, project, type):
    """
    显示导出数据窗口
    """
    if nc:
        if not project:
            set_head_notice('Please select project!', 'warning')
        elif not type:
            set_head_notice('Please select export type!', 'warning')
        else:
            maskData.set_export_project(project, type)
            fileSelecter.open_export_box()

@callback(
    Input('mv-checkbox-mask', 'checked'),
    Input('mv-checkbox-contour', 'checked'),
    State('mv-store-figure-type', 'data')
)
def switch_mask_contour(checked_mask, checked_contour, figure_type):
    """
    开关mask与contour
    """
    patch = Patch()
    patch['data'][1]['visible'] = checked_mask
    patch['data'][2]['visible'] = checked_contour
    if figure_type['left']:
        set_props('mv-graph-left', dict(figure=patch))
    if figure_type['right']:
        set_props('mv-graph-right', dict(figure=patch))

@callback(
    Input("mv-graph-right", "relayoutData"),
    State("mv-graph-right", "figure"),
    State("mv-store-figure-type", "data"),
)
def update_relayout(rightLayout, rightFig, figure_type):
    """
    同步左图状态
    """
    if rightFig is None or rightLayout is None or figure_type['left']!=figure_type['right']:
        raise PreventUpdate
    patch = Patch()
    patch['layout']['xaxis'] = rightFig['layout']['xaxis']
    patch['layout']['yaxis'] = rightFig['layout']['yaxis']
    set_props('mv-graph-left', dict(figure=patch))

@callback(
    Input("mv-graph-left", "relayoutData"),
    State("mv-graph-left", "figure"),
    State("mv-store-figure-type", "data"),
)
def update_relayout(leftLayout, leftFig, figure_type):
    """
    同步右图状态
    """
    if leftFig is None or leftLayout is None or not figure_type['contrast'] or figure_type['left']!=figure_type['right']:
        return
    patch = Patch()
    patch['layout']['xaxis'] = leftFig['layout']['xaxis']
    patch['layout']['yaxis'] = leftFig['layout']['yaxis']
    set_props('mv-graph-right', dict(figure=patch))

@callback(
    Output('mv-graph-right', 'figure', allow_duplicate=True),
    Output('mv-store-figure-type', 'data', allow_duplicate=True),
    Input('mv-select-contrast', 'value'),
    State('mv-select-project', 'value'),
    State('mv-select-slice', 'value'),
    State('mv-checkbox-mask', 'checked'),
    State('mv-checkbox-contour', 'checked'),
    running=[
        (Output('mv-select-contrast', 'disabled'), True, False),
        (Output('mv-select-slice', 'disabled'), True, False),
    ],
    prevent_initial_call=True,
)
def right_figure_change(figureType, project, slice, showMask, showContour):
    """
    切换右图
    """
    if figureType and project and slice:
        patch = Patch()
        if figureType=='overlay':
            patch['right'] = False
        else:
            patch['right'] = True
        graph = maskData.get_graph(project, slice, figureType, showMask=showMask, showContour=showContour)
        return graph, patch
    return no_update, no_update

@callback(
    Output('mv-graph-left', 'figure', allow_duplicate=True),
    Output('mv-store-figure-type', 'data', allow_duplicate=True),
    Input('mv-select-figure', 'value'),
    State('mv-select-project', 'value'),
    State('mv-select-slice', 'value'),
    State('mv-checkbox-mask', 'checked'),
    State('mv-checkbox-contour', 'checked'),
    running=[
        (Output('mv-select-figure', 'disabled'), True, False),
        (Output('mv-select-slice', 'disabled'), True, False),
    ],
    prevent_initial_call=True,
)
def left_figure_change(figureType, project, slice, showMask, showContour):
    """
    切换左图
    """
    if figureType and project and slice:
        patch = Patch()
        if figureType=='overlay':
            patch['left'] = False
        else:
            patch['left'] = True
        graph = maskData.get_graph(project, slice, figureType, showMask=showMask, showContour=showContour)
        return graph, patch
    return no_update, no_update

@callback(
    Input('mv-select-slice', 'value'),
    State('mv-select-project', 'value'),
    State('mv-select-figure', 'value'),
    State('mv-select-contrast', 'value'),
    State('mv-checkbox-mask', 'checked'),
    State('mv-checkbox-contour', 'checked'),
    State('mv-store-figure-type', 'data'),
    running=[
        (Output('mv-select-figure', 'disabled'), True, False),
        (Output('mv-select-contrast', 'disabled'), True, False),
        (Output('mv-select-slice', 'disabled'), True, False),
    ]
)
def slice_change(slice, project, figType, contrastType, showMask, showContour, figure_type):
    """
    切换切片时，更新图像
    """
    if project and slice:
        figTypes = set(maskData.get_figure_types(project, slice))
        if figType in figTypes:
            graph = maskData.get_graph(project, slice, figType, showMask=showMask, showContour=showContour)
            set_props('mv-graph-left', dict(figure=graph))
        else:
            set_props('mv-select-figure', dict(value=None))
        if contrastType in figTypes:
            if figure_type['contrast']:
                graph = maskData.get_graph(project, slice, contrastType, showMask=showMask, showContour=showContour)
                set_props('mv-graph-right', dict(figure=graph))
        else:
            set_props('mv-select-contrast', dict(value=None))

@callback(
    Output('mv-select-slice', 'value'),
    Output('mv-select-figure', 'value'),
    Input('mv-select-project', 'value'),
    State('mv-store-figure-type', 'data')
)
def reset_select_value(project, figure_type):
    """
    重置选择框值
    """
    if project:
        if figure_type['contrast']:
            set_props('mv-select-contrast', dict(value=None))
        return None, None
    return no_update, no_update

@callback(
    Output('mv-contrast-figure-space', 'style', allow_duplicate=True),
    Output('mv-splitter', 'items', allow_duplicate=True),
    Output('mv-store-figure-type', 'data', allow_duplicate=True),
    Input('mv-button-remove-contrast', 'nClicks'),   
    State('mv-splitter', 'items'),
    prevent_initial_call=True
)
def remove_contrast_space(nClicks, items): 
    """
    移除对比图像
    """
    if nClicks:
        items.pop()
        patch = Patch()
        patch['right'] = False
        patch['contrast'] = False
        return {'width':'100%', 'display':'none'}, items, patch
    return {'width':'100%'}, no_update, no_update


@callback(
    Output('mv-contrast-figure-space', 'style'),
    Output('mv-splitter', 'items'),
    Output('mv-store-figure-type', 'data'),
    Output('mv-select-contrast', 'value'),
    Input('mv-button-add-contrast', 'nClicks'),  
    State('mv-splitter', 'items'),
    State('mv-select-contrast', 'value'),
    State('mv-store-figure-type', 'data')
)
def add_contrast_space(nClicks, items, contrast_type, figure_type):
    """
    添加对比图像
    """
    if nClicks and not figure_type['contrast']:
        items.append(contrast_graph)
        patch = Patch()
        patch['contrast'] = True
        return {'width':'100%'}, items, patch, contrast_type
    raise PreventUpdate

@callback(
    Output('mv-select-contrast', 'options'),
    Input('mv-select-contrast-tooltip', 'open'),
    State('mv-select-project', 'value'),
    State('mv-select-slice', 'value'),
    prevent_initial_call=True
)
def update_contrast_options(open, project, slice):
    """
    更新对比图像类型列表
    """
    if open and project and slice:
        figure_types = maskData.get_figure_types(project, slice)
        return figure_types
    return no_update

@callback(
    Output('mv-select-figure', 'options'),
    Input('mv-select-figure-tooltip', 'open'),
    State('mv-select-project', 'value'),
    State('mv-select-slice', 'value'),
    prevent_initial_call=True
)
def update_figure_options(open, project, slice):
    """
    更新图像类型列表
    """
    if open and project and slice:
        figure_types = maskData.get_figure_types(project, slice)
        return figure_types
    return no_update

@callback(
    Output('mv-select-slice', 'options'),
    Input('mv-select-slice-tooltip', 'open'),
    State('mv-select-project', 'value'),
    prevent_initial_call=True
)
def update_slice_options(open, project):
    """
    更新切片列表
    """
    if open and project:
        slices = maskData.get_project_slices(project)
        return slices
    return no_update

@callback(
    Output('mv-select-project', 'options'),
    Input('mv-select-project-tooltip', 'open'),
    prevent_initial_call=True
)
def update_project_options(open):
    """
    更新项目列表
    """
    if open:
        projects = maskData.get_exist_projects()
        return projects
    return no_update

@callback(
    Output('mv-select-export-project', 'options'),
    Input('mv-select-export-project-tooltip', 'open'),
    prevent_initial_call=True
)
def update_export_project_options(open):
    """
    更新导出项目列表
    """
    if open:
        projects = maskData.get_exist_projects()
        return projects
    return no_update

@callback(
    Output('mv-select-export-type', 'options'),
    Input('mv-select-export-type-tooltip', 'open'),
    State('mv-select-export-project', 'value'),
    prevent_initial_call=True
)
def update_export_type_options(open, project):
    """
    更新导出结果选项
    """
    if open and project:
        types = set()
        for slice in maskData.get_project_slices(project):
            curTypes = maskData.get_figure_types(project, slice)
            curTypes = set(curTypes)
            if len(types) > 0:
                types = types.intersection(curTypes)
            else:
                types = types.union(curTypes)
        types.discard('overlay')
        types = list(types)
        types.sort()
        return types
    return no_update