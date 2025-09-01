from dash import set_props
from dataManager.workspace import *
from dataManager.maskviewer_d import maskData

graph_visible_style = {'width': '100%', 'height':'90%', 'margin':'auto'}
def update_graph_with_type(slice, graph, showMask, showContour):
    """
    graphtype改变，更改图像，配准，分割或扩展
    """
    if graph=='registration':
        set_props('mv-checkbox-mask', dict(disabled=True))
        set_props('mv-checkbox-contour', dict(disabled=True))
    else:
        set_props('mv-checkbox-mask', dict(disabled=False))
        set_props('mv-checkbox-contour', dict(disabled=False))
    if slice:
        if graph=='registration':
            before, aligned = maskData.get_registration_figure()
            set_props('mv-graph-left', dict(figure=before, style=graph_visible_style))
            set_props('mv-graph-right', dict(figure=aligned, style=graph_visible_style))
        else:
            cellpose_fig, watershed_fig = maskData.get_contrast_mask_contour_figure(showMask=showMask, showContour=showContour)
            set_props('mv-graph-left', dict(figure=watershed_fig, style=graph_visible_style))
            set_props('mv-graph-right', dict(figure=cellpose_fig, style=graph_visible_style))

def update_graph_with_slice(taskname, slice, graph, showMask, showContour):
    """
    slice改变， 更改图像，配准，分割或扩展
    """
    if taskname and slice:
        maskData.read_slice_adata(taskname, slice)
        if graph=='registration':
            before, aligned = maskData.get_registration_figure()
            set_props('mv-graph-left', dict(figure=before, style=graph_visible_style))
            set_props('mv-graph-right', dict(figure=aligned, style=graph_visible_style))
        else:
            cellpose_fig, watershed_fig = maskData.get_contrast_mask_contour_figure(showMask=showMask, showContour=showContour)
            set_props('mv-graph-left', dict(figure=watershed_fig, style=graph_visible_style))
            set_props('mv-graph-right', dict(figure=cellpose_fig, style=graph_visible_style))
    
def update_taskname(taskname):
    """
    task改变，更新store， 切片列表
    """
    set_props('maskviewer-store-taskname', dict(data=taskname))
    set_props('maskviewer-select-slice', dict(value=None))
def restore_initial_data(lastTaskName):
    """
    恢复网页初始数据
    """
    if lastTaskName:
        set_props('maskviewer-select-taskname', dict(value=lastTaskName))