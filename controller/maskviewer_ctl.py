from dash import set_props
from dataManager.workspace import *
from dataManager.maskviewer_d import maskData


def update_taskname(taskname):
    """
    task改变，更新store， 切片列表
    """
    set_props('maskviewer-store-taskname', dict(data=taskname))
    set_props('maskviewer-select-slice', dict(value=None))
def update_graph(graph):
    """
    更改图像，配准，分割或扩展
    """
    if graph=='registration':
        set_props('mv-checkbox-mask', dict(disabled=True))
        set_props('mv-checkbox-contour', dict(disabled=True))
    else:
        set_props('mv-checkbox-mask', dict(disabled=False))
        set_props('mv-checkbox-contour', dict(disabled=False))
def restore_initial_data(lastTaskName):
    """
    恢复网页初始数据
    """
    if lastTaskName:
        set_props('maskviewer-select-taskname', dict(value=lastTaskName))