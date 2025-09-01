from controller.expansion_ctl import *
from dash import callback, Input, Output, State, no_update, Patch
from dash.exceptions import PreventUpdate
from pages.components.fileSelecter import fileSelecter
from controller.auth import *
from dataManager.segmentation_d import segData



@callback(
    Output('expandtask-select-taskname', 'options'),
    Input('expandtask-select-taskname-tooltip', 'open'),
    prevent_initial_call=True
)
def update_tasklist_options(open):
    """
    获取最新分割任务列表
    """
    if open:
        tasks = list(segData.get_exist_tasks())
        tasks.sort()
        return tasks
    return no_update

@callback(
    Input('expand-button-submitTask', 'nClicks'),
    State('expand-input-patchsize', 'value')
)
def submit_task(nClicks, patchsize):
    """
    提交胞域扩增任务
    """
    print(patchsize)

@callback(
    Output('expand-input-patchsize', 'disabled'),
    Output('expand-input-patchsize', 'value'),
    Input('expand-select-mode', 'value')
)
def update_patchsize_disabled(mode):
    """
    根据分割模式选择是否禁用patchSize输入框
    """
    if mode.startswith('p'):
        return False, 300
    else:
        return True, 0