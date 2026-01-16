from controller.segmentation_ctl import *
from dash import callback, Input, Output, State, no_update
from pages.components.fileSelecter import fileSelecter
from controller.auth import verify_modify_permission
from websocket.websocket import ms

@callback(
    Input('seg-button-showBug', 'nClicks'),
    State('seg-select-project', 'value'),
    prevent_initial_call=True
)
def show_bug(nc, project):
    """
    显示运行时异常
    """
    if nc:
        raise_runtime_bug(project)

@callback(
    Output('seg-select-project', 'value', allow_duplicate=True),
    Input('seg-refresh-status', 'children'),
    State('seg-select-project', 'value'),
    prevent_initial_call=True
)
def refresh_segmentation_status(_, projectname):
    """
    刷新分割项目状态
    """
    return projectname

@callback(
    Input('seg-start-project', 'nClicks'),
    State('seg-select-project', 'value'),
    prevent_initial_call=True
)
def start_segment_project(nc, projectname):
    """
    启动分割项目
    """
    if nc and verify_modify_permission():
        if projectname:
            start_segmentation_project(projectname)
        else:
            set_head_notice('Please select a project to start !', type='warning')

@callback(
    Input('seg-init-restore', 'n_intervals'),
    State('seg-store-project', 'data')
)
def restore_previous_project(_, projectName):
    """
    恢复初始数据
    """
    restore_init_data(projectName)

@callback(
    Input('seg-select-project', 'value'),
    State('seg-store-project', 'data')
)
def update_table_project_data(projectName, previousProject):
    """
    根据下拉列表选项更新项目数据
    """
    if projectName:
        ms.clinetSend(ms._updateParams_c2s, params=dict(project=projectName))
        update_table_with_project(projectName)
        if projectName != previousProject:
            set_props('seg-store-project', {'data': projectName})

@callback(
    Output('seg-select-project', 'options'),
    Input('seg-select-project-tooltip', 'open'),
    prevent_initial_call=True
)
def update_projectlist_options(open):
    """
    更新项目下拉列表选项
    """
    if open:
        projects = segData.get_exist_projects()
        return projects
    return no_update

@callback(
    Input('seg-delete-project-confirm', 'confirmCounts'),
    State('seg-select-project', 'value'),
    running=[
        (Output('seg-delete-project', 'loading'), True, False),
        (Output('seg-start-project', 'disabled'), True, False),
        (Output('seg-select-project', 'disabled'), True, False)
    ],
    prevent_initial_call=True
)
def delete_project(nc, projectName):
    """
    删除项目
    """
    if nc and projectName and verify_modify_permission():
        project_running = segData.has_running_project(projectName)
        if project_running:
            set_head_notice(f'{projectName} is running, please wait for it to complete !', type='warning')
            return
        delete_project_from_disk(projectName)

@callback(
    Input('seg-button-submitProject', 'nClicks'),
    State('seg-input-projectName', 'value'),
    State('seg-project-filename', 'type'),
    State('seg-select-modelType', 'value'),
    State('seg-input-diameter', 'value'),
    State('seg-input-batchsize', 'value'),
    State('seg-checkbox-useGPU', 'checked'),
    prevent_initial_call=True
)
def submit_project(nc, projectName, fileStatus, modelType, diameter, batchsize, useGPU):
    """
    检查提交的项目，并将项目持久化到本地
    """
    if nc:
        process_submited_project(projectName, fileStatus, modelType, diameter, batchsize, useGPU)

@callback(
    Output('seg-project-filename', 'children', allow_duplicate='True'),
    Output('seg-project-filename', 'type', allow_duplicate='True'),
    Input('seg-dragger-upload-project', 'lastUploadTaskRecord'),
    prevent_initial_call=True
)
def update_upload_status(lastRecord):
    """
    监听项目文件上传状态
    """
    if lastRecord['taskStatus']=='success':
        return lastRecord['fileName'], 'success'
    else:
        set_head_notice(lastRecord['fileName']+' upload failed, please check file format!', type='error')
        return 'No file', 'secondary'

@callback(
    Input('seg-button-project', 'nClicks'),
    prevent_initial_call=True
)
def open_upload_box(nc):
    """
    打开上传项目文件窗口
    """
    if nc:
        fileSelecter.open_import_box()

@callback(
    Output('seg-button-submitProject', 'disabled'),
    Output('seg-input-projectName', 'status'),
    Input('seg-input-projectName', 'value')
)
def validate_project_name(project_name):
    """
    验证项目名称是否为空
    """
    if not project_name:
        return True, 'error'
    return False, None

@callback(
    Output('seg-modal-newproject', 'visible'),
    Output('seg-project-filename', 'children'),
    Output('seg-project-filename', 'type'),
    Output('seg-input-projectName', 'value'),
    Input('seg-button-newproject', 'nClicks'),
    prevent_initial_call=True
)
def show_newproject_modal(nClicks):
    """
    显示创建项目对话框
    """
    if nClicks and verify_modify_permission():
        return True, 'No file', 'secondary', None
    return no_update, no_update, no_update, no_update
