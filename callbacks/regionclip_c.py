from controller.regionclip_ctl import *
from dash import callback, Input, Output, State, no_update, Patch
from pages.components.fileSelecter import fileSelecter
from dash.exceptions import PreventUpdate

@callback(
    Output('clip-stain-image', 'src', allow_duplicate=True),
    Output('clip-gem-image', 'src', allow_duplicate=True),
    Input('clip-select-slice', 'value'),
    Input('clip-select-clipname', 'value'),
    State('clip-select-taskname', 'value'),
    prevent_initial_call=True
)
def update_clipped_images(slicename, clipName, taskName):
    """
    基于切片和裁剪名称更新裁剪图像
    """
    if not taskName or not slicename or not clipName:
        return no_update, no_update
    stain_path = clipData.get_task_clipName_stain_path(taskName, slicename, clipName)
    gem_image_path = clipData.get_task_clipName_gemImage_path(taskName, slicename, clipName)
    if stain_path is None:
        stain_src = no_update
        set_props('clip-stain-image', dict(style={'visibility':'hidden'}))
    else:
        stain_src = get_image_base64(stain_path)
        set_props('clip-stain-image', dict(style={'visibility':'visible'}))
    if gem_image_path is None:
        gem_src = no_update
        set_props('clip-gem-image', dict(style={'visibility':'hidden'}))
    else:
        gem_src = get_image_base64(gem_image_path)
        set_props('clip-gem-image', dict(style={'visibility':'visible'}))
    return stain_src, gem_src

@callback(
    Output('clip-stain-image', 'src'),
    Output('clip-gem-image', 'src'),
    Input('clip-button-startClip', 'nClicks'),
    State('clip-graph-original', 'relayoutData'),
    State('clip-select-taskname', 'value'),
    State('clip-select-slice', 'value'),
    State('clip-select-clipname', 'value'),
    running=[
        (Output('clip-button-startClip', 'loading'), True, False)
    ],
)
def start_clip(nc, relayoutData, taskName, slicename, clipName):
    """
    开始裁剪图像
    """
    if not nc:
        return no_update, no_update
    if not verify_modify_permission():
        return no_update, no_update
    if not taskName:
        set_head_notice('Please select a task first!', type='warning')
        return no_update, no_update
    if not slicename:
        set_head_notice('Please select a slice first!', type='warning')
        return no_update, no_update
    if not clipName:
        set_head_notice('Please select a clip name first!', type='warning')
        return no_update, no_update
    x_start = int(relayoutData.get('xaxis.range[0]', 0))
    x_end = int(relayoutData.get('xaxis.range[1]', 0))
    y_start = int(relayoutData.get('yaxis.range[1]', 0))
    y_end = int(relayoutData.get('yaxis.range[0]', 0))
    if (x_start >= x_end or y_start >= y_end):
        set_head_notice('Invalid clip region selected!', type='error')
        return no_update, no_update
    stain_src, gem_src = get_clipped_images(taskName, slicename, clipName, x_start, x_end, y_start, y_end)
    set_props('clip-stain-image', dict(style=None))
    set_props('clip-gem-image', dict(style=None))
    return stain_src, gem_src

@callback(
    Input('clip-delete-task-confirm', 'confirmCounts'),
    State('clip-select-taskname', 'value'),
    running=[
        (Output('clip-delete-task', 'loading'), True, False),
        (Output('clip-select-taskname', 'disabled'), True, False)
    ],
    prevent_initial_call=True
)
def delete_task(nc, taskName):
    """
    删除任务
    """
    if nc and taskName and verify_modify_permission():
        delete_task_from_disk(taskName)
        set_props('clip-stain-image', dict(style={'visibility':'hidden'}))
        set_props('clip-gem-image', dict(style={'visibility':'hidden'}))
        set_props('clip-graph-original', dict(style={'visibility':'hidden', 'height':'calc(90vh - 50px)', 'width':'100%'}))
        set_props('clip-select-slice', dict(value=None))
        set_props('clip-select-clipname', dict(value=None))
        set_props('clip-select-taskname', dict(value=None))
@callback(
    Output('clip-graph-original', 'figure'),
    Output('clip-graph-original', 'style'),
    Input('clip-select-slice', 'value'),
    State('clip-select-taskname', 'value'),
    running=[
        (Output('clip-text-slice', 'children'), 'Loading...', 'Select Slice')
    ]
)
def update_orifig_clipnames(slicename, taskname):
    """
    基于切片更新原始图像以及裁剪名称列表
    """
    if slicename:
        stain_fig = get_slice_figure(taskname, slicename)
        showStyle = {'visibility':True, 'height':'calc(90vh - 50px)', 'width':'100%'}
        return stain_fig, showStyle
    return no_update, no_update
@callback(
    Output('clip-select-slice', 'value'),
    Output('clip-select-slice', 'options'),
    Input('clip-select-taskname', 'value')
)
def update_slicelist_options(taskname):
    """
    基于任务更新切片列表
    """
    if taskname:
        slices = list(clipData.get_task_slices(taskname))
        slices.sort()
        slice = slices[0] if slices else None
        return slice, slices
    return no_update, no_update

@callback(
    Output('clip-select-clipname', 'options', allow_duplicate=True),
    Input('clip-select-clipname-tooltip', 'open'),
    State('clip-select-taskname', 'value'),
    prevent_initial_call=True
)
def update_tasklist_options(open, taskName):
    """
    更新任务裁剪项目列表
    """
    if not taskName:
        return no_update
    if open:
        clips = list(clipData.get_taskclips(taskName))
        clips.sort()
        return clips
    return no_update

@callback(
    Output('clip-select-taskname', 'options'),
    Input('clip-select-taskname-tooltip', 'open'),
    prevent_initial_call=True
)
def update_tasklist_options(open):
    """
    更新任务列表
    """
    if open:
        tasks = list(clipData.get_exist_tasks())
        tasks.sort()
        return tasks
    return no_update

@callback(
    Input('clip-button-submitTaskList', 'nClicks'),
    State('clip-input-taskname', 'value'),
    State('clip-tasklist-filename', 'type'),
    State('main-title-username', 'children'),
    prevent_initial_call=True
)
def submit_tasklist(nc, taskName, fileStatus, username):
    """
    检查提交的任务列表，并将任务持久化到本地
    """
    if nc:
        process_submited_tasklist(taskName, fileStatus, username)

@callback(
    Output('clip-tasklist-filename', 'children', allow_duplicate='True'),
    Output('clip-tasklist-filename', 'type', allow_duplicate='True'),
    Input('clip-dragger-upload', 'lastUploadTaskRecord'),
    prevent_initial_call=True
)
def upload_status(lastRecord):
    """
    监听文件上传状态
    """
    if lastRecord['taskStatus']=='success':
        return lastRecord['fileName'], 'success'
    else:
        set_head_notice(lastRecord['fileName']+' upload failed, please check file format!', type='error')
        return 'No file', 'secondary'

@callback(
    Input('clip-button-importTaskList', 'nClicks'),
    prevent_initial_call=True  
)
def open_import_box(nc):
    """
    打开导入任务列表文件窗口
    """
    if nc:
        fileSelecter.open_import_box()

@callback(
    Output('clip-select-clipname', 'value', allow_duplicate=True),
    Output('clip-select-clipname', 'options', allow_duplicate=True),
    Input('clip-delete-clipName-confirm', 'confirmCounts'),
    State('clip-select-taskname', 'value'),
    State('clip-select-clipname', 'value'),
    running=[
        (Output('clip-delete-clipName', 'loading'), True, False)
    ],
    prevent_initial_call=True
)
def delete_clip(nc, taskName, clipName):
    """
    删除裁剪任务
    """
    if not nc:
        return no_update, no_update
    if not verify_modify_permission():
        return no_update, no_update
    if not taskName:
        set_head_notice('Please select a task first!', type='warning')
        return no_update, no_update
    if not clipName:
        set_head_notice('Please select a clip name to delete!', type='warning')
        return no_update, no_update
    clipData.delete_taskclip(taskName, clipName)
    set_head_notice(clipName+' has been removed from your disk !', type='success')
    set_props('clip-stain-image', dict(style={'visibility':'hidden'}))
    set_props('clip-gem-image', dict(style={'visibility':'hidden'}))
    clips = clipData.get_taskclips(taskName)
    return None, clips
@callback(
    Output('clip-select-clipname', 'value'),
    Output('clip-select-clipname', 'options'),
    Output('clip-modal-inputClipName', 'visible', allow_duplicate=True),
    Input('clip-button-submitClipName', 'nClicks'),
    State('clip-input-clipName', 'value'),
    State('clip-select-taskname', 'value'),
    prevent_initial_call=True
)
def submit_clip_name(nc, clipName, taskName):
    """
    提交分片名称
    """
    if not nc:
        return no_update, no_update, no_update
    checked = True
    if not taskName:
        checked = False
        set_head_notice('Task Name cannot be empty', type='error')
    if not clipName:
        checked = False
        set_head_notice('Clip Name cannot be empty', type='error')
    if not checked:
        return no_update, no_update, no_update
    clipData.create_taskclip(taskName, clipName)
    clips = clipData.get_taskclips(taskName)
    set_head_notice(clipName+' is added successfully!', type='success')
    return clipName, clips, False

@callback(
    Output('clip-modal-inputClipName', 'visible'),
    Input('clip-add-clipName', 'nClicks'),
    prevent_initial_call=True
)
def show_setPatchName_modal(nClicks):
    """
    显示设置分块名称对话框
    """
    if nClicks and verify_modify_permission():
        return True
    return no_update

@callback(
    Output('clip-modal-newtask', 'visible'),
    Output('clip-tasklist-filename', 'children'),
    Output('clip-tasklist-filename', 'type'),
    Input('clip-button-newtask', 'nClicks'),
    prevent_initial_call=True
)
def show_newtask_modal(nClicks):
    """
    显示创建任务对话框
    """
    if nClicks and verify_modify_permission():
        return True, 'No file', 'secondary'
    return no_update, no_update, no_update