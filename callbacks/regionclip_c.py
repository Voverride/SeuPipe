from controller.regionclip_ctl import *
from dash import callback, Input, Output, State, no_update, Patch
from pages.components.fileSelecter import fileSelecter
from dash.exceptions import PreventUpdate

import plotly.express as px
import imageio.v3 as iio
stain_path = '/data1/share/stereo-seq/2nd_array/ssDNA_registered_split/1-2_embryo_rectangle/embryo_1-2_E7.75_z16.tif'
img = iio.imread(stain_path)
stain_fig = px.imshow(
    img,
    color_continuous_scale='gray',
    binary_format="jpeg",
    binary_compression_level=5
)
stain_fig.update_layout(
    coloraxis_showscale=False,
    xaxis_visible=False,
    yaxis_visible=False
)

from PIL import Image
import numpy as np
combined_path = '/home/zhouyb/pyProject/SeuPipe/assets/combined.png'
img_pil = Image.open(combined_path)
img = np.array(img_pil)

zmin = img.min()
zmax = img.max()

combined_fig = px.imshow(
    img,
    color_continuous_scale='gray' if img.ndim == 2 else None,
    binary_format='png',
    zmin=zmin,
    zmax=zmax,
)

combined_fig.update_layout(
    xaxis_visible=False,
    yaxis_visible=False,
    coloraxis_showscale=False,
)
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

@callback(
    Output('clip-graph-original', 'figure'),
    Output('clip-image-gemOverlay', 'src'),
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
        stain_fig, overlay_fig = get_slice_figure(taskname, slicename)
        print("=============================================================================")
        # print(overlay_fig)
        return stain_fig, overlay_fig
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
    clips = list(clipData.get_taskclips(taskName))
    clips.sort()
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