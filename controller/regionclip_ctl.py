from controller.auth import *
from dataManager.regionclip_d import clipData
from utils.commonfuc import get_current_date
from dash import no_update
import plotly.express as px
import imageio.v3 as iio
import pandas as pd
from skimage.transform import resize
import numpy as np
import re
import matplotlib.pyplot as plt
import base64
from io import BytesIO

def get_gem_on_stain_to_base64(stain_path, gem_path, alpha=0.5):
    """
    将 GEM 数据叠加到染色图像上，并返回 data:image/png;base64,... 格式的字符串
    """
    df = pd.read_csv(gem_path, sep='\t', comment='#')
    if not {'x', 'y'}.issubset(df.columns):
        raise ValueError("GEM 文件缺少 x 或 y 列")

    x_raw = df['x'].values
    y_raw = df['y'].values
    x_norm = (x_raw - x_raw.min()).astype(int)
    y_norm = (y_raw - y_raw.min()).astype(int)

    max_x = x_norm.max()
    max_y = y_norm.max()
    
    mtx = np.zeros((max_x + 1, max_y + 1), dtype=np.float32)
    valid_mask = (x_norm >= 0) & (y_norm >= 0)
    mtx[x_norm[valid_mask], y_norm[valid_mask]] = 1.0

    gem_img = mtx

    stain_img = iio.imread(stain_path)
    if stain_img.ndim == 2:
        stain_rgb = np.stack([stain_img, stain_img, stain_img], axis=-1)
    else:
        stain_rgb = stain_img.copy()
    stain_h, stain_w = stain_rgb.shape[:2]

    gem_resized = resize(gem_img, (stain_h, stain_w), order=0, anti_aliasing=False, preserve_range=True)

    overlay = np.zeros_like(stain_rgb, dtype=np.float32)
    overlay[:, :, 0] = gem_resized  # Red
    overlay[:, :, 1] = gem_resized  # Green

    blended = (1 - alpha) * stain_rgb.astype(np.float32) + alpha * (overlay * 255)
    blended = np.clip(blended, 0, 255).astype(np.uint8)

    plt.figure(figsize=(10, 10))
    plt.imshow(blended)
    plt.axis('off')
    plt.tight_layout(pad=0)

    buffer = BytesIO()
    plt.savefig(buffer, format='png', bbox_inches='tight', pad_inches=0)
    buffer.seek(0)
    image_png = buffer.getvalue()
    graph = base64.b64encode(image_png).decode('utf-8')
    plt.close()

    return f'data:image/png;base64,{graph}'


def get_slice_figure(taskname, slicename):
    """
    获取指定任务和切片的图像数据
    """
    taskInfo = clipData.get_slice_info(taskname, slicename)
    stain_fig = no_update
    overlay_fig = no_update
    if taskInfo is not None:
        img_path = taskInfo.get('image', None)
        gem_path = taskInfo.get('gem', None)
        if img_path is not None:
            img = iio.imread(img_path)
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
        if gem_path is not None:
            overlay_fig = get_gem_on_stain_to_base64(img_path, gem_path)
    return stain_fig, overlay_fig

def delete_task_from_disk(taskName):
    """
    从磁盘删除任务
    """
    clipData.delete_task(taskName)
    set_head_notice(taskName+' has been removed from your disk !', type='success')
    set_props('clip-select-taskname', dict(value=None))
    set_props('clip-select-slice', dict(value=None, options=[]))
    set_props('clip-graph-original', dict(figure={}, style={'visibility': 'hidden'}))
def process_submited_tasklist(taskName, fileStatus, username):
    """
    处理提交的任务列表，并将任务持久化到本地
    """
    if taskName is None or taskName.strip()=='':
        set_head_notice('Task Name cannot be empty', type='error')
        return False
    taskList = clipData.get_exist_tasks()
    if taskName in taskList:
        set_head_notice('Task Name already exists', type='warning')
        return False
    if fileStatus!='success':
        set_head_notice('Please upload your file first', type='warning')
        return False
    
    clipData.set_temptask_metadata(taskName,{
        'Creator': username,
        'Date': get_current_date()
    })

    clipData.save_temptask()
    set_head_notice('Task '+taskName+' created successfully!', type='success')
    set_props('clip-modal-newtask', dict(visible=False))
    set_props('clip-select-taskname', dict(options=list(taskList), value=taskName))

def parse_regionclip_tasklist(lines:list):
    """
    解析切片裁剪任务列表
    """
    # data=[
    #     {
    #         'z': 1,
    #         'image': 'image2',
    #         'gem': 'gem1',
    #         'registration': {'status': 'success', 'text': 'success'},
    #         'segmentation': {'status': 'success', 'text': 'success'},
    #     }
    # ],
    try:
        data = []
        for line in lines:
            z, image, gem, *extra = re.split(r'[,]+', line.strip())
            data.append({
                'z': int(z),
                'image': image.strip(),
                'gem': gem.strip(),
            })
        clipData.set_temptask_data(data)
    except Exception as e:
        clipData.reset_temptask_data()
        raise e
    
def read_regionclip_file(fpath:str):
    """
    读取裁剪任务列表文件
    """
    success = True
    try:
        with open(fpath, 'r') as f:
            parse_regionclip_tasklist(f.readlines())
    except Exception as e:
        success = False

    filename = os.path.basename(fpath)
    if success:
        set_props('clip-tasklist-filename', dict(type='success', children=filename))
        set_head_notice(filename+' import successfully!', type='success')
    else:
        set_props('clip-tasklist-filename', dict(type='secondary', children='No file'))
        set_head_notice(filename+' import failed, please check your file format', type='error')
    return success