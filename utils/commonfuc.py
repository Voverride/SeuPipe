import pandas as pd
from PIL import Image
import plotly.express as px
import io
import base64
import numpy as np
from io import BytesIO
from datetime import datetime
import socket
import json
import torch
import tensorflow as tf
import gc
import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import find_contours
from scipy import sparse
from PIL import Image
import matplotlib.colors as mcolors
import io
from utils.colors import *
from dataManager.workspace import *
import plotly.graph_objects as go
import random
import pickle
import cv2

def get_mask_contour_figure(stain, mask, contour, title=None, showmask=False, showcontour=True, mask_opacity=0.7):
    """
    获取stain， mask，contour， 叠加交互图
    """
    if stain.ndim == 2:
        stain_normalized = (stain / np.max(stain)) * 255
        stain_rgb = np.stack([stain_normalized]*3, axis=-1)
    else:
        stain_rgb = stain[..., :3]
    
    fig = get_imgfig_withplotly(stain_rgb, title=title)

    fig.add_trace(go.Image(
        source=f"data:image/png;base64,{array_to_base64(mask)}",
        opacity=mask_opacity,
        visible=showmask,
    ))

    fig.add_trace(go.Image(
        source=f"data:image/png;base64,{array_to_base64(contour)}",
        visible=showcontour,
    ))
    return fig

def get_current_date():
    """
    获取当前日期
    """
    now = datetime.now()
    formatted_time = now.strftime("%Y-%m-%d %H:%M")
    return formatted_time

def compress_image(image, target_height=1024, target_width=1024):
    """
    将单通道图像压缩到指定尺寸，并返回区域映射函数。

    参数:
        image (np.ndarray): 原始单通道图像，shape (H, W)
        target_height (int): 目标高度
        target_width (int): 目标宽度

    返回:
        compressed_img (np.ndarray): 压缩后的图像 (target_height, target_width)
        map_region_to_original (function): 
            输入:
                top_left = (i1, j1)  # 压缩图中区域左上角（包含）
                bottom_right = (i2, j2)  # 压缩图中区域右下角（包含）
            输出:
                (orig_top_left, orig_bottom_right)
                其中 orig_top_left = (orig_i1, orig_j1)
                      orig_bottom_right = (orig_i2, orig_j2)
                所有坐标均为整数，且在原始图像范围内，闭区间（包含端点）
    """
    assert len(image.shape) == 2, "图像必须是单通道（2D）"
    orig_h, orig_w = image.shape

    if orig_h <= target_height and orig_w <= target_width:
        return image, lambda i1, j1, i2, j2: (i1, j1, i2, j2)
    compressed_img = cv2.resize(image, (target_width, target_height), interpolation=cv2.INTER_AREA)

    scale_h = orig_h / target_height
    scale_w = orig_w / target_width

    def map_region_to_original(i1, j1, i2, j2):
        """
        将压缩图中的闭区间矩形区域映射回原图的闭区间区域。

        参数:
            (i1, j1) —— 压缩图左上角，包含
            (i2, j2) —— 压缩图右下角，包含

        返回:
            orig_i1, orig_j1, orig_i2, orig_j2
            所有坐标为整数，闭区间，不越界。
        """
        if not (0 <= i1 <= i2 < target_height and 0 <= j1 <= j2 < target_width):
            raise ValueError(
                f"无效区域：top_left=({i1},{j1}), bottom_right=({i2},{j2})，"
                f"图像尺寸为 {target_height}x{target_width}"
            )

        orig_i_start_f = i1 * scale_h
        orig_i_end_f   = (i2 + 1) * scale_h
        orig_j_start_f = j1 * scale_w
        orig_j_end_f   = (j2 + 1) * scale_w

        orig_i1 = int(np.floor(orig_i_start_f))
        orig_i2 = int(np.ceil(orig_i_end_f) - 1)
        orig_j1 = int(np.floor(orig_j_start_f))
        orig_j2 = int(np.ceil(orig_j_end_f) - 1)

        orig_i1 = max(0, orig_i1)
        orig_j1 = max(0, orig_j1)
        orig_i2 = min(orig_h - 1, orig_i2)
        orig_j2 = min(orig_w - 1, orig_j2)

        if orig_i1 > orig_i2:
            orig_i2 = orig_i1
        if orig_j1 > orig_j2:
            orig_j2 = orig_j1

        return orig_i1, orig_j1, orig_i2, orig_j2

    return compressed_img, map_region_to_original

def write_pkl(obj, filepath):
    """
    将对象写入到指定的 .pkl 文件中。

    :param obj: 要写入的 Python 对象
    :param filepath: 目标 .pkl 文件的路径（字符串）
    """
    with open(filepath, 'wb') as file:
        pickle.dump(obj, file)

def read_pkl(filepath):
    """
    从指定的 .pkl 文件中读取对象。

    :param filepath: 源 .pkl 文件的路径（字符串）
    :return: 从文件中读取的 Python 对象
    """
    with open(filepath, 'rb') as file:
        obj = pickle.load(file)
    return obj

def generate_cell_masks(mask_array, scale=255, hex_colors=None, background=np.nan):
    """
    输入: 
        mask_array (numpy.ndarray or scipy.sparse.csr_matrix): 细胞编号矩阵（0=背景，非零值=细胞编号）
        scale (int): 输出RGB值的缩放因子（默认255）
        hex_colors (list): 16进制颜色列表，如 ['#FF0000', '#00FF00', '#0000FF'],默认None
        background: 背景填充值(默认NaN)
    输出:
        colored_rgb (numpy.ndarray): RGB格式的细胞区域图（shape: h,w,3），背景为NaN
        contour_rgb (numpy.ndarray): RGB格式的细胞轮廓图（shape: h,w,3），背景为NaN
    """
    if not sparse.issparse(mask_array):
        mask_array = sparse.csr_matrix(mask_array.astype(np.uint16))
    
    h, w = mask_array.shape
    colored_rgb = np.full((h, w, 3), background, dtype=np.float32)
    contour_rgb = np.full((h, w, 3), background, dtype=np.float32)
    
    cell_ids = np.unique(mask_array.data)
    cell_ids = cell_ids[cell_ids != 0]
    
    if hex_colors == None:
        colors = plt.cm.tab20(np.linspace(0, 1, len(cell_ids)))
        colors = colors[:, :3] 
    else:
        colors = np.array([mcolors.hex2color(c) for c in hex_colors])
    
    for i, cell_id in enumerate(cell_ids):
        cell_mask = (mask_array == cell_id)
        dense_mask = cell_mask.toarray()
        
        color = colors[i % len(colors)]
        
        rows, cols = cell_mask.nonzero()
        colored_rgb[rows, cols] = color
        
        contours = find_contours(dense_mask, level=0.5)
        for contour in contours:
            for y, x in contour.astype(int):
                if 0 <= y < h and 0 <= x < w:
                    contour_rgb[y, x] = color
    
    if scale != 1:
        colored_rgb = colored_rgb * scale
        contour_rgb = contour_rgb * scale
    
    return colored_rgb, contour_rgb

def generate_cell_masks_rgba(mask_array, scale=255, hex_colors=primaryColors):
    """
    输入: 
        mask_array (numpy.ndarray or scipy.sparse.csr_matrix): 细胞编号矩阵（0=背景，非零值=细胞编号）
        scale (int): 输出RGB和Alpha值的缩放因子（默认255）
        hex_colors (list): 16进制颜色列表，如 ['#FF0000', '#00FF00', '#0000FF'],默认None
    输出:
        colored_rgba (numpy.ndarray): RGBA格式的细胞区域图（shape: h,w,4）
        contour_rgba (numpy.ndarray): RGBA格式的细胞轮廓图（shape: h,w,4）
    """
    if not sparse.issparse(mask_array):
        mask_array = sparse.csr_matrix(mask_array.astype(np.uint16))
    
    h, w = mask_array.shape
    colored_rgba = np.zeros((h, w, 4), dtype=np.float32)
    contour_rgba = np.zeros((h, w, 4), dtype=np.float32)
    
    cell_ids = np.unique(mask_array.data)
    cell_ids = cell_ids[cell_ids != 0]
    
    if hex_colors is None:
        colors = plt.cm.tab20(np.linspace(0, 1, len(cell_ids)))
        colors = colors[:, :3]
    else:
        colors = np.array([mcolors.hex2color(c) for c in hex_colors])
    
    for cell_id in cell_ids:
        cell_mask = (mask_array == cell_id)
        dense_mask = cell_mask.toarray()
        color = random.choice(colors)
        
        rows, cols = cell_mask.nonzero()
        colored_rgba[rows, cols, :3] = color
        colored_rgba[rows, cols, 3] = 1.0
        
        contours = find_contours(dense_mask, level=0.5)
        for contour in contours:
            for y, x in contour.astype(int):
                if 0 <= y < h and 0 <= x < w:
                    contour_rgba[y, x, :3] = color
                    contour_rgba[y, x, 3] = 1.0
    
    if scale != 1:
        colored_rgba[..., :3] *= scale
        colored_rgba[..., 3] *= scale
        contour_rgba[..., :3] *= scale
        contour_rgba[..., 3] *= scale
    
    return colored_rgba, contour_rgba

def generate_expression_mask(expr_data, scale=255, expr_cmap='cividis', background=np.nan):
    """
    输入: 
        expr_data: 基因表达矩阵(稀疏或密集格式)
        scale (int): RGB值缩放因子(默认255)
        expr_cmap (str): 基因表达热图的colormap名称
        background: 背景填充值(默认NaN)
    输出:
        expr_rgb (numpy.ndarray): 基因表达热图(shape: h,w,3), 背景为指定值
    """
    if sparse.issparse(expr_data):
        expr_data = expr_data.toarray()
    
    h, w = expr_data.shape
    expr_rgb = np.full((h, w, 3), background, dtype=np.float32)
    
    nonzero_mask = expr_data > 0
    if np.any(nonzero_mask):
        vmin = np.quantile(expr_data[nonzero_mask], 0.05)
        vmax = np.quantile(expr_data[nonzero_mask], 0.95)
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        
        cmap = plt.get_cmap(expr_cmap)
        expr_rgb[nonzero_mask] = cmap(norm(expr_data[nonzero_mask]))[..., :3]
    
    if scale != 1:
        expr_rgb = expr_rgb * scale
    
    return expr_rgb

def release_memory():
    """
    释放显存和内存
    """
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    try:
        tf.keras.backend.clear_session()
    except:
        pass
    gc.collect()

def read_json(file_path):
    """
    读取 JSON 文件并返回数据
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    return data
def write_json(file_path, data):
    """
    将数据写入 JSON 文件
    """
    with open(file_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)
def get_local_ip():
    """
    获取本机网络 IP 地址
    """
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as s:
            s.connect(("8.8.8.8", 80))
            ip = s.getsockname()[0]
        return ip
    except Exception:
        return "127.0.0.1"
def array_to_base64(arr):
    """
    将 NumPy 数组转换为 Base64 编码的 PNG 图片
    """
    pil_img = Image.fromarray(arr)
    buf = BytesIO()
    pil_img.save(buf, format='PNG')
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()

def get_imgfig_withplotly(img_array, title=None):
    """
    基于图像矩阵绘制交互图形
    """
    fig = px.imshow(img_array)
    fig.update_layout(
        autosize=True,
        title=dict(
            text=title,
            x=0.5,
        ),
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
    )
    return fig

def array_to_base64(img_array):
    """将 NumPy 数组转换为 Base64 编码的 PNG 图片"""
    pil_img = Image.fromarray(img_array.astype(np.uint8))
    buffer = io.BytesIO()
    pil_img.save(buffer, format="PNG")
    return base64.b64encode(buffer.getvalue()).decode('utf-8')

def is_all_numeric(series)->bool:
    """
    判断dataframe中某一列的数据是否为纯数值
    """
    if pd.api.types.is_numeric_dtype(series):
        return True
    converted = pd.to_numeric(series, errors='coerce')
    original_na = series.isna().sum()
    new_na = converted.isna().sum()
    return new_na == original_na