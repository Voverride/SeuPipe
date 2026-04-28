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
import dill

def convert_mtx_to_base64_image(mtx):
    """
    将矩阵转换为 base64 编码的 PNG 图像
    输出严格正方形，图像居中，四周留白均匀且透明
    """
    
    height, width = mtx.shape[:2]
    max_dim = max(height, width)
    img_width_ratio = width / max_dim
    img_height_ratio = height / max_dim
    fig, ax = plt.subplots(figsize=(10, 10), dpi=100)
    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)
    ax.imshow(mtx, cmap='gray', aspect='auto')
    ax.axis('off')
    left = (1 - img_width_ratio) / 2
    bottom = (1 - img_height_ratio) / 2
    ax.set_position([left, bottom, img_width_ratio, img_height_ratio])
    buffer = BytesIO()
    plt.savefig(buffer, format='png', pad_inches=0, 
                facecolor='none', edgecolor='none', dpi=100, 
                transparent=True)
    buffer.seek(0)
    image_png = buffer.getvalue()
    graph = base64.b64encode(image_png).decode('utf-8')
    plt.close(fig)
    return f'data:image/png;base64,{graph}'

def write_func(func, path):
    """
    将函数写入文件
    """
    with open(path, 'wb') as f:
        dill.dump(func, f)

def read_func(path):
    """
    从文件中读取函数
    """
    with open(path, 'rb') as f:
        func = dill.load(f)
    return func

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

def compress_image(image, max_dimension=1024):
    """
    将单通道图像等比例压缩，确保最长边不超过指定尺寸。
    保持原始宽高比，不拉伸变形。

    参数:
        image (np.ndarray): 原始单通道图像，shape (H, W)
        max_dimension (int): 目标最大边长，默认 1024

    返回:
        compressed_img (np.ndarray): 压缩后的图像
        map_region_to_original (function): 
            输入：压缩图中的闭区间矩形区域 (i1, j1, i2, j2)
            输出：原图中的闭区间矩形区域 (orig_i1, orig_j1, orig_i2, orig_j2)
    """
    import cv2
    import numpy as np
    
    assert len(image.shape) == 2, "图像必须是单通道（2D）"
    orig_h, orig_w = image.shape
    
    if orig_h <= max_dimension and orig_w <= max_dimension:
        return image, lambda i1, j1, i2, j2: (i1, j1, i2, j2)
    
    scale = max_dimension / max(orig_h, orig_w)
    target_height = int(orig_h * scale)
    target_width = int(orig_w * scale)
    
    compressed_img = cv2.resize(image, (target_width, target_height), interpolation=cv2.INTER_AREA)
    
    scale_h = orig_h / target_height
    scale_w = orig_w / target_width
    
    def map_region_to_original(i1, j1, i2, j2):
        """
        将压缩图中的闭区间矩形区域映射回原图的闭区间区域。
        如果输入范围越界，会自动裁剪到图像边界。
        """
        i1 = max(0, min(i1, target_height - 1))
        i2 = max(0, min(i2, target_height - 1))
        j1 = max(0, min(j1, target_width - 1))
        j2 = max(0, min(j2, target_width - 1))
        
        if i1 > i2:
            i1, i2 = i2, i1
        if j1 > j2:
            j1, j2 = j2, j1
        
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
    def json_serializer(obj):
        """自定义JSON序列化器，处理NumPy和pandas数据类型"""
        if isinstance(obj, (np.integer, np.int8, np.int16, np.int32, np.int64,
                           np.uint8, np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.floating, np.float16, np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.bool_):
            return bool(obj)
        elif hasattr(obj, 'to_json'):  # 处理pandas Series/DataFrame
            return json.loads(obj.to_json())
        elif hasattr(obj, 'to_dict'):  # 处理其他可转换为字典的对象
            return obj.to_dict()
        elif isinstance(obj, (set, tuple)):  # 处理集合和元组
            return list(obj)
        else:
            raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")
    with open(file_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, default=json_serializer, ensure_ascii=False, indent=4)
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