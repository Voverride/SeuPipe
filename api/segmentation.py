import warnings
from numba import errors
import tensorflow as tf
import torch
import matplotlib.pyplot as plt
from tensorflow.python.ops.numpy_ops import np_config
gpus = tf.config.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
    except RuntimeError as e:
        print(e)
np_config.enable_numpy_behavior()
warnings.filterwarnings("ignore", category=errors.NumbaWarning)
import spateo as st
from cellpose import models
import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import find_contours
from scipy import sparse
from PIL import Image
import plotly.graph_objects as go
import matplotlib.colors as mcolors
import io
import pickle
from utils.colors import *
import traceback
from dataManager.workspace import *
def run_segmentation_task(taskInfo):
    """
    启动分割任务
    """
    try:
        taskName = taskInfo['taskname']
        metadata = taskInfo['metadata']
        model = metadata['model']
        diameter = None
        if metadata['diameter']>0:
            diameter = metadata['diameter']
        batchsize = metadata['batchsize']
        gpu = metadata['GPU']
        data = taskInfo['data'] 
        current = 0
        size = 2*len(data)
        for idx, item in enumerate(data):
            try:
                image_file = item['image']
                bin_file = item['gem']
                reg = item['registration']
                set_progress_status(reg, 'processing')
                adata = refine_alignment(image_file, bin_file, cmap='cividis', dpi=300)
            except Exception as e:
                set_progress_status(reg, 'error')
                taskInfo['exception'] = traceback.format_exc()
                return
            set_progress_status(reg, 'success')
            current += 1
            metadata['progress'] = round(current/size, 2)
            
            try:
                seg = item['segmentation']
                set_progress_status(seg, 'processing')
                adata = segment_cells_cellpose(adata, modeltype=model, batch_size=batchsize, diameter=diameter, gpu=gpu, hexcolors=primaryColors)
                adata = segment_cells_watershed(adata, hexcolors=primaryColors)
                z = item['z']
                root_path = get_segmentation_workspace()
                slice_folder = os.path.join(root_path, taskName, 'slices')
                check_path(slice_folder)
                adata.write_h5ad(os.path.join(slice_folder, f'z_{z}.h5ad'))

            except Exception as e:
                set_progress_status(seg, 'error')
                taskInfo['exception'] = traceback.format_exc()
                return
            
            set_progress_status(seg, 'success')
            current += 1
            metadata['progress'] = round(current/size, 2)  
            taskInfo_copy = taskInfo.copy()
            taskInfo_copy['thread'] = None
            with open(os.path.join(root_path, taskName, 'tasklist.pkl'), 'wb') as f:
                pickle.dump(taskInfo_copy, f)

    except Exception as e:
        taskInfo['exception'] = traceback.format_exc()
        return

def set_progress_status(progress:dict, status):
    """
    设置进展状态
    """
    progress['status'] = status
    if status == 'success':
        progress['text'] = 'success'
    elif status == 'processing':
        progress['text'] = 'running'
    elif status == 'warning':
        progress['text'] = 'waiting'
    elif status == 'error':
        progress['text'] = 'failed'

def refine_alignment(image_file, bin_file, cmap='cividis', dpi=300):
    """
    对齐图像和基因表达，并保存图片，返回对齐后结果
    """
    adata = st.io.read_bgi_agg(bin_file, image_file)
    adata.layers['unspliced'] = adata.X
    before = adata.layers['stain'].copy()
    st.cs.refine_alignment(adata, mode='rigid', transform_layers=['stain'])
    fig1, ax1 = plt.subplots(dpi=dpi)
    ax1.imshow(before, cmap=cmap)
    st.pl.imshow(adata, 'unspliced', ax=ax1, alpha=0.6, cmap='Reds', vmax=10, use_scale=False, save_show_or_return='return')
    ax1.set_title('')
    ax1.axis('off')
    before_aligned = pltfig_to_array(fig1, dpi=dpi)
    adata.uns['before'] = before_aligned
    plt.close()
    fig2, ax2 = plt.subplots()
    ax2.imshow(adata.layers['stain'], cmap=cmap)
    st.pl.imshow(adata, 'unspliced', ax=ax2, alpha=0.6, cmap='Reds', vmax=10, use_scale=False, save_show_or_return='return')
    ax2.set_title('')
    ax2.axis('off')
    after_aligned = pltfig_to_array(fig2, dpi=dpi)
    adata.uns['aligned'] = after_aligned
    plt.close()
    return adata

def pltfig_to_array(fig, dpi=300):
    """
    将plt图像转成numpy矩阵
    """
    with io.BytesIO() as buf:
        fig.savefig(buf, format='png', bbox_inches='tight', pad_inches=0, dpi=dpi)
        buf.seek(0)
        pil_image = Image.open(buf)
        arr = np.array(pil_image)
        pil_image.close()
    return arr

def segment_cells_cellpose(adata, modeltype='cyto3', diameter=None, batch_size=8, layer='stain', output_layer='cellpose_mask', gpu=True, hexcolors=None):
    """
    基于cellpose进行细胞分割
    """
    use_gpu = False
    if gpu:
        try:
            use_gpu = torch.cuda.is_available()
        except Exception as e:
            use_gpu = False
    try:
        model = models.Cellpose(model_type=modeltype, gpu=use_gpu)
        masks, flows, styles, diams = model.eval(adata.layers[layer], batch_size=batch_size, diameter=diameter)
        adata.layers[output_layer] = sparse.csr_matrix(masks.astype(np.uint16))
    except Exception as e:
        raise e
    finally:
        del model
        if use_gpu:
            torch.cuda.empty_cache()
    mask_array = adata.layers[output_layer]
    mask, contour = generate_cell_masks(mask_array, hex_colors=hexcolors)
    adata.uns['cellpose_mask'] = mask
    adata.uns['cellpose_contour'] = contour
    return adata

def segment_cells_watershed(adata, layer='stain', output_layer='watershed_mask', hexcolors=None):
    """
    基于分水岭算法分割细胞
    """
    st.cs.mask_nuclei_from_stain(adata, otsu_classes = 4, otsu_index=1)
    st.cs.find_peaks_from_mask(adata, layer, 7)
    st.cs.watershed(adata, layer, 5, out_layer=output_layer)
    try:
        del adata.layers['stain_mask']
        del adata.layers['stain_distances']
        del adata.layers['stain_markers']
    except Exception as e:
        pass
    adata.layers[output_layer] = sparse.csr_matrix(adata.layers[output_layer].astype(np.uint16))
    mask_array = adata.layers[output_layer]
    mask, contour = generate_cell_masks(mask_array, hex_colors=hexcolors)
    adata.uns['watershed_mask'] = mask
    adata.uns['watershed_contour'] = contour
    return adata

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