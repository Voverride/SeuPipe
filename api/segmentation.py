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
from scipy import sparse
from PIL import Image
import io
import pickle
from utils.colors import *
import traceback
from dataManager.workspace import *
from utils.commonfuc import *
from dataManager.segmentation_d import segData
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
                adata_path = segData.get_seg_adata_path(taskName, z)
                write_result(adata, taskName, z)
                adata.write_h5ad(adata_path)

            except Exception as e:
                set_progress_status(seg, 'error')
                taskInfo['exception'] = traceback.format_exc()
                return
            
            set_progress_status(seg, 'success')
            current += 1
            metadata['progress'] = round(current/size, 2)  
            taskInfo_copy = taskInfo.copy()
            taskInfo_copy['thread'] = None
            task_folder = segData.get_seg_task_folder(taskName)
            with open(os.path.join(task_folder, 'tasklist.pkl'), 'wb') as f:
                pickle.dump(taskInfo_copy, f)

    except Exception as e:
        taskInfo['exception'] = traceback.format_exc()
        return
    
def write_result(ad, taskName, zIndex):
    """
    写入分割结果
    """
    before_path = segData.get_registration_before_path(taskName, zIndex)
    aligned_path = segData.get_registration_aligned_path(taskName, zIndex)
    watershed_mask_figure_path = segData.get_seg_watershed_mask_figure_path(taskName, zIndex)
    cellpose_mask_figure_path = segData.get_seg_cellpose_mask_figure_path(taskName, zIndex)
    watershed_contour_figure_path = segData.get_seg_watershed_contour_figure_path(taskName, zIndex)
    cellpose_contour_figure_path = segData.get_seg_cellpose_contour_figure_path(taskName, zIndex)
    stain_figure_path = segData.get_seg_stain_figure_path(taskName, zIndex)

    write_pkl(ad.layers['stain'], stain_figure_path)
    write_pkl(ad.uns['before'], before_path)
    write_pkl(ad.uns['aligned'], aligned_path)
    write_pkl(ad.uns['watershed_mask'], watershed_mask_figure_path)
    write_pkl(ad.uns['cellpose_mask'], cellpose_mask_figure_path)
    write_pkl(ad.uns['watershed_contour'], watershed_contour_figure_path)
    write_pkl(ad.uns['cellpose_contour'], cellpose_contour_figure_path)
    del ad.uns

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
    mask, contour = generate_cell_masks_rgba(mask_array, hex_colors=hexcolors)
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
    mask, contour = generate_cell_masks_rgba(mask_array, hex_colors=hexcolors)
    adata.uns['watershed_mask'] = mask
    adata.uns['watershed_contour'] = contour
    return adata