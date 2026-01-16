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
from utils.colors import *
import traceback
from dataManager.workspace import *
from utils.commonfuc import *
from dataManager.segmentation_d import segData
from dataManager.maskviewer_d import maskData
import traceback
from utils.typing import Status

def run_segmentation_project(observed_info):
    """
    运行分割项目
    """
    observed_info['running'] = True
    try:
        seg_step = 'segmentation'
        post_step = 'postprocess'
        watershed_layer = 'watershed_mask'
        cellpose_layer = 'cellpose_mask'
        project = observed_info['name']
        model = observed_info['model']
        diameter = None
        if observed_info['diameter']>0:
            diameter = observed_info['diameter']
        batchsize = observed_info['batchsize']
        gpu = observed_info['gpu']
        slices = observed_info['slices']
        current = 0
        size = 2*len(slices)
        for slice in slices:
            img_file = slice['img']
            bin_file = slice['gem']
            set_slice_status(observed_info, slice, seg_step, Status.PROCESSING, notify=True)
            adata = st.io.read_bgi_agg(bin_file, img_file)
            adata = segment_cells_cellpose(adata, output_layer=cellpose_layer, modeltype=model, batch_size=batchsize, diameter=diameter, gpu=gpu)
            adata = segment_cells_watershed(adata, output_layer=watershed_layer)
            z = slice['z']
            result_path = segData.get_result_path(project, z)
            adata.write_h5ad(result_path)
            set_slice_status(observed_info, slice, seg_step, Status.SUCCESS)
            set_slice_status(observed_info, slice, post_step, Status.PROCESSING)
            current += 1
            observed_info['progress'] = round(current/size, 2)

            fig, ax = plt.subplots(dpi=300)
            ax.imshow(adata.layers['stain'], cmap='cividis')
            adata.layers['unspliced'] = adata.X
            st.pl.imshow(adata, 'unspliced', ax=ax, alpha=0.6, cmap='Reds', vmax=10, use_scale=False, save_show_or_return='return')
            ax.set_title('')
            ax.axis('off')
            overlay = pltfig_to_array(fig)
            overlay_fig = get_imgfig_withplotly(overlay)
            overlay_path = maskData.get_overlay_path(project, z)
            write_pkl(overlay_fig, overlay_path)
            cellpose_mask_array = adata.layers[cellpose_layer]
            watershed_mask_array = adata.layers[watershed_layer]
            cellpose_mask, cellpose_contour = generate_cell_masks_rgba(cellpose_mask_array)
            watershed_mask, watershed_contour = generate_cell_masks_rgba(watershed_mask_array)
            stain = adata.layers['stain']
            stain_path = segData.get_stain_path(project, z)
            stain_fig = px.imshow(stain, color_continuous_scale='gray')
            write_pkl(stain_fig, stain_path)
            cellpose_fig = get_mask_contour_figure(stain, cellpose_mask, cellpose_contour)
            cellpose_path = maskData.get_cellpose_path(project, z)
            write_pkl(cellpose_fig, cellpose_path)
            watershed_fig = get_mask_contour_figure(stain, watershed_mask, watershed_contour)
            watershed_path = maskData.get_watershed_path(project, z)
            write_pkl(watershed_fig, watershed_path)
            del adata, stain, overlay, cellpose_mask_array, watershed_mask_array, cellpose_contour, watershed_contour, cellpose_mask, watershed_mask, watershed_fig, cellpose_fig
            set_slice_status(observed_info, slice, post_step, Status.SUCCESS)
            current += 1
            observed_info['progress'] = round(current/size, 2)
    except Exception as e:
        for slice in slices:
            if slice[seg_step]['status']==Status.PROCESSING:
                set_slice_status(observed_info, slice, seg_step, Status.ERROR)
                break
            if slice[post_step]['status']==Status.PROCESSING:
                set_slice_status(observed_info, slice, post_step, Status.ERROR)
                break
        observed_info['exception'] = traceback.format_exc()
    finally:
        if observed_info['exception'] is None:
            observed_info['running'] = False
        project_info = dict(observed_info)
        project_info['running'] = False
        segData.write_metadata_to_disk(project, project_info)
        segData.remove_running_project(project)

def segment_cells_cellpose(adata, modeltype='cyto3', diameter=None, batch_size=8, layer='stain', output_layer='cellpose_mask', gpu=True):
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
    return adata

def segment_cells_watershed(adata, layer='stain', output_layer='watershed_mask'):
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
    return adata

def set_slice_status(observed_info, slice, step: str, status: Status, notify=False):
    """
    设置进度状态
    """
    step_status = slice[step]
    status_text = status
    if status == Status.PROCESSING:
        status_text = 'running'
    if status == Status.ERROR:
        status_text = 'failed'
    step_status['status'] = status
    step_status['text'] = status_text
    if notify:
        observed_info['running'] = True

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