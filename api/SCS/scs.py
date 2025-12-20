import math
import numpy as np
from api.SCS import preprocessing, transformer, postprocessing
import os
from utils.typing import Status
from utils.commonfuc import write_json
from dataManager.expansion_d import expData
from websocket.websocket import ws
import traceback

def set_progress(taskInfo, progress: float):
    """
    设置任务进度
    """
    taskInfo['progress'] = progress
def set_slice_status(taskInfo, slice, stepName, status, text, notify=False):
    """
    设置切片处理状态
    """
    slice[stepName] = {'status': status, 'text': text}
    if notify:
        taskInfo['running'] = True

def complete_expansion_task(taskInfo, taskName):
    """
    完成胞域扩增任务
    """
    taskInfo['running'] = False
    write_json(expData.get_expTaskMetadata_path(taskName), dict(taskInfo))
    expData.delete_running_task(taskName, ws.notifyUpdateExpansionStatus)

def stop_expansion_task(taskInfo, taskName):
    """
    异常终止胞域扩增任务
    """
    for slice in taskInfo['slices']:
        steps = ['preprocess', 'train', 'postprocess', 'patchprocess']
        stopMarker = False
        for step in steps:
            if slice[step]['status'] == Status.PROCESSING:
                stopMarker = True
                set_slice_status(taskInfo, slice, 'patchprocess', Status.ERROR, slice['patchprocess']['text'])
                if step != 'patchprocess':
                    set_slice_status(taskInfo, slice, step, Status.ERROR, 'failed')
                break
        if stopMarker:
            break
    taskInfo['running'] = False
    taskInfo['exception'] = traceback.format_exc()
    write_json(expData.get_expTaskMetadata_path(taskName), dict(taskInfo))
    expData.delete_running_task(taskName, ws.notifyUpdateExpansionStatus)

def segment_cells(adata, bin_file, project_path, taskInfo, slice, layer_name='watershed_mask', prealigned=False, align=None, patch_size=0, bin_size=3, n_neighbor=50, epochs=100, r_estimate=15, val_ratio=0.0625):
    """
    Parameters:
        bin_file - string, tsv file for detected RNAs
        image_file - string, staining image of the tissue (.tiff)
        prealigned - boolean, if the staining image is prealigned with the sequencing spots coordinates, default False
        align - string, alignment method used to align the input staining image and sequencing spots ('rigid', 'non-rigid' or 'None'), default None
        patch_size - int, length and width (spots) of patches, if greater than zero, the input section will be cut into patches and processed patch by patch, if zero, the section will be process as a whole, default 0
        bin_size - int, the length and width (spots) of regions that will be merged as a bin, default 3
        n_neighbor - int, the number of nearest neighbors who will be considered when make predictions for one spot in the transformer model, default 50
        epochs - int, the training epochs of the transformer model, default 100
        r_estimate - int, the estimated radius (spots) of cells, used to calculate the priors for transformer predictions, default 15
        val_ratio - float, the fraction of the patch set aside for validation, default 0.0625 (1/4 height x 1/4 width)
    """
    tmp_dir = os.path.join(project_path, 'tmp')
    result_dir = os.path.join(project_path, 'result')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    if patch_size == 0:
        set_slice_status(taskInfo, slice, 'preprocess', Status.PROCESSING, 'running')
        set_slice_status(taskInfo, slice, 'patchprocess', Status.PROCESSING, '1/1', notify=True)
        preprocessing.preprocess(adata, bin_file, tmp_dir, layer_name, prealigned, align, 0, 0, patch_size, bin_size, n_neighbor)
        set_slice_status(taskInfo, slice, 'preprocess', Status.SUCCESS, 'success')
        set_slice_status(taskInfo, slice, 'train', Status.PROCESSING, 'running', notify=True)
        transformer.train(project_path, 0, 0, patch_size, epochs, val_ratio)
        set_slice_status(taskInfo, slice, 'train', Status.SUCCESS, 'success')
        set_slice_status(taskInfo, slice, 'postprocess', Status.PROCESSING, 'running', notify=True)
        postprocessing.postprocess(project_path, layer_name, 0, 0, patch_size, bin_size, r_estimate)
        set_slice_status(taskInfo, slice, 'postprocess', Status.SUCCESS, 'success')
    else:
        rmax, cmax = adata.shape
        n_patches = math.ceil(rmax / patch_size) * math.ceil(cmax / patch_size)
        current_patch = 0
        for startr in range(0, rmax, patch_size):
            for startc in range(0, cmax, patch_size):
                current_patch += 1
                set_slice_status(taskInfo, slice, 'patchprocess', Status.PROCESSING, f'{current_patch}/{n_patches}')
                set_slice_status(taskInfo, slice, 'preprocess', Status.PROCESSING, 'running')
                set_slice_status(taskInfo, slice, 'train', Status.WARNING, 'waiting')
                set_slice_status(taskInfo, slice, 'postprocess', Status.WARNING, 'waiting', notify=True)
                preprocessing.preprocess(adata, bin_file, tmp_dir, layer_name, prealigned, align, startr, startc, patch_size, bin_size, n_neighbor)
                set_slice_status(taskInfo, slice, 'preprocess', Status.SUCCESS, 'success')
                set_slice_status(taskInfo, slice, 'train', Status.PROCESSING, 'running', notify=True)
                transformer.train(project_path, startr, startc, patch_size, epochs, val_ratio)
                set_slice_status(taskInfo, slice, 'train', Status.SUCCESS, 'success')
                set_slice_status(taskInfo, slice, 'postprocess', Status.PROCESSING, 'running', notify=True)
                postprocessing.postprocess(project_path, layer_name, startr, startc, patch_size, bin_size, r_estimate)
                set_slice_status(taskInfo, slice, 'postprocess', Status.SUCCESS, 'success')