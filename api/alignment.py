import time
import traceback
from dataManager.alignment_d import alidata
from utils.typing import StepStatus
import paste as pst
from paste2 import PASTE2, projection
import pandas as pd
import numpy as np
import scanpy as sc
import torch
import ot

COUNTS_KEY = 'SeuPipe_Alignment_Counts'

def set_steps_status(observed_metadata, step, status, notify=False):
    """
    设置对齐任务步骤状态
    """
    observed_metadata['steps'][step] = status
    if notify:
        observed_metadata['running'] = True

def run_ali_project(model, device, observed_metadata):
    """
    运行对齐项目
    """
    if model == 'paste2':
        paste2(observed_metadata)
    else:
        use_gpu = False
        if device == 'GPU':
            use_gpu = True
        paste1(observed_metadata, use_gpu=use_gpu)

def paste1(observed_metadata, use_gpu=False, use_rep=None):
    """
    paste1算法
    """
    try:
        cur_step = 'preprocess' 
        set_steps_status(observed_metadata, cur_step, StepStatus.PROCESS, notify=True)
        adata_path = alidata.get_initialdata_path(observed_metadata['projectName'])
        adata = sc.read_h5ad(adata_path)
        x_field = observed_metadata['x']
        y_field = observed_metadata['y']
        z_field = observed_metadata['z']
        ensure_numeric_fields(adata, x_field, y_field, z_field)
        field, _ = validate_exp_data(adata)
        if use_rep==None and field!=None:
            use_rep = field
        adata.obsm['spatial'] = adata.obs[[x_field, y_field]].values
        z = list(adata.obs[z_field].unique())
        z.sort(key=lambda num:float(num))
        slices = [adata[adata.obs[z_field]==slice] for slice in z]
        backend = ot.backend.TorchBackend()
        pis = []
        set_steps_status(observed_metadata, cur_step, StepStatus.FINISH)
        cur_step = 'alignment'
        set_steps_status(observed_metadata, cur_step, StepStatus.PROCESS, notify=True)
        for i in range(len(slices)-1):
            pis.append(pst.pairwise_align(slices[i], slices[i+1], alpha=0.3, use_gpu = use_gpu, backend=backend, use_rep=use_rep, gpu_verbose=False))
            percent = int((i+2)/len(slices)*100)
            set_steps_status(observed_metadata, 'percent', percent, notify=True)
        set_steps_status(observed_metadata, cur_step, StepStatus.FINISH)
        cur_step = 'postprocess'
        set_steps_status(observed_metadata, cur_step, StepStatus.PROCESS, notify=True)
        new_slices = pst.stack_slices_pairwise(slices, pis)
        for slice in new_slices:
            index = slice.obs_names
            x_aligned = slice.obsm['spatial'][:, 0]
            y_aligned = slice.obsm['spatial'][:, 1]
            adata.obs.loc[index, x_field] = x_aligned
            adata.obs.loc[index, y_field] = y_aligned
        alidata.update_coordinate(observed_metadata['projectName'], adata.obs)
        set_steps_status(observed_metadata, cur_step, StepStatus.FINISH, notify=True)
        set_steps_status(observed_metadata, cur_step, StepStatus.FINISH, notify=True)
    except:
        set_steps_status(observed_metadata, cur_step, StepStatus.ERROR)
        observed_metadata['exception'] = traceback.format_exc()
    finally:
        time.sleep(1)
        observed_metadata['running'] = False
        project = observed_metadata['projectName']
        alidata.update_project_info(project, dict(observed_metadata))
        alidata.remove_alignment_task(project)
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

def paste2(observed_metadata, use_rep=None):
    """
    paste2算法
    """
    try:
        cur_step = 'preprocess' 
        set_steps_status(observed_metadata, cur_step, StepStatus.PROCESS, notify=True)
        adata_path = alidata.get_initialdata_path(observed_metadata['projectName'])
        adata = sc.read_h5ad(adata_path)
        x_field = observed_metadata['x']
        y_field = observed_metadata['y']
        z_field = observed_metadata['z']
        ensure_numeric_fields(adata, x_field, y_field, z_field)
        field, _ = validate_exp_data(adata)
        if use_rep==None and field!=None:
            use_rep = field
        adata.obsm['spatial'] = adata.obs[[x_field, y_field]].values
        z = list(adata.obs[z_field].unique())
        z.sort(key=lambda num:float(num))
        slices = [adata[adata.obs[z_field]==slice] for slice in z]
        pis = []
        set_steps_status(observed_metadata, cur_step, StepStatus.FINISH)
        cur_step = 'alignment'
        set_steps_status(observed_metadata, cur_step, StepStatus.PROCESS, notify=True)
        for i in range(len(slices)-1):
            pis.append(PASTE2.partial_pairwise_align(slices[i], slices[i+1], s=0.7, use_rep=use_rep, verbose=False))
            percent = int((i+2)/len(slices)*100)
            set_steps_status(observed_metadata, 'percent', percent, notify=True)
        time.sleep(1)
        set_steps_status(observed_metadata, cur_step, StepStatus.FINISH)
        cur_step = 'postprocess'
        set_steps_status(observed_metadata, cur_step, StepStatus.PROCESS, notify=True)
        new_slices = projection.partial_stack_slices_pairwise(slices, pis)
        for slice in new_slices:
            index = slice.obs_names
            x_aligned = slice.obsm['spatial'][:, 0]
            y_aligned = slice.obsm['spatial'][:, 1]
            adata.obs.loc[index, x_field] = x_aligned
            adata.obs.loc[index, y_field] = y_aligned
        alidata.update_coordinate(observed_metadata['projectName'], adata.obs)
        set_steps_status(observed_metadata, cur_step, StepStatus.FINISH, notify=True)
    except:
        set_steps_status(observed_metadata, cur_step, StepStatus.ERROR)
        observed_metadata['exception'] = traceback.format_exc()
    finally:
        time.sleep(1)
        observed_metadata['running'] = False
        project = observed_metadata['projectName']
        alidata.update_project_info(project, dict(observed_metadata))
        alidata.remove_alignment_task(project)
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

def validate_exp_data(adata):
    """
    验证数据是否正确
    """
    field = COUNTS_KEY
    X_backend = adata.X
    flag = True
    if flag and adata.raw is not None and adata.raw.X is not None:
        try:
            adata.X = adata.raw.X.copy()
            flag = False
        except Exception as e:
            flag = True
    elif flag and 'counts' in adata.layers:
        try:
            adata.X = adata.layers["counts"].copy()
            flag = False
        except Exception as e:
            flag = True
    elif flag:
        adata.X = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4) 
    min_exp = np.min(adata.X)
    adata.obsm[field] = adata.X
    adata.X = X_backend
    if min_exp>=0:
        return field, True
    return None, False

def ensure_numeric_fields(adata, x_field, y_field, z_field):
    """
    确保字段为数值类型
    """
    if not pd.api.types.is_numeric_dtype(adata.obs[x_field]):
        adata.obs[x_field] = pd.to_numeric(adata.obs[x_field], errors='raise')
    if not pd.api.types.is_numeric_dtype(adata.obs[y_field]):
        adata.obs[y_field] = pd.to_numeric(adata.obs[y_field], errors='raise')
    if not pd.api.types.is_numeric_dtype(adata.obs[z_field]):
        adata.obs[z_field] = pd.to_numeric(adata.obs[z_field], errors='raise')