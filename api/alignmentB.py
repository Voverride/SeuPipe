import paste as pst
from paste2 import PASTE2, projection
import pandas as pd
from anndata import AnnData
import numpy as np
from typing import Tuple
from dataManager.alignmentB_d import alidata
import scanpy as sc
import torch
import ot

COUNTS_KEY = 'SeuPipe_Alignment_Counts'
def paste1(adata:AnnData, x_field:str, y_field:str, z_field:str, use_gpu:bool=False, use_rep:str=None, alistatus:dict=None)->bool:
    """
    paste1算法
    """
    steps = None
    if alistatus:
        steps = alistatus['steps']
    try:
        set_steps_status(steps, 1, 'running', True)
        ensure_numeric_fields(adata, x_field, y_field, z_field)
        field, status = validate_exp_data(adata)
    except Exception as e:
        set_steps_status(steps, 1, 'failed', True)
        if alistatus:
            alistatus['exception'] = e
            return False
        else:
            raise e
        
    set_steps_status(steps, 1, 'complete', True)

    if use_rep==None and field!=None:
        use_rep = field
    try:
        set_steps_status(steps, 2, 'running', True)
        adata.obsm['spatial'] = adata.obs[[x_field, y_field]].values
        z = list(adata.obs[z_field].unique())
        z.sort(key=lambda num:int(num))
        slices = [adata[adata.obs[z_field]==slice] for slice in z]
    except Exception as e:
        set_steps_status(steps, 2, 'failed', True)
        if alistatus:
            alistatus['exception'] = e
            return False
        else:
            raise e
    set_steps_status(steps, 2, 'complete', True)
    backend = ot.backend.TorchBackend() 
    pis = []
    try:
        set_steps_status(steps, 3, 'running', True)
        for i in range(len(slices)-1):
            pis.append(pst.pairwise_align(slices[i], slices[i+1], use_gpu = use_gpu, backend=backend, use_rep=use_rep, gpu_verbose=False))
            set_steps_status(steps, 3, 'percent', int((i+2)/len(slices)*100))
    except Exception as e:
        set_steps_status(steps, 3, 'failed', True)
        if alistatus:
            alistatus['exception'] = e
            return False
        else:
            raise e
    set_steps_status(steps, 3, 'complete', True)
    try:
        set_steps_status(steps, 4, 'running', True)
        new_slices = pst.stack_slices_pairwise(slices, pis)
    except Exception as e:
        set_steps_status(steps, 4, 'failed', True)
        if alistatus:
            alistatus['exception'] = e
            return False
        else:
            raise e
    set_steps_status(steps, 4, 'complete', True)
    try:
        set_steps_status(steps, 5, 'running', True)
        for slice in new_slices:
            index = slice.obs_names
            x_aligned = slice.obsm['spatial'][:, 0]
            y_aligned = slice.obsm['spatial'][:, 1]
            adata.obs.loc[index, 'x_aligned'] = x_aligned
            adata.obs.loc[index, 'y_aligned'] = y_aligned
        if field!=None:
            del adata.obsm[field]
        adata.obsm['X_spatial_registered'] = np.array(adata.obs[['x_aligned', 'y_aligned', z_field]], dtype=np.float64)
    except Exception as e:
        set_steps_status(steps, 5, 'failed', True)
        if alistatus:
            alistatus['exception'] = e
            return False
        else:
            raise e
    set_steps_status(steps, 5, 'complete', True)
    alidata.set_x_aligned_field('x_aligned')
    alidata.set_y_aligned_field('y_aligned')
    alidata.plot_3d_scatter()
    if COUNTS_KEY in adata.obsm:
        del adata.obsm[COUNTS_KEY]
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return True

def paste2(adata:AnnData, x_field:str, y_field:str, z_field:str, use_rep:str=None, alistatus:dict=None)->bool:
    """
    paste2算法
    """
    steps = None
    if alistatus:
        steps = alistatus['steps']
    try:
        set_steps_status(steps, 1, 'running', True)
        ensure_numeric_fields(adata, x_field, y_field, z_field)
    except Exception as e:
        set_steps_status(steps, 1, 'failed', True)
        if alistatus:
            alistatus['exception'] = e
            return False
        else:
            raise e
    field, status = validate_exp_data(adata)
    if not status:
        set_steps_status(steps, 1, 'experror', True)
        set_steps_status(steps, 1, 'failed', True)
        return False
    set_steps_status(steps, 1, 'complete', True)
    if use_rep==None and field!=None:
        use_rep = field
    try:
        set_steps_status(steps, 2, 'running', True)
        adata.obsm['spatial'] = adata.obs[[x_field, y_field]].values
        z = list(adata.obs[z_field].unique())
        z.sort(key=lambda num:int(num))
        slices = [adata[adata.obs[z_field]==slice] for slice in z]
    except Exception as e:
        set_steps_status(steps, 2, 'failed', True)
        if alistatus:
            alistatus['exception'] = e
            return False
        else:
            raise e
    set_steps_status(steps, 2, 'complete', True)
    pis = []
    try:
        set_steps_status(steps, 3, 'running', True)
        for i in range(len(slices)-1):
            pis.append(PASTE2.partial_pairwise_align(slices[i], slices[i+1], s=0.7, use_rep=use_rep, verbose=False))
            set_steps_status(steps, 3, 'percent', int((i+2)/len(slices)*100))
    except Exception as e:
        set_steps_status(steps, 3, 'failed', True)
        if alistatus:
            alistatus['exception'] = e
            return False
        else:
            raise e
    set_steps_status(steps, 3, 'complete', True)
    try:
        set_steps_status(steps, 4, 'running', True)
        new_slices = projection.partial_stack_slices_pairwise(slices, pis)
    except Exception as e:
        set_steps_status(steps, 4, 'failed', True)
        if alistatus:
            alistatus['exception'] = e
            return False
        else:
            raise e
    set_steps_status(steps, 4, 'complete', True)
    try:
        set_steps_status(steps, 5, 'running', True)
        for slice in new_slices:
            index = slice.obs_names
            x_aligned = slice.obsm['spatial'][:, 0]
            y_aligned = slice.obsm['spatial'][:, 1]
            adata.obs.loc[index, 'x_aligned'] = x_aligned
            adata.obs.loc[index, 'y_aligned'] = y_aligned
        if field!=None:
            del adata.obsm[field]
        adata.obsm['X_spatial_registered'] = np.array(adata.obs[['x_aligned', 'y_aligned', z_field]], dtype=np.float64)
    except Exception as e:
        set_steps_status(steps, 5, 'failed', True)
        if alistatus:
            alistatus['exception'] = e
            return False
        else:
            raise e
    set_steps_status(steps, 5, 'complete', True)
    alidata.set_x_aligned_field('x_aligned')
    alidata.set_y_aligned_field('y_aligned')
    alidata.plot_3d_scatter()
    if COUNTS_KEY in adata.obsm:
        del adata.obsm[COUNTS_KEY]
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return True

def validate_exp_data(adata:AnnData)->Tuple[str, bool]:
    """
    Validate the experimental data to ensure that all values are non-negative.

    Args:
        adata (AnnData): The annotated data object.

    Returns:
        tuple: A tuple containing:
            - field (str): The name of the field ('counts_for_paste') if data is valid.
            - bool: True if the data is valid (non-negative), False otherwise.
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
    elif  flag:
        adata.X = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4) 
    min_exp = np.min(adata.X)
    adata.obsm[field] = adata.X
    adata.X = X_backend
    if min_exp>=0:
        return field, True
    return None, False

def ensure_numeric_fields(adata:AnnData, x_field:str, y_field:str, z_field:str):
    """
    Ensure that the specified fields (x_field, y_field, z_field) in the adata object are numeric.

    Args:
        adata (AnnData): The annotated data object.
        x_field (str): The name of the field to check for numeric values in the x-axis.
        y_field (str): The name of the field to check for numeric values in the y-axis.
        z_field (str): The name of the field to check for numeric values in the z-axis.
    """
    if not pd.api.types.is_numeric_dtype(adata.obs[x_field]):
        adata.obs[x_field] = pd.to_numeric(adata.obs[x_field], errors='raise')
    if not pd.api.types.is_numeric_dtype(adata.obs[y_field]):
        adata.obs[y_field] = pd.to_numeric(adata.obs[y_field], errors='raise')
    if not pd.api.types.is_numeric_dtype(adata.obs[z_field]):
        adata.obs[z_field] = pd.to_numeric(adata.obs[z_field], errors='raise')
def set_steps_status(steps:dict, step:int, field:str, status:bool)->None:
    if steps is None:
        return
    if field=='complete' or field=='failed':
        steps[step]['running'] = False
    steps[step][field] = status