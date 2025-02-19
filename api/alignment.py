import paste as pst
from paste2 import PASTE2, projection
import pandas as pd
from anndata import AnnData
import numpy as np
from typing import Tuple
import ot

def paste1(adata:AnnData, x_field:str, y_field:str, z_field:str, use_gpu:bool=False, use_rep:str=None, steps:dict=None)->bool:
    """
    Align spatial slices of data using the PASTE1 algorithm. This function takes in a single AnnData object,
    performs alignment of slices based on spatial coordinates, and updates the AnnData object with the aligned
    spatial coordinates.

    Args:
        adata (AnnData): The annotated data object containing the spatial coordinates and other experimental data.
        x_field (str): The field in `adata.obs` that represents the x-coordinate of the spatial data.
        y_field (str): The field in `adata.obs` that represents the y-coordinate of the spatial data.
        z_field (str): The field in `adata.obs` that represents the slice or z-coordinate for each slice of data.
        use_gpu (bool, optional): Whether to use GPU for computation in alignment. Defaults to False.
        use_rep (str, optional): The representation field for alignment, specifying how the data should be transformed for matching. Defaults to None.
        steps (dict, optional): An object to track and log progress status through the various stages of alignment. If provided, status updates are logged.

    Returns:
        bool: True if alignment was successful, False otherwise.

    Raises:
        Exception: If any of the steps in the process fail, an exception is raised and progress is logged.
    """
    try:
        ensure_numeric_fields(adata, x_field, y_field, z_field)
    except Exception as e:
        set_steps_status(steps, 1, 'failed', True)
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
        adata.obsm['spatial'] = adata.obs[[x_field, y_field]].values
        z = list(adata.obs[z_field].unique())
        z.sort(key=lambda num:int(num))
        slices = [adata[adata.obs[z_field]==slice] for slice in z]
    except Exception as e:
        set_steps_status(steps, 2, 'failed', True)
        raise e
    set_steps_status(steps, 2, 'complete', True)
    backend = ot.backend.TorchBackend() 
    pis = []
    try:
        for i in range(len(slices)-1):
            pis.append(pst.pairwise_align(slices[i], slices[i+1], use_gpu = use_gpu, backend=backend, use_rep=use_rep, gpu_verbose=False))
            set_steps_status(steps, 3, 'percent', int((i+2)/len(slices)*100))
    except Exception as e:
        set_steps_status(steps, 3, 'failed', True)
        raise e
    set_steps_status(steps, 3, 'complete', True)
    try:
        new_slices = pst.stack_slices_pairwise(slices, pis)
    except Exception as e:
        set_steps_status(steps, 4, 'failed', True)
        raise e
    set_steps_status(steps, 4, 'complete', True)
    try:
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
        raise e
    set_steps_status(steps, 5, 'complete', True)
    return True

def paste2(adata:AnnData, x_field:str, y_field:str, z_field:str, use_rep:str=None, steps:dict=None)->bool:
    """
    Align slices using the paste2 algorithm, which is a method for registering slices based on spatial coordinates
    and a specified representation field. It aligns slices in 2D space and stacks them based on pairwise alignments.

    Args:
        adata (AnnData): The annotated data object containing spatial and experimental data.
        x_field (str): The field in `adata.obs` that represents the x-coordinate for spatial data.
        y_field (str): The field in `adata.obs` that represents the y-coordinate for spatial data.
        z_field (str): The field in `adata.obs` that corresponds to the z-coordinate or slice identifier.
        use_rep (str, optional): The representation field for alignment. If not provided, defaults to None. If provided, it will be used to compute the pairwise alignment between slices.
        steps (str, optional): A tracking object to report progress and status at each step. Defaults to None.

    Returns:
        bool: True if alignment was successful, False otherwise.

    Steps:
        1. Ensure that the fields for x, y, and z coordinates are numeric.
        2. Validate the experimental data in `adata` for correctness.
        3. Create spatial coordinates from the provided x and y fields.
        4. Sort slices by the z-field and split the data into individual slices.
        5. Perform pairwise alignment between consecutive slices.
        6. Stack the slices based on the pairwise alignments.
        7. Align the coordinates and store the results in `adata`.
    """
    try:
        ensure_numeric_fields(adata, x_field, y_field, z_field)
    except Exception as e:
        set_steps_status(steps, 1, 'failed', True)
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
        adata.obsm['spatial'] = adata.obs[[x_field, y_field]].values
        z = list(adata.obs[z_field].unique())
        z.sort(key=lambda num:int(num))
        slices = [adata[adata.obs[z_field]==slice] for slice in z]
    except Exception as e:
        set_steps_status(steps, 2, 'failed', True)
        raise e
    set_steps_status(steps, 2, 'complete', True)
    pis = []
    try:
        for i in range(len(slices)-1):
            pis.append(PASTE2.partial_pairwise_align(slices[i], slices[i+1], s=0.7, use_rep=use_rep, verbose=False))
            set_steps_status(steps, 3, 'percent', int((i+2)/len(slices)*100))
    except Exception as e:
        set_steps_status(steps, 3, 'failed', True)
        raise e
    set_steps_status(steps, 3, 'complete', True)
    try:
        new_slices = projection.partial_stack_slices_pairwise(slices, pis)
    except Exception as e:
        set_steps_status(steps, 4, 'failed', True)
        raise e
    set_steps_status(steps, 4, 'complete', True)
    try:
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
        raise e
    set_steps_status(steps, 5, 'complete', True)
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
    field = 'counts_for_paste'
    min_exp = np.min(adata.X)
    if min_exp>=0:
        return None, True
    if adata.raw!=None:
        min_exp = np.min(adata.raw.X)
        if min_exp>=0:
            adata.obsm[field] = adata.raw.X
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
    steps[step][field] = status