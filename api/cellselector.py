import numpy as np
import scanpy as sc
import pandas as pd
import anndata as ad
import time
import os
import traceback
from utils.typing import StepIcon
from utils.commonfuc import write_json
from dataManager.cellselector_d import selData
from dataManager.segmentation_d import segData
from controller.maskviewer_ctl import generate_layer_adata

def run_cluster(project_name, observed_metadata):
    """
    运行聚类分析
    """
    try:
      set_step_status(observed_metadata, 'step1', StepIcon.PROCESSING)
      observed_metadata['running'] = True
      project_info = selData.get_project_info(project_name)
      seg_project_name = project_info.get('project_name')
      resolution = project_info.get('resolution', 0.8)
      iteration = project_info.get('iteration', 0)
      if iteration == 0:
        iteration = -1
      mask_label = project_info.get('result_field')
      projectInfo = segData.get_project_info(seg_project_name)
      set_step_status(observed_metadata, 'step1', StepIcon.SUCCESS)
      set_step_status(observed_metadata, 'step2', StepIcon.PROCESSING, notify=True)
      adatas = []
      size = len(projectInfo['slices'])
      cnt = 0
      for slice in projectInfo['slices']:
          cnt += 1
          z = slice['z']
          bin_path = slice['gem']
          result_path = segData.get_result_path(seg_project_name, z)
          segResult = sc.read_h5ad(result_path)
          mask_matrix = segResult.layers[mask_label]
          df = pd.read_csv(bin_path, sep='\t', comment='#')
          adata = generate_layer_adata(z, df, mask_matrix, centered=False)
          adatas.append(adata)
          set_step_status(observed_metadata, 'percent', int(cnt/size*100), notify=True)
      merged_adata = ad.concat(
          adatas,
          join='outer',
          fill_value=0,
          index_unique=None,
          merge='first'
      )
      del adatas
      set_step_status(observed_metadata, 'step2', StepIcon.SUCCESS)
      set_step_status(observed_metadata, 'step3', StepIcon.PROCESSING, notify=True)
      sc.pp.normalize_total(merged_adata, target_sum=1e4)
      sc.pp.log1p(merged_adata)
      sc.pp.highly_variable_genes(merged_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
      merged_adata_hvg = merged_adata[:, merged_adata.var.highly_variable].copy()
      del merged_adata
      sc.pp.scale(merged_adata_hvg, max_value=10)
      sc.tl.pca(merged_adata_hvg, svd_solver='arpack')
      n_pcs = auto_select_pcs(merged_adata_hvg)
      n_neighbors = adaptive_neighbors_selection(merged_adata_hvg, n_pcs)
      sc.pp.neighbors(merged_adata_hvg, n_neighbors=n_neighbors, n_pcs=n_pcs)
      set_step_status(observed_metadata, 'step3', StepIcon.SUCCESS)
      set_step_status(observed_metadata, 'step4', StepIcon.PROCESSING, notify=True)
      sc.tl.umap(merged_adata_hvg)
      sc.tl.leiden(merged_adata_hvg, resolution=resolution, key_added='cluster', flavor='igraph', n_iterations=iteration)
      df = merged_adata_hvg.obs.copy()[['x', 'y', 'z', 'cluster']]
      del merged_adata_hvg
      df['selected'] = True
      set_step_status(observed_metadata, 'step4', StepIcon.SUCCESS)
      set_step_status(observed_metadata, 'step5', StepIcon.PROCESSING, notify=True)
      slices = set(df['z'].tolist())
      for slice in slices:
          df_slice = df[df['z']==slice]
          position_path = selData.get_position_path(project_name, slice)
          df_slice.to_csv(position_path, sep='\t')
      set_step_status(observed_metadata, 'step5', StepIcon.SUCCESS, notify=True)
    except Exception as _:
        time.sleep(2)
        for i in range(1, 6):
          step_name = f'step{i}'
          if observed_metadata['steps'][step_name] == StepIcon.PROCESSING:
            set_step_status(observed_metadata, step_name, StepIcon.ERROR)
        observed_metadata['exception'] = traceback.format_exc()
    finally:
        time.sleep(2)
        observed_metadata['running'] = False
        metadata = dict(observed_metadata)
        project_folder = selData.get_project_folder(project_name)
        write_json(os.path.join(project_folder, 'metadata.json'), metadata)
        selData.remove_clustering(project_name)
        

def set_step_status(observed_metadata, step, status, notify=False):
    """
    设置步骤状态
    """
    observed_metadata['steps'][step] = status
    if notify:
        observed_metadata['running'] = True


def auto_select_pcs(adata, min_pcs=10, max_pcs=50, threshold=0.9):
    """
    自动选择主成分数
    """
    try:
        variance_ratio = adata.uns['pca']['variance_ratio']
        cumulative_variance = variance_ratio.cumsum()
        
        threshold_indices = np.where(cumulative_variance >= threshold)[0]
        
        if len(threshold_indices) > 0:
            pcs_by_threshold = threshold_indices[0] + 1
        else:
            threshold_85 = np.where(cumulative_variance >= 0.85)[0]
            if len(threshold_85) > 0:
                pcs_by_threshold = threshold_85[0] + 1
            else:
                pcs_by_threshold = int(len(cumulative_variance) * 0.8)
        
        n_components = len(variance_ratio)
        variance_change = np.diff(variance_ratio)
        
        if len(variance_change) > 0:
            elbow_point = None
            for i in range(1, min(50, len(variance_change))):
                if variance_change[i] < variance_change[0] * 0.1:
                    elbow_point = i + 1
                    break
            
            if elbow_point is None:
                elbow_point = int(n_components * 0.5)
        else:
            elbow_point = int(n_components * 0.5)
        
        combined_pcs = int((pcs_by_threshold + elbow_point) / 2)
        
        selected_pcs = max(min_pcs, min(combined_pcs, max_pcs, n_components))
        return selected_pcs
        
    except Exception as _:
        return min_pcs

def adaptive_neighbors_selection(adata, n_pcs):
    """
    根据主成分数自适应选择邻居数
    更多PCs → 更多邻居
    """
    base_neighbors = 15
    if n_pcs < 20:
        n_neighbors = base_neighbors - 5
    elif n_pcs < 30:
        n_neighbors = base_neighbors
    elif n_pcs < 40:
        n_neighbors = base_neighbors + 5
    else:
        n_neighbors = base_neighbors + 10
    
    n_cells = adata.n_obs
    if n_cells < 1000:
        n_neighbors = max(5, n_neighbors - 5)
    elif n_cells > 10000:
        n_neighbors = min(50, n_neighbors + 10)
    
    if n_neighbors % 2 == 0:
        n_neighbors += 1
    
    return n_neighbors