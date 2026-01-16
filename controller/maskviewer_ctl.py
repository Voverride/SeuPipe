from dataManager.workspace import *
from dataManager.maskviewer_d import maskData
from dataManager.segmentation_d import segData
from controller.notice import set_head_notice
import feffery_antd_components as fac
import scanpy as sc
from dash import dcc
import anndata as ad
from scipy.sparse import csr_matrix
import numpy as np
import pandas as pd

graph_visible_style = {'width': '100%', 'height':'90%', 'margin':'auto'}

contrast_graph = {
    'children': fac.AntdCenter(
        dcc.Graph(
            id="mv-graph-right", 
            config={'displaylogo':False}, 
            style={'display':'block', 'height':'90vh', 'width':'100%', 'visibility':'visible'}
        ),
    ),
    'collapsible': True,
}

def generate_layer_adata(z, df, mask, centered=True):
    """
    生成单层anndata
    """
    countField = 'MIDCount' if 'MIDCount' in df.columns else 'MIDCounts'
    max_row, max_col = mask.shape
    df['cellID'] = 0
    valid_mask = (df['x'] >= 0) & (df['x'] < max_row) & (df['y'] >= 0) & (df['y'] < max_col)
    if valid_mask.any():
        valid_x = df.loc[valid_mask, 'x'].values
        valid_y = df.loc[valid_mask, 'y'].values
        df.loc[valid_mask, 'cellID'] = np.asarray(mask[valid_x, valid_y]).flatten()
    df = df[df['cellID']!=0]
    expr_df = df.groupby(['cellID', 'geneID'], observed=True)[countField].sum().reset_index()
    
    cells = sorted(expr_df['cellID'].unique())
    genes = sorted(expr_df['geneID'].unique())
    
    cell_to_idx = {cell: i for i, cell in enumerate(cells)}
    gene_to_idx = {gene: j for j, gene in enumerate(genes)}
    
    row_indices = expr_df['cellID'].map(cell_to_idx).values
    col_indices = expr_df['geneID'].map(gene_to_idx).values
    data = expr_df[countField].values
    
    X = csr_matrix((data, (row_indices, col_indices)), shape=(len(cells), len(genes)))
    
    center_df = df.groupby('cellID', observed=True)[['x', 'y']].mean().reindex(cells)

    new_cell_ids = [f"z{z}_{cell}" for cell in cells]

    x_abs = center_df['x'].values
    y_abs = center_df['y'].values
    mean_x = x_abs.mean()
    mean_y = y_abs.mean()
    
    if centered:
        x_centered = x_abs - mean_x
        y_centered = y_abs - mean_y
    else:
        x_centered = x_abs
        y_centered = y_abs

    adata = ad.AnnData(
        X = X,
        obs = pd.DataFrame(
            {
                'x': y_centered,
                'y': x_centered,
                'z': z,
                'image_x': y_abs,
                'image_y': x_abs
            },
            index = new_cell_ids
        ),
        var = pd.DataFrame(index = genes)
    )

    adata.obs['cellID'] = new_cell_ids
    
    return adata
    

def export_maskviewer_file(path):
    """
    导出整合后的数据
    """
    export_project = maskData.get_export_project()
    project = export_project.get('project')
    type = export_project.get('type')
    projectInfo = segData.get_project_info(project)
    adatas = []
    for slice in projectInfo['slices']:
        z = slice['z']
        bin_path = slice['gem']
        result_path = segData.get_result_path(project, z)
        mask_label = f'{type}_mask'
        segResult = sc.read_h5ad(result_path)
        mask_matrix = segResult.layers[mask_label]
        df = pd.read_csv(bin_path, sep='\t', comment='#')
        adata = generate_layer_adata(z, df, mask_matrix)
        adatas.append(adata)
    merged_adata = ad.concat(
        adatas,
        join='outer',
        fill_value=0,
        index_unique=None,
        merge='first'
    )
    merged_adata.layers['counts'] = merged_adata.X.copy()
    sc.pp.normalize_total(merged_adata, target_sum=1e4)
    sc.pp.log1p(merged_adata)
    merged_adata.write_h5ad(path)
    set_head_notice('Export successfully!', 'success')
    return True