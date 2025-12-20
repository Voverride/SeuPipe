from dash import set_props
from dataManager.workspace import *
from dataManager.maskviewer_d import maskData
from dataManager.segmentation_d import segData
from controller.notice import set_head_notice
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix
import numpy as np
import pandas as pd

graph_visible_style = {'width': '100%', 'height':'90%', 'margin':'auto'}


def generate_layer_adata(z, df, mask):
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
    x_centered = x_abs - mean_x
    y_centered = y_abs - mean_y

    adata = ad.AnnData(
        X = X,
        obs = pd.DataFrame(
            {
                'x': -y_centered,
                'y': x_centered,
                'z': z,
                'image_x': x_abs,
                'image_y': y_abs
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
    export_task = maskData.get_export_task()
    taskname = export_task.get('taskname')
    type = export_task.get('type')
    taskInfo = segData.read_taskinfo(taskname)
    adatas = []
    for slice in taskInfo['data']:
        z = slice['z']
        bin_path = slice['gem']
        result_path = segData.get_seg_adata_path(taskname, z)
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
    merged_adata.write(path)
    set_head_notice('Export successfully!', 'success')
    return True

def update_graph_with_type(taskname, slice, graph, showMask, showContour, leftGraph, rightGraph):
    """
    graphtype改变，更改图像，配准，分割或扩展
    """
    if graph=='registration':
        set_props('mv-checkbox-mask', dict(disabled=True))
        set_props('mv-checkbox-contour', dict(disabled=True))
        set_props('maskviewer-left-graph', dict(disabled=True))
        set_props('maskviewer-right-graph', dict(disabled=True))
    if slice:
        if graph=='registration':
            before, aligned = maskData.get_registration_figure(taskname, slice)
            set_props('mv-graph-left', dict(figure=before, style=graph_visible_style))
            set_props('mv-graph-right', dict(figure=aligned, style=graph_visible_style))
        else:
            set_props('maskviewer-left-graph', dict(disabled=False))
            set_props('maskviewer-right-graph', dict(disabled=False))
            set_props('mv-checkbox-mask', dict(disabled=False))
            set_props('mv-checkbox-contour', dict(disabled=False))
            left_graph = maskData.get_graph(taskname, slice, leftGraph, showMask=showMask, showContour=showContour)
            right_graph = maskData.get_graph(taskname, slice, rightGraph, showMask=showMask, showContour=showContour)
            set_props('mv-graph-left', dict(figure=left_graph, style=graph_visible_style))
            set_props('mv-graph-right', dict(figure=right_graph, style=graph_visible_style))

def update_graph_with_slice(taskname, slice, graphType, showMask, showContour, leftGraph, rightGraph):
    """
    slice改变， 更改图像，配准，分割或扩展
    """
    if taskname and slice:
        if graphType=='registration':
            before, aligned = maskData.get_registration_figure(taskname, slice)
            set_props('mv-graph-left', dict(figure=before, style=graph_visible_style))
            set_props('mv-graph-right', dict(figure=aligned, style=graph_visible_style))
        else:
            left_graph = maskData.get_graph(taskname, slice, leftGraph, showMask=showMask, showContour=showContour)
            right_graph = maskData.get_graph(taskname, slice, rightGraph, showMask=showMask, showContour=showContour)
            set_props('mv-graph-left', dict(figure=left_graph, style=graph_visible_style))
            set_props('mv-graph-right', dict(figure=right_graph, style=graph_visible_style))
            set_props('maskviewer-left-graph', dict(disabled=False))
            set_props('maskviewer-right-graph', dict(disabled=False))
            set_props('mv-checkbox-mask', dict(disabled=False))
            set_props('mv-checkbox-contour', dict(disabled=False))
    
def update_taskname(taskname):
    """
    task改变，更新store， 切片列表
    """
    set_props('maskviewer-store-taskname', dict(data=taskname))
    set_props('maskviewer-select-slice', dict(value=None))
def restore_initial_data(lastTaskName):
    """
    恢复网页初始数据
    """
    if lastTaskName:
        set_props('maskviewer-select-taskname', dict(value=lastTaskName))