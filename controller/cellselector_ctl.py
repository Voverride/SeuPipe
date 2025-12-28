from dataManager.workspace import *
from dataManager.cellselector_d import selData
from dataManager.segmentation_d import segData
from utils.observer import observer
from websocket.websocket import ws
import plotly.express as px
import pandas as pd
import scanpy as sc
import anndata as ad
from controller.maskviewer_ctl import generate_layer_adata
from controller.notice import *


def export_cellselector_file(path):
    """
    导出筛选后的数据
    """
    export_project = selData.get_export_project()
    project_info = selData.get_project_info(export_project)
    project_name = project_info.get('project_name')
    mask_label = project_info.get('result_field')
    taskInfo = segData.read_taskinfo(project_name)
    adatas = []
    for slice in taskInfo['data']:
        z = slice['z']
        bin_path = slice['gem']
        position = selData.get_cell_position(export_project, z)
        if position is None:
            set_head_notice(f'Slice {z} segmentation is in progress. Please wait...', type='warning')
            return False
        selected_index = position[position['selected']].index
        if selected_index.empty:
            continue
        result_path = segData.get_seg_adata_path(project_name, z)
        segResult = sc.read_h5ad(result_path)
        mask_matrix = segResult.layers[mask_label]
        df = pd.read_csv(bin_path, sep='\t', comment='#')
        adata = generate_layer_adata(z, df, mask_matrix)
        common_index = selected_index.intersection(adata.obs.index)
        if common_index.empty:
            continue
        adata_sub = adata[common_index, :]
        adatas.append(adata_sub)
    if len(adatas)==0:
        set_head_notice('No cell selected', 'warning')
        return False
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
def select_cell(project, slice, pointIndex: list, operation):
    """
    筛选细胞
    """
    isWritting = selData.is_slice_writting(project, slice)
    if isWritting:
        set_head_notice('Slice is busy. Please wait a moment', type='warning')
        return
    selData.set_slice_writting(project, slice)
    signal = observer.observe({'complete':False}, ws.notifyUpdateCellSelectorStatus, project=project, slice=slice)
    position = selData.get_cell_position(project, slice)
    if operation == 'retain':
        position['selected'] = position.index.isin(pointIndex)
    elif operation == 'remove':
        position.loc[pointIndex, 'selected'] = False
    else:
        position['selected'] = True
    signal['complete'] = True
    observer.disobserve(signal, ws.notifyUpdateCellSelectorStatus)
    del signal
    position_path = selData.get_position_path(project, slice)
    position.to_csv(position_path, sep='\t')
    selData.remove_slice_writting(project, slice)
def get_slice_graph(project, slice):
    """
    更新切片图像
    """
    position = selData.get_cell_position(project, slice)
    if position is None:
        seg_mask = selData.get_seg_mask(project, slice)
        if seg_mask is None:
            set_head_notice(f'Slice {slice} segmentation is in progress. Please wait...', type='warning')
            return None
        position = selData.generate_cell_position(project, slice, seg_mask)

    img = selData.get_stain_image(project, slice)
    fig = px.imshow(img, color_continuous_scale='gray')
    position = position[position['selected']]
    fig.update_layout(
        coloraxis_showscale=False,
        xaxis_visible=False,
        yaxis_visible=False,
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        margin=dict(l=0, r=0, t=0, b=0)
    )
    fig.add_scatter(
        x=position['x'],
        y=position['y'],
        mode='markers',
        customdata=position.index,
        marker=dict(color='#f6b894', opacity=1, size=7),
        hoverinfo='skip'
    )
    return fig

def get_result_options(project: str):
    """
    获取分割结果选项
    """
    slices_folder = segData.get_seg_slice_folder(project)
    slices = [slice for slice in os.listdir(slices_folder) if os.path.isdir(os.path.join(slices_folder, slice))]
    sel_projects = set(selData.get_exist_projects())
    result = set()
    for slice in slices:
        slice_index = slice.split('_')[1]
        figure_folder = segData.get_seg_figure_folder(project, slice_index)
        figures = [fig for fig in os.listdir(figure_folder) if 'mask' in fig]
        for figure in figures:
            result.add(figure.split('_mask')[0])
    result = [type for type in result if f'{project}_{type}' not in sel_projects]
    result.sort()
    return result
        
def get_seg_projects():
    """
    获取分割项目列表
    """
    segProjects = segData.get_exist_tasks()
    return segProjects
