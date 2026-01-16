from dataManager.workspace import *
from api.cellselector import *
from dataManager.maskviewer_d import maskData
import threading
from utils.colors import spot_colors
from dash import set_props
from utils.observer import observer
from websocket.websocket import ws
from controller.notice import *

def show_bug_info(project):
    """
    显示聚类异常信息
    """
    project_info = selData.get_project_info(project)
    exception = project_info.get('exception', None)
    if exception is not None:
        set_aside_notice('Runtime Error', exception, 'error')

def update_cluster_status(project):
    """
    更新细胞筛选项目的聚类状态
    """
    project_info = selData.get_project_info(project)
    exception = project_info.get('exception', None)
    if selData.is_clustering(project) or exception:
        if exception:
            set_props('sel-cluster-showBug', dict(style={'backgroundColor':'#bb5548'}))
        else:
            set_props('sel-cluster-showBug', dict(style={'display':'none'}))
        creator = project_info.get('creator', None)
        set_props('sel-cluster-creator', dict(children=creator))
        steps = project_info.get('steps', {})
        set_props('sel-modal-clusteringStatus', dict(visible=True))
        for step, status in steps.items():
            if step.startswith('percent'):
                set_props(f'sel-cluster-{step}', dict(percent=status)) 
            else:
                set_props(f'sel-cluster-{step}', dict(icon=status))        
        return True
    return False

def create_project(project, result, resolution, iteration):
    """
    创建细胞筛选项目
    """
    project_name, metadata = selData.create_project(project, result, resolution, iteration)
    if resolution is not None:
        observered_metadata = selData.add_clustering(project_name, metadata)
        threading.Thread(target=run_cluster, args=(project_name, observered_metadata)).start()
    return project_name

def export_cellselector_file(path):
    """
    导出筛选后的数据  
    """
    export_project = selData.get_export_project()
    project_info = selData.get_project_info(export_project)
    project_name = project_info.get('project_name')
    mask_label = project_info.get('result_field')
    projectInfo = segData.get_project_info(project_name)
    adatas = []
    for slice in projectInfo['slices']:
        z = slice['z']
        bin_path = slice['gem']
        position = selData.get_cell_position(export_project, z)
        if position is None:
            set_head_notice(f'Slice {z} segmentation is in progress. Please wait...', type='warning')
            return False
        selected_index = position[position['selected']].index
        if selected_index.empty:
            continue
        result_path = segData.get_result_path(project_name, z)
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

def update_image_style(fig, show_image, show_scatter, clusters, scatter_color='#defcff', scatter_size=10, line_color='black', line_width=1):
    """
    是否显示图像和散点图, 并更新散点图样式
    """
    fig['data'][0]['visible'] = show_image
    for i in range(1, clusters+1):
        fig['data'][i]['visible'] = show_scatter
    if show_scatter:
        for i in range(1, clusters+1):
            if clusters == 1:
                fig['data'][i]['marker']['color'] = scatter_color
            fig['data'][i]['marker']['size'] = scatter_size
            fig['data'][i]['marker']['line']['color'] = line_color
            fig['data'][i]['marker']['line']['width'] = line_width
    if show_scatter and not show_image:
        fig['layout']['plot_bgcolor'] = 'rgba(0, 0, 0, 1)'
    else:
        fig['layout']['plot_bgcolor'] = 'rgba(0, 0, 0, 0)'
    return fig

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

    seg_project = project.rsplit('_', 1)[0]
    fig = selData.get_stain_fig(seg_project, slice)
    position = position[position['selected']]
    clusters = position['cluster'].unique()
    for i, cluster in enumerate(clusters):
        cluster_data = position[position['cluster'] == cluster]
        color = spot_colors[i % len(spot_colors)]
        fig.add_scatter(
            x=cluster_data['x'],
            y=cluster_data['y'],
            mode='markers',
            customdata=cluster_data.index,
            marker=dict(
                color=color,
                opacity=1,
                size=10,
                line_color='black',
                line_width=1
            ),
            hovertemplate='Cell: %{customdata}<br>Cluster: ' + str(cluster) + '<extra></extra>',
        )

    fig.update_layout(
        coloraxis_showscale=False,
        xaxis_visible=False,
        yaxis_visible=False,
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        legend=dict(visible=False),
        margin=dict(l=0, r=0, t=0, b=0)
    )

    return fig, len(clusters)

def get_result_options(project: str):
    """
    获取分割结果选项
    """
    result = None
    for slice in maskData.get_project_slices(project):
        curTypes = maskData.get_figure_types(project, slice)
        curTypes = set(curTypes)
        if result is None:
            result = curTypes
        else:
            result = result.intersection(curTypes)
    sel_projects = set(selData.get_exist_projects())
    result.discard('overlay')
    result = [type for type in result if f'{project}_{type}' not in sel_projects]
    result.sort()
    return result
        
def get_seg_projects():
    """
    获取分割项目列表
    """
    segProjects = segData.get_exist_projects()
    return segProjects
