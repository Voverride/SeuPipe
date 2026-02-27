from dash import set_props, Patch
import scanpy as sc
from api.annotation import *
from utils.typing import StepStatus
import threading
import shutil
from controller.notice import set_aside_notice
from dataManager.annotation_d import annData
from controller.notice import set_head_notice

def export_annotation_file(path):
    """
    导出注释文件
    """
    project = annData.get_export_project()
    if project is None:
        set_head_notice('Timeout export, please click export button again!', type='warning')
        return False
    querydata_path = annData.get_querydata_path(project)
    shutil.copy(querydata_path, path)
    set_head_notice('export successfully!', type='success')
    return True

def get_slicer_result(project_name, z_min, z_max, spot_size, border_width, border_color):
    """
    获取切片结果图
    """
    query_path = annData.get_querydata_path(project_name)
    querydata = sc.read_h5ad(query_path, backed='r')
    querydf = querydata.obs.copy()
    querydata.file.close()
    metadata = annData.get_project_info(project_name)
    x = metadata.get('x', None)
    y = metadata.get('y', None)
    z = metadata.get('z', None)
    label = metadata.get('label', None)
    querydf = querydf[(querydf[z] >= z_min) & (querydf[z] <= z_max)]
    fig = plot_3d_scatter(querydf, x, y, z, label, spot_size, border_width, border_color)
    return fig

def show_bug_info(project):
    """
    显示聚类异常信息
    """
    project_info = annData.get_project_info(project)
    exception = project_info.get('exception', None)
    if exception is not None:
        set_aside_notice('Runtime Error', exception, 'error')

def start_ann_project(projectname):
    """
    开始注释项目
    """
    observed_metadata = annData.create_annotation_task(projectname)
    thread = threading.Thread(target=scvi_annotation, args=(observed_metadata,))
    thread.start()

def update_project_metadata(projectname, spot_size, border_width, border_color):
    """
    更新项目元数据
    """
    patch = Patch()
    if not projectname:
        set_props('ann-bug-panel', dict(style={'display':'none'}))
        set_props('ann-start-project', dict(disabled=False))
        set_props('ann-delete-project', dict(disabled=False))
        patch[0]['status'] = StepStatus.WAIT
        patch[1]['status'] = StepStatus.WAIT
        patch[2]['status'] = StepStatus.WAIT
        set_props('ann-annotation-steps', dict(steps=patch))
        set_props('ann-table-metadata', dict(data=[{}]))
        set_props('ann-table-refquery', dict(data=[{}]))
        set_props('ann-graph-result', dict(figure=None))
        set_props('ann-graph-heatmap', dict(figure=None))
        return
    metadata = annData.get_project_info(projectname)
    exception = metadata.get('exception', None)
    if exception:
        set_props('ann-bug-panel', dict(style=None))
    else:
        set_props('ann-bug-panel', dict(style={'display':'none'}))
    metadata['rm_mt'] = str(metadata['rm_mt'])
    metadata['rm_ribo'] = str(metadata['rm_ribo'])
    metadata['rm_hb'] = str(metadata['rm_hb'])
    metadata['use_hvg'] = str(metadata['use_hvg'])
    steps = metadata.get('steps', {})
    preprocess = steps.get('preprocess', StepStatus.WAIT)
    training = steps.get('training', StepStatus.WAIT)
    postprocess = steps.get('postprocess', StepStatus.WAIT)
    percent = steps.get('percent', 0)
    patch[0]['status'] = preprocess
    patch[1]['status'] = training
    patch[2]['status'] = postprocess
    table_metadata = [metadata]
    activeKey = 'ProjectInfo'
    heatmapFig = None
    resultFig = None
    if annData.is_running(projectname):
        set_props('ann-start-project', dict(disabled=True))
        set_props('ann-delete-project', dict(disabled=True))
    else:
        diff_gene_fig = annData.get_diffgene_fig(projectname)
        result_fig = annData.get_annotation_fig(projectname)
        if diff_gene_fig is not None:
            activeKey = 'Results'
            heatmapFig = diff_gene_fig
        if result_fig is not None:
            activeKey = 'Results'
            resultFig = result_fig
            classes = annData.get_classes(projectname)
            if classes is not None:
                for i in range(classes):
                    resultFig['data'][i]['marker']['size'] = spot_size
                    resultFig['data'][i]['marker']['line']['width'] = border_width
                    resultFig['data'][i]['marker']['line']['color'] = border_color
        set_props('ann-start-project', dict(disabled=False))
        set_props('ann-delete-project', dict(disabled=False))
    z_min, z_max = annData.get_z_range(projectname)
    if z_min is not None and z_max is not None:
        set_props('ann-slider-slicer', dict(min=z_min, max=z_max, value=[z_min, z_max]))
    set_props('ann-content-tabs', dict(activeKey=activeKey))
    set_props('ann-graph-heatmap', dict(figure=heatmapFig))
    set_props('ann-graph-result', dict(figure=resultFig))
    set_props('ann-annotation-steps', dict(steps=patch, percent=percent))
    set_props('ann-table-metadata', dict(data=table_metadata))
    set_props('ann-table-refquery', dict(data=table_metadata))

def create_project(refdata, querydata, label, x, y, z, rm_mt, rm_ribo, rm_hb, use_hvg, n_layers, n_hiddens, n_latent, epochs, batch_size, dropout, projectname):
    """
    创建项目
    """
    exist_projects = annData.get_exist_projects()
    if projectname in exist_projects:
        set_head_notice(f'Project Name {projectname} Already Exists, Please Input Another Name', type='warning')
        return
    annData.create_project(refdata, querydata, label, x, y, z, rm_mt, rm_ribo, rm_hb, use_hvg, n_layers, n_hiddens, n_latent, epochs, batch_size, dropout, projectname)
    set_head_notice('Created Successfully', type='success')
    set_props('ann-select-project', dict(value=projectname))
    set_props('ann-modal-newproject', dict(visible=False))

def read_annotask_querydata(path:str)->bool:
    """
    读取查询数据
    """
    adata = sc.read_h5ad(path, backed='r')
    obs_fields = adata.obs.columns.tolist()
    obs_fields.sort()
    adata.file.close()
    set_props('ann-store-querydata', dict(data=path))
    set_props('ann-select-x', {'options':obs_fields, 'value':None})
    set_props('ann-select-y', {'options':obs_fields, 'value':None})
    set_props('ann-select-z', {'options':obs_fields, 'value':None})
    set_head_notice('import successfully!', 'success')
    return True

def read_annotask_refdata(path:str)->bool:
    """
    读取参考数据
    """
    adata = sc.read_h5ad(path, backed='r')
    obs_fields = adata.obs.columns.tolist()
    obs_fields.sort()
    adata.file.close()
    set_props('ann-store-refdata', dict(data=path))
    set_props('ann-select-label', {'options':obs_fields, 'value':None})
    set_head_notice('import successfully!', 'success')
    return True

def open_new_project_modal():
    """
    显示新建项目面板
    """
    set_props('ann-store-refdata', dict(data=None))
    set_props('ann-store-querydata', dict(data=None))
    set_props('ann-modal-newproject', dict(visible=True))
    set_props('ann-select-label', dict(value=None, options=[]))
    set_props('ann-select-x', dict(value=None, options=[]))
    set_props('ann-select-y', dict(value=None, options=[]))
    set_props('ann-select-z', dict(value=None, options=[]))
    set_props('ann-projectname', dict(value=None, status='error'))