import scanpy as sc
from dataManager.alignment_d import alidata
from controller.notice import set_head_notice, set_aside_notice
from utils.typing import StepStatus
from api.alignment import *
from controller.auth import verify_modify_permission
import pandas as pd
import threading
import feffery_antd_components as fac
from typing import Tuple, List
import numpy as np
from dash import Patch, set_props
import os

def update_reference_slice(project, slice, actSlice, color, icon):
    """
    更新参考切片显示配置
    """
    patch = Patch()
    znums = alidata.get_znums(project)
    for i, znum in enumerate(znums):
        if znum == slice:
            patch['data'][i]['marker']['color'] = color
            patch['data'][i]['visible'] = True
        elif znum != actSlice:
            patch['data'][i]['visible'] = False
    set_props('ali-graph-left', dict(figure=patch))
    if icon=='antd-minus':
        set_props('ali-graph-right', dict(figure=patch))

def update_active_slice(project, slice, refSlice, color, icon):
    """
    更新当前切片显示配置
    """
    patch = Patch()
    znums = alidata.get_znums(project)
    for i, znum in enumerate(znums):
        if znum == slice:
            patch['data'][i]['marker']['color'] = color
            patch['data'][i]['visible'] = True
            if i > 0:
                set_props('ma-selecter-referenceSlice', dict(value=znums[i-1]))
        elif znum != refSlice:
            patch['data'][i]['visible'] = False
    set_props('ali-graph-left', dict(figure=patch))
    if icon=='antd-minus':
        set_props('ali-graph-right', dict(figure=patch))

def slicer_layers(zmin, zmax, znums, icon):
    """
    更新切片配置
    """
    patch = Patch()
    for i, znum in enumerate(znums):
        if zmin <= znum <= zmax:
            patch['data'][i]['visible'] = True
        else:
            patch['data'][i]['visible'] = False
    set_props('ali-graph-left', dict(figure=patch))
    if icon=='antd-minus':
        set_props('ali-graph-right', dict(figure=patch))


def update_figure_spot(figure, spotSize, borderWidth, borderColor):
    """
    更新散点配置
    """
    for trace in figure['data']:
        trace['marker']['size'] = spotSize
        trace['marker']['line']['width'] = borderWidth
        trace['marker']['line']['color'] = borderColor

def update_spot_properties(spotSize, borderWidth, borderColor, znum):
    """
    更新点属性
    """
    patch = Patch()
    for i in range(znum):
        patch['data'][i]['marker']['size'] = spotSize
        patch['data'][i]['marker']['line']['width'] = borderWidth
        patch['data'][i]['marker']['line']['color'] = borderColor
    set_props('ali-graph-left', dict(figure=patch))
    set_props('ali-graph-right', dict(figure=patch))

def get_contrast_figure(project):
    """
    获取对比图像
    """
    figure = alidata.get_initialfig(project)
    if figure is None:
        project_info = alidata.get_project_info(project)
        x = project_info['x']
        y = project_info['y']
        z = project_info['z']
        initial_data = alidata.get_initialdata_path(project)
        adata = sc.read_h5ad(initial_data, backed='r')
        df = adata.obs.copy()
        adata.file.close()
        figure = alidata.update_initialfig(project, df, x, y, z)
    return figure

def show_bug_info(project):
    """
    显示对齐异常信息
    """
    project_info = alidata.get_project_info(project)
    exception = project_info.get('exception', None)
    if exception is not None:
        set_aside_notice('Runtime Error', exception, 'error')

def start_alignment_project(projectname, model, device):
    """
    开始对齐项目
    """
    observed_metadata = alidata.create_alignment_task(projectname)
    thread = threading.Thread(target=run_ali_project, args=(model, device, observed_metadata))
    thread.start()


def update_project_metadata(project, spotSize, borderWidth, borderColor, userTriggered=False):
    """
    更新项目元数据
    """
    patch = Patch()
    if not project:
        set_props('ali-bug-panel', dict(style={'display':'none'}))
        set_props('ali-start-project', dict(disabled=False))
        set_props('ali-delete-project', dict(disabled=False))
        set_props('ali-table-metadata', dict(data=[{}]))
        set_props('ali-graph-left', dict(figure=None))
        set_props('ali-select-project', dict(disabled=False))
        patch[0]['status'] = StepStatus.WAIT
        patch[1]['status'] = StepStatus.WAIT
        patch[2]['status'] = StepStatus.WAIT
        set_props('ali-alignment-steps', dict(steps=patch))
        return       
    isRunning = alidata.is_running(project) 
    metadata = alidata.get_project_info(project)
    if not isRunning or userTriggered:
        table_metadata = [metadata]
        set_props('ali-table-metadata', dict(data=table_metadata))
        exception = metadata.get('exception', None)
        if exception:
            set_props('ali-bug-panel', dict(style=None))
        else:
            set_props('ali-bug-panel', dict(style={'display':'none'}))
        resultfig = alidata.get_resultfig(project)
        if resultfig is not None:
            update_figure_spot(resultfig, spotSize, borderWidth, borderColor)
            scale = resultfig['layout']['scene']['aspectratio']['x']
            scalePatch = Patch()
            scalePatch['resultScale'] = scale
            set_props('ali-store-figureScale', scalePatch)
            set_props('ali-graph-left', dict(figure=resultfig))
        zmin = metadata['zmin']
        zmax = metadata['zmax']
        znum = metadata['znum']
        znums = metadata['znums']
        zfield = metadata['z']
        props = dict()
        props['min'] = zmin
        props['max'] = zmax
        props['value'] = [zmin, zmax]
        props['tooltipPrefix'] = f'{zfield}: '
        set_props('ma-slider-z', props)
        set_props('ma-slider-znum', dict(data=znum))
        set_props('ma-slider-znums', dict(data=znums))
        set_props('ma-selecter-activeSlice', dict(options=znums))
        set_props('ma-selecter-referenceSlice', dict(options=znums))
        znumsReversed = [{'slices': z} for z in znums[::-1]]
        set_props('ma-table-SyncSlices', dict(data=znumsReversed))
    steps = metadata.get('steps', {})
    set_props('ali-progress', dict(percent=steps.get('percent', 0)))
    steps = metadata.get('steps', {})
    preprocess = steps.get('preprocess', StepStatus.WAIT)
    alignment = steps.get('alignment', StepStatus.WAIT)
    postprocess = steps.get('postprocess', StepStatus.WAIT)
    percent = steps.get('percent', 0)
    patch[0]['status'] = preprocess
    patch[1]['status'] = alignment
    patch[2]['status'] = postprocess
    set_props('ali-alignment-steps', dict(steps=patch, percent=percent))
    if isRunning:
        set_props('ali-start-project', dict(disabled=True))
        set_props('ali-delete-project', dict(disabled=True))
    else:
        set_props('ali-start-project', dict(disabled=False))
        set_props('ali-delete-project', dict(disabled=False))

def export_alignment_file(path: str) -> bool:
    """
    导出对齐数据文件
    """
    try:
        # 假设已经有处理好的adata对象
        # adata.write_h5ad(path)
        set_head_notice('Export successfully!', 'success')
        return True
    except Exception as e:
        set_head_notice('Export failed!', 'error')
        return False


def read_alignment_file(path: str) -> bool:
    """
    读取对齐数据文件
    """
    adata = sc.read_h5ad(path, backed='r')
    obs_fields = list(adata.obs.columns)
    obs_fields.sort()
    adata.file.close()
    set_props('ali-select-x', {'options': obs_fields, 'value': None})
    set_props('ali-select-y', {'options': obs_fields, 'value': None})
    set_props('ali-select-z', {'options': obs_fields, 'value': None})
    set_props('ali-store-importdata', {'data': path})
    set_head_notice('Import successfully!', 'success')
    return True


def open_new_project_modal():
    """
    显示新建项目面板，并清空所有状态
    """
    set_props('ali-modal-newproject', dict(visible=True))
    set_props('ali-select-x', dict(value=None, options=[]))
    set_props('ali-select-y', dict(value=None, options=[]))
    set_props('ali-select-z', dict(value=None, options=[]))
    set_props('ali-store-importdata', dict(data=None))
    set_props('ali-projectname', dict(value=None, status='error'))

def create_project(importdata, x, y, z, projectname):
    """
    创建对齐项目
    """
    exist_projects = alidata.get_exist_projects()
    if projectname in exist_projects:
        set_head_notice(f'Project Name {projectname} Already Exists, Please Input Another Name', type='warning')
        return
    
    alidata.create_project(importdata, x, y, z, projectname)
    set_head_notice('Project Created Successfully', type='success')
    set_props('ali-select-project', dict(value=projectname))
    set_props('ali-modal-newproject', dict(visible=False))

def update_project_list():
    """
    更新项目列表
    """
    project_list = alidata.get_exist_projects()
    return project_list
