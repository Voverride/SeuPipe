import scanpy as sc
from dataManager.alignment_d import alidata
from controller.notice import set_head_notice, set_aside_notice
from utils.typing import StepStatus
from api.alignment import *
from utils.colors import get_color_map, get_scale_colors
import threading
import gc
import feffery_antd_components as fac
import numpy as np
from dash import Patch, set_props, html

def set_color_bygene(patch, project, gene, colorType='C1'):
    """
    设置切片颜色基于基因表达
    """
    if colorType=='C1':
        colorType='red'
    elif colorType=='C2':
        colorType='green'
    elif colorType=='C3':
        colorType='blue'
    initPath = alidata.get_initialdata_path(project)
    initData = sc.read_h5ad(initPath, backed='r')
    geneData = initData[:, gene].to_memory()
    initData.file.close()
    del initData
    gc.collect()
    normcountField = alidata.get_normcount_field()
    geneCounts = geneData[:, gene].layers[normcountField]
    minValue = float(np.min(geneCounts))
    maxValue = float(np.max(geneCounts))
    znums = alidata.get_znums(project)
    initFig = alidata.get_initialfig(project)
    for idx in range(len(znums)):
        obsIndex = initFig['data'][idx]['customdata'].flatten()
        geneExp = geneData[obsIndex, gene].layers[normcountField]
        colorList = get_scale_colors(geneExp, minValue, maxValue, colorType=colorType)
        patch['data'][idx]['marker']['color'] = colorList

def get_field_legends(project, field, startIndex, colorType='COLORS_60', fetchSize=30):
    """
    获取字段图例
    """
    coordinate = alidata.get_coordinate(project)
    labels = set(coordinate[field].unique())
    colorMap = get_color_map(labels, type=colorType)
    legend_items = list(labels)
    legend_items.sort()
    legendColor = {}
    length = len(legend_items)
    for i in range(startIndex, startIndex+fetchSize):
        if i<length:
            label = legend_items[i]
            legendColor[label] = colorMap[label]
        else:
            break
    hasLegend = True
    if startIndex+fetchSize>=length:
        hasLegend = False
    legends = [get_legend_item(legendColor[key], key) for key in legendColor]
    return legends, hasLegend

def get_legend_item(color, label):
    """
    获取图例项
    """
    return fac.AntdTooltip(
        html.Div(
            children=[
                html.Div(
                    style={
                        'width': '10px',
                        'height': '10px',
                        'borderRadius': '50%',
                        'backgroundColor': color,
                        'marginRight': '10px',
                    }
                ),
                html.Div(
                    label, 
                    style={
                        'marginBottom': '4px',
                        'whiteSpace': 'nowrap',
                        'overflow': 'hidden', 
                        'textOverflow': 'ellipsis', 
                        'maxWidth': '150px', 
                    }
                ),
            ],
            style={
                'display': 'flex',
                'alignItems': 'center',
            }
        ), 
        title=fac.AntdText(label, style={'color':'white'}),
        color=color
    )

def set_color_byfield(patch, project, field, colorType='COLORS_60'):
    """
    基于字段分配颜色
    """
    coordinate = alidata.get_coordinate(project)
    labels = set(coordinate[field].unique())
    colorMap = get_color_map(labels, type=colorType)
    znums = alidata.get_znums(project)
    initFig = alidata.get_initialfig(project)
    for idx in range(len(znums)):
        obsIndex = initFig['data'][idx]['customdata'].flatten()
        fieldValues = coordinate.loc[obsIndex, field]
        patch['data'][idx]['marker']['color'] = [colorMap[val] for val in fieldValues]

def get_sync_slices(project, selectedRowKeys)->set:
    """
    获取同步切片
    """
    syncSlices = set()
    znum = alidata.get_znums(project)
    znum = znum[::-1]
    for key in selectedRowKeys:
        syncSlices.add(znum[int(key)])
    return syncSlices

def update_relayoutfig(patch, relayoutData):
    """
    更新视图
    """
    if relayoutData is None:
        return
    if 'scene.camera' in relayoutData:
        camera = relayoutData['scene.camera']
        if 'projection' in camera:
            patch['layout']['scene']['projection']['type'] = camera['projection']['type']
        if 'center' in camera:
            patch['layout']['scene']['camera']['center'] = camera['center']
        if 'eye' in camera:
            patch['layout']['scene']['camera']['eye'] = camera['eye']
        if 'up' in camera:
            patch['layout']['scene']['camera']['up'] = camera['up']
    if 'scene.aspectratio' in relayoutData:
        patch['layout']['scene']['aspectmode'] = 'manual'
        patch['layout']['scene']['aspectratio'] = relayoutData['scene.aspectratio']

def get_transformed_coord(project, transMtx, slices):
    """
    获取转换坐标
    """
    metadata = alidata.get_project_info(project)
    x = metadata['x']
    y = metadata['y']
    z = metadata['z']
    coordinate = alidata.get_coordinate(project)
    operations = dict()
    for slice in slices:
        idx = alidata.get_slice_index(project, slice)
        obsIndex = coordinate[coordinate[z]==slice].index
        points_ori = np.array(coordinate.loc[obsIndex, [x, y]])
        points_new = transform_points(points_ori,transMtx)
        coordinate.loc[obsIndex, [x, y]] = points_new
        operations[idx] = {
            'x': list(points_new[:, 0]),
            'y': list(points_new[:, 1]),
        }
    alidata.update_manual_adjust_status(project, operations)

def transform_points(points:np.ndarray, coordTransMtx:list)->np.ndarray:
    """
    转换坐标
    """
    dx, dy, angle_deg = coordTransMtx
    center = np.mean(points, axis=0)    
    points_centered = points - center
    angle_rad = np.radians(angle_deg)
    rotation_matrix = np.array([
        [np.cos(angle_rad), -np.sin(angle_rad)],
        [np.sin(angle_rad),  np.cos(angle_rad)]
    ])
    rotated_points = points_centered @ rotation_matrix
    transformed_points = rotated_points + center + np.array([dx, dy])
    return transformed_points


def update_reference_slice(project, slice, actSlice, color, icon, colorMode, colorField, colorGene):
    """
    更新参考切片显示配置
    """
    patch = Patch()
    znums = alidata.get_znums(project)
    for i, znum in enumerate(znums):
        if znum == slice:
            if colorMode=='Field' and not colorField or colorMode=='Gene' and not colorGene or colorMode=='Custom':
                patch['data'][i]['marker']['color'] = color
            patch['data'][i]['visible'] = True
        elif znum != actSlice:
            patch['data'][i]['visible'] = 'legendonly'
    set_props('ali-graph-left', dict(figure=patch))
    if icon=='antd-minus':
        set_props('ali-graph-right', dict(figure=patch))

def update_active_slice(project, slice, refSlice, color, icon, colorMode, colorField, colorGene):
    """
    更新当前切片显示配置
    """
    patch = Patch()
    znums = alidata.get_znums(project)
    for i, znum in enumerate(znums):
        if znum == slice:
            if colorMode=='Field' and not colorField or colorMode=='Gene' and not colorGene or colorMode=='Custom':
                patch['data'][i]['marker']['color'] = color
            patch['data'][i]['visible'] = True
            if i > 0:
                set_props('ma-selecter-referenceSlice', dict(value=znums[i-1]))
        elif znum != refSlice:
            patch['data'][i]['visible'] = 'legendonly'
    # patch['layout']['scene']['camera']['eye'] = {'x': 1.2149787566408422e-05, 'y': -0.006897060070669372, 'z': 2.1650525237080886}
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
            patch['data'][i]['visible'] = 'legendonly'
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
            set_props('ali-store-figureScale', dict(data=scalePatch))
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
        set_props('ma-selecter-activeSlice', dict(options=znums, value=None))
        set_props('ma-selecter-referenceSlice', dict(options=znums, value=None))
        znumsReversed = [{'slices': z} for z in znums[::-1]]
        set_props('ma-table-SyncSlices', dict(data=znumsReversed))
        fields = alidata.get_fields(project)
        genelist = alidata.get_gene_list(project)
        set_props('ma-selecter-colorField', dict(options=fields, value=None))
        set_props('ma-selecter-colorGene', dict(options=genelist, value=None))
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

def export_alignment_file(path):
    """
    导出对齐数据文件
    """
    project = alidata.get_export_project()
    if project is None:
        set_head_notice('Timeout export, please click export button again!', type='warning')
        return False
    initial_path = alidata.get_initialdata_path(project)
    initialData = sc.read_h5ad(initial_path)
    coordinate = alidata.get_coordinate(project)
    metadata = alidata.get_project_info(project)
    x = metadata['x']
    y = metadata['y']
    initialData.obs[f'{x}_aligned'] = coordinate[x].values
    initialData.obs[f'{y}_aligned'] = coordinate[y].values
    normcountField = alidata.get_normcount_field()
    if normcountField in initialData.layers:
        del initialData.layers[normcountField]
    initialData.write_h5ad(path)
    set_head_notice('export successfully!', type='success')
    return True

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
