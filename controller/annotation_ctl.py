from dataManager.annotation_d import annData
from controller.notice import set_head_notice
from dash import set_props
import scanpy as sc
import threading
from utils.commonfuc import *
import plotly.express as px
from api.annotation import *
import numpy as np
from dataManager.workspace import get_annotation_workspace
import plotly.graph_objects as go

result_fig_style = {'width': '100%', 'height':'90%', 'margin':'auto'}
heatmap_fig_style = {'width': '100%', 'height':'90%', 'margin':'auto'}
queryumap_fig_style = {'width': '100%', 'aspectRatio': '3/2', 'minWidth':'30vw', 'maxWidth':'50vw', 'margin':'auto'}
def export_annotation_file(path:str)->bool:
    """
    导出注释结果数据
    """
    access_path = path
    if not os.path.exists(path):
        access_path = os.path.dirname(path)
    if not os.access(access_path, os.W_OK):
        set_head_notice('You have no permission to write data to this folder !', 'warning')
        return False
    status = False
    try:
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        annData.get_querydata().write_h5ad(path, compression='gzip')
    except Exception as e:
        set_head_notice('export failed!', 'error')
        raise e
    else:
        status = True
        set_head_notice('export successfully!', 'success')
    return status

def check_result_data():
    """
    检查导出数据是否存在
    """
    querydata = annData.get_querydata()
    if querydata is None:
        set_head_notice('There was no query data exist, Please importing previsously !', 'warning')
        return False
    return True

def run_scvi(rm_mt, rm_ribo, rm_hb, use_hvg, n_layers, n_hiddens, n_latent, epochs, batch_size, dropout, usrname):
    """
    执行scvi细胞注释程序
    """
    status = check_required_data()
    if not status:
        return False
    refdata = annData.get_refdata()
    querydata = annData.get_querydata()
    label_field = annData.get_labelfield()
    annData.reset_anntstatus()
    annstatus = annData.get_anntstatus()
    thread = threading.Thread(
        target=scvi_annotation, 
        args=(refdata, querydata, label_field, rm_mt, rm_ribo, rm_hb, use_hvg, n_layers, n_hiddens, n_latent, epochs, batch_size, dropout, annData)
    )
    annData.set_annstatus_props(usrname, thread, rm_mt, rm_hb, rm_ribo, use_hvg, epochs, batch_size, n_layers, dropout, n_hiddens, n_latent)
    thread.start()
    reset_annstatus_progress(usrname)
    return True

# def run_tosica(rm_mt, rm_ribo, rm_hb, use_hvg, epoch, depth, batchsize, lr, gmt, usrname):
#     """
#     执行tosica细胞注释程序
#     """
#     status = check_required_data(gmt)
#     if not status:
#         return False
#     refdata = annData.get_refdata()
#     querydata = annData.get_querydata()
#     label_field = annData.get_labelfield()
#     project = annData.get_project_name()
#     annotation_path = get_annotation_workspace()
#     annData.reset_anntstatus()
#     annstatus = annData.get_anntstatus()
#     thread = threading.Thread(
#         target=tosica, 
#         args=(refdata, querydata, label_field, project, annotation_path, rm_mt, rm_ribo, rm_hb, use_hvg, epoch, depth, batchsize, lr, gmt), 
#         kwargs={'annstatus': annstatus}
#     )
#     annData.set_annstatus_props(usrname, thread, rm_mt, rm_hb, rm_ribo, use_hvg, epoch, batchsize, gmt, depth, lr)
#     thread.start()
#     reset_annstatus_progress(usrname)
#     evaluate_fig = get_evaluate_fig(annstatus)
#     set_props('annotation-graph-evaluate', dict(figure=evaluate_fig))
#     return True

def reset_annstatus_progress(usrname)->None:
    """
    启动训练任务时设置页面组件状态
    """
    set_props('annotask-timeline-creator', dict(children=usrname))
    set_props('annotask-button-showBug', dict(style={'backgroundColor':'#bb5548', 'display':'none'}))

    for i in range(1, 4):
        set_step_status(i, 'schedule', 0)
    set_disabled_status(True)


def update_annotation_progress():
    """
    轮回查询注释进度
    """
    annotstatus = annData.get_anntstatus()
    steps = annotstatus['steps']
    thread = annotstatus['thread']
    success = False
    for i in range(1, 4):
        if steps[i]['failed']:
            set_disabled_status(False)
            set_step_status(i, 'failed')
            set_props('annotation-interval', {'disabled':True})
            e = annotstatus['exception']
            if e:
                set_props('annotask-button-showBug', dict(style={'backgroundColor':'#bb5548'}))  
            return
        elif steps[i]['complete']:
            set_step_status(i, 'complete', 100)
            if i==2:
                set_props('annot-percent', dict(percent=100))
            if i==3:
                success = True
        else:
            if thread is not None and thread.is_alive():
                set_step_status(i, 'running', steps[2]['percent'])
            return
    if success:
        set_props('annotation-interval', {'disabled':True})
        set_disabled_status(False)
        result_fig = annData.get_resultfig()
        queryumap_fig = annData.get_queryumap()
        heatmap_fig = annData.get_heatmap()
        set_props('annotation-graph-result', dict(figure=result_fig))
        set_props('annotation-graph-result', dict(style=result_fig_style))
        set_props('annotation-graph-queryumap', dict(figure=queryumap_fig))
        set_props('annotation-graph-queryumap', dict(style=queryumap_fig_style))
        set_props('annotation-graph-heatmap', dict(figure=heatmap_fig))
        set_props('annotation-graph-heatmap', dict(style=heatmap_fig_style))
        

# def update_annotation_progress(epoch=None):
#     """
#     轮回查询注释进度
#     """
#     annotstatus = annData.get_anntstatus()
#     steps = annotstatus['steps']
#     curEpoch = len(annotstatus['evaluate']['trainAcc'])
#     evaluate_fig = None
#     if epoch is not None and curEpoch!= epoch:
#         evaluate_fig = get_evaluate_fig(annotstatus)
#         set_props('annotask-epoch', dict(data=curEpoch))
#     thread = annotstatus['thread']
#     success = False
#     for i in range(1, 4):
#         if steps[i]['failed']:
#             set_disabled_status(False)
#             set_step_status(i, 'failed')
#             set_props('annotation-interval', {'disabled':True})
#             e = annotstatus['exception']
#             if e:
#                 set_props('annotask-button-showBug', dict(style={'backgroundColor':'#bb5548'}))  
#             return
#         elif steps[i]['complete']:
#             set_step_status(i, 'complete', 100)
#             if i==2 and evaluate_fig is not None:
#                 set_props('annotation-graph-evaluate', dict(figure=evaluate_fig))
#             if i==3:
#                 success = True
#         else:
#             if thread is not None and thread.is_alive():
#                 set_step_status(i, 'running', steps[2]['percent'])
#                 if i==2 and evaluate_fig is not None:
#                     set_props('annotation-graph-evaluate', dict(figure=evaluate_fig))
#             return
#     if success and epoch is not None:
#         set_props('annotation-interval', {'disabled':True})
#         set_disabled_status(False)
#         resultfig = plot_3d_scatter()
#         annData.set_resultfig(resultfig)
#         refumap = plot_2d_umap(ref=True)
#         queryumap = plot_2d_umap(ref=False)
#         annData.set_refumap(refumap)
#         annData.set_queryumap(queryumap)
#         set_props('annotation-graph-result', dict(figure=resultfig))
#         set_props('annotation-graph-refumap', dict(figure=refumap))
#         set_props('annotation-graph-queryumap', dict(figure=queryumap))
def set_step_status(step:int, status:str, percent:int=None)->None:
    """
    设置每一步状态
    """
    id = 'annot-step'+str(step)
    icon_map = {
        'schedule': 'md-schedule',
        'running': 'antd-loading',
        'complete': 'fc-ok',
        'failed': 'fc-high-priority'
    }
    icon = icon_map[status]
    set_props(id, dict(icon=icon))
    if step==2 and percent is not None:
        set_props('annot-percent', dict(percent=percent))
def set_disabled_status(disabled:bool)->None:
    """
    设置参数控制和按钮是否为禁用状态
    """
    set_props('anntask-button-train', dict(disabled=disabled))
    set_props('annotask-button-refdata', dict(disabled=disabled))
    set_props('annotask-button-querydata', dict(disabled=disabled))
    set_props('annotask-select-label', dict(disabled=disabled))
    set_props('annotask-select-x', dict(disabled=disabled))
    set_props('annotask-select-y', dict(disabled=disabled))
    set_props('annotask-select-z', dict(disabled=disabled))

    
def check_required_data()->bool:
    """
    执行注释前检查数据完整性
    """
    refdata = annData.get_refdata()
    querydata = annData.get_querydata()
    if refdata is None:
        set_head_notice('Reference data is empty !', 'warning')
        return False
    if querydata is None:
        set_head_notice('Query data is empty !', 'warning')
        return False
    label = annData.get_labelfield()
    if label is None:
        set_head_notice('label field is empty !', 'warning')
        return False
    
    xfield = annData.get_xfield()
    if xfield is None:
        set_head_notice('x field is empty !', 'warning')
        return False
    
    yfield = annData.get_yfield()
    if yfield is None:
        set_head_notice('y field is empty !', 'warning')
        return False
    
    zfield = annData.get_zfield()
    if zfield is None:
        set_head_notice('z field is empty !', 'warning')
        return False
    return True

def get_evaluate_fig(anntstatus, train_color='#5383c3', valid_color='#d0826c', grid_color='#e5e4e6', line_width=2, width=None, height=None):
    """
    绘制训练验证曲线
    """ 
    evaluate = anntstatus['evaluate']
    epoch = anntstatus['epoch']
    train_acc = evaluate['trainAcc']
    valid_acc = evaluate['valAcc']
    train_loss = evaluate['trainLoss']
    valid_loss = evaluate['valLoss']
    epochs = list(range(epoch))
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=epochs,
        y=train_acc,
        name='Train Acc',
        line=dict(color=train_color, width=line_width),
        yaxis='y2',
    ))

    fig.add_trace(go.Scatter(
        x=epochs,
        y=valid_acc,
        name='Valid Acc',
        line=dict(color=valid_color, width=line_width),
        yaxis='y2',
    ))

    fig.add_trace(go.Scatter(
        x=epochs,
        y=train_loss,
        name='Train Loss',
        line=dict(color=train_color, width=line_width, dash='dot'),
        yaxis='y1',
    ))

    fig.add_trace(go.Scatter(
        x=epochs,
        y=valid_loss,
        name='Valid Loss',
        line=dict(color=valid_color, width=line_width, dash='dot'),
        yaxis='y1',
    ))

    fig.update_layout(
        title=dict(
            text='Training & Validation Metrics',
            x=0.5,
        ),
        xaxis=dict(
            title='Epochs',
            gridcolor=grid_color,
            dtick=1,
        ),
        yaxis=dict(
            title='Loss',
        ),
        yaxis2=dict(
            title='Accuracy',
            gridcolor=grid_color,
            overlaying='y',
            side='right',
            tickformat=".0%",
            range=[0, 1],
        ),
        legend=dict(
            x=1.1,
            y=0.9,
        ),
        hovermode='x unified',
        plot_bgcolor='white',
        width=width,
        height=height,
    )
    return fig
def restore_initial_data()->None:
    """
    恢复网页数据
    """
    ref_fields = annData.get_refdata_fields()
    query_fields = annData.get_querydata_fields()
    set_props('annotask-select-label', {'options':ref_fields, 'value':annData.get_labelfield()})
    set_props('annotask-select-x', {'options':query_fields, 'value':annData.get_xfield()})
    set_props('annotask-select-y', {'options':query_fields, 'value':annData.get_yfield()})
    set_props('annotask-select-z', {'options':query_fields, 'value':annData.get_zfield()})
    zfield_min = annData.get_zfield_min()
    zfield_max = annData.get_zfield_max()
    if zfield_min is not None and zfield_max is not None:
        props = dict()
        props['min'] = zfield_min
        props['max'] = zfield_max
        props['value'] = [zfield_min, zfield_max]
        props['tooltipPrefix'] = f'{annData.get_zfield()}: '
        set_props('annotation-slider-slicer', props)

    anntstatus = annData.get_anntstatus()
    set_props('annotask-remove-mt', dict(checked=anntstatus['rmMt']))
    set_props('annotask-remove-ribo', dict(checked=anntstatus['rmRibo']))
    set_props('annotask-remove-hb', dict(checked=anntstatus['rmHb']))
    set_props('annotask-use-hvg', dict(checked=anntstatus['useHvg']))
    set_props('annotask-batchsize', dict(value=anntstatus['batch_size']))
    set_props('annotask-nlatent', dict(value=anntstatus['n_latent']))
    set_props('annotask-dropout', dict(value=anntstatus['dropout']))
    set_props('annotask-nlayers', dict(value=anntstatus['n_layers']))
    set_props('annotask-epochs', dict(value=anntstatus['epochs']))
    set_props('annotask-nhiddens', dict(value=anntstatus['n_hiddens']))
    set_props('annot-percent', dict(children=anntstatus['steps'][2]['percent']))

    set_props('annotask-timeline-creator', dict(children=anntstatus['creator']))
    thread = anntstatus['thread']
    if thread is not None and thread.is_alive():
        set_props('annotation-interval', dict(disabled=False))
        set_disabled_status(True)

    update_annotation_progress()

    result_fig = annData.get_resultfig()
    if result_fig is not None:
        set_props('annotation-graph-result', dict(figure=result_fig))
        set_props('annotation-graph-result', dict(style=result_fig_style))

    queryumap = annData.get_queryumap()
    if queryumap is not None:
        set_props('annotation-graph-queryumap', dict(figure=queryumap))
        set_props('annotation-graph-queryumap', dict(style=queryumap_fig_style))    

    heatmap_fig = annData.get_heatmap()
    if heatmap_fig is not None:
        set_props('annotation-graph-heatmap', dict(figure=heatmap_fig))
        set_props('annotation-graph-heatmap', dict(style=heatmap_fig_style))

def check_queryfield_type(field)->bool:
    """
    检查查询字段对应的坐标字段数据类型
    """
    adata = annData.get_querydata()
    if adata is None:
        return False
    status = is_all_numeric(adata.obs[field])
    if status:
        return True
    set_head_notice(field+' is not numeric!', type='warning')
    return False

def read_annotask_querydata(path:str)->bool:
    """
    读取查询数据
    """
    status = False
    annData.reset_queryprops()
    try:
        adata = sc.read_h5ad(path)
        annData.set_querydata(adata)
        annData.set_queryname(os.path.basename(path))
        obs_fields = annData.get_querydata_fields()
    except Exception as e:
        set_head_notice('import failed!', 'error')
        raise e
    else:
        status = True
        set_props('annotask-select-x', {'options':obs_fields, 'value':None})
        set_props('annotask-select-y', {'options':obs_fields, 'value':None})
        set_props('annotask-select-z', {'options':obs_fields, 'value':None})
        props = dict()
        props['min'] = 0
        props['max'] = 0
        props['value'] =[]
        props['tooltipPrefix'] = ''
        set_props('annotation-slider-slicer', props)
        set_head_notice('import successfully!', 'success')
    return status
def read_annotask_refdata(path:str)->bool:
    """
    读取参考数据
    """
    annData.reset_refprops()
    status = False
    try:
        adata = sc.read_h5ad(path)
        annData.set_refdata(adata)
        annData.set_refname(os.path.basename(path))
        obs_fields = annData.get_refdata_fields()
    except Exception as e:
        set_head_notice('import failed!', 'error')
        raise e
    else:
        status = True
        set_props('annotask-select-label', {'options':obs_fields, 'value':None})
        set_head_notice('import successfully!', 'success')
    return status


def debug_run_tosica():
    ad1 = sc.read_h5ad('/data1/zhouyb/public/data/scMouseEmbryo/embryo_sc_mouse_E7.75_wy.h5ad')
    ad2 = sc.read_h5ad('/data1/zhouyb/public/data/stMouseEmbryo/embryo_1-2-E7.75_min400_Ann_HC0.5.h5ad')
    # ad1 = sc.read_h5ad('/data1/zhouyb/public/data/stMouseEmbryo/embryo_1-2-E7.75_min400_Ann_HC0.5.h5ad')
    # ad2 = sc.read_h5ad('/data1/zhouyb/public/data/stMouseEmbryo/embryo_1-2-E7.75_scs_aligned.h5ad')
    # ad2 = ad1
    # ad2 = sc.read_h5ad('/data1/zhouyb/public/data/BSAlignmentTestData/spatpy_scvi_annot.h5ad')
    annData.set_refdata(ad1)
    annData.set_querydata(ad2)

    annData.set_labelfield('celltype')
    annData.set_xfield('x')
    annData.set_yfield('y')
    annData.set_zfield('z')


    # annData.set_refname('embryo_1-2-E7.75_min400_Ann_HC0.5.h5ad')
    # annData.set_queryname('embryo_1-2-E7.75_scs_aligned.h5ad')

    # annstatus = annData.get_anntstatus()
    # annstatus['steps'] = {
    #     1:{'running':False, 'complete':True, 'failed':False},
    #     2:{'running':False, 'complete':True, 'failed':False, 'percent':100},
    #     3:{'running':False, 'complete':True, 'failed':False},
    # }

    # annstatus['creator'] = 'zhouyb'

    # resultfig = plot_3d_scatter()
    # annData.set_resultfig(resultfig)
    # ad1.obs[['umapX', 'umapY']] = ad1.obs[['x', 'y']]
    # ad2.obs[['umapX', 'umapY']] = ad2.obs[['x', 'y']]
    # refumap = plot_2d_umap(ref=True)
    # annData.set_refumap(refumap)
    # queryumap = plot_2d_umap(ref=False)
    # annData.set_queryumap(queryumap)
    
# debug_run_tosica()