import anndata
from lightning.pytorch.callbacks import Callback
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import VotingClassifier
from sklearn.neighbors import KNeighborsClassifier
import scanpy as sc
from utils.colors import *
from utils.typing import StepStatus
from dataManager.annotation_d import annData
import traceback
import time
import os
import torch
import numpy as np
import plotly.express as px
import pandas as pd
import scvi

BATCH_KEY = 'SeuPipe_Batch'
DEGCOUNT_KEY = 'SeuPipe_DegCount'

def scvi_annotation(observed_metadata):
    """
    执行细胞注释程序
    """
    try:
        cur_step = 'preprocess'
        set_steps_status(observed_metadata, 'preprocess', StepStatus.PROCESS, notify=True)
        project_name = observed_metadata['project_name']
        rm_mt = observed_metadata['rm_mt']
        rm_ribo = observed_metadata['rm_ribo']
        rm_hb = observed_metadata['rm_hb']
        use_hvg = observed_metadata['use_hvg']
        x = observed_metadata['x']
        y = observed_metadata['y']
        z = observed_metadata['z']
        label_field = observed_metadata['label']
        ref_path = annData.get_refdata_path(project_name)
        query_path = annData.get_querydata_path(project_name)
        refdata = sc.read_h5ad(ref_path)
        querydata = sc.read_h5ad(query_path)
        adata_combined = scvi_preprocessing(refdata, querydata, label_field, rm_mt, rm_ribo, rm_hb, use_hvg)
        cur_step = 'training'
        set_steps_status(observed_metadata, 'preprocess', StepStatus.FINISH)
        set_steps_status(observed_metadata, 'training', StepStatus.PROCESS, notify=True)
        n_layers = observed_metadata['n_layers']
        n_hiddens = observed_metadata['n_hiddens']
        n_latent = observed_metadata['n_latent']
        epochs = observed_metadata['epochs']
        batch_size = observed_metadata['batch_size']
        dropout = observed_metadata['dropout']
        adata_combined = get_scvi_latent(adata_combined, n_layers, n_hiddens, n_latent, epochs, batch_size, dropout, observed_metadata)
        cur_step = 'postprocess'
        set_steps_status(observed_metadata, 'training', StepStatus.FINISH)
        set_steps_status(observed_metadata, 'percent', 0)
        set_steps_status(observed_metadata, 'postprocess', StepStatus.PROCESS, notify=True)
        label_field = observed_metadata['label']
        adata_result = annotation_with_scvi_latend(adata_combined, label_field)
        query_index = querydata.obs.index
        querydata.obs[label_field] = ''
        querydata.obs.loc[query_index, label_field] = adata_result.obs.loc[query_index, label_field]
        classes = len(set(querydata.obs[label_field].unique()))
        z_min = querydata.obs[z].min()
        z_max = querydata.obs[z].max()
        if isinstance(z_min, (np.int64, np.int32, np.int16, np.int8)):
            z_min = int(z_min)
        else:
            z_min = float(z_min)
        if isinstance(z_max, (np.int64, np.int32, np.int16, np.int8)):
            z_max = int(z_max)
        else:
            z_max = float(z_max)
        annData.set_parameter(project_name, dict(classes=classes, z_min=z_min, z_max=z_max))
        heatmap_fig = get_diffgene_heatmap(querydata, label_field)
        result_fig = plot_3d_scatter(querydata, x, y, z, label_field)
        annData.set_diffgene_fig(project_name, heatmap_fig)
        annData.set_annotation_fig(project_name, result_fig)
        querydata.write_h5ad(query_path)
        set_steps_status(observed_metadata, 'postprocess', StepStatus.FINISH)
    except Exception as e:
        set_steps_status(observed_metadata, cur_step, StepStatus.ERROR)
        observed_metadata['exception'] = traceback.format_exc()
    finally:
        time.sleep(1)
        observed_metadata['running'] = False
        project = observed_metadata['project_name']
        annData.update_project_info(project, dict(observed_metadata))
        annData.remove_annotation_task(project)
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

def scvi_preprocessing(adata_ref, adata_query, field, rm_mt, rm_ribo, rm_hb, use_hvg):
    """
    scvi数据预处理与整合
    """
    cache_raw(adata_ref)
    cache_raw(adata_query)
    adata_ref_tmp = sc.AnnData(
        X=adata_ref.X.copy(),
        obs=adata_ref.obs.copy(),
        var=adata_ref.var.copy(),
        layers=adata_ref.layers.copy()
    )
    adata_query_tmp = sc.AnnData(
        X=adata_query.X.copy(),
        obs=adata_query.obs.copy(),
        var=adata_query.var.copy(),
        layers=adata_query.layers.copy()
    )
    adata_ref = adata_ref_tmp
    adata_query = adata_query_tmp
    adata_ref = filter_genes(adata_ref, remove_mt=rm_mt, remove_ribo=rm_ribo, remove_hb=rm_hb).copy()
    adata_query = filter_genes(adata_query, remove_mt=rm_mt, remove_ribo=rm_ribo, remove_hb=rm_hb).copy()
    
    common_genes = adata_ref.var_names.intersection(adata_query.var_names)
    adata_ref = adata_ref[:, common_genes].copy()
    adata_query = adata_query[:, common_genes].copy()

    if use_hvg:
        adata_marker = filter_marker_genes_adaptive(adata_ref, field)

        if adata_marker is not None:
            adata_ref = adata_marker
            common_genes = adata_ref.var_names.intersection(adata_query.var_names)
            adata_ref = adata_ref[:, common_genes].copy()
            adata_query = adata_query[:, common_genes].copy()

    adata_ref.obs[BATCH_KEY] = "ref"
    adata_query.obs[BATCH_KEY] = "query"

    adata_combined = anndata.concat(
        [adata_ref, adata_query],
        axis=0,                  
        join="outer",
        merge="unique",       
        label=BATCH_KEY,
        keys=["ref", "query"]
    )
    # 这里因界面修改为use_marker，代码逻辑需要调整为use_marker为False时才进行hvg筛选，为了避免修改数据库仍然使用use_hvg变量
    if not use_hvg:
        sc.pp.highly_variable_genes(
            adata_combined,
            flavor="seurat_v3",
            n_top_genes=2000,
            layer="counts",
            batch_key=BATCH_KEY,
            subset=True
        )
    return adata_combined
def get_scvi_latent(adata_combined, n_layers, n_hiddens, n_latent, epochs, batch_size, dropout, observed_metadata):
    """
    训练scvi模型,提取潜在特征
    """
    if torch.cuda.is_available():
        torch.set_float32_matmul_precision('medium')
        torch.backends.cudnn.benchmark = True
    
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['NUMEXPR_NUM_THREADS'] = '1'

    scvi.model.SCVI.setup_anndata(
        adata_combined,
        batch_key=BATCH_KEY,
        layer="counts"
    )

    model = scvi.model.SCVI(adata_combined, n_hidden=n_hiddens, n_latent=n_latent, n_layers=n_layers, dropout_rate=dropout)
    
    def get_optimal_print_interval(total_batches, target_updates=10):
        """
        根据总batch数自动计算合适的print_interval
        target_updates: 期望的进度更新次数
        """
        if total_batches <= target_updates:
            return 1
        return max(1, total_batches // target_updates)

    class CustomProgressPrinter(Callback):
        def __init__(self, target_updates=10):
            self.target_updates = target_updates
            self.print_interval = None
            self.previous_progress = 0
            self.min_update_interval = 2

        def on_train_epoch_start(self, trainer, pl_module):
            total_batches = trainer.num_training_batches * trainer.max_epochs
            self.print_interval = get_optimal_print_interval(total_batches, self.target_updates)
            epoch = str(trainer.current_epoch+1)+'/'+str(trainer.max_epochs)
        
        def on_train_batch_end(self, trainer, pl_module, outputs, batch, batch_idx):
            if batch_idx % self.print_interval == 0:
                batches = trainer.num_training_batches
                pre_epoch = trainer.current_epoch
                all_epoch = trainer.max_epochs
                progress = (batch_idx + 1 + pre_epoch*batches) / (batches * all_epoch) * 100
                if progress - self.previous_progress >= self.min_update_interval:
                    self.previous_progress = progress
                    set_steps_status(observed_metadata, 'percent', int(progress), notify=True)

    model.train(
        max_epochs=epochs,
        batch_size=batch_size,
        callbacks=[CustomProgressPrinter()],
        enable_progress_bar=False,
        accelerator='gpu' if torch.cuda.is_available() else 'cpu',
    )

    latent_combined = model.get_latent_representation(adata_combined)
    normalized_exp = model.get_normalized_expression(adata_combined)
    adata_combined.obsm["X_scVI"] = latent_combined
    adata_combined.layers["X_normalized_scVI"] = normalized_exp.values
    return adata_combined

def plot_3d_scatter(querydata, x, y, z, label, marker_size:int=5, boarder_width:int=1, boarder_color:str='#0d0015'):
    """
    绘制三维散点图
    """
    if isinstance(querydata, pd.DataFrame):
        df = querydata
    else:
        df = querydata.obs.copy()
    fields = set(df[label].unique())
    cmap = get_color_map(fields, type='COLORS_60')
    df = df.sort_values(by=label)
    fig = px.scatter_3d(
        df, x=x, y=y, z=z, 
        color=label, 
        color_discrete_map=cmap,
    )

    fig.update_traces(
        marker=dict(
            size=marker_size,
            line=dict(color=boarder_color, width=boarder_width)
        )
    )
    fig.update_layout(
        autosize=True,
        title=dict(
            text='3D Annotation Visualization',
            x=0.5,
        ),
        legend=dict(
            itemsizing='constant',
            title=None
        ),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        scene=dict(
            xaxis=dict(visible=False, showgrid=True, showbackground=False, showticklabels=False, title=''), 
            yaxis=dict(visible=False, showgrid=True, showbackground=False, showticklabels=False, title=''), 
            zaxis=dict(visible=False, showgrid=False, showbackground=False, showticklabels=False, title=''),
            camera=dict(projection=dict(type='orthographic'))
        ),
    )
    return fig

def annotation_with_scvi_latend(adata, label_field):
    """
    对查询数据进行注释
    """
    latent_ref = adata[adata.obs[BATCH_KEY] == 'ref'].obsm["X_scVI"]
    labels_ref = adata[adata.obs[BATCH_KEY] == 'ref'].obs[label_field]
    latent_query = adata[adata.obs[BATCH_KEY] == 'query'].obsm["X_scVI"]
    adata_query = adata[adata.obs[BATCH_KEY] == 'query'].copy()
    estimators = [
        ("knn", KNeighborsClassifier(n_neighbors=15, weights="distance", metric="cosine")),
        ("svm", SVC(probability=True, class_weight="balanced")),
        ("rf", RandomForestClassifier(n_estimators=200, max_depth=20)),
        ("mlp", MLPClassifier(hidden_layer_sizes=(128, 64), activation="relu", alpha=0.001, batch_size=128, early_stopping=True))
    ]
    ensemble = VotingClassifier(estimators, voting="soft")
    ensemble.fit(latent_ref, labels_ref)
    predicted_labels = ensemble.predict(latent_query)
    adata_query.obs[label_field] = predicted_labels
    return adata_query

def get_diffgene_heatmap(adata_ori, label_field, top_n=5, min_logfc=0.25, min_pval=0.05):
    """
    计算差异表达基因并生成热图
    """
    import warnings
    warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)
    
    adata = adata_ori.copy()
    
    if 'counts' in adata.layers:
        adata.X = adata.layers['counts'].copy()
    else:
        print("Warning: 'counts' layer not found, using adata.X")
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    celltype_counts = adata.obs[label_field].value_counts()
    valid_celltypes = set(celltype_counts[celltype_counts >= 2].index.tolist())
    adata_filtered = adata[adata.obs[label_field].isin(valid_celltypes)].copy()
    
    if len(valid_celltypes) == 0:
        raise ValueError("Error: No valid cell types with >= 2 cells")
    
    sc.tl.rank_genes_groups(
        adata_filtered,
        groupby=label_field,
        method='wilcoxon',
        use_raw=False,
        key_added='markers',
        n_genes=adata_filtered.shape[1]
    )
    
    result_df = sc.get.rank_genes_groups_df(adata_filtered, group=None, key='markers')
    
    selected_genes = []
    deg_info = {}
    
    for cell_type in adata_filtered.obs[label_field].cat.categories:
        type_res = result_df[result_df['group'] == cell_type].copy()
        
        significant = type_res[
            (type_res['pvals_adj'] < min_pval) & 
            (type_res['logfoldchanges'] > min_logfc)
        ]
        
        top_genes = significant.head(top_n)['names'].tolist()
        selected_genes.extend(top_genes)
        
        deg_info[cell_type] = {
            'n_significant': len(significant),
            'n_selected': len(top_genes),
            'genes': top_genes
        }
    
    unique_genes = list(dict.fromkeys(selected_genes))
    
    valid_genes = [g for g in unique_genes if g in adata_filtered.var_names]
    
    if len(valid_genes) < 5:
        raise ValueError(f"Too few valid genes ({len(valid_genes)}). Try relaxing thresholds.")
    
    celltypes = sorted(adata_filtered.obs[label_field].cat.categories.tolist())
    
    heatmap_data = []
    for gene in valid_genes:
        row = []
        for cell_type in celltypes:
            mask = (adata_filtered.obs[label_field] == cell_type)
            exp_values = adata_filtered[mask, gene].X.toarray().flatten() if hasattr(adata_filtered[mask, gene].X, 'toarray') else adata_filtered[mask, gene].X.flatten()
            mean_exp = np.mean(exp_values)
            row.append(mean_exp)
        heatmap_data.append(row)
    
    heatmap_data = np.array(heatmap_data)
    
    def minmax_scale_safe(row):
        row_min = np.min(row)
        row_max = np.max(row)
        if row_max - row_min < 1e-10:
            return np.zeros_like(row)
        return (row - row_min) / (row_max - row_min)
    
    def zscore_scale_safe(row):
        """
        Z-score 标准化：z = (x - mean) / std
        添加除零保护：如果标准差接近 0，返回全 0
        """
        row_mean = np.mean(row)
        row_std = np.std(row)
        if row_std < 1e-10:
            return np.zeros_like(row)
        return (row - row_mean) / row_std
    
    # minmax_data = np.apply_along_axis(minmax_scale_safe, axis=1, arr=heatmap_data)
    zscore_data = np.apply_along_axis(zscore_scale_safe, axis=1, arr=heatmap_data)
    
    fig = px.imshow(
        zscore_data,
        labels=dict(x="Cell Type", y="Gene", color="Z-score"),
        x=celltypes,
        y=valid_genes,
        color_continuous_scale='Cividis',
        aspect='auto'
    )
    
    fig.update_layout(
        autosize=True,
        title=dict(
            text='Differential Gene Expression Heatmap',
            x=0.5,
            font=dict(size=16)
        ),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        coloraxis_colorbar=dict(
            title='Z-score',
            titleside='right',
            titlefont=dict(size=12)
        ),
        xaxis_title=None,
        yaxis_title=None,
    )
    
    fig.update_xaxes(
        tickangle=45, 
        tickfont=dict(size=10)
    )
    
    fig.update_yaxes(
        tickfont=dict(size=9)
    )
    return fig

# def get_diffgene_heatmap(adata_ori, label_field):
#     """
#     返回基因表达热图
#     """
#     import warnings
#     warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)
#     adata = adata_ori.copy()
#     cache_raw(adata)
#     adata.layers[DEGCOUNT_KEY] = adata.layers['counts'].copy()
#     sc.pp.normalize_total(adata, target_sum=1e4, layer=DEGCOUNT_KEY)
#     sc.pp.log1p(adata, layer=DEGCOUNT_KEY)
#     celltype_counts = adata.obs[label_field].value_counts()
#     valid_celltypes = set(celltype_counts[celltype_counts >= 2].index.tolist())
#     adata_filtered = adata[adata.obs[label_field].isin(valid_celltypes)].copy()
#     sc.tl.rank_genes_groups(
#         adata_filtered,
#         groupby=label_field,
#         layer=DEGCOUNT_KEY, 
#         method='wilcoxon',
#         n_genes=adata_filtered.shape[1],
#         use_raw=False
#     )
#     adata_ori.uns['rank_genes_groups'] = adata_filtered.uns['rank_genes_groups'].copy()
#     deg_results = adata_filtered.uns['rank_genes_groups']
#     celltypes = sorted(deg_results['names'].dtype.names)
#     genes = []
#     for ct in celltypes:
#         top3_genes = adata_filtered.uns['rank_genes_groups']['names'][ct][:3]
#         genes.extend(top3_genes)

#     def get_mean_expression(celltype, gene):
#         adtmp = adata_filtered[adata_filtered.obs[label_field] == celltype, gene]
#         mean_exp = np.mean(adtmp.layers[DEGCOUNT_KEY])
#         return mean_exp

#     heatmap_data = [[get_mean_expression(cell, gene) for cell in celltypes]for gene in genes]
#     def minmax_scale(row):
#         return (row - np.min(row)) / (np.max(row) - np.min(row))
    
#     minmax_data = np.apply_along_axis(minmax_scale, axis=1, arr=heatmap_data)

#     fig = px.imshow(
#         minmax_data,
#         labels=dict(x="celltype", y="gene", color="mean_expression"),
#         x=celltypes,
#         y=genes,
#         color_continuous_scale='Cividis',
#         aspect='auto'
#     )
#     fig.update_layout(
#         height=800,
#         width=800,
#         autosize=True,
#         title=dict(
#             text='Differential Gene Heatmap',
#             x=0.5,
#         ),
#         paper_bgcolor='rgba(0,0,0,0)',
#         plot_bgcolor='rgba(0,0,0,0)',
#         coloraxis_colorbar_title_text='',
#         xaxis_title=None,
#         yaxis_title=None,
#     )
#     fig.update_xaxes(
#         tickangle=45, 
#         tickfont=dict(size=10)
#     )
#     fig.update_yaxes(
#         tickfont=dict(size=10)
#     )
#     return fig

def set_steps_status(observed_metadata, step, status, notify=False):
    """
    设置注释任务步骤状态
    """
    observed_metadata['steps'][step] = status
    if notify:
        observed_metadata['running'] = True

def cache_raw(adata):
    """
    缓存原始计数
    """
    if 'counts' in adata.layers:
        return
    if adata.raw is not None and adata.raw.X is not None:
        adata.layers["counts"] = adata.raw.X.copy()
    else:
        adata.layers["counts"] = adata.X.copy()

def filter_genes(adata, remove_mt=True, remove_ribo=True, remove_hb=True):
    """
    Filter genes based on their names.
    Args:
        adata: AnnData object
        remove_mt: bool, default True
        remove_ribo: bool, default True
        remove_hb: bool, default True
    Returns:
        adata: AnnData object
    """
    if not remove_mt and not remove_ribo and not remove_hb:
        return adata
    lower_var_names = adata.var_names.str.lower().str
    condition = True
    if remove_mt:
        condition = condition & ~lower_var_names.startswith("mt-")
    if remove_ribo:
        condition = condition & ~lower_var_names.startswith(("rps", "rpl"))
    if remove_hb:
        condition = condition & ~lower_var_names.contains("^hb[^(p)]")
    adata = adata[:, condition]
    return adata

def filter_marker_genes_adaptive(
    adata, 
    field='celltype', 
    min_logfc=0.25, 
    min_pval=0.05,
    max_markers_per_type=50,
    min_markers_per_type=10,
    min_cells_per_type=2
):
    """
    自适应筛选 marker 基因的函数，基于统计显著性动态确定数量。
    """
    adata.obs[field] = adata.obs[field].astype('category')

    celltype_counts = adata.obs[field].value_counts()
    valid_celltypes = celltype_counts[celltype_counts >= min_cells_per_type].index.tolist()
    if len(valid_celltypes) == 0:
        raise ValueError(f"Error: No cell types with >= {min_cells_per_type} cells")
    adata_filtered = adata[adata.obs[field].isin(valid_celltypes)].copy()
    
    adata_filtered.obs[field] = adata_filtered.obs[field].cat.remove_unused_categories()

    adata_filtered.X = adata_filtered.layers['counts'].copy()
    
    sc.pp.normalize_total(adata_filtered, target_sum=1e4)
    sc.pp.log1p(adata_filtered)
    
    sc.tl.rank_genes_groups(
        adata_filtered, 
        groupby=field, 
        method='wilcoxon', 
        use_raw=False,
        key_added='markers',
        n_genes=adata_filtered.shape[1]
    )
    
    result_df = sc.get.rank_genes_groups_df(adata_filtered, group=None, key='markers')
    selected_genes = set()
    
    for cell_type in adata_filtered.obs[field].cat.categories:
        type_res = result_df[result_df['group'] == cell_type].copy()
        
        type_res = type_res.sort_values(by=['pvals_adj', 'logfoldchanges'], 
                                         ascending=[True, False])
        
        significant = type_res[
            (type_res['pvals_adj'] < min_pval) & 
            (type_res['logfoldchanges'] > min_logfc)
        ]
        
        n_selected = len(significant)
        n_selected = max(min_markers_per_type, min(n_selected, max_markers_per_type))
        
        top_genes = significant.head(n_selected)['names'].tolist()
        selected_genes.update(top_genes)
    
    valid_genes = [g for g in selected_genes if g in adata_filtered.var_names]
    
    if len(valid_genes) == 0:
        return None
    
    adata_marker = adata_filtered[:, valid_genes].copy()
    return adata_marker