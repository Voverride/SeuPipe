# from .TOSICA import *
import anndata
from lightning.pytorch.callbacks import Callback
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import VotingClassifier
from sklearn.neighbors import KNeighborsClassifier
import scanpy as sc
import torch
import numpy as np
import plotly.express as px
import pandas as pd
import scvi
import os

def scvi_annotation(refdata, querydata, label_field, rm_mt, rm_ribo, rm_hb, use_hvg, n_layers, n_hiddens, n_latent, epochs, batch_size, dropout, annData):
    """
    执行细胞注释程序
    """
    annstatus = annData.get_anntstatus()
    steps = annstatus['steps']
    try:
        adata_combined = scvi_preprocessing(refdata, querydata, rm_mt, rm_ribo, rm_hb, use_hvg)
    except Exception as e:
        set_steps_status(steps, 1, 'failed', True)
        if annstatus:
            annstatus['exception'] = e
            return False
        else:
            raise e
    set_steps_status(steps, 1, 'complete', True)
    try:   
        adata_combined = get_scvi_latent(adata_combined, n_layers, n_hiddens, n_latent, epochs, batch_size, dropout, annstatus)
    except Exception as e:
        set_steps_status(steps, 2, 'failed', True)
        if annstatus:
            annstatus['exception'] = e
            return False
        else:
            raise e
    set_steps_status(steps, 2, 'complete', True)
    try:
        adata_result = annotation_with_scvi_latend(adata_combined, label_field)
        query_index = querydata.obs.index
        querydata.obs.loc[query_index, label_field] = adata_result.obs.loc[query_index, label_field]
        adata_result.obs[['umapX', 'umapY']] = adata_result.obsm['X_umap']
        querydata.obs.loc[query_index, 'umapX'] = adata_result.obs.loc[query_index, 'umapX']
        querydata.obs.loc[query_index, 'umapY'] = adata_result.obs.loc[query_index, 'umapY']
        queryumap_fig = plot_2d_umap(annData, ref=False)
        annData.set_queryumap(queryumap_fig)
        heatmap_fig = get_diffgene_heatmap(querydata, label_field)
        annData.set_heatmap(heatmap_fig)
        result_fig = plot_3d_scatter(annData)
        annData.set_resultfig(result_fig)
        del adata_result
    except Exception as e:
        set_steps_status(steps, 3, 'failed', True)
        if annstatus:
            annstatus['exception'] = e
            return False
        else:
            raise e
    set_steps_status(steps, 3, 'complete', True)
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return True

def scvi_preprocessing(adata_ref, adata_query, rm_mt, rm_ribo, rm_hb, use_hvg):
    """
    scvi数据预处理与整合
    """
    cache_raw(adata_ref)
    cache_raw(adata_query)
    common_genes = adata_ref.var_names.intersection(adata_query.var_names)
    adata_ref = adata_ref[:, common_genes]
    adata_query = adata_query[:, common_genes]
    adata_ref = filter_genes(adata_ref, remove_mt=rm_mt, remove_ribo=rm_ribo, remove_hb=rm_hb).copy()
    adata_query = filter_genes(adata_query, remove_mt=rm_mt, remove_ribo=rm_ribo, remove_hb=rm_hb).copy()
    adata_ref.obs["spatpy_batch"] = "ref"
    adata_query.obs["spatpy_batch"] = "query"
    adata_combined = anndata.concat(
        [adata_ref, adata_query],
        axis=0,                  
        join="outer",
        merge="unique",       
        label="spatpy_batch",
        keys=["ref", "query"]
    )
    if use_hvg:
        sc.pp.highly_variable_genes(
            adata_combined,
            flavor="seurat_v3",
            n_top_genes=2000,
            layer="counts",
            batch_key="spatpy_batch",
            subset=True
        )
    return adata_combined
def get_scvi_latent(adata_combined, n_layers, n_hiddens, n_latent, epochs, batch_size, dropout, annstatus):
    """
    训练scvi模型,提取潜在特征
    """
    steps = annstatus['steps']
    scvi.model.SCVI.setup_anndata(
        adata_combined,
        batch_key="spatpy_batch",
        layer="counts"
    )

    model = scvi.model.SCVI(adata_combined, n_hidden=n_hiddens, n_latent=n_latent, n_layers=n_layers, dropout_rate=dropout)
    
    class CustomProgressPrinter(Callback):
        def __init__(self, print_interval=1):
            self.print_interval = print_interval

        def on_train_epoch_start(self, trainer, pl_module):
            steps[2]['epoch'] = str(trainer.current_epoch+1)+'/'+str(trainer.max_epochs)
        
        def on_train_batch_end(self, trainer, pl_module, outputs, batch, batch_idx):
            if batch_idx % self.print_interval == 0:
                batches = trainer.num_training_batches
                pre_epoch = trainer.current_epoch
                all_epoch = trainer.max_epochs
                progress = (batch_idx + 1 + pre_epoch*batches) / (batches * all_epoch) * 100
                steps[2]['percent'] = int(progress)

    model.train(
        max_epochs=epochs,
        batch_size=batch_size,
        callbacks=[CustomProgressPrinter(print_interval=1)],
        enable_progress_bar=False
    )

    latent_combined = model.get_latent_representation(adata_combined)
    normalized_exp = model.get_normalized_expression(adata_combined)
    adata_combined.obsm["X_scVI"] = latent_combined
    adata_combined.layers["X_normalized_scVI"] = normalized_exp.values
    return adata_combined

def plot_3d_scatter(annData, min_z=None, max_z=None, marker_size:int=5, boarder_width:int=1, boarder_color:str='#0d0015', grid_color:str = '#5F9EA0'):
    """
        绘制三维散点图
    """
    df = annData.get_querydata().obs.copy()
    x = annData.get_xfield()
    y = annData.get_yfield()
    z = annData.get_zfield()
    if min_z is not None and max_z is not None:
        df = df[(df[z]>=min_z) & (df[z]<=max_z)]
    label = annData.get_labelfield()
    cmap = annData.get_cmap()
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
            xaxis=dict(visible=False, showgrid=True, showbackground=False, showticklabels=False, gridcolor=grid_color, tickfont=dict(color=grid_color), title=''), 
            yaxis=dict(visible=False, showgrid=True, showbackground=False, showticklabels=False, gridcolor=grid_color, tickfont=dict(color=grid_color), title=''), 
            zaxis=dict(visible=False, showgrid=False, showbackground=False, showticklabels=False, gridcolor=grid_color),
            camera=dict(projection=dict(type='orthographic'))
        ),
    )
    return fig

def plot_2d_umap(annData, ref:bool=False, marker_size:int=5, boarder_width:int=0, boarder_color:str='#0d0015'):
    """
    绘制二维UMAP图
    """
    if ref:
        df = annData.get_refdata().obs.copy()
        title = 'Reference UMAP'
    else:
        df = annData.get_querydata().obs.copy()
        title = 'Predicted UMAP'
    label = annData.get_labelfield()
    cmap = annData.get_cmap()
    df = df.sort_values(by=label)
    fig = px.scatter(df, x='umapX', y='umapY', color=label, color_discrete_map=cmap)
    fig.update_traces(
        marker=dict(
            size=marker_size,
            line=dict(color=boarder_color, width=boarder_width)
        )
    )
    fig.update_layout(
        autosize=True,
        title=dict(
            text=title,
            x=0.5,
        ),
        legend=dict(
            itemsizing='constant',
            title=None
        ),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(
            visible=False,
            showticklabels=False,
            showgrid=False,
        ),
        yaxis=dict(
            visible=False,
            showticklabels=False,
            showgrid=False,
        ),
    )
    return fig

def annotation_with_scvi_latend(adata, label_field):
    """
    对查询数据进行注释
    """
    latent_ref = adata[adata.obs["spatpy_batch"] == 'ref'].obsm["X_scVI"]
    labels_ref = adata[adata.obs["spatpy_batch"] == 'ref'].obs[label_field]
    latent_query = adata[adata.obs["spatpy_batch"] == 'query'].obsm["X_scVI"]
    adata_query = adata[adata.obs["spatpy_batch"] == 'query'].copy()
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

    sc.pp.neighbors(adata_query, n_neighbors=15, use_rep='X_scVI')
    sc.tl.umap(adata_query, min_dist=0.3, spread=1.0, random_state=42)
    return adata_query

def tosica(ref_adata, query_adata, label_field, project, annotation_path, rm_mt, rm_ribo, rm_hb, use_hvg, epoch, depth, batchsize, lr, gmt, annstatus):
    """
    执行细胞注释程序
    """
    return False
    steps=None
    evaluate=None
    if annstatus is not None:
        steps = annstatus['steps']
        evaluate = annstatus['evaluate']
    
    refdata = ref_adata.copy()
    querydata = query_adata.copy()
    try:
        cache_raw(refdata)
        cache_raw(querydata)
        refdata = filter_genes(refdata, remove_mt=rm_mt, remove_hb=rm_hb, remove_ribo=rm_ribo)
        if use_hvg:
            common_genes = refdata.var_names.intersection(querydata.var_names)
            refdata = refdata[:, common_genes].copy()
            use_raw(refdata)
            normalize(refdata)
            sc.pp.highly_variable_genes(refdata)
            refdata = refdata[:, refdata.var['highly_variable']]
        common_genes = refdata.var_names.intersection(querydata.var_names)
        refdata = refdata[:, common_genes].copy()
        querydata = querydata[:, common_genes].copy()
        use_raw(refdata)
        use_raw(querydata)
        normalize(refdata)
        normalize(querydata)
    except Exception as e:
        set_steps_status(steps, 1, 'failed', True)
        if annstatus:
            annstatus['exception'] = e
            return False
        else:
            raise e
    set_steps_status(steps, 1, 'complete', True)
    try:    
        project_path = os.path.join(annotation_path, project)
        train(refdata, gmt_path=gmt, label_name=label_field, epochs=epoch, project=project_path, batch_size=batchsize, depth=depth, lr=lr, train_steps=steps[2], evaluate=evaluate)
    except Exception as e:
        set_steps_status(steps, 2, 'failed', True)
        if annstatus:
            annstatus['exception'] = e
            return False
        else:
            raise e
    set_steps_status(steps, 2, 'complete', True)
    try:
        model_weight_path = os.path.join(project_path, 'model.pth')
        prediction = pre(querydata, model_weight_path = model_weight_path, project=project_path)
        query_index = query_adata.obs.index
        ref_index = ref_adata.obs.index
        query_adata.obs.loc[query_index, label_field] = prediction.obs.loc[query_index, 'Prediction']
        query_adata.obs.loc[query_index, 'Probability'] = prediction.obs.loc[query_index, 'Probability']
        umap(refdata)
        umap(querydata)
        refdata.obs[['umapX', 'umapY']] = refdata.obsm['X_umap']
        querydata.obs[['umapX', 'umapY']] = querydata.obsm['X_umap']
        ref_adata.obs.loc[ref_index, 'umapX'] = refdata.obs.loc[ref_index, 'umapX']
        ref_adata.obs.loc[ref_index, 'umapY'] = refdata.obs.loc[ref_index, 'umapY']
        query_adata.obs.loc[query_index, 'umapX'] = querydata.obs.loc[query_index, 'umapX']
        query_adata.obs.loc[query_index, 'umapY'] = querydata.obs.loc[query_index, 'umapY']

    except Exception as e:
        set_steps_status(steps, 3, 'failed', True)
        if annstatus:
            annstatus['exception'] = e
            return False
        else:
            raise e
    set_steps_status(steps, 3, 'complete', True)
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

def get_diffgene_heatmap(adata_ori, label_field):
    """
    返回基因表达热图
    """
    import warnings
    warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)
    adata = adata_ori.copy()
    cache_raw(adata)
    adata.layers['spatpy_deg_count'] = adata.layers['counts'].copy()
    sc.pp.normalize_total(adata, target_sum=1e4, layer='spatpy_deg_count')
    sc.pp.log1p(adata, layer='spatpy_deg_count')
    celltype_counts = adata.obs[label_field].value_counts()
    valid_celltypes = celltype_counts[celltype_counts >= 2].index.tolist()
    adata_filtered = adata[adata.obs[label_field].isin(valid_celltypes)].copy()
    sc.tl.rank_genes_groups(
        adata_filtered,
        groupby=label_field,
        layer='spatpy_deg_count', 
        method='t-test',
        n_genes=2000,
        use_raw=False
    )
    adata_ori.uns['rank_genes_groups'] = adata_filtered.uns['rank_genes_groups'].copy()
    deg_results = adata_filtered.uns['rank_genes_groups']
    celltypes = sorted(deg_results['names'].dtype.names)
    genes = []
    for ct in celltypes:
        top3_genes = adata_filtered.uns['rank_genes_groups']['names'][ct][:3]
        genes.extend(top3_genes)

    def get_mean_expression(celltype, gene):
        adtmp = adata_filtered[adata_filtered.obs[label_field] == celltype, gene]
        mean_exp = np.mean(adtmp.layers['spatpy_deg_count'])
        return mean_exp

    heatmap_data = [[get_mean_expression(cell, gene) for cell in celltypes]for gene in genes]
    def minmax_scale(row):
        return (row - np.min(row)) / (np.max(row) - np.min(row))
    minmax_data = np.apply_along_axis(minmax_scale, axis=1, arr=heatmap_data)
    fig = px.imshow(
        minmax_data,
        labels=dict(x="celltype", y="gene", color="mean_expression"),
        x=celltypes,
        y=genes,
        color_continuous_scale='Cividis',
        aspect='auto'
    )
    fig.update_layout(
        height=800,
        width=800,
        autosize=True,
        title=dict(
            text='Differential Gene Heatmap',
            x=0.5,
        ),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        coloraxis_colorbar_title_text='',
        xaxis_title=None,
        yaxis_title=None,
    )
    fig.update_xaxes(
        tickangle=45, 
        tickfont=dict(size=10)
    )
    fig.update_yaxes(
        tickfont=dict(size=10)
    )
    return fig

def umap(ad):
    sc.pp.pca(ad)
    sc.pp.neighbors(ad)
    sc.tl.umap(ad)
    sc.tl.leiden(ad)

def set_steps_status(steps:dict, step:int, field:str, status:bool)->None:
    if steps is None:
        return
    steps[step][field] = status  
  
def use_raw(adata):
    """
    使用原始计数
    """
    if 'counts' in adata.layers:
        adata.X = adata.layers["counts"].copy()

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

def normalize(adata):
    """
    Normalize the data
    Args:
        adata: AnnData object
    Returns:
        None
    """
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

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