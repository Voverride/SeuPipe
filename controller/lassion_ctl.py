from dataManager.lassion_d import lasData
from dash import set_props
import numpy as np
from utils.commonfuc import *

def select_cell(pointIndex, dataset):
    """
    筛选细胞
    """
    df = lasData.get_dataset_data(dataset)
    if pointIndex is not None:
        pointIndex = np.array(pointIndex)
        df['selected'] = False
        df['selected'] = df.index.isin(set(pointIndex.flatten()))
    set_props('lasso-download-exportData', dict(data={
        'filename': f'{dataset}.csv',
        'content': df.to_csv()
    }))

def update_relayout(patch, relayoutData):
    """
    更新视图
    """
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