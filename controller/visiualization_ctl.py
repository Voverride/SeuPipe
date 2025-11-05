from dataManager.visiualization_d import visData
from utils.commonfuc import *


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