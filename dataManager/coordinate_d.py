import pandas as pd
from dataManager.workspace import *
from utils.commonfuc import *
import cachetools
import os
import time


ali_coordinate = {}
ali_coordinateTime = {}
ali_initfig = cachetools.TTLCache(maxsize=200, ttl=3600)

def get_figure_folder(project):
    """
    获取对齐项目可视化文件夹路径
    """
    project_folder = get_project_folder(project)
    figure_folder = os.path.join(project_folder, 'figure')
    check_path(figure_folder)
    return figure_folder

def get_initialfig_path(project):
    """
    获取初始数据可视化结果路径
    """
    figure_folder = get_figure_folder(project)
    initialfig_path = os.path.join(figure_folder, 'initialfig.pkl')
    return initialfig_path

def get_initialfig(project):
    """
    获取初始数据可视化结果
    """
    if project in ali_initfig:
        initialfig = ali_initfig[project]
    else:
        initialfig_path = get_initialfig_path(project)
        if os.path.exists(initialfig_path):
            try:
                initialfig = read_pkl(initialfig_path)
            except:
                initialfig = None
        else:
            initialfig = None
    ali_initfig[project] = initialfig
    return initialfig

def get_project_folder(project_name):
    """
    获取对齐项目文件夹路径
    """
    ali_path = get_alignment_workspace()
    project_path = os.path.join(ali_path, project_name)
    check_path(project_path)
    return project_path

def get_coordinate_folder(project):
    """
    获取对齐项目坐标文件夹路径
    """
    project_folder = get_project_folder(project)
    coordinate_folder = os.path.join(project_folder, 'coordinate')
    check_path(coordinate_folder)
    return coordinate_folder

def get_coordinate_path(project):
    """
    获取对齐项目坐标路径
    """
    coordinate_folder = get_coordinate_folder(project)
    coordinate_path = os.path.join(coordinate_folder, 'coordinate.csv')
    return coordinate_path

def get_coordinate(project):
    """
    获取对齐项目坐标
    """
    if project in ali_coordinate:
        coordinate = ali_coordinate[project]
    else:
        coordinate_path = get_coordinate_path(project)
        coordinate = pd.read_csv(coordinate_path, index_col=0)
    ali_coordinate[project] = coordinate
    ali_coordinateTime[project] = time.time()
    return coordinate
