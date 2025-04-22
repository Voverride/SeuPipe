from tinydb import TinyDB
from tinydb import Query
import os

directory = 'database/'
if not os.path.exists(directory):
    os.makedirs(directory)
    
db = TinyDB('database/db.json')
table_name = 'workspace'
table = db.table(table_name)
if table_name not in db.tables():
    table.insert({'workspace':None})
def set_workpase(workspace_path:str)->None:
    """
    Set workspace path
    Args:
        workspace_path(str): workspace path
    Returns:
        None
    """
    workspace_path = os.path.join(workspace_path, 'SeuPipeWorkspace')
    if not os.path.exists(workspace_path):
        os.makedirs(workspace_path)

    table.update({'workspace':workspace_path})
    
def get_workspace()->str:
    """
    Get the current workspace path.

    Returns:
        str: The current workspace path.
    """
    workspace = table.all()[0]['workspace']
    if workspace is not None:
        check_path(workspace)
    return workspace
def get_annotation_workspace()->str:
    path = os.path.join(get_workspace(), 'Annotation')
    check_path(path)
    return path

def get_segmentation_workspace()->str:
    path = os.path.join(get_workspace(), 'Segmentation')
    check_path(path)
    return path

def get_cache_workspace()->str:
    path = os.path.join(get_workspace(), 'Cache')
    check_path(path)
    return path

def check_path(path):
    if not os.path.exists(path):
        os.makedirs(path)