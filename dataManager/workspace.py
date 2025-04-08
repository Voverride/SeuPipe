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
        annotation_path = os.path.join(workspace_path, 'Annotation')
        os.makedirs(annotation_path)

    table.update({'workspace':workspace_path})
    
def get_workspace()->str:
    """
    Get the current workspace path.

    Returns:
        str: The current workspace path.
    """
    workspace = table.all()[0]['workspace']
    if workspace is not None:
        if not os.path.exists(workspace):
            os.makedirs(workspace)
    return workspace

def get_annotation_workspace()->str:
    return os.path.join(get_workspace(), 'Annotation')