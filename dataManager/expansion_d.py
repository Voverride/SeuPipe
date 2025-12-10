from dataManager.workspace import *
from controller.auth import get_host
from dataManager.users import search_user
from utils.commonfuc import *
from dash import Patch
from dataManager.segmentation_d import segData
import shutil
import pickle

class ExpansionData:
    def __init__(self):
        self._runningTask = {} # {key:taskname, value:taskinfo}
        self._userVirtualTask = {} # {key:username, value:{key:taskname, value:taskinfo}}
        self._temptask = {} #{key:username, value:taskinfo}
        # {
        #     'taskname': None,
        #     'thread':None,
        #     'exception':None,
        #     'metadata':None,
        #     {
        #         'creator': 'zhouyb',
        #         'model': 'cyto3',
        #         'diameter': 10,
        #         'batchsize': 8,
        #         'GPU': 'True',
        #         'progress': 0
        #     }
        #     'data': None
        #     [
        #         {
        #             'image': 'image2',
        #             'gem': 'gem1',
        #             'z': 1,
        #             'registration': {'status': 'success', 'text': 'success'},
        #             'segmentation': {'status': 'success', 'text': 'success'},
        #         }
        #     ],
        # }

    def create_running_task(self, taskInfo):
        """
        创建正在运行的任务
        """
        self._runningTask[taskInfo['taskname']] = taskInfo
    def has_running_task(self, taskName:str):
        """
        检查当前任务是否正在运行
        """
        return taskName in self._runningTask
    
    def get_request_usrname(self):
        """
        获取请求的用户名
        """
        host = get_host()
        user = search_user(usrhost=host)
        if len(user)==0:
            return None
        username = user[0]['usrname']
        if username not in self._temptask:
            self._temptask[username] = {}
        return username

    def delete_task(self, taskName):
        """
        从磁盘删除任务
        """
        task_folder = os.path.join(get_expansion_workspace(), taskName)
        if os.path.exists(task_folder):
            shutil.rmtree(task_folder)

    def read_taskinfo(self, taskName:str):
        """
        读取用户创建的任务清单
        """
        try:
            root_path = get_expansion_workspace()
            json_path = os.path.join(root_path, taskName, 'metadata.json')
        except Exception as e:
            return None
        if os.path.exists(json_path):
            return read_json(json_path)
        else:
            return None

    def save_expansion_task(self, userName, taskName, mode, patchsize, binsize, epochs, diameter, neighbor):
        """
        保存用户创建的任务清单
        """
        exp_workspace = get_expansion_workspace()
        task_folder = os.path.join(exp_workspace, taskName)
        check_path(task_folder)
        segTask = segData.read_taskinfo(taskName)
        slices = [
            {
                'z': task['z'],
                'registration': task['registration'],
                'segmentation': task['segmentation'],
                'preprocess': {'status': 'warning', 'text': 'waiting'},
                'train': {'status': 'warning', 'text': 'waiting'},
                'postprocess': {'status': 'warning', 'text': 'waiting'},
            }
            for task in segTask['data']
        ]
        metadata = {
            'creator': userName,
            'mode': mode,
            'patchSize': patchsize,
            'binSize': binsize,
            'epochs': epochs,
            'diameter': diameter,
            'neighbors': neighbor,
            'progress': 0,
            'slices': slices
        }
        output_path = os.path.join(task_folder, 'metadata.json')
        write_json(output_path, metadata)

    def get_exist_tasks(self):
        """
        获取已创建的任务列表
        """
        exp_path = get_expansion_workspace()
        tasks = os.listdir(exp_path)
        tasks.sort()
        return [task for task in tasks if os.path.isdir(os.path.join(exp_path, task))]

expData = ExpansionData()