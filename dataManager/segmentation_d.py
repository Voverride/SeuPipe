from dataManager.workspace import *
from controller.auth import get_host
from dataManager.users import search_user
from dash import Patch
import shutil
import pickle

class SegmentationData:
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

    def get_diff(self, taskName:str):
        """
        获取虚拟dom状态与当前执行任务状态的差异
        """
        patch_metadata = None
        patch_data = None
        error = None
        def different(virtualTask, runningTask):
            pmetadata = None
            pdata = None
            error_message = None
            virtual_pct = virtualTask['metadata']['progress']
            running_pct = runningTask['metadata']['progress']
            if virtual_pct != running_pct:
                virtualTask['metadata']['progress'] = running_pct
                pmetadata = Patch()
                pmetadata[0]['progress'] = running_pct
            if runningTask['exception'] is not None:
                virtualTask['exception'] = runningTask['exception']
                error_message = runningTask['exception']
            virtual_data = virtualTask['data']
            running_data = runningTask['data']
            for idx, running_item in enumerate(running_data):
                virReg = virtual_data[idx]['registration']
                runReg = running_item['registration']
                virSeg = virtual_data[idx]['segmentation']
                runSeg = running_item['segmentation']
                if virReg['status'] != runReg['status']:
                    virReg['status'] = runReg['status']
                    virReg['text'] = runReg['text']
                    if pdata is None:
                        pdata = Patch()
                    pdata[idx]['registration'] = virReg
                if virSeg['status'] != runSeg['status']:
                    virSeg['status'] = runSeg['status']
                    virSeg['text'] = runSeg['text']
                    if pdata is None:
                        pdata = Patch()
                    pdata[idx]['segmentation'] = virSeg
            return pmetadata, pdata, error_message

        if taskName in self._runningTask:
            runningTask = self._runningTask[taskName]
            userVirtualTask = self.get_userVirtualTask(taskName)
            patch_metadata, patch_data, error = different(userVirtualTask, runningTask)
            
            if runningTask['thread'] is not None and not runningTask['thread'].is_alive():
                runningTask['thread'] = None
                seg_workspace = get_segmentation_workspace()
                task_folder = os.path.join(seg_workspace, taskName)
                check_path(task_folder)
                output_path = os.path.join(task_folder, 'tasklist.pkl')
                with open(output_path, 'wb') as f:
                    pickle.dump(runningTask, f)
                del self._runningTask[taskName]
        
        elif self.has_userVirtualTask(taskName):
            runningTask = self.read_taskinfo(taskName)
            userVirtualTask = self.get_userVirtualTask(taskName)
            patch_metadata, patch_data, error = different(userVirtualTask, runningTask)
            self.remove_userVirtualTask(taskName)
        
        return patch_metadata, patch_data, error

    def has_running_task(self, taskName:str):
        """
        检查当前任务是否正在运行
        """
        return taskName in self._runningTask
    def has_userVirtualTask(self, taskName:str):
        """
        检查当前用户是否存在尚未更新完的虚拟dom状态
        """
        username = self.get_request_usrname()
        if username not in self._userVirtualTask:
            return False
        if taskName not in self._userVirtualTask[username]:
            return False
        return True
    
    def remove_userVirtualTask(self, taskName:str):
        """
        删除用户界面虚拟dom状态
        """
        username = self.get_request_usrname()
        if username in self._userVirtualTask:
            if taskName in self._userVirtualTask[username]:
                del self._userVirtualTask[username][taskName]
    def get_userVirtualTask(self, taskName:str):
        """
        获取用户界面虚拟dom状态
        """
        username = self.get_request_usrname()
        if username not in self._userVirtualTask:
            self._userVirtualTask[username] = {}
        if taskName not in self._userVirtualTask[username]:
            self._userVirtualTask[username][taskName] = self.read_taskinfo(taskName)
        return self._userVirtualTask[username][taskName]

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
        task_folder = os.path.join(get_segmentation_workspace(), taskName)
        if os.path.exists(task_folder):
            shutil.rmtree(task_folder)

    def read_taskinfo(self, taskName:str):
        """
        读取用户创建的任务清单
        """
        try:
            root_path = get_segmentation_workspace()
            pkl_path = os.path.join(root_path, taskName, 'tasklist.pkl')
        except Exception as e:
            return None
        if os.path.exists(pkl_path):
            with open(pkl_path, 'rb') as f:
               tasklist = pickle.load(f)
               return tasklist
        else:
            return None
        
    def reset_taskinfo(self, taskInfo:str):
        """
        重置任务状态
        """
        taskInfo['thread']=None
        taskInfo['exception']=None
        taskInfo['metadata']['progress']=0
        for item in taskInfo['data']:
            item['registration']['status']='warning'
            item['segmentation']['status']='warning'
            item['registration']['text']='waiting'
            item['segmentation']['text']='waiting'

        
    def save_taskinfo(self, taskInfo:dict):
        """
        保存任务任务状态
        """
        seg_workspace = get_segmentation_workspace()
        taskName = taskInfo['taskname']
        task_folder = os.path.join(seg_workspace, taskName)
        check_path(task_folder)
        pkl_path = os.path.join(task_folder, 'tasklist.pkl')
        with open(pkl_path, 'wb') as f:
            pickle.dump(taskInfo, f)

    def save_temptask(self):
        """
        保存用户创建的任务清单
        """
        username = self.get_request_usrname()
        taskname = self._temptask[username]['taskname']
        seg_workspace = get_segmentation_workspace()
        task_folder = os.path.join(seg_workspace, taskname)
        check_path(task_folder)
        output_path = os.path.join(task_folder, 'tasklist.pkl')
        with open(output_path, 'wb') as f:
            pickle.dump(self._temptask[username], f)
    def reset_temptask_data(self):
        """
        重置用户创建的任务清单
        """
        username = self.get_request_usrname()
        self._temptask[username]['data'] = None

    def set_temptask_data(self, data:list):
        """
        设置用户创建的任务清单
        """
        username = self.get_request_usrname()
        self._temptask[username]['data'] = data

    def set_temptask_metadata(self, taskname:str, metadate:dict):
        """
        设置进度和元数据
        """
        username = self.get_request_usrname()
        self._temptask[username]['taskname'] = taskname
        self._temptask[username]['thread'] = None
        self._temptask[username]['exception'] = None
        self._temptask[username]['metadata'] = metadate
    def get_exist_tasks(self):
        """
        获取已创建的任务列表
        """
        seg_path = get_segmentation_workspace()
        tasks = os.listdir(seg_path)
        return {task for task in tasks if os.path.isdir(os.path.join(seg_path, task))}

segData = SegmentationData()