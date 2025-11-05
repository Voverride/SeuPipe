from dataManager.workspace import *
from controller.auth import get_host
from dataManager.users import search_user
import shutil
import pickle

class RegionClipData:
    def __init__(self):
        self._temptask = {} #{key:username, value:taskinfo}
        taskinfo = {
            'taskname': 'taskname',
            'metadata': {'Creator', 'Date'},
            'data': {
                'z': {'z', 'img', 'gem'}
            },
        }

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
    
    def get_slice_info(self, taskName, slicename):
        """
        获取任务信息
        """
        taskInfo = self.read_taskinfo(taskName)
        if taskInfo is None:
            return None
        sliceInfo = taskInfo.get('data', {}).get(slicename)
        return sliceInfo
    
    def get_task_slices(self, taskName):
        """
        获取任务切片列表
        """
        taskInfo = self.read_taskinfo(taskName)
        if taskInfo is None:
            return []
        slices = list(taskInfo['data'].keys())
        return slices
    
    def get_taskclips_folder(self, taskName):
        """
        获取任务裁剪项目目录
        """
        if not taskName:
            return None
        task_folder = os.path.join(get_regionclip_workspace(), taskName)
        clips_folder = os.path.join(task_folder, 'Clips')
        check_path(clips_folder)
        return clips_folder

    def get_taskclips(self, taskName):
        """
        获取当前任务下的所有裁剪项目
        """
        if not taskName:
            return None
        clips_folder = self.get_taskclips_folder(taskName)
        if not clips_folder:
            return None
        clips = os.listdir(clips_folder)
        return clips

    def create_taskclip(self, taskName, clipName):
        """
        在当前任务下创建裁剪项目
        """
        if not taskName or not clipName:
            return None
        clips_folder = self.get_taskclips_folder(taskName)
        if not clips_folder:
            return None
        clip_path = os.path.join(clips_folder, clipName)
        if not os.path.exists(clip_path):
            os.makedirs(clip_path)

    def delete_task(self, taskName):
        """
        从磁盘删除任务
        """
        if not taskName:
            return
        task_folder = os.path.join(get_regionclip_workspace(), taskName)
        if os.path.exists(task_folder):
            shutil.rmtree(task_folder)

    def read_taskinfo(self, taskName:str):
        """
        读取用户创建的任务清单
        """
        try:
            root_path = get_regionclip_workspace()
            pkl_path = os.path.join(root_path, taskName, 'tasklist.pkl')
        except Exception as e:
            return None
        if os.path.exists(pkl_path):
            with open(pkl_path, 'rb') as f:
               tasklist = pickle.load(f)
               return tasklist
        else:
            return None

    def save_temptask(self):
        """
        保存用户创建的任务清单
        """
        username = self.get_request_usrname()
        taskname = self._temptask[username]['taskname']
        clip_workspace = get_regionclip_workspace()
        task_folder = os.path.join(clip_workspace, taskname)
        check_path(task_folder)
        output_path = os.path.join(task_folder, 'tasklist.pkl')
        with open(output_path, 'wb') as f:
            pickle.dump(self._temptask[username], f)
        del self._temptask[username]

    def set_temptask_metadata(self, taskname:str, metadate:dict):
        """
        设置进度和元数据
        """
        username = self.get_request_usrname()
        self._temptask[username]['taskname'] = taskname
        self._temptask[username]['metadata'] = metadate
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
        data_dict = {d['z']:d for d in data}
        self._temptask[username]['data'] = data_dict

    def get_exist_tasks(self):
        """
        获取已创建的任务列表
        """
        clip_path = get_regionclip_workspace()
        tasks = os.listdir(clip_path)
        return {task for task in tasks if os.path.isdir(os.path.join(clip_path, task))}

clipData = RegionClipData()