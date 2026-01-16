from dataManager.workspace import *
from controller.auth import get_host
from dataManager.users import search_user
from utils.commonfuc import *
from dataManager.segmentation_d import segData
from dataManager.maskviewer_d import maskData
from utils.typing import Status
import shutil
from utils.observer import observer

class ExpansionData:
    def __init__(self):
        self._runningTask = {} # {key:taskname, value:taskinfo}

    def get_running_task(self, taskName):
        """
        获取正在运行的任务
        """
        if taskName in self._runningTask:
            return self._runningTask[taskName]
        return None
    def create_running_task(self, taskName, callback):
        """
        创建正在运行的任务
        """
        taskInfo = self.read_taskinfo(taskName)
        taskInfo['running'] = True
        taskInfo['exception'] = None
        slices = taskInfo['slices']
        patchprocess_text = '0/0'
        cnt = 0
        for slice in slices:
            expansion_path = maskData.get_expansion_path(taskName, slice['z'])
            status = Status.WARNING
            text = 'waiting'
            patchprocess_text = '1/1'
            if os.path.exists(expansion_path):
                status = Status.SUCCESS
                text = 'success'
                cnt+=1
            slice['preprocess']['status'] = status
            slice['preprocess']['text'] = text
            slice['train']['status'] = status
            slice['train']['text'] = text
            slice['postprocess']['status'] = status
            slice['postprocess']['text'] = text
            slice['patchprocess']['status'] = status
            slice['patchprocess']['text'] = patchprocess_text
        
        taskInfo['progress'] = cnt/len(slices)

        boserved_task = observer.observe(taskInfo, callback)
        self._runningTask[taskName] = boserved_task

    def delete_running_task(self, taskName, callback):
        """
        删除正在运行的任务
        """
        if taskName in self._runningTask:
            runningTask = self._runningTask[taskName]
            observer.disobserve(runningTask, callback)
            del self._runningTask[taskName]
    
    def is_task_running(self, taskName):
        """
        判断任务是否正在运行
        """
        if taskName in self._runningTask:
            taskInfo = self._runningTask[taskName]
            return taskInfo.get('running', False)
        return False
    
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
        slices = maskData.get_project_slices(taskName)
        for slice in slices:
            expansion_result_path = maskData.get_expansion_path(taskName, slice)
            if os.path.exists(expansion_result_path):
                os.remove(expansion_result_path)

    def read_taskinfo(self, taskName:str):
        """
        读取用户创建的任务清单
        """
        if taskName in self._runningTask:
            taskInfo = dict(self._runningTask[taskName])
            if taskInfo:
                return taskInfo
        try:
            root_path = get_expansion_workspace()
            json_path = os.path.join(root_path, taskName, 'metadata.json')
        except Exception as e:
            return None
        if os.path.exists(json_path):
            return read_json(json_path)
        else:
            return None
        
    def get_expTask_folder(self, taskName):
        """
        获取扩增任务文件夹
        """
        task_folder = os.path.join(get_expansion_workspace(), taskName)
        check_path(task_folder)
        return task_folder
    
    def get_expTaskSlices_folder(self, taskName):
        """
        获取扩增任务切片文件夹
        """
        task_slices_folder = os.path.join(get_expansion_workspace(), taskName, 'slices')
        check_path(task_slices_folder)
        return task_slices_folder
    
    def get_expTaskSlices_zfolder(self, taskName, zIndex):
        """
        获取扩增任务切片SCS输出指定z切片文件夹
        """
        task_zfolder = os.path.join(get_expansion_workspace(), taskName, 'slices', f'z_{zIndex}')
        check_path(task_zfolder)
        return task_zfolder
    
    def get_expTask_result_folder(self, taskName, zIndex):
        """
        获取扩增任务结果文件夹
        """
        task_result_folder = os.path.join(get_expansion_workspace(), taskName, 'slices', f'z_{zIndex}', 'result', 'cell_mask')
        check_path(task_result_folder)
        return task_result_folder
    
    def get_expTaskMetadata_path(self, taskName):
        """
        获取扩增任务元数据文件路径
        """
        task_metadata_path = os.path.join(get_expansion_workspace(), taskName, 'metadata.json')
        return task_metadata_path
        
    def save_expansion_task(self, userName, taskName, mode, patchsize, binsize, epochs, diameter, neighbor):
        """
        保存用户创建的任务清单
        """
        exp_workspace = get_expansion_workspace()
        task_folder = os.path.join(exp_workspace, taskName)
        check_path(task_folder)
        slices_folder = os.path.join(task_folder, 'slices')
        check_path(slices_folder)
        segTask = segData.get_project_info(taskName)
        slices = [
            {
                'z': task['z'],
                'img': task['img'],
                'gem': task['gem'],
                'segmentation': task['segmentation'],
                'seg_postprocess': task['postprocess'],
                'preprocess': {'status': 'warning', 'text': 'waiting'},
                'train': {'status': 'warning', 'text': 'waiting'},
                'postprocess': {'status': 'warning', 'text': 'waiting'},
                'patchprocess': {'status': 'warning', 'text': '0/0'},
            }
            for task in segTask['slices']
        ]
        metadata = {
            'running': False,
            'taskName': taskName,
            'creator': userName,
            'mode': mode,
            'patchSize': patchsize,
            'binSize': binsize,
            'epochs': epochs,
            'diameter': diameter,
            'neighbors': neighbor,
            'progress': 0,
            'slices': slices,
            'exception': None
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