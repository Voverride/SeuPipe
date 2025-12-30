from dataManager.workspace import *
from controller.auth import get_host
from dataManager.users import search_user
import shutil
from websocket.websocket import ws
import numpy as np
import pandas as pd
import imageio.v3 as iio
import pickle
from utils.observer import observer
import cachetools
from functools import partial
from utils.commonfuc import compress_image, write_json, read_json

class RegionClipData:
    def __init__(self):
        self._temptask = {} #{key:username, value:taskinfo}
        # taskinfo = {
        #     'taskname': 'taskname',
        #     'metadata': {'Creator', 'Date'},
        #     'data': {
        #         'z': {'z', 'image', 'gem'}
        #     },
        # }
        self._gemCache = cachetools.TTLCache(maxsize=10, ttl=1800)
        self._imgCache = cachetools.TTLCache(maxsize=10, ttl=1800)
        self._cordMapFunc = {}
        self._runningTask = {}
        self._max_threads = 5
        self._threads = {}
        # threads = {
        #     'username': {
        #         'active': {
        #             'key': 'thread'
        #         },
        #         'pending': {
        #             'key': 'thread'
        #         }
        #     }
        # }
        self._notifyfunc = 'notifyfunc'
        
    def add_thread(self, key, thread):
        """
        添加线程
        """
        username = self.get_request_usrname()
        if username not in self._threads:
            self._threads[username] = {
                'active': {},
                'pending': {}
            }
        active = self._threads[username]['active']
        pending = self._threads[username]['pending']
        if len(active) < self._max_threads:
            thread.start()
            active[key] = thread
        else:
            pending[key] = thread

    def add_running_task(self, taskName, slice, clipName):
        """
        添加运行中的任务
        """
        key = f'{taskName}_{slice}_{clipName}'
        username = self.get_request_usrname()
        status = self.read_taskclip_status(taskName, clipName, slice)
        status['operator'] = username
        status['running'] = True
        status[self._notifyfunc] = ws.notifyUpdateRegionClipStatus
        observed_status = observer.observe(
            status, 
            status[self._notifyfunc], 
            taskName=taskName,
            slice=slice,
            clipName=clipName
        )
        self._runningTask[key] = observed_status
        self._runningTask[key]['exception'] = None
    
    def remove_running_task(self, taskName, slice, clipName):
        """
        删除运行中的任务
        """
        key = f'{taskName}_{slice}_{clipName}'
        if key not in self._runningTask:
            return
        observed_status = self._runningTask[key]
        username = observed_status['operator']
        notify_func = observed_status[self._notifyfunc]
        active = self._threads[username]['active']
        pending = self._threads[username]['pending']
        observer.disobserve(observed_status, notify_func)
        del self._runningTask[key]
        if key in active:
            del active[key]
            if len(active) < self._max_threads:
                for key, thread in pending.items():
                    if key in active:
                        continue
                    thread.start()
                    active[key] = thread
                    del pending[key]
                    if len(active) >= self._max_threads:
                        break

    def has_running_task(self, taskName, slice, clipName):
        """
        判断任务是否正在运行
        """
        key = f'{taskName}_{slice}_{clipName}'
        if key not in self._runningTask:
            return False
        return self._runningTask[key]['running']

    def set_cord_map_func(self, taskName, sliceName, func):
        """
        设置用户获取坐标映射函数
        """
        key = f'{taskName}_{sliceName}'
        self._cordMapFunc[key] = func

    def get_cord_map_func(self, taskName, slicename):
        """
        用户获取坐标映射函数
        """
        key = f'{taskName}_{slicename}'
        return self._cordMapFunc.get(key, None)
    
    def get_imgCache(self, taskName, slicename):
        """
        获取任务图像缓存
        """
        img_key = f'{taskName}_{slicename}'
        if img_key in self._imgCache:
            img = self._imgCache[img_key]
            # 重置缓存计时
            del self._imgCache[img_key]
        else:
            taskInfo = self.get_slice_info(taskName, slicename)
            if taskInfo is None:
                return None
            img_path = taskInfo.get('image', None)
            if img_path is None:
                return None
            img = iio.imread(img_path)
        self._imgCache[img_key] = img
        return img

    def get_gemCache(self, taskName, slicename, cache=False):
        """
        获取任务gem缓存
        """
        gem_key = f'{taskName}_{slicename}'
        if cache and gem_key in self._gemCache:
            gem = self._gemCache[gem_key]
            # 重置缓存计时
            del self._gemCache[gem_key]
        else:
            taskInfo = self.get_slice_info(taskName, slicename)
            if taskInfo is None:
                return None
            gem_path = taskInfo.get('gem', None)
            if gem_path is None:
                return None
            gem = pd.read_csv(gem_path, sep='\t', comment='#')
            countField = 'MIDCounts' if 'MIDCounts' in gem else 'MIDCount'
            gem = gem[['geneID', 'x', 'y', countField]].rename(columns={countField: 'MIDCount'})
        if cache:
            self._gemCache[gem_key] = gem
        return gem

    def get_stain_img(self, taskName, slicename):
        """
        获取任务指定切片的染色压缩图像
        """
        img = self.get_imgCache(taskName, slicename)
        compressd_img, cord_map_func = compress_image(img)
        self.set_cord_map_func(taskName, slicename, cord_map_func)
        return compressd_img
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
        slices.sort()
        return slices
    def get_output_imgName(self, taskName, sliceName, clipName):
        """
        获取输出图像文件名称
        """
        imgName = f'{taskName}_z{sliceName}_{clipName}.tif'
        return imgName
    
    def get_output_gemName(self, taskName, sliceName, clipName):
        """
        获取输出图像文件名称
        """
        imgName = f'{taskName}_z{sliceName}_{clipName}.gem'
        return imgName
    def get_output_gemImgName(self, taskName, sliceName, clipName):
        """
        获取输出图像文件名称
        """
        imgName = f'{taskName}_z{sliceName}_{clipName}.png'
        return imgName
    def get_task_clipName_gemImage_path(self, taskName, sliceName, clipName):
        """
        获取任务裁剪项目下基因信息缩略图路径
        """
        if not taskName or not clipName:
            return None
        clip_gem_folder = self.get_task_clipName_gem_folder(taskName, clipName)
        gemImgName = self.get_output_gemImgName(taskName, sliceName, clipName)
        gem_image_path = os.path.join(clip_gem_folder, gemImgName)
        if not os.path.exists(gem_image_path):
            return None
        return gem_image_path
    
    def get_task_clipName_gem_path(self, taskName, sliceName, clipName):
        """
        获取任务裁剪项目下基因信息路径
        """
        if not taskName or not clipName:
            return None
        clip_gem_folder = self.get_task_clipName_gem_folder(taskName, clipName)
        gemName = self.get_output_gemName(taskName, sliceName, clipName)
        gem_path = os.path.join(clip_gem_folder, gemName)
        if not os.path.exists(gem_path):
            return None
        return gem_path
    
    def get_task_clipName_stain_path(self, taskName, sliceName, clipName):
        """
        获取任务裁剪项目下染色图像路径
        """
        if not taskName or not clipName:
            return None
        clip_stain_folder = self.get_task_clipName_stain_folder(taskName, clipName)
        imgName = self.get_output_imgName(taskName, sliceName, clipName)
        stain_path = os.path.join(clip_stain_folder, imgName)
        if not os.path.exists(stain_path):
            return None
        return stain_path
    
    def get_task_clipName_gem_folder(self, taskName, clipName):
        """
        获取任务裁剪项目下基因信息目录
        """
        if not taskName or not clipName:
            return None
        clip_folder = self.get_task_clipName_folder(taskName, clipName)
        clip_gem_folder = os.path.join(clip_folder, 'GEM')
        check_path(clip_gem_folder)
        return clip_gem_folder

    def get_task_clipName_stain_folder(self, taskName, clipName):
        """
        获取任务裁剪项目下染色图像目录
        """
        if not taskName or not clipName:
            return None
        clip_folder = self.get_task_clipName_folder(taskName, clipName)
        clip_stain_folder = os.path.join(clip_folder, 'Stain')
        check_path(clip_stain_folder)
        return clip_stain_folder
    def get_task_clipName_folder(self, taskName, clipName):
        """
        获取任务裁剪项目目录
        """
        if not taskName or not clipName:
            return None
        clips_folder = self.get_taskclips_folder(taskName)
        clip_folder = os.path.join(clips_folder, str(clipName))
        check_path(clip_folder)
        return clip_folder    
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
        clips.sort()
        return clips
    
    def delete_taskclip(self, taskName, clipName):
        """
        在当前任务下删除裁剪项目
        """
        if not taskName or not clipName:
            return
        clips_folder = self.get_taskclips_folder(taskName)
        if not clips_folder:
            return
        clip_path = os.path.join(clips_folder, clipName)
        if os.path.exists(clip_path):
            shutil.rmtree(clip_path)

    def write_taskclip_status(self, taskName, clipName, slice, status):
        """
        写入任务裁剪项目切片状态
        """
        status_folder = self.get_task_clipName_status_folder(taskName, clipName)
        status_path = os.path.join(status_folder, f'z_{slice}.json')
        write_json(os.path.join(status_path), status)

    def get_running_task_status(self, taskName, clipName, slice):
        """
        获取正在运行中的任务状态
        """
        key = f'{taskName}_{slice}_{clipName}'
        if key in self._runningTask:
            return self._runningTask[key]
        return None

    def read_taskclip_status(self, taskName, clipName, slice):
        """
        读取任务裁剪项目切片状态
        """
        status = self.get_running_task_status(taskName, clipName, slice)
        if status is not None:
            status = dict(status)
            if self._notifyfunc in status:
                del status[self._notifyfunc]
            return status
        status_folder = self.get_task_clipName_status_folder(taskName, clipName)
        status_path = os.path.join(status_folder, f'z_{slice}.json')
        if not os.path.exists(status_path):
            return None
        return read_json(status_path)

    def get_task_clipName_status_folder(self, taskName, clipName):
        """
        获取任务裁剪项目状态目录
        """
        if not taskName or not clipName:
            return None
        clip_status_folder = os.path.join(self.get_task_clipName_folder(taskName, clipName), 'Status')
        check_path(clip_status_folder)
        return clip_status_folder

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
        slices = self.get_task_slices(taskName)
        init_status = {
            'running': False,
            'exception': None
        }
        for slice in slices:
            self.write_taskclip_status(taskName, clipName, slice, init_status)
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
        设置元数据
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