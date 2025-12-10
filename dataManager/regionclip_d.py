from dataManager.workspace import *
from controller.auth import get_host
from dataManager.users import search_user
import shutil
import numpy as np
import pandas as pd
import imageio.v3 as iio
import pickle

class RegionClipData:
    def __init__(self):
        self._temptask = {} #{key:username, value:taskinfo}
        taskinfo = {
            'taskname': 'taskname',
            'metadata': {'Creator', 'Date'},
            'data': {
                'z': {'z', 'image', 'gem'}
            },
        }

    def clip_image(self, taskName, slicename, clipName, x_start, x_end, y_start, y_end):
        """
        裁剪染色图像和gem, 生成裁剪后gem热图缩略图, 并保存到指定目录
        """
        slice_info = self.get_slice_info(taskName, slicename)
        image_path = slice_info.get('image')
        gem_path = slice_info.get('gem')
        if not image_path or not gem_path:
            return

        img = iio.imread(image_path)
        gem = pd.read_csv(gem_path, sep='\t', comment='#')
        gem['x'] = gem['x'] - gem['x'].min()
        gem['y'] = gem['y'] - gem['y'].min()
        clipped_img = img[y_start:y_end, x_start:x_end]
        clipped_gem = gem[
            gem['y'].between(x_start, x_end) & 
            gem['x'].between(y_start, y_end)
        ].copy()

        clipped_image_folder = self.get_task_clipName_stain_folder(taskName, clipName)
        clipped_gem_folder = self.get_task_clipName_gem_folder(taskName, clipName)
        clipped_image_path = os.path.join(clipped_image_folder, f'{taskName}_{slicename}_{clipName}.tif')
        clipped_gem_path = os.path.join(clipped_gem_folder, f'{taskName}_{slicename}_{clipName}.gem')
        clipped_gem_image_path = os.path.join(clipped_gem_folder, f'{taskName}_{slicename}_{clipName}.png')

        iio.imwrite(clipped_image_path, clipped_img)
        clipped_gem.to_csv(clipped_gem_path, sep='\t', index=False)

        if clipped_gem.empty:
            mtx = np.zeros((100, 100), dtype=np.uint8)
        else:
            grouped = clipped_gem.groupby(['x', 'y'])['MIDCounts'].sum().reset_index()

            x_raw = grouped['x'].values
            y_raw = grouped['y'].values
            counts = grouped['MIDCounts'].values

            x_norm = (x_raw - x_raw.min()).astype(int)
            y_norm = (y_raw - y_raw.min()).astype(int)

            max_x = x_norm.max()
            max_y = y_norm.max()

            mtx = np.zeros((max_x + 1, max_y + 1), dtype=np.float32)
            mtx[x_norm, y_norm] = counts

            mtx = np.log1p(mtx)

            max_val = mtx.max()
            if max_val > 0:
                mtx = (mtx / max_val * 255).astype(np.uint8)
            else:
                mtx = mtx.astype(np.uint8)

        iio.imwrite(clipped_gem_image_path, mtx, cmap='hot')
        
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
    
    def get_task_clipName_gemImage_path(self, taskName, sliceName, clipName):
        """
        获取任务裁剪项目下基因信息缩略图路径
        """
        if not taskName or not clipName:
            return None
        clip_gem_folder = self.get_task_clipName_gem_folder(taskName, clipName)
        gem_image_path = os.path.join(clip_gem_folder, f'{taskName}_{sliceName}_{clipName}.png')
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
        gem_path = os.path.join(clip_gem_folder, f'{taskName}_{sliceName}_{clipName}.gem')
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
        stain_path = os.path.join(clip_stain_folder, f'{taskName}_{sliceName}_{clipName}.tif')
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