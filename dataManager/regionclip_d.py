from dataManager.workspace import *
from controller.auth import get_host
from dataManager.users import search_user
import shutil
import numpy as np
import pandas as pd
import imageio.v3 as iio
import pickle
import cachetools
from utils.commonfuc import compress_image, write_json

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
        self._gemCache = cachetools.TTLCache(maxsize=10, ttl=18000)
        self._imgCache = cachetools.TTLCache(maxsize=10, ttl=18000)
        self._cordMapFunc = {}

    def set_cord_map_func(self, func):
        """
        设置用户获取坐标映射函数
        """
        username = self.get_request_usrname()
        self._cordMapFunc[username] = func

    def get_cord_map_func(self):
        """
        用户获取坐标映射函数
        """
        username = self.get_request_usrname()
        return self._cordMapFunc.get(username, None)
    
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

    def get_gemCache(self, taskName, slicename):
        """
        获取任务gem缓存
        """
        gem_key = f'{taskName}_{slicename}'
        if gem_key in self._gemCache:
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
        self._gemCache[gem_key] = gem
        return gem

    def get_stain_img(self, taskName, slicename):
        """
        获取任务指定切片的染色压缩图像
        """
        img = self.get_imgCache(taskName, slicename)
        compressd_img, cord_map_func = compress_image(img)
        self.set_cord_map_func(cord_map_func)
        return compressd_img

    def clip_image(self, taskName, slicename, clipName, row_start, row_end, col_start, col_end):
        """
        裁剪染色图像和gem, 生成裁剪后gem热图缩略图, 并保存到指定目录
        """
        img = self.get_imgCache(taskName, slicename)
        gem = self.get_gemCache(taskName, slicename)
        if img is None or gem is None:
            raise Exception('Invalid image or gem')
        gem['x'] = gem['x'] - gem['x'].min()
        gem['y'] = gem['y'] - gem['y'].min()
        cordMapFunc = self.get_cord_map_func()
        mapped_row_start, mapped_col_start, mapped_row_end, mapped_col_end = cordMapFunc(row_start, col_start, row_end, col_end)
        clipped_img = img[mapped_row_start:mapped_row_end+1, mapped_col_start:mapped_col_end+1]
        clipped_gem = gem[
            gem['x'].between(mapped_row_start, mapped_row_end) & 
            gem['y'].between(mapped_col_start, mapped_col_end)
        ].copy()

        clipped_image_folder = self.get_task_clipName_stain_folder(taskName, clipName)
        clipped_gem_folder = self.get_task_clipName_gem_folder(taskName, clipName)
        clipped_image_path = os.path.join(clipped_image_folder, f'{taskName}_z{slicename}_{clipName}.tif')
        clipped_gem_path = os.path.join(clipped_gem_folder, f'{taskName}_z{slicename}_{clipName}.gem')
        clipped_gem_image_path = os.path.join(clipped_gem_folder, f'{taskName}_z{slicename}_{clipName}.png')
        clip_bounds_path = os.path.join(clipped_image_folder, f'{taskName}_z{slicename}_{clipName}_clip_bounds.json')
        taskInfo = self.get_slice_info(taskName, slicename)
        img_path = taskInfo.get('image', None)
        gem_path = taskInfo.get('gem', None)

        clip_bounds = {
            'original_image': img_path,
            'original_gem': gem_path,
            'clipped_image': clipped_image_path,
            'clipped_gem': clipped_gem_path,
            'clip_bounds': {
                'row_start': mapped_row_start,
                'row_end': mapped_row_end,
                'col_start': mapped_col_start,
                'col_end': mapped_col_end
            }
        }
        write_json(clip_bounds_path, clip_bounds)

        clipped_gem['x'] = clipped_gem['x'] - mapped_row_start
        clipped_gem['y'] = clipped_gem['y'] - mapped_col_start

        iio.imwrite(clipped_image_path, clipped_img)
        clipped_gem.to_csv(clipped_gem_path, sep='\t', index=False)

        if clipped_gem.empty:
            mtx = np.zeros((100, 100), dtype=np.uint8)
        else:
            countField = 'MIDCounts' if 'MIDCounts' in clipped_gem else 'MIDCount'
            grouped = clipped_gem.groupby(['x', 'y'])[countField].sum().reset_index()

            x_raw = grouped['x'].values
            y_raw = grouped['y'].values
            counts = grouped[countField].values

            max_row, max_col = clipped_img.shape

            mtx = np.zeros((max_row, max_col), dtype=np.float32)
            mtx[x_raw, y_raw] = counts

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
        gem_image_path = os.path.join(clip_gem_folder, f'{taskName}_z{sliceName}_{clipName}.png')
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
        gem_path = os.path.join(clip_gem_folder, f'{taskName}_z{sliceName}_{clipName}.gem')
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
        stain_path = os.path.join(clip_stain_folder, f'{taskName}_z{sliceName}_{clipName}.tif')
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