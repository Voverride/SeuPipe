from dataManager.workspace import *
from utils.commonfuc import *
from utils.typing import StepStatus
from websocket.websocket import ws
from utils.observer import observer
from controller.auth import get_request_usrname
from sklearn.preprocessing import normalize
import atexit
import time
import threading
import scanpy as sc
import cachetools
import shutil
import bisect

class AlignmentData:
    def __init__(self, cleanup_interval=3600, ttl=3600):
        self._alignment = {}
        self._ismodify = {}
        self._normCountField = 'seupipe_alicounts'
        self._projects = cachetools.TTLCache(maxsize=200, ttl=3600)
        self._export_project = cachetools.TTLCache(maxsize=200, ttl=1800)
        self._genelist = cachetools.TTLCache(maxsize=200, ttl=3600)
        self._fields = cachetools.TTLCache(maxsize=200, ttl=3600)
        self._initfig = cachetools.TTLCache(maxsize=200, ttl=3600)
        self._coordinate = {}
        self._coordinateTime = {}
        self._cleanup_interval = cleanup_interval
        self._ttl = ttl
        self._stop_cleanup = False
        self._cleanup_thread = None
        if os.environ.get("WERKZEUG_RUN_MAIN") is None:
            self.start_cleanup_thread()
            atexit.register(self.save_all_coordinates)

    def get_normcount_field(self):
        """
        获取标准化计数字段名
        """
        return self._normCountField

    def get_fields(self, project):
        """
        获取项目所有字段
        """
        fields = self._fields.get(project, None)
        if fields is None:
            coordinate = self._coordinate.get(project, None)
            if coordinate is None:
                return []
            fields = coordinate.columns.tolist()
            fields.sort()
            self._fields[project] = fields
        return fields
    
    def get_gene_list(self, project):
        """
        获取项目所有基因
        """
        gene_list = self._genelist.get(project, None)
        if gene_list is None:
            initial_path = self.get_initialdata_path(project)
            initialData = sc.read_h5ad(initial_path, backed='r')
            gene_list = initialData.var.index.tolist()
            initialData.file.close()
            gene_list.sort()
            self._genelist[project] = gene_list
        return gene_list

    def set_export_project(self, projectname):
        """
        设置用户导出项目
        """
        username = get_request_usrname()
        self._export_project[username] = projectname
    
    def get_export_project(self):
        """
        获取用户导出项目
        """
        username = get_request_usrname()
        return self._export_project.get(username, None)

    def start_cleanup_thread(self):
        """
        启动后台清理线程
        """
        def cleanup_loop():
            while not self._stop_cleanup:
                self._cleanup_expired()
                time.sleep(self._cleanup_interval)
        
        self._cleanup_thread = threading.Thread(target=cleanup_loop, daemon=True)
        self._cleanup_thread.start()
    
    def _cleanup_expired(self):
        """
        清理过期数据
        """
        current_time = time.time()
        expired_keys = []
        
        for project, store_time in list(self._coordinateTime.items()):
            age = current_time - store_time
            if age > self._ttl:
                expired_keys.append(project)
        
        for project in expired_keys:
            coordinate = self._coordinate[project]
            self.update_coordinate(project, coordinate)
            del self._coordinate[project]
            del self._coordinateTime[project] 

    def save_all_coordinates(self):
        """
        保存所有坐标到磁盘
        """
        for key in list(self._coordinate.keys()):
            if key in self._coordinate:
                coordinate = self._coordinate[key]
                self.update_coordinate(key, coordinate)
                
    def is_project_modify(self, project):
        """
        判断项目是否正在被手动修改
        """
        return self._ismodify.get(project, None)
    
    def set_project_modify(self, project):
        """
        设置项目正在被手动修改
        """
        username = get_request_usrname()
        self._ismodify[project] = username

    def unset_project_modify(self, project):
        """
        取消项目正在被手动修改
        """
        if project in self._ismodify:
            del self._ismodify[project]

    def get_coord_trans_mtx(self, dx=None, dy=None, reg=None):
        """
        获取坐标转换矩阵
        """
        if not dx and not dy and not reg:
            return None
        if dx:
            return [dx, 0, 0]
        if dy:
            return [0, dy, 0]
        if reg:
            return [0, 0, reg]
        return None

    def get_syncslice_index(self, project, slice):
        """
        获取同步切片索引
        """
        znums = self.get_znums(project)
        idx = bisect.bisect_left(znums, slice)
        return len(znums) - idx -1
    
    def get_slice_index(self, project, slice):
        """
        获取同步切片索引
        """
        znums = self.get_znums(project)
        idx = bisect.bisect_left(znums, slice)
        return idx
    
    def get_znums(self, projectname):
        """
        获取切片列表
        """
        metadata = self.get_project_info(projectname)
        znums = metadata.get('znums', [])
        return znums

    def create_alignment_task(self, projectname):
        """
        创建对齐任务
        """
        metadata = self.get_project_info(projectname)
        metadata['steps'] = {
            'preprocess': StepStatus.WAIT,
            'alignment': StepStatus.WAIT,
            'postprocess': StepStatus.WAIT,
            'percent': 0,
        }
        metadata['exception'] = None
        observed_metadata = observer.observe(metadata, ws.notifyUpdateAlignmentStatus, project=projectname)
        self._alignment[projectname] = observed_metadata
        return observed_metadata
    
    def update_manual_adjust_status(self, project, operations):
        """
        更新手动调整任务状态
        """
        signal = dict()
        username = get_request_usrname()
        observed_signal = observer.observe(signal, ws.notifyUpdateManualAdjustStatus, project=project, operations=operations, username=username)
        observed_signal['update'] = True
        observer.disobserve(observed_signal, ws.notifyUpdateManualAdjustStatus)
        del signal, observed_signal

    def remove_alignment_task(self, projectname):
        """
        移除对齐任务
        """
        if projectname in self._alignment:
            observer.disobserve(self._alignment[projectname], ws.notifyUpdateAlignmentStatus)
            del self._alignment[projectname]
    
    def is_running(self, projectname):
        """
        检查是否存在对齐任务
        """
        if projectname not in self._alignment:
            return False
        metadata = self.get_project_info(projectname)
        return metadata.get('running', False)
    def get_exist_projects(self):
        """
        获取存在的对齐项目列表
        """
        ali_path = get_alignment_workspace()
        projects = os.listdir(ali_path)
        projects = [p for p in projects if os.path.isdir(os.path.join(ali_path, p))]
        projects.sort()
        return projects

    def get_project_info(self, project_name: str):
        """
        获取对齐项目信息
        """
        project_info = dict(self._alignment.get(project_name, {}))
        if project_info: 
            return project_info
        elif project_name in self._projects:
            project_info = self._projects[project_name]
            del self._projects[project_name]
        else:
            metadata_path = self.get_metadata_path(project_name)
            project_info = read_json(metadata_path)
        self._projects[project_name] = project_info
        return project_info

    def get_metadata_path(self, project_name: str):
        """
        获取对齐项目元数据路径
        """
        project_folder = self.get_project_folder(project_name)
        metadata_path = os.path.join(project_folder, 'metadata.json')
        return metadata_path

    def get_project_folder(self, project_name: str):
        """
        获取对齐项目文件夹路径
        """
        ali_path = get_alignment_workspace()
        project_path = os.path.join(ali_path, project_name)
        check_path(project_path)
        return project_path
    
    def delete_project(self, project: str):
        """
        删除对齐选项目
        """
        project_folder = self.get_project_folder(project)
        if project in self._projects:
            del self._projects[project]
        if project in self._coordinate:
            del self._coordinate[project]
        if project in self._coordinateTime:
            del self._coordinateTime[project]
        if project in self._initfig:
            del self._initfig[project]
        if project in self._genelist:
            del self._genelist[project]
        if project in self._fields:
            del self._fields[project]
        if os.path.exists(project_folder):
            shutil.rmtree(project_folder)

    def get_initial_folder(self, project):
        """
        获取初始数据文件夹
        """
        project_folder = self.get_project_folder(project)
        initial_folder = os.path.join(project_folder, 'initial')
        check_path(initial_folder)
        return initial_folder
    
    def get_initialdata_path(self, project):
        """
        获取原始数据路径
        """
        initial_folder = self.get_initial_folder(project)
        initial_path = os.path.join(initial_folder, 'data.h5ad')
        return initial_path
    
    def get_coordinate_folder(self, project):
        """
        获取对齐项目坐标文件夹路径
        """
        project_folder = self.get_project_folder(project)
        coordinate_folder = os.path.join(project_folder, 'coordinate')
        check_path(coordinate_folder)
        return coordinate_folder
    
    def get_coordinate_path(self, project):
        """
        获取对齐项目坐标路径
        """
        coordinate_folder = self.get_coordinate_folder(project)
        coordinate_path = os.path.join(coordinate_folder, 'coordinate.csv')
        return coordinate_path
    
    def get_coordinate(self, project):
        """
        获取对齐项目坐标
        """
        if project in self._coordinate:
            coordinate = self._coordinate[project]
        else:
            coordinate_path = self.get_coordinate_path(project)
            coordinate = pd.read_csv(coordinate_path, index_col=0)
        self._coordinate[project] = coordinate
        self._coordinateTime[project] = time.time()
        return coordinate
    
    def update_coordinate(self, project, coordinate):
        """
        更新对齐项目坐标
        """
        self._coordinate[project] = coordinate
        self._coordinateTime[project] = time.time()
        coordinate_path = self.get_coordinate_path(project)
        coordinate.to_csv(coordinate_path)
    
    def get_figure_folder(self, project):
        """
        获取对齐项目可视化文件夹路径
        """
        project_folder = self.get_project_folder(project)
        figure_folder = os.path.join(project_folder, 'figure')
        check_path(figure_folder)
        return figure_folder
    
    def get_resultfig(self, project):
        """
        获取对齐项目可视化结果
        """
        metadata = self.get_project_info(project)
        x = metadata['x']
        y = metadata['y']
        z = metadata['z']
        coordinate = self.get_coordinate(project)
        resultfig = self.plot_3d_scatter(project, coordinate, x, y, z)
        return resultfig
    
    def update_initialfig(self, project, coordinate, x, y, z):
        """
        更新初始数据可视化结果
        """
        initialfig_path = self.get_initialfig_path(project)
        initialfig = self.plot_3d_scatter(project, coordinate, x, y, z, initial=True)
        write_pkl(initialfig, initialfig_path)
        return initialfig
    
    def get_initialfig_path(self, project):
        """
        获取初始数据可视化结果路径
        """
        figure_folder = self.get_figure_folder(project)
        initialfig_path = os.path.join(figure_folder, 'initialfig.pkl')
        return initialfig_path
    
    def get_initialfig(self, project):
        """
        获取初始数据可视化结果
        """
        if project in self._initfig:
            initialfig = self._initfig[project]
        else:
            initialfig_path = self.get_initialfig_path(project)
            if os.path.exists(initialfig_path):
                try:
                    initialfig = read_pkl(initialfig_path)
                except:
                    initialfig = None
            else:
                initialfig = None
        self._initfig[project] = initialfig
        return initialfig
    
    def update_project_info(self, project_name: str, metadata: dict):
        """
        更新对齐项目信息
        """
        metadata_path = self.get_metadata_path(project_name)
        write_json(metadata_path, metadata)

    def create_project(self, importdata, x, y, z, projectname):
        """
        创建对齐项目
        """
        initialdata_path = self.get_initialdata_path(projectname)
        data = sc.read_h5ad(importdata)
        self.detect_and_normalize_counts(data)
        coordinate = data.obs.copy()
        data.write_h5ad(initialdata_path)
        del data
        self.update_coordinate(projectname, coordinate)
        alidata.update_initialfig(projectname, coordinate, x, y, z)
        zmin = coordinate[z].min()
        zmax = coordinate[z].max()
        znums = list(set(coordinate[z].unique()))
        znums.sort()
        znum = len(znums)
        x_min = np.min(coordinate[x])
        x_max = np.max(coordinate[x])
        y_min = np.min(coordinate[y])
        y_max = np.max(coordinate[y])
        min_val = min(x_min, y_min)
        max_val = max(x_max, y_max)
        dist = max_val-min_val
        max_val += dist/2
        min_val -= dist/2
        iniRange = max_val-min_val
        x_scale = iniRange/(x_max-x_min)
        y_scale = iniRange/(y_max-y_min)
        iniScale = (x_scale+y_scale)/2
        metadata = {
            'creator': get_request_usrname(),
            'projectName': projectname,
            'date': get_current_date(),
            'exception': None,
            'running': False,
            'initialData': importdata,
            'x': x,
            'y': y,
            'z': z,
            'zmin': zmin,
            'zmax': zmax,
            'znum': znum,
            'znums': znums,
            'iniRange': iniRange,
            'iniScale': iniScale,
            'steps': {
                'preprocess': StepStatus.WAIT,
                'alignment': StepStatus.WAIT,
                'postprocess': StepStatus.WAIT,
                'percent': 0,
            }
        }
        metadata_path = self.get_metadata_path(projectname)
        write_json(metadata_path, metadata)

    def detect_and_normalize_counts(self, adata, log1p=True):
        """
        检测原始计数并进行标准化
        """
        def is_raw_counts(matrix):
            """
            判断矩阵是否为原始计数
            """
            if matrix is None:
                return False

            if sparse.issparse(matrix):
                data_sample = matrix.data[:min(10000, len(matrix.data))]
                min_val = matrix.min()
                max_val = matrix.max()
            else:
                if matrix.size > 10000:
                    data_sample = matrix.ravel()[np.random.choice(matrix.size, 10000, replace=False)]
                else:
                    data_sample = matrix.ravel()
                min_val = np.min(data_sample)
                max_val = np.max(data_sample)
            
            if min_val < 0:
                return False
            
            if len(data_sample) > 0:
                non_integer_ratio = np.sum(data_sample != np.floor(data_sample)) / len(data_sample)
                if non_integer_ratio > 0.1:
                    return False
            
            if max_val < 100 and log1p:
                return False
            if max_val > 0 and max_val < 50:
                median_val = np.median(data_sample[data_sample > 0]) if len(data_sample[data_sample > 0]) > 0 else 0
                if median_val < 10 and non_integer_ratio > 0.05:
                    return False
            return True
        
        def normalize_counts(matrix, log1p=True, target_sum=1e4):
            """
            对原始计数进行标准化
            """
            if sparse.issparse(matrix):
                normalized = matrix.astype(np.float32)
            else:
                normalized = matrix.astype(np.float32)
            
            if sparse.issparse(normalized):
                normalized = normalize(normalized, norm='l1', axis=1) * target_sum
            else:
                row_sums = normalized.sum(axis=1, keepdims=True)
                row_sums[row_sums == 0] = 1
                normalized = normalized / row_sums * target_sum
            
            if log1p:
                normalized = np.log1p(normalized)
            
            return normalized
                
        print("🔍 检测原始计数...")
        
        raw_matrix = getattr(adata.raw, 'X', None) if adata.raw is not None else None
        if raw_matrix is not None and is_raw_counts(raw_matrix):
            print("📊 使用 raw.X 作为原始计数")
            raw_dense = raw_matrix.toarray() if sparse.issparse(raw_matrix) else raw_matrix
            adata.layers[self._normCountField] = normalize_counts(raw_dense, log1p=log1p)
            return True
        
        counts_layer = adata.layers.get('counts', None)
        if counts_layer is not None and is_raw_counts(counts_layer):
            print("📊 使用 layers['counts'] 作为原始计数")
            counts_dense = counts_layer.toarray() if sparse.issparse(counts_layer) else counts_layer
            adata.layers[self._normCountField] = normalize_counts(counts_dense, log1p=log1p)
            return True
        
        if is_raw_counts(adata.X):
            print("📊 使用 X 作为原始计数")
            x_dense = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
            adata.layers[self._normCountField] = normalize_counts(x_dense, log1p=log1p)
            return True
        
        print("⚠️ 未检测到原始计数，所有数据可能已标准化")
        print("📊 复制 X 到 layers['" + self._normCountField + "']")
        
        adata.layers[self._normCountField] = adata.X.copy()
        
        return False

    def get_graph_view(self, actSlice:float, refSlice:float)->dict:
        """
        获取图像视角
        """
        eye = {
            'bottom': {'x': -3.1294811103329714e-20, 'y': -1.3257190114690003e-16, 'z': -2.1650635094610964},
            'top': {'x': 1.2149787566408422e-05, 'y': -0.006897060070669372, 'z': 2.1650525237080886},
            'center': {'x': 1.25, 'y': 1.25, 'z': 1.25}
        }
        if not actSlice or not refSlice:
            return eye['center']
        if float(actSlice)>=float(refSlice):
            return eye['top']
        return eye['bottom']

    def plot_3d_scatter(self, project, coordinate, x, y, z, initial=False, marker_size:int=5, boarder_width:int=1, boarder_color:str='#0d0015', grid_color:str = '#5F9EA0'):
        """
        绘制3D散点图
        """
        df = coordinate.copy()
        df = df.sort_values(by=z)
        x_min = np.min(df[x])
        x_max = np.max(df[x])
        y_min = np.min(df[y])
        y_max = np.max(df[y])
        min_val = min(x_min, y_min)
        max_val = max(x_max, y_max)
        dist = max_val-min_val
        max_val += dist/2
        min_val -= dist/2
        dist = max_val-min_val
        if initial:
            step = dist/20
            x_scale = dist/(x_max-x_min)
            y_scale = dist/(y_max-y_min)
            scale = (x_scale+y_scale)/2
        else:
            metadata = self.get_project_info(project)
            iniRange = metadata['iniRange']
            iniScale = metadata['iniScale']
            multy_factor = dist/iniRange
            step = dist/20/multy_factor
            scale = multy_factor*iniScale

        df[z] = df[z].astype(str)
        df['index'] = df.index
        cmap = get_color_map(set(df[z]), type='COLORS_60')
        fig = px.scatter_3d(
            df, x=x, y=y, z=z, 
            color=z, 
            custom_data='index',
            color_discrete_map=cmap,
        )
        fig.update_traces(
            marker=dict(
                size=marker_size,
                line=dict(color=boarder_color, width=boarder_width)
            )
        )
        fig.update_layout(
            autosize=True,
            uirevision=True,
            legend=dict(
                itemsizing='constant',
                traceorder = 'reversed',
                title=None
            ),
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            scene=dict(
                aspectratio=dict(x=scale, y=scale, z=1),
                xaxis=dict(visible=False, showgrid=True, showbackground=False, showticklabels=False, gridcolor=grid_color, tickfont=dict(color=grid_color), title='', dtick=step, range=[min_val, max_val]), 
                yaxis=dict(visible=False, showgrid=True, showbackground=False, showticklabels=False, gridcolor=grid_color, tickfont=dict(color=grid_color), title='', dtick=step, range=[min_val, max_val]), 
                zaxis=dict(visible=False, showgrid=False, showbackground=False, showticklabels=False, gridcolor=grid_color),
                camera=dict(projection=dict(type='orthographic'))
            ),
        )
        return fig

alidata = AlignmentData()