from api.SCS.scs import *
from dataManager.segmentation_d import segData
import scanpy as sc
import re
from utils.commonfuc import *

def run_scs(taskName):
    """
    运行胞域扩增任务
    """
    taskInfo = expData.get_running_task(taskName)
    taskInfo['running'] = True
    patchSize = taskInfo['patchSize']
    binSize = taskInfo['binSize']
    epochs = taskInfo['epochs']
    diameter = taskInfo['diameter']
    neighbors = taskInfo['neighbors']
    taskName = taskInfo['taskName']
    slices = taskInfo['slices']
    total_slices = len(slices)
    try:
        for i, slice in enumerate(slices, start=1):
            bin_file = slice['gem']
            z_index = slice['z']
            adata_path = segData.get_seg_adata_path(taskName, z_index)
            project_path = expData.get_expTaskSlices_zfolder(taskName, z_index)
            adata = sc.read_h5ad(adata_path)
            segment_cells(adata, bin_file, project_path, taskInfo, slice, layer_name='watershed_mask', patch_size=patchSize, bin_size=binSize, n_neighbor=neighbors, epochs=epochs, r_estimate=diameter)
            release_memory()
            rows, cols = adata.layers['stain'].shape
            result_folder = expData.get_expTask_result_folder(taskName, z_index)
            expansion_mask = process_cell_masks(result_folder, rows, cols)
            adata.layers['expansion_mask'] = expansion_mask
            mask_fig, contour_fig = generate_cell_masks_rgba(expansion_mask, hex_colors=None)
            expansion_mask_path = segData.get_seg_expansion_mask_figure_path(taskName, z_index)
            expansion_contour_path = segData.get_seg_expansion_contour_figure_path(taskName, z_index)
            write_pkl(mask_fig, expansion_mask_path)
            write_pkl(contour_fig, expansion_contour_path)
            adata.write_h5ad(adata_path)
            set_slice_status(taskInfo, slice, 'patchprocess', Status.SUCCESS, slice['patchprocess']['text'])
            percent = i / total_slices
            set_progress(taskInfo, percent)
        release_memory()
    except Exception as e:
        stop_expansion_task(taskInfo, taskName)
    complete_expansion_task(taskInfo, taskName)

def process_cell_masks(folder_path, original_height, original_width):
    """
    处理细胞分割结果文件，生成完整的细胞ID矩阵
    
    Args:
        folder_path: 包含分割结果文件的文件夹路径
        original_height: 原图高度
        original_width: 原图宽度
    
    Returns:
        numpy.ndarray: 完整的细胞ID矩阵，未识别区域为0
    """
    result_matrix = np.zeros((original_height, original_width), dtype=int)
    
    files = [f for f in os.listdir(folder_path) if f.startswith('spot2cell_') and f.endswith('.txt')]
    
    next_cell_id = 1
    
    for file_name in files:
        match = re.search(r'spot2cell_(\d+):(\d+):(\d+):(\d+)\.txt', file_name)
        if not match:
            continue
            
        start_row, start_col, block_height, block_width = map(int, match.groups())
        
        file_path = os.path.join(folder_path, file_name)
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        id_mapping = {}
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
                
            parts = line.split('\t')
            if len(parts) != 2:
                continue
                
            coord_part, original_id = parts[0], int(parts[1])
            row_idx, col_idx = map(int, coord_part.split(':'))
            
            global_row = start_row + row_idx
            global_col = start_col + col_idx
            
            if global_row < original_height and global_col < original_width:
                if original_id not in id_mapping:
                    id_mapping[original_id] = next_cell_id
                    next_cell_id += 1
                
                result_matrix[global_row, global_col] = id_mapping[original_id]
    
    return result_matrix