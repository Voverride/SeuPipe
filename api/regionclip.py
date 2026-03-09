import os
import imageio.v3 as iio
import numpy as np
import traceback
from utils.commonfuc import write_json, compress_image, write_pkl, convert_mtx_to_base64_image
from dataManager.regionclip_d import clipData
def clip_image(taskName, slicename, clipName, row_start, row_end, col_start, col_end):
    """
    裁剪染色图像和gem, 生成裁剪后gem热图缩略图, 并保存到指定目录
    """
    observed_status = clipData.get_running_task_status(taskName, clipName, slicename)
    try:
        img = clipData.get_imgCache(taskName, slicename)
        gem = clipData.get_gemCache(taskName, slicename)
        if img is None or gem is None:
            raise Exception('Invalid image or gem')
        cordMapFunc = clipData.get_cord_map_func(taskName, slicename)
        mapped_row_start, mapped_col_start, mapped_row_end, mapped_col_end = cordMapFunc(row_start, col_start, row_end, col_end)
        clipped_img = img[mapped_row_start:mapped_row_end+1, mapped_col_start:mapped_col_end+1]
        clipped_gem = gem[
            gem['x'].between(mapped_row_start, mapped_row_end) & 
            gem['y'].between(mapped_col_start, mapped_col_end)
        ].copy()

        clipped_image_folder = clipData.get_task_clipName_stain_folder(taskName, clipName)
        clipped_gem_folder = clipData.get_task_clipName_gem_folder(taskName, clipName)
        clipped_image_path = os.path.join(clipped_image_folder, f'{taskName}_z{slicename}_{clipName}.tif')
        clipped_gem_path = os.path.join(clipped_gem_folder, f'{taskName}_z{slicename}_{clipName}.gem')
        clipped_gem_image_path = os.path.join(clipped_gem_folder, f'{taskName}_z{slicename}_{clipName}.png')
        clip_bounds_path = os.path.join(clipped_image_folder, f'{taskName}_z{slicename}_{clipName}_clip_bounds.json')
        taskInfo = clipData.get_slice_info(taskName, slicename)
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

        stainfig_path = clipData.get_task_clipName_stainBase64_path(taskName, slicename, clipName)
        compressed_stain, _ = compress_image(clipped_img)
        stain_base64 = convert_mtx_to_base64_image(compressed_stain)
        write_pkl(stain_base64, stainfig_path)

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
        gemfig_path = clipData.get_task_clipName_gemBase64_path(taskName, slicename, clipName)
        compressed_mtx, _ = compress_image(mtx)
        gem_base64 = convert_mtx_to_base64_image(compressed_mtx)
        write_pkl(gem_base64, gemfig_path)

    except Exception as e:
        observed_status['exception'] = traceback.format_exc()
    finally:
        observed_status['running'] = False
    status = dict(observed_status)
    if clipData._notifyfunc in status:
        del status[clipData._notifyfunc]
    clipData.write_taskclip_status(taskName, clipName, slicename, status)
    clipData.remove_running_task(taskName, slicename, clipName)
