import ctypes
import subprocess
from dataManager.alignment_d import alidata
from dataManager.annotation_d import annData
def get_available_gpus():
    """
    使用 nvidia-smi 获取当前系统中可用的 GPU 列表。
    返回 GPU 编号列表（从 0 开始）。
    """
    try:
        result = subprocess.run(
            ['nvidia-smi', '--query-gpu=index', '--format=csv,noheader'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        if result.returncode == 0:
            gpus = [int(idx.strip()) for idx in result.stdout.splitlines()]
            return gpus
        else:
            return []
    except FileNotFoundError:
        return []
    
def reset_gpu_memory():
    annstatus = annData.get_anntstatus()
    alistatus = alidata.get_alistatus()
    ann_thread = annstatus['thread']
    ali_thread = alistatus['thread']
    ann_alive = ann_thread is not None and ann_thread.is_alive()
    ali_alive = ali_thread is not None and ali_thread.is_alive()
    if ann_alive or ali_alive:
        return
    try:
        cu = ctypes.CDLL('libcudart.so')
        gpus = get_available_gpus()
        for gpu in gpus:
            cu.cudaSetDevice(gpu)
            cu.cudaDeviceReset()
        import torch
        torch.cuda.set_device(0)
    except Exception as e:
        pass