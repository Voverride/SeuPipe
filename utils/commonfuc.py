import pandas as pd
from PIL import Image
import plotly.express as px
import io
import base64
import numpy as np
from io import BytesIO
from datetime import datetime

def array_to_base64(arr):
    """
    将 NumPy 数组转换为 Base64 编码的 PNG 图片
    """
    pil_img = Image.fromarray(arr)
    buf = BytesIO()
    pil_img.save(buf, format='PNG')
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()

def get_current_date():
    """
    获取当前日期
    """
    now = datetime.now()
    formatted_time = now.strftime("%Y/%m/%d %H:%M")
    return formatted_time
def get_imgfig_withplotly(img_array, title=None):
    """
    基于图像矩阵绘制交互图形
    """
    fig = px.imshow(img_array)
    fig.update_layout(
        autosize=True,
        title=dict(
            text=title,
            x=0.5,
        ),
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
    )
    return fig

def array_to_base64(img_array):
    """将 NumPy 数组转换为 Base64 编码的 PNG 图片"""
    pil_img = Image.fromarray(img_array.astype(np.uint8))
    buffer = io.BytesIO()
    pil_img.save(buffer, format="PNG")
    return base64.b64encode(buffer.getvalue()).decode('utf-8')

def is_all_numeric(series)->bool:
    """
    判断dataframe中某一列的数据是否为纯数值
    """
    if pd.api.types.is_numeric_dtype(series):
        return True
    converted = pd.to_numeric(series, errors='coerce')
    original_na = series.isna().sum()
    new_na = converted.isna().sum()
    return new_na == original_na