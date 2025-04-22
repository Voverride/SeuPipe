import pandas as pd
import traceback

# def construct_exception(exc):
#     """
#     构造异常的详细文本描述
#     """
#     # exc_type = type(exc)
#     # exc_value = exc
#     # exc_traceback = exc.__traceback__
    
#     # full_traceback = traceback.format_exception(exc_type, exc_value, exc_traceback)
    
#     # error_description = {
#     #     'Exception Type': exc_type.__name__,
#     #     'Error Message': str(exc_value),
#     #     'Traceback': ''.join(full_traceback)
#     # }
#     # return error_description
#     # return traceback.format_exc()

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