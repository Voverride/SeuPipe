import numpy as np

def hex_to_rgb(hex_color):
    """将16进制颜色转换为RGB元组"""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

def rgb_to_hex(rgb):
    """将RGB元组转换为16进制颜色"""
    return '#' + ''.join(f'{int(c):02x}' for c in rgb)

def lighten_color(hex_color, factor=0.5):
    """将颜色变为浅色，factor表示变亮的程度"""
    rgb = hex_to_rgb(hex_color)
    lightened_rgb = tuple(min(255, int(c + (255 - c) * factor)) for c in rgb)
    return rgb_to_hex(lightened_rgb)

def lighten_colors(hex_colors):
    """给定16进制颜色列表，生成对应的浅色列表"""
    return [lighten_color(color) for color in hex_colors]

def get_color_map(labels, type=None):
    """
        生成一个颜色映射，根据不同的 `type` 选择颜色。
        参数:
        - labels (list): 标签列表，颜色将根据标签分配。
        - type (str): 颜色类型。可选值为：
            -  None，不进行颜色映射（默认值）。
            - 'primary'： 主要颜色
            - 'ui'： UI 软件界面色系
            - 'binary'：二值颜色
            - 'tab20': tab20 色盘
            - 'set3': set3 色盘
            - 'paired': Paired色盘
            - 'viridis': viridis色盘
            - 'plasma' : plasma色盘
            - 'COLORS_60' : COLORS_60色盘
        返回:
        - dict: 标签到颜色的映射字典。
    """
    labels = list(labels)
    labels.sort()
    if type=='primary':
        colorList = primaryColors
    elif type=='ui':
        colorList = uiColors
    elif type=='binary':
        colorList = binaryColors
    elif type=='tab20':
        colorList = tab20Colors
    elif type=='set3':
        colorList = set3Colors
    elif type=='paired':
        colorList = pairedColors
    elif type=='viridis':
        colorList = viridisColors
    elif type=='plasma':
        colorList = plasmaColors
    elif type=='COLORS_60':
        colorList = COLORS_60
    else:
        return None
    colorMap = {label:colorList[i%len(colorList)] for i, label in enumerate(labels)}
    return colorMap

def get_scale_colors(nums:list, minValue, maxValue, colorType='blue', rangeDist=[0, 0.25, 0.75, 1]):
    """
        对于给定数组返回对应连续变化颜色值
        nums: 数据列表
        colorType: 可选值 red, green, blue, 默认值 blue
        rangeDist: 对应颜色的变化区间, 默认值[0, 0.25, 0.75, 1]
    """
    colors = blueScaleColor
    if colorType=='red':
        colors = redScaleColor
    if colorType=='green':
        colors = greenScaleColor
    colorMap = {val:color for val, color in zip(rangeDist, colors)}
    return [get_color_for_value(val, minValue, maxValue, colorMap) for val in nums]

def get_color_for_value(value, min_value, max_value, color_map):
    """
        获取任意数值在min和max之间对应的渐变色
    """
    normalized_value = (value - min_value) / (max_value - min_value + 0.01)
    
    sorted_keys = sorted(color_map.keys())
    
    lower_key = max([key for key in sorted_keys if key <= normalized_value], default=0)
    upper_key = min([key for key in sorted_keys if key > normalized_value], default=1)
    
    def rgb(hex):
        return np.array([int(hex[i:i+2], 16) / 255 for i in (1, 3, 5)])
    
    def hex(rgb):
        return '#' + ''.join([f'{int(c * 255):02x}' for c in rgb])

    lower_color = rgb(color_map[lower_key])
    upper_color = rgb(color_map[upper_key])
    
    t = (normalized_value - lower_key) / (upper_key - lower_key) if upper_key != lower_key else 0
    interpolated_color = lower_color + t * (upper_color - lower_color)

    return hex(interpolated_color)

# 渐变色

redScaleColor = ['#dcdddd', '#a3a3a2','#d7003a', '#640125']
blueScaleColor = ['#dcdddd', '#a3a3a2', '#1e50a2', '#17184b']
greenScaleColor = ['#dcdddd', '#a3a3a2', '#028760', '#005243']


# 主要 色盘：--------------------------------------------------
primaryColors = ['#8870ad','#647a4f','#C72228','#B51D8D','#139992','#FBBE92','#EF5A9D','#ff891c','#005579','#8EC792','#CDE088','#3F84AA','#cc7818','#65A83E','#C3C388','#532C8A','#f79083','#9e6762','#f9decf','#EF4E22','#8DB5CE','#354E23','#c9a997','#FACB12','#F397C0','#c19f70','#0F4A9C','#DFCDE4','#635547','#C594BF','#DABE99','#989898','#1a1a1a','#f6bfcb','#7f6874']

# ui主题 色盘：--------------------------------------------------
uiColors=['#698aab', '#867ba9', '#a58f86', '#5F9EA0', '#ca8269']

# 二值 色盘：--------------------------------------------------
binaryColors=['#867ba9', '#ca8269']

# tab20 色盘：--------------------------------------------------
tab20Colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#9f9f9f', '#ff9f0e', '#e58c01', '#3d5f3e', '#9c47a0', '#c0e517', '#595ed5', '#7267d0', '#19c5bb', '#b6b6b6']

# Set3 色盘：--------------------------------------------------
set3Colors = ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5', '#ffed6f']

# Paired 色盘：--------------------------------------------------
pairedColors = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928']

# viridis 色盘：--------------------------------------------------
viridisColors = ['#440154', '#482576', '#3e4a89', '#31688e', '#26828e', '#1f9d8a', '#6cce5b', '#b4de2c', '#fee725']

# plasma 色盘：--------------------------------------------------
plasmaColors = ['#0d0887', '#46039f', '#7201a8', '#9c179e', '#d61c61', '#f19c4c', '#f8c54a', '#f8f8ff']

# COLORS_60 色盘 ------------------------------------------------
COLORS_60 = [
    "#c0301d", "#b4cab9", "#008000", "#1166FF", "#00CED1", "#49784e", "#DEB887", "#008080", "#FFA500", "#483D8B", 
    "#FF69B4", "#AFEEEE", "#D2B48C", "#FFDAB9", "#DA70D6", "#FFA07A", "#DC143C", "#00008B", "#9932CC", "#16677b",
    "#800080", "#8D54ad", "#4ec2b3", "#e98ee3", "#CD5C5C", "#F0E68C", "#D2691E", "#6495ED", "#800000", "#808000", 
    "#fbde7e", "#2E8B57", "#D8BFD8", "#2F4F4F", "#FF8C00", "#8Bcc00", "#d1c681", "#482a44", "#e3c67d", "#406441", "#697e66", "#533158",
    "#008B8B", "#4B0082", "#5F9EA0", "#BDB76B", "#000000", "#FF6347", "#708090", "#FF7F50", "#B8860B", 
    "#006400", "#556B2F", "#E9967A", "#8FBC8F", "#94338a", "#696969",
    "#00BFFF", "#1E90FF", "#B22222", "#228B22", "#A52A2A",
    "#20B2AA", "#a23919", "#9b7d24", "#443645", "#39607d", "#e6dbcb", "#355b6e", "#423e55", "#b0b173", "#376d78", 
    "#2d5aa2", "#486343", "#552141", "#4682B4"]

lightened_colors = lighten_colors(primaryColors)