from dash import Input, Output, callback, State
from setting import setting
from dash.exceptions import PreventUpdate
import pages.regionclip_p as regionclip
import pages.cellselector_p as cellselector
import pages.segmentation_p as segmentation
import pages.expansion_p as expansion
import pages.maskviewer_p as maskviewer
import pages.annotation_p as annotation
import pages.alignment_p as alignment
# import pages.visiualization_p as visualization
# import pages.lassion_p as lassion

# import pages.filtering_p as filtering

# 注意 fileSelecter 模块的注释要改过来
from pages.components.fileSelecter import fileSelecter
import pages.passcode_p as passcode
from websocket.message import ms
from dataManager.workspace import *
from controller.auth import *
from dash import Patch


menu = [
    {'title':setting.regionClip, 'icon':'fc-repair', 'page':regionclip},
    {'title': setting.segmentation, 'icon':'fc-radar-chart', 'page':segmentation},
    {'title':setting.expansion, 'icon':'fc-mind-map', 'page':expansion},
    {'title': setting.maskViewer, 'icon':'fc-data-sheet', 'page':maskviewer},
    {'title':setting.cellSelector, 'icon':'fc-serial-tasks', 'page':cellselector},
    # {'title':setting.filtering, 'icon':'fc-multiple-inputs', 'page':filtering},
    {'title':setting.annotation, 'icon':'fc-view-details', 'page':annotation},
    {'title':setting.alignment, 'icon':'fc-workflow', 'page':alignment},
    # {'title':setting.visualization, 'icon':'fc-scatter-chart', 'page':visualization},
    # {'title':setting.lassion, 'icon':'fc-briefcase', 'page':lassion},
]

menu.append({'title':setting.passcode, 'icon':'fc-tree-structure', 'page':passcode})

@callback(
    Output('main-menu-item', 'menuItems'),
    Input('main-title-username', 'children'),
)
def set_authorized_views(usrname):
    if not usrname:
        raise PreventUpdate
    legnth = len(menu)
    if usrname!=admin:
        legnth-=1
    menuItems = [
        {
            'component': menu[i]['title'],
            'props': {
                'key': i,
                'title': menu[i]['title'],
                'icon': menu[i]['icon'],
            },
        } for i in range(legnth)
    ]
    return menuItems

@callback(
    Input('main-button-logout', 'confirmCounts'),
    State('userid', 'data'),
)
def confirm_logout(click, userid):
    if click:
        logout(userid)
        ms.clinetSend(ms._logout_c2s)
@callback(
    Output('main-title-header', 'children'), 
    Input('main-menu-item', 'currentKey'),
    Input('main-title-username', 'children'),
)
def update_header_title_by_menuItem(currentKey, usrname):
    if not usrname:
        raise PreventUpdate
    key = get_key(currentKey, usrname)
    menuItem = menu[key]['title']
    ms.clinetSend(ms._updateMenuItem_c2s, menuItem=menuItem)
    return menuItem

@callback(
    Output('main-sider-control', 'children'), 
    Input('main-menu-item', 'currentKey'),
    Input('main-title-username', 'children'),
)
def update_sider_control_by_menuItem(currentKey, usrname):
    if not usrname:
        raise PreventUpdate
    key = get_key(currentKey, usrname)
    return menu[key]['page'].control_panel

@callback(
    Output('main-center-content', 'children'), 
    Input('main-menu-item', 'currentKey'),
    Input('main-title-username', 'children'),
)
def update_center_content_by_menuItem(currentKey, usrname):
    if not usrname:
        raise PreventUpdate
    key = get_key(currentKey, usrname)
    return menu[key]['page'].content_panel

def get_key(currentKey, usrname)->int:
    length = len(menu)
    if usrname!=admin:
        length-=1
    key = int(currentKey)%length
    return key

@callback(
    Output('main-error-area', 'style'),
    Input('main-error-area-close', 'nClicks')
)
def close_error_area(nc):
    if nc:
        pat = Patch()
        pat['display'] = 'none'
        return pat
    raise PreventUpdate

@callback(
    Input('SeuPipe', 'style'),
)
def select_workspace(style):
    """
    登录成功后弹出选择工作目录窗口
    """
    if style is None:
        workspace = get_workspace()
        if workspace is None:
            fileSelecter.open_workspace_box()

@callback(
    Input('ws', 'message')
)
def verify_user_state(message):
    """
    刷新页面时，检查用户登录状态
    """
    ms.clientParse(message)