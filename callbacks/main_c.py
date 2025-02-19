from dash import Input, Output, callback, State
from dash.exceptions import PreventUpdate
import pages.segmentation_p as segmentation
import pages.annotation_p as annotation
import pages.alignment_p as alignment
import pages.passcode_p as passcode
from controller.auth import *
from dash import Patch

menu = [
    {'title':'Segmentation', 'icon':'fc-mind-map', 'page':segmentation},
    {'title':'Annotation', 'icon':'fc-view-details', 'page':annotation},
    {'title':'Alignment', 'icon':'fc-workflow', 'page':alignment},
]

menu.append({'title':'Passcode', 'icon':'fc-tree-structure', 'page':passcode})

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
    Output('main-refresh', 'href', allow_duplicate=True),
    Input('main-button-logout', 'confirmCounts'),
    State('userid', 'data'),
    prevent_initial_call=True
)
def confirm_logout(click, userid):
    if click:
        logout(userid)
        return '/'
    raise PreventUpdate

@callback(
    Input('auth-interval', 'n_intervals')
)
def verify_usrhost(n):
    verify_host()

@callback(
    Output('main-title-header', 'children'), 
    Input('main-menu-item', 'currentKey'),
    Input('main-title-username', 'children'),
)
def update_header_title_by_menuItem(currentKey, usrname):
    if not usrname:
        raise PreventUpdate
    key = get_key(currentKey, usrname)
    return menu[key]['title']

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