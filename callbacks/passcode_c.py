from dash import Input, Output, callback, State, no_update, Patch
from dash.exceptions import PreventUpdate
import feffery_antd_components as fac
from controller.passcode_ctl import *


@callback(
    Output('main-notice-area', 'children', allow_duplicate=True),
    Output('passcode-table', 'data', allow_duplicate=True),
    Output('passcode-table', 'selectedRowKeys'),
    Input('delete-passcode', 'nClicks'),
    State('passcode-table', 'selectedRows'),
    prevent_initial_call=True
)
def delete_passcode(nClicks, rows):
    if nClicks and rows:
        ids = [usr['id'] for usr in rows]
        keys = [int(usr['key']) for usr in rows]
        keys.sort(reverse=True)
        data = Patch()
        for key in keys:
            del data[key]
        selectedRowKeys = []
        notice = fac.AntdMessage(content='Remove Successfully !', type='success')
        remove_passcode(ids)
        return notice, data, selectedRowKeys
    raise PreventUpdate

@callback(
    Output('main-notice-area', 'children', allow_duplicate=True),
    Output('create-passcode', 'visible', allow_duplicate=True),
    Output('passcode-table', 'data', allow_duplicate=True),
    Input('create-submit', 'nClicks'),
    State('input-passcode', 'value'),
    State('select-permission', 'value'),
    State('select-status', 'value'),
    prevent_initial_call=True
)
def submit_new_passcode(nClicks, passcode, permission, status):
    if nClicks:
        notice = no_update
        visible = no_update
        data = no_update
        if not passcode:
            notice = fac.AntdMessage(content='Failed to create, The passcode cannot be empty !', type='error')
        else:
            res = create_passcode(passcode, permission, status)
            if res:
                visible = False
                data = Patch()
                data.append(res)
                notice = fac.AntdMessage(content='Create Successfully !', type='success')
            else:
                notice = fac.AntdMessage(content='Failed to create, This passcode was exists in database !', type='error')
        return notice, visible, data
    raise PreventUpdate

@callback(
    Output('create-passcode', 'visible', allow_duplicate=True),
    Input('create-cancel', 'nClicks'),
    prevent_initial_call=True
)
def close_create_box(nClicks):
    if nClicks:
        return False
    raise PreventUpdate

@callback(
    Output('create-passcode', 'visible', allow_duplicate=True),
    Input('new-passcode', 'nClicks'),
    prevent_initial_call=True
)
def show_create_box(nClicks):
    if nClicks:
        return True
    raise PreventUpdate

@callback(
    Output('main-notice-area', 'children', allow_duplicate=True),
    Input('passcode-table', 'recentlySelectRow'),
    Input('passcode-table', 'recentlySelectDataIndex'),
    Input('passcode-table', 'recentlySelectValue'),
    prevent_initial_call=True
)
def modify_permission_status(row, index, value):
    if not value:
        raise PreventUpdate
    id = row['id']
    update_permission_status(id, index, value)
    return fac.AntdMessage(content='Update Successfully !', type='success')

@callback(
    Output('main-notice-area', 'children', allow_duplicate=True),
    Output('passcode-table', 'data', allow_duplicate=True),
    Input('passcode-table', 'recentlyChangedRow'),
    Input('passcode-table', 'recentlyChangedColumn'),
    prevent_initial_call=True
)
def modify_passcode(row, col):
    if not col:
        raise PreventUpdate
    id = row['id']
    new_passcode = row['Passcode']
    key = int(row['key'])
    notice = no_update
    patch = no_update

    if new_passcode=='':
        notice = fac.AntdMessage(content='Failed to update, The passcode cannot be empty !', type='error')
        patch = Patch()
        old_passcode = get_old_passcode(id)
        patch[key]['Passcode'] = old_passcode
    else:
        notice = fac.AntdMessage(content='Update Successfully !', type='success')
        update_passcode(id, new_passcode)
    return notice, patch

@callback(
    Output('passcode-table', 'data'),
    Input('main-title-header', 'children'),
)
def initial_user_table(header):
    if header=='Passcode':
        data = get_all_passcode()
        return data
    raise PreventUpdate