from dash import html, callback, Input, State, Output, no_update
from dash.exceptions import PreventUpdate
import feffery_antd_components as fac
from controller.auth import verify_user
import os
import getpass

username = getpass.getuser()

login_box = fac.AntdModal(
    [
        fac.AntdTabs(
            items=[
                {
                    'key': 'Account',
                    'label': 'Account',
                    'forceRender':True,
                    'children': html.Div(
                        fac.AntdCenter(
                            fac.AntdSpace(
                                [
                                    fac.AntdInput(
                                        id='login-username',
                                        value=username,
                                        placeholder='Username for Linux',
                                        prefix=fac.AntdIcon(icon='antd-user'),
                                        disabled=True,
                                        style={'height':'35px'}
                                    ),
                                    fac.AntdInput(
                                        id='login-password',
                                        placeholder='Password for Linux',
                                        prefix=fac.AntdIcon(icon='antd-key'),
                                        style={'height':'35px'},
                                        mode='password',
                                    ),
                                ],
                                style={'width':'70%'},
                                direction='vertical',
                                size='middle'
                            ),
                            style={'height':150},
                        )
                    )
                },
                {
                    'key': 'Passcode',
                    'label': 'Passcode',
                    'forceRender':True,
                    'children': html.Div(
                        fac.AntdCenter(
                            fac.AntdSpace(
                                [
                                    fac.AntdInput(
                                        id='login-passcode',
                                        placeholder='Enter Passcode',
                                        prefix=fac.AntdIcon(icon='antd-idcard'),
                                        style={'height':'35px'}
                                    ),
                                ],
                                style={'width':'70%'},
                                direction='vertical',
                                size='middle'
                            ),
                            style={'height':150},
                        )
                    )
                }
            ],
            id='login-tabs',
            activeKey='Account'
        ),
        html.Div(
            fac.AntdSpace(
                [
                    fac.AntdButton('Reset', id='login-reset', type='primary', style={'backgroundColor':'#ab6953'}),
                    fac.AntdButton('Login', id='login-submit', type='primary', style={'backgroundColor':'#507ea4'})
                ],
                size='middle'
            ),
            style={'text-align':'right'}
        )
    ],
    id='login-box', 
    visible=False,
    closable=False,
    keyboard=False,
    maskClosable=False,
    forceRender=True,
    destroyOnClose=False,
    maskStyle={'backgroundColor':'#5F9EA0'},
    title='User Login'
)

@callback(
    Output('login-password', 'value'),
    Output('login-passcode', 'value'),
    Input('login-reset', 'nClicks'),
    State('login-tabs', 'activeKey')
)
def reset_input(nClicks, key):
    if nClicks:
        if key=='Account':
            return None, no_update
        else:
            return no_update, None
    raise PreventUpdate

@callback(
    Output('main-notice-area', 'children', allow_duplicate=True),
    Input('login-submit', 'nClicks'),
    State('login-tabs', 'activeKey'),
    State('login-password', 'value'),
    State('login-passcode', 'value'),
    running=[
        (Output("login-submit", "disabled"), True, False),
    ],
    prevent_initial_call=True
)
def usr_login(nClicks, key, password, passcode):
    if nClicks:
        notice = no_update
        if key=='Account':
            usrname = os.getlogin()
            access = verify_user(usrname, password)
            if access:
                notice = fac.AntdMessage(content='Login Successfully. Welcome !', type='success')  
            else:
                notice = fac.AntdMessage(content='Incorrect password. Please try again.', type='error')
        else:
            access = verify_user(passcode)
            if access:
                notice = fac.AntdMessage(content='Login Successfully. Welcome !', type='success')  
            else:
                notice = fac.AntdMessage(content='Invalid passcode. Please try again.', type='error')
        return notice
    raise PreventUpdate