from dash import html, dcc
import feffery_antd_components as fac
from callbacks.main_c import *
from .components.login import login_box
from .components.fileSelecter import fileSelecter

main_layout = html.Div([
    login_box,
    html.Div(id='main-notice-area'), 
    fac.AntdAlert(
        message='Error',
        description=html.Div(
            [
                fac.AntdIcon(icon='fc-rules'),
                fac.AntdText(' Info:', strong=True),
                fac.AntdParagraph(id = 'main-error-area-info'),
                fac.AntdIcon(icon='fc-upload'),
                fac.AntdText(' Input:', strong=True),
                fac.AntdParagraph(id='main-error-area-input'),
                fac.AntdIcon(icon='fc-share'),
                fac.AntdText(' Traceback:', strong=True),
                fac.AntdParagraph(id='main-error-area-traceback'),
            ],
            style={
                'maxHeight': '89vh',
                'overflowY': 'auto',
                'width':'16vw'
            },
        ),
        type='error',
        id='main-error-area',
        action=fac.AntdButton('╳', size='small', type='text', id='main-error-area-close'),
        showIcon=True,
        style={'display': 'none'}
    ),
    dcc.Store(id="key-pressed-events"),
    dcc.Store(id="userid"),
    dcc.Location(id='main-refresh', refresh=True),
    dcc.Interval(id="init-restore", interval=1, max_intervals=1),
    dcc.Interval(
        id='auth-interval',
        interval=1000,
    ),
    fac.AntdSpin(
        html.Div(id='main-loading-area'), 
        text='loading', 
        fullscreen=True,
    ),
    fileSelecter.get_box(),
    fac.AntdLayout(
        [
            fac.AntdSider(
                [
                    fac.AntdButton(
                        id='main-button-trigger',
                        icon=fac.AntdIcon(
                            id='main-icon-menuItem',
                            icon='antd-arrow-left',
                            style={'fontSize': '14px'},
                        ),
                        shape='circle',
                        type='text',
                        style={
                            'position': 'absolute',
                            'zIndex': 1,
                            'top': 25,
                            'right': -13,
                            'boxShadow': 'rgb(0 0 0 / 10%) 0px 4px 10px 0px',
                            'background': 'white',
                        },
                    ),
                    fac.AntdMenu(
                        id='main-menu-item',
                        persistence=True,
                        defaultSelectedKey='0',
                        menuItems=[],
                        mode='inline',
                        style={'height': '100%', 'overflow': 'hidden auto'},
                    ),
                ],
                id='main-sider-collapse',
                collapsible=True,
                collapsedWidth=60,
                trigger=None,
                style={'position': 'relative'},
            ),
            
            fac.AntdLayout(
                [
                    fac.AntdHeader(
                        [
                            fac.AntdTitle(
                                '', level=4, 
                                id='main-title-header',
                                style={'color': 'white', 'margin': '0'}
                            ),
                            fac.AntdSpace(
                                [
                                    fac.AntdText(
                                        children=None, 
                                        id='main-title-username',
                                        style={'color': 'white', 'margin': '0', 'fontSize':'18px'}
                                    ),
                                    fac.AntdPopconfirm(
                                        fac.AntdButton(
                                            type='link',
                                            icon=fac.AntdIcon(icon='md-power-settings-new'),
                                            iconPosition='end',
                                            style={'color': 'white', 'margin': '0', 'fontSize':'18px'}
                                        ),
                                        id='main-button-logout',
                                        locale='en-us',
                                        arrow='hide',
                                        okText='yes',
                                        placement='bottomLeft',
                                        title='Confirm Logout?'
                                    )
                                ],
                                style={'position': 'absolute', 'right': '50px'}
                            )
                        ],
                        style={
                            'display': 'flex',
                            'alignItems': 'center',
                            'backgroundColor': '#5F9EA0',
                            'height':'7vh'
                        },
                    ),
                    fac.AntdLayout(
                        [
                            fac.AntdSider(
                                id='main-sider-control',
                                width=250,
                                style={'backgroundColor': 'white'},
                            ),
                            fac.AntdLayout(
                                fac.AntdContent(
                                    fac.AntdCenter(
                                        id='main-center-content',
                                        style={'height': '100%', 'display':'flex', 'width':'100%', 'backgroundColor':'white'}
                                    ),
                                    style={'backgroundColor': 'white', 'display': 'flex', 'alignItems': 'center', 'justifyContent': 'center', 'minHeight': '93vh'},
                                )
                            ),
                        ]
                    ),
                ],
                id='main-layout-right'
            )
        ],  style={'height': '100vh'}
    )
], id='SeuPipe', style={'display':'none'})