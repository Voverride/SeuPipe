from dash import html
from dash import dcc
import feffery_antd_components as fac
from callbacks.maskviewer_c import *

control_panel = html.Div(
    fac.AntdCard(
        fac.AntdSpace(
            [
                dcc.Interval(id="init-restore-maskviewer", interval=1, max_intervals=1),
                dcc.Store(id='maskviewer-store-taskname', storage_type='local'),
                fac.AntdSpace(
                    [
                        fac.AntdText('Select Task'),
                        fac.AntdTooltip(
                            fac.AntdSelect(
                                id='maskviewer-select-taskname',
                                debounceWait=300,
                                locale='en-us',
                                allowClear=False,
                                style={'width':'100%'}
                            ),
                            id='maskviewer-select-taskname-tooltip',
                            title=fac.AntdText('select task'),
                            open=False,
                            color='white'
                        ),
                    ],
                    direction='vertical',
                    size='small',
                    style={'width':'100%'}
                ),
                fac.AntdSpace(
                    [
                        fac.AntdText('Select Slice'),
                        fac.AntdTooltip(
                            fac.AntdSelect(
                                id='maskviewer-select-slice',
                                debounceWait=300,
                                locale='en-us',
                                allowClear=False,
                                style={'width':'100%'}
                            ),
                            id='maskviewer-select-slice-tooltip',
                            title=fac.AntdText('select slice'),
                            open=False,
                            color='white'
                        ),
                    ],
                    direction='vertical',
                    size='small',
                    style={'width':'100%'}
                ),
                fac.AntdSpace(
                    [
                        fac.AntdText('Select Graph'),
                        fac.AntdTooltip(
                            fac.AntdSelect(
                                id='maskviewer-select-graph',
                                options=['registration', 'segmentation'],
                                value='registration',
                                debounceWait=300,
                                locale='en-us',
                                allowClear=False,
                                style={'width':'100%'}
                            ),
                            id='maskviewer-select-graph-tooltip',
                            open=False,
                            color='white'
                        ),
                    ],
                    direction='vertical',
                    size='small',
                    style={'width':'100%'}
                ),
                fac.AntdCenter(
                    fac.AntdSpace(
                        [
                            fac.AntdCheckbox(id='mv-checkbox-mask', label='mask', checked=True, disabled=True),
                            fac.AntdCheckbox(id='mv-checkbox-contour', label='contour', disabled=True),
                        ],
                        size='large'
                    ),
                )
            ],
            size='large',
            direction='vertical',
            style={'width':'100%'},
        ),
        headStyle={'display': 'none'},
        style={'height':'93vh', 'maxHeight': '93vh','overflowY': 'auto'},
    ), 
    id='segmentation-control_panel',
)

content_panel = html.Div(
    [
        fac.AntdCenter(
            fac.AntdSplitter(
                items=[
                    {
                        'children': fac.AntdCenter(
                            dcc.Graph(
                                id="mv-graph-left", 
                                config={'displaylogo':False}, 
                                style={'display':'block', 'height':'90vh', 'width':'100%', 'visibility':'hidden'}
                            ),
                        ),
                        'collapsible': True,
                    },
                    {
                        'children': fac.AntdCenter(
                            dcc.Graph(
                                id="mv-graph-right", 
                                config={'displaylogo':False}, 
                                style={'display':'block', 'height':'90vh', 'width':'100%', 'visibility':'hidden'}
                            ),
                        ),
                        'collapsible': True,
                    }    
                ],
            ),
        )
    ],
    id='mv-content_panel',
    style={'width':'calc(100vw - 480px)', 'height':'100%'},
)