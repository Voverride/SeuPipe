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
                        fac.AntdText('Select Project'),
                        fac.AntdTooltip(
                            fac.AntdSelect(
                                id='maskviewer-select-taskname',
                                debounceWait=300,
                                locale='en-us',
                                allowClear=False,
                                style={'width':'100%'}
                            ),
                            id='maskviewer-select-taskname-tooltip',
                            title=fac.AntdText('select project'),
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
                        fac.AntdText('Select Slice', id='mv-select-slice-label'),
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
                        fac.AntdText('Select Graph', id='mv-select-graph-label'),
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
                fac.AntdSpace(
                    [
                        fac.AntdText('Left Graph', id='mv-left-graph-label'),
                        fac.AntdTooltip(
                            fac.AntdSelect(
                                id='maskviewer-left-graph',
                                options=[],
                                value='watershed',
                                debounceWait=300,
                                locale='en-us',
                                allowClear=False,
                                disabled=True,
                                style={'width':'100%'}
                            ),
                            id='maskviewer-left-graph-tooltip',
                            title=fac.AntdText('select left graph'),
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
                        fac.AntdText('Right Graph', id='mv-right-graph-label'),
                        fac.AntdTooltip(
                            fac.AntdSelect(
                                id='maskviewer-right-graph',
                                options=[],
                                value='cellpose',
                                debounceWait=300,
                                locale='en-us',
                                allowClear=False,
                                disabled=True,
                                style={'width':'100%'}
                            ),
                            id='maskviewer-right-graph-tooltip',
                            title=fac.AntdText('select right graph'),
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
                            fac.AntdCheckbox(id='mv-checkbox-mask', label='mask', checked=False, disabled=True),
                            fac.AntdCheckbox(id='mv-checkbox-contour', label='contour', checked=True, disabled=True),
                        ],
                        size='large'
                    ),
                ),
                fac.AntdDivider(),
                fac.AntdSpace(
                    [
                        fac.AntdText('Select Export Type', id='mv-select-export-type-label'),
                        fac.AntdTooltip(
                            fac.AntdSelect(
                                id='maskviewer-select-export-type',
                                options=[],
                                debounceWait=300,
                                locale='en-us',
                                allowClear=False,
                                style={'width':'100%'}
                            ),
                            id='maskviewer-select-export-type-tooltip',
                            title=fac.AntdText('select export type'),
                            open=False,
                            color='white'
                        ),
                    ],
                    direction='vertical',
                    size='small',
                    style={'width':'100%'}
                ),
                fac.AntdButton(
                    'Export Data', 
                    type='primary',
                    id='maskviewer-button-exportData',
                    block=True,
                    icon=fac.AntdIcon(icon='antd-save'),
                    style={'backgroundColor':'#ca8269'}
                )
            ],
            size='large',
            direction='vertical',
            style={'width':'100%'},
        ),
        styles={'header': {'display': 'none'}},
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