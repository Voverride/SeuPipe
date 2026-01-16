from dash import html
from dash import dcc
import feffery_antd_components as fac
from callbacks.maskviewer_c import *

control_panel = html.Div(
    fac.AntdCard(
        fac.AntdSpace(
            [
                fac.AntdSpace(
                    [
                        fac.AntdText('Select Project'),
                        fac.AntdTooltip(
                            fac.AntdSelect(
                                id='mv-select-project',
                                debounceWait=300,
                                locale='en-us',
                                allowClear=False,
                                style={'width':'100%'}
                            ),
                            id='mv-select-project-tooltip',
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
                                id='mv-select-slice',
                                debounceWait=300,
                                locale='en-us',
                                allowClear=False,
                                style={'width':'100%'}
                            ),
                            id='mv-select-slice-tooltip',
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
                        fac.AntdText('Select Figure', id='mv-select-figure-label'),
                        fac.AntdSpace(
                            [
                                fac.AntdTooltip(
                                    fac.AntdSelect(
                                        id='mv-select-figure',
                                        options=[],
                                        debounceWait=300,
                                        locale='en-us',
                                        allowClear=False,
                                        style={'width':'160px'}
                                    ),
                                    id='mv-select-figure-tooltip',
                                    title=fac.AntdText('select figure'),
                                    open=False,
                                    color='white'
                                ),
                                fac.AntdTooltip(
                                    fac.AntdButton(
                                        id='mv-button-add-contrast',
                                        size='small',
                                        icon=fac.AntdIcon(icon='antd-plus', style={'color': "#BDBDBD"}),
                                        style={'width':'30px', 'height':'30px'}
                                    ),
                                    title=fac.AntdText('add contrast figure'),
                                    color='white'
                                )
                            ],
                            size='small',
                            style={'width':'100%'}
                        ),
                    ],
                    direction='vertical',
                    size='small',
                    style={'width':'100%'}
                ),
                fac.AntdSpace(
                    [
                        fac.AntdText('Contrast Figure', id='mv-contrast-figure-label'),
                        fac.AntdSpace(
                            [
                                fac.AntdTooltip(
                                    fac.AntdSelect(
                                        id='mv-select-contrast',
                                        options=[],
                                        debounceWait=300,
                                        locale='en-us',
                                        allowClear=False,
                                        style={'width':'160px'}
                                    ),
                                    id='mv-select-contrast-tooltip',
                                    title=fac.AntdText('select contrast figure'),
                                    open=False,
                                    color='white'
                                ),
                                fac.AntdTooltip(
                                    fac.AntdButton(
                                        id='mv-button-remove-contrast',
                                        size='small',
                                        icon=fac.AntdIcon(icon='antd-minus', style={'color': "#BDBDBD"}),
                                        style={'width':'30px', 'height':'30px'}
                                    ),
                                    title=fac.AntdText('remove contrast figure'),
                                    color='white'
                                )
                            ],
                            size='small',
                            style={'width':'100%'}
                        ),
                    ],
                    id='mv-contrast-figure-space',
                    direction='vertical',
                    size='small',
                    style={'width':'100%', 'display':'none'}
                ),
                fac.AntdCenter(
                    fac.AntdSpace(
                        [
                            fac.AntdCheckbox(id='mv-checkbox-mask', label='mask', checked=False),
                            fac.AntdCheckbox(id='mv-checkbox-contour', label='contour', checked=True),
                        ],
                        size='large'
                    ),
                ),
                fac.AntdDivider(),
                fac.AntdSpace(
                    [
                        fac.AntdText('Select Export Project'),
                        fac.AntdTooltip(
                            fac.AntdSelect(
                                id='mv-select-export-project',
                                debounceWait=300,
                                locale='en-us',
                                allowClear=False,
                                style={'width':'100%'}
                            ),
                            id='mv-select-export-project-tooltip',
                            title=fac.AntdText('select export project'),
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
                        fac.AntdText('Select Export Type', id='mv-select-export-type-label'),
                        fac.AntdTooltip(
                            fac.AntdSelect(
                                id='mv-select-export-type',
                                options=[],
                                debounceWait=300,
                                locale='en-us',
                                allowClear=False,
                                style={'width':'100%'}
                            ),
                            id='mv-select-export-type-tooltip',
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
                    id='mv-button-exportData',
                    block=True,
                    icon=fac.AntdIcon(icon='antd-save'),
                    style={'backgroundColor':'#ca8269'}
                ),
                dcc.Store(id='mv-store-figure-type', data={'left':False, 'right':False, 'contrast':False}),
            ],
            size='large',
            direction='vertical',
            style={'width':'100%'},
        ),
        styles={'header': {'display': 'none'}},
        style={'height':'93vh', 'maxHeight': '93vh','overflowY': 'auto'},
    ), 
    id='mv-control_panel',
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
                                style={'display':'block', 'height':'90vh', 'width':'100%', 'visibility':'visible'}
                            ),
                        ),
                        'collapsible': True,
                    } 
                ],
                id='mv-splitter',
            ),
        )
    ],
    id='mv-content_panel',
    style={'width':'calc(100vw - 480px)', 'height':'100%'},
)