from dash import html, dcc
import feffery_antd_components as fac
from callbacks.visiualization_c import *

control_panel = html.Div(
    fac.AntdCard(
        fac.AntdSpace(
            [
                fac.AntdDivider('Dataset Selection', innerTextOrientation='left', id='visiualization-control-dataset-divider'),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        options=datasets,
                        value=initial_dataset,
                        id='visiualization-select-dataset',
                        allowClear=False,
                        locale='en-us',
                        style={'width': '100%'},
                    ),
                    id='visiualization-select-dataset-tooltip',
                    title=fac.AntdText('select dataset'), color='white'
                ),
                fac.AntdDivider('Spot Properties', innerTextOrientation='left'),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        options=[1,2,3,4,5,6,7,8,9,10],
                        value=5,
                        id='visiualization-select-spotSize',
                        allowClear=False,
                        locale='en-us',
                        style={'width': '100%'},
                    ),
                    title=fac.AntdText('Spot Size'), color='white'
                ),
                fac.AntdSpace(
                    [
                        fac.AntdTooltip(
                            fac.AntdSelect(
                                options=[0,1],
                                allowClear=False,
                                value=1,
                                id='visiualization-select-borderWidth',
                                locale='en-us',
                                style={'width': '152px'},
                            ), 
                            title=fac.AntdText('Border Width'), color='white'
                        ),
                        fac.AntdTooltip(
                            fac.AntdColorPicker(
                                disabledAlpha=False,
                                locale='en-us',
                                presets=[
                                    {'colors': ['#0d0015', '#474a4d', '#1c305c', '#640125'], 'label': 'presets'}
                                ],
                                value='#1c305c',
                                id='visiualization-colorPicker-boarderColor',
                            ),
                            title=fac.AntdText('Border Color'), color='white'
                        ),
                    ],
                    size='middle',
                    style={'width':'100%'}
                ),
            ],
            size='middle',
            direction='vertical',
            style={'width':'100%'},
        ),
        headStyle={'display': 'none'},
        style={'height':'93vh', 'maxHeight': '93vh','overflowY': 'auto'},
    ), 
    id='visiualization-control_panel',
)

content_panel = html.Div(
    [
        fac.AntdCenter(
            dcc.Graph(
                config={'displaylogo':False}, 
                id='visiualization-graph-result', 
                style={'height':'100%', 'width': '100%', 'display':'none'},
            ),
            style={'height': '100%', 'width': '100%'}
        ),
        dcc.Store(id='visiualization-store-relayoutData', data={}),
    ],
    id='visiualization-content-panel',
    style={'width':'calc(100vw - 480px)', 'height':'100%', 'position':'relative'}
)