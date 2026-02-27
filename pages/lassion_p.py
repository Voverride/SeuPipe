from dash import html, dcc
import feffery_antd_components as fac
from callbacks.lassion_c import *

control_panel = html.Div(
    fac.AntdCard(
        fac.AntdSpace(
            [
                fac.AntdDivider('Dataset Selection', innerTextOrientation='left', id='lasso-control-dataset-divider'),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        options=datasets,
                        value=initial_dataset,
                        id='lasso-select-dataset',
                        allowClear=False,
                        locale='en-us',
                        style={'width': '100%'},
                    ),
                    id='lasso-select-dataset-tooltip',
                    title=fac.AntdText('select dataset'), color='white'
                ),
                fac.AntdDivider('Spot Properties', innerTextOrientation='left'),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        options=[1,2,3,4,5,6,7,8,9,10],
                        value=5,
                        id='lasso-select-spotSize',
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
                                id='lasso-select-borderWidth',
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
                                id='lasso-colorPicker-boarderColor',
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
        styles={'header': {'display': 'none'}},
        style={'height':'93vh', 'maxHeight': '93vh','overflowY': 'auto'},
    ), 
    id='lasso-control_panel',
)

content_panel = html.Div(
    fac.AntdSpace(
        [
            fac.AntdSpace(
                [
                    dcc.Download(id="lasso-download-exportData"),
                    fac.AntdButton(
                        'Download CSV', 
                        type='primary',
                        id='lasso-button-download',
                        icon=fac.AntdIcon(icon='antd-download'),
                        style={'backgroundColor':"#698aab"}
                    )
                ],
                size='middle',
                style={'marginTop':'24.5px'},
            ),
            fac.AntdCenter(
                dcc.Graph(
                    id="lasso-graph-result", 
                    config={'displaylogo':False}, 
                    style={'height':'calc(90vh - 50px)', 'width':'100%'}
                ),
            )
        ],
        direction='vertical',
        size='middle',
        style={'width':'calc(100vw - 480px)', 'height':'100%'}
    ),
    id='lasso-content-panel',
    style={'width':'100%', 'height':'100%', 'display':'flex', 'justifyContent':'center'}
)