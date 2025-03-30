from callbacks.alignment_c import *
from dash import html
from .components.graphSetting import GraphSettingDrawer
import feffery_antd_components as fac
from dash import dcc

control_panel = html.Div(
    fac.AntdCard(
        fac.AntdSpace(
            [
                dcc.Interval(id="init-restore-alignment", interval=1, max_intervals=1),
                dcc.Interval(interval=1000, disabled=True, id='alignment-interval'),
                dcc.Interval(id="alignment-event-loop", interval=1000),
                fac.AntdButton(
                    'Import Data', 
                    type='primary',
                    id='alignment-button-importData',
                    icon=fac.AntdIcon(icon='antd-cloud-upload'),
                    style={'backgroundColor':'#698aab', 'width': '100%'}
                ),
                fac.AntdDivider(),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        placeholder='select x field',
                        optionFilterMode='case-insensitive',
                        allowClear=False,
                        options=[],
                        id='alignment-select-x',
                        locale='en-us',
                        style={'width': '100%'},
                    ),
                    title=fac.AntdText('select x field'), color='white'
                ),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        placeholder='select y field',
                        optionFilterMode='case-insensitive',
                        allowClear=False,
                        options=[],
                        id='alignment-select-y',
                        locale='en-us',
                        style={'width': '100%'},
                    ),
                    title=fac.AntdText('select y field'), color='white'
                ),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        placeholder='select z field',
                        optionFilterMode='case-insensitive',
                        allowClear=False,
                        options=[],
                        id='alignment-select-z',
                        locale='en-us',
                        style={'width': '100%'},
                    ),
                    title=fac.AntdText('select z field'), color='white'
                ),
                fac.AntdButton(
                    'Plot Figure', 
                    type='primary',
                    id='alignment-button-plotFig',
                    block=True,
                    icon=fac.AntdIcon(icon='antd-bar-chart'),
                    style={'backgroundColor':'#867ba9'}
                ),
                fac.AntdDivider(),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        placeholder='select model',
                        optionFilterMode='case-insensitive',
                        allowClear=False,
                        value='paste1',
                        options=['paste1', 'paste2'],
                        id='alignment-select-model',
                        locale='en-us',
                        style={'width': '100%'},
                    ),
                    title=fac.AntdText('select model'), color='white'
                ),
                fac.AntdCenter(
                    fac.AntdRadioGroup(
                        options=[
                            {
                                'label': 'CPU',
                                'value': 'CPU'
                            },
                            {
                                'label': 'GPU',
                                'value': 'GPU'
                            }      
                        ],
                        defaultValue='CPU',
                        block=True,
                        id='alignment-radio-device',
                        style={'width': '100%'},
                    ),
                ),
                fac.AntdSpace(
                    [
                        fac.AntdButton(
                            'Align Slices', 
                            type='primary',
                            id='alignment-button-alignSlices',
                            icon=fac.AntdIcon(icon='pi-stack'),
                            style={'backgroundColor':'#5F9EA0', 'width':'160px'}
                        ),
                        fac.AntdTooltip(
                            fac.AntdButton(
                                type='primary',
                                id='alignment-button-showProgress',
                                icon=fac.AntdIcon(icon='antd-monitor'),
                                style={'backgroundColor':'#5F9EA0'}
                            ),
                            title=fac.AntdText('show alignment progress'), color='white'
                        ),
                    ]
                ),
                fac.AntdButton(
                    'Graph Setting', 
                    type='primary',
                    id='alignment-button-graphSetting',
                    block=True,
                    icon=fac.AntdIcon(icon='antd-setting'),
                    style={'backgroundColor':'#a58f86'}
                ),
                fac.AntdDivider(),
                fac.AntdButton(
                    'Export Data', 
                    type='primary',
                    id='alignment-button-exportData',
                    block=True,
                    icon=fac.AntdIcon(icon='antd-save'),
                    style={'backgroundColor':'#ca8269'}
                )
            ],
            size='middle',
            direction='vertical',
            style={'width':'100%'},
        ),
        headStyle={'display': 'none'},
        style={'height':'93vh', 'maxHeight': '93vh','overflowY': 'auto'},
    ), 
    id='alignment-control_panel',
)

content_panel = html.Div(
    [
        GraphSettingDrawer,
        fac.AntdPopupCard(
            [
                fac.AntdSpace(
                    [
                        fac.AntdTimeline(
                            items=[
                                {
                                    'content': html.Div(
                                        [
                                            html.Div(
                                                [
                                                    fac.AntdText('Creator:', strong=True, style={'marginRight':'6px'}),
                                                    fac.AntdText(type='success', id='alignment-timeline-creator'),  
                                                ],
                                                style={'marginRight':'20px'}
                                            ),
                                            fac.AntdButton(
                                                'Show Bug',
                                                id='alignment-button-showBug',
                                                size='small', 
                                                type='primary', 
                                                icon=fac.AntdIcon(icon='antd-bug'), 
                                                style={'backgroundColor':'#bb5548', 'display':'none'}
                                            )           
                                        ],
                                        style={'width':'100%', 'display': 'flex'}
                                    ),
                                    'icon':fac.AntdAvatar(size='small'),
                                },
                                {
                                    'content': 'Data Validation',
                                    'color':'gray',
                                    'icon': fac.AntdIcon(icon='md-schedule', id='step1')
                                },
                                {
                                    'content': 'Data Preparation',
                                    'color':'gray',
                                    'icon': fac.AntdIcon(icon='md-schedule', id='step2')
                                },
                                {
                                    'content': fac.AntdSpace(
                                        [
                                            'Slice Alignment',
                                            fac.AntdProgress(percent=0, id='alignment-alipercent', style={'width': '200px'}),
                                        ],
                                        size='middle'
                                    ),
                                    'color':'gray',
                                    'icon': fac.AntdIcon(icon='md-schedule', id='step3')
                                },
                                {
                                    'content': 'Stack Slices',
                                    'color':'gray',
                                    'icon': fac.AntdIcon(icon='md-schedule', id='step4')
                                },
                                {
                                    'content': 'Update Data',
                                    'color':'gray',
                                    'icon': fac.AntdIcon(icon='md-schedule', id='step5')
                                }
                            ],
                            id = 'alignment-timeline-alistatus'
                        )
                    ],
                    direction='vertical',
                    size='middle',
                    style={'width': '100%', 'marginTop':'30px'},
                ),
                fac.AntdParagraph(
                    [
                        fac.AntdIcon(icon='antd-notification-two-tone'),
                        ' The aligned coordinates (x and y) are stored in ',
                         fac.AntdText('adata.obs', code=True, type='success'), 'under the columns ',
                         fac.AntdText('x_aligned', code=True, type='success'), ' and ', fac.AntdText('y_aligned', code=True, type='success'), 
                         '. The final spatially registered coordinates are saved in ',
                         fac.AntdText("adata.obsm['X_spatial_registered']", code=True, type='success')
                    ],
                    id='alignment-paragraph-fieldnotice',
                    style={'marginTop':'-20px', 'display':'none'}
                )
            ],
            id='alignment-popupcard-alignTask',
            title='Alignment Progress',
            width='500px',
            closable=True,
            destroyOnClose=False,
            draggable=True,
            visible=False,
        ),
        fac.AntdCenter(
            fac.AntdSplitter(
                items=[
                    {
                        'children': fac.AntdCenter(
                            dcc.Graph(
                                id="alignment-graph-aligned", 
                                config={'displaylogo':False}, 
                                style={'display':'block', 'height':'90vh', 'width':'100%', 'visibility':'hidden'}
                            ),
                        ),
                        'collapsible': True,
                    },
                    {
                        'children': fac.AntdCenter(
                            dcc.Graph(
                                id="alignment-graph-origion", 
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
    id='alignment-content_panel',
    style={'width':'95%', 'height':'100%'},
)