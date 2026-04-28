from callbacks.alignment_c import *
from dash import html
from .components.manualAdjust import ManualAdjustDrawer
import feffery_antd_components as fac
from dash import dcc

new_project_modal = fac.AntdModal(
    fac.AntdSpace(
        [
            fac.AntdDivider(),
            fac.AntdSpace(
                [
                    fac.AntdButton(
                        'Import Data', 
                        id='ali-button-importdata',
                        type='primary', 
                        icon=fac.AntdIcon(icon='antd-cloud-download'), 
                        style={'backgroundColor':'#5383c3', 'width':'130px'}
                    ),
                    dcc.Store(id='ali-store-importdata', data=None),
                    fac.AntdTooltip(
                        fac.AntdSelect(
                            placeholder='x field',
                            optionFilterMode='case-insensitive',
                            allowClear=False,
                            options=[],
                            style={'width':'100px'},
                            id='ali-select-x',
                            locale='en-us',
                        ),
                        title=fac.AntdText('select x field'), color='white'
                    ),
                    fac.AntdTooltip(
                        fac.AntdSelect(
                            placeholder='y field',
                            optionFilterMode='case-insensitive',
                            allowClear=False,
                            options=[],
                            style={'width':'100px'},
                            id='ali-select-y',
                            locale='en-us',
                        ),
                        title=fac.AntdText('select y field'), color='white'
                    ),
                    fac.AntdTooltip(
                        fac.AntdSelect(
                            placeholder='z field',
                            optionFilterMode='case-insensitive',
                            allowClear=False,
                            options=[],
                            style={'width':'100px'},
                            id='ali-select-z',
                            locale='en-us',
                        ),
                        title=fac.AntdText('select z field'), color='white'
                    )
                ],
                size='middle'
            ),
            fac.AntdDivider(),
            fac.AntdSpace(
                [
                    fac.AntdInput(
                        placeholder='input project name',
                        id='ali-projectname',
                        autoComplete="off",
                        status='error',
                        style={'width':'330px'}
                    ),
                    fac.AntdButton(
                        'Submit', 
                        id='ali-button-submit',
                        type='primary', 
                        icon=fac.AntdIcon(icon='md-launch'), 
                        style={'backgroundColor': '#d0826c', 'width': '130px'}
                    ),
                ],
                size='middle',
                style={'marginBottom':'15px'}
            ),
        ],
        size='middle',
        direction='vertical',
        style={'width':'100%'},
    ),
    id='ali-modal-newproject', 
    title='Create New Project',
    mask=False,
    width=530,
    maskClosable=False,
    visible=False
)

control_panel = html.Div(
    fac.AntdCard(
        fac.AntdSpace(
            [
                fac.AntdButton(
                    'Create New Project', 
                    type='primary',
                    id='ali-button-newproject',
                    icon=fac.AntdIcon(icon='antd-plus'),
                    style={'backgroundColor':'#698aab', 'width': '100%'}
                ),
                fac.AntdDivider(),
                fac.AntdText('Select Project', id='ali-text-project'),
                fac.AntdSpace(
                    [
                        fac.AntdTooltip(
                            fac.AntdSelect(
                                id='ali-select-project',
                                placeholder='Select Project',
                                debounceWait=300,
                                locale='en-us',
                                allowClear=False,
                                style={'width':'160px'}
                            ),
                            id='ali-select-project-tooltip',
                            open=False,
                            title=fac.AntdText('Select Project'), 
                            color='white'
                        ),
                        fac.AntdTooltip(
                            fac.AntdButton(
                                id='ali-button-add-contrast',
                                size='small',
                                icon=fac.AntdIcon(id='ali-icon-add-contrast', icon='antd-plus', style={'color': "#BDBDBD"}),
                                style={'width':'30px', 'height':'30px'}
                            ),
                            title=fac.AntdText('Add Contrast Figure', id='ali-text-add-contrast'), 
                            color='white'
                        )
                    ]
                ),
                fac.AntdDivider(),
                fac.AntdText('Select Model', id='ali-text-model'),
                fac.AntdSelect(
                    optionFilterMode='case-insensitive',
                    allowClear=False,
                    value='paste2',
                    options=['paste2', 'paste1'],
                    id='ali-select-model',
                    locale='en-us',
                    style={'width': '100%', 'marginBottom':'15px'},
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
                        id='ali-radio-device',
                        style={'width': '100%'},
                    ),
                ),
                fac.AntdDivider(),
                fac.AntdButton(
                    'Manual Adjust', 
                    type='primary',
                    id='ali-button-manualAdjust',
                    block=True,
                    icon=fac.AntdIcon(icon='bi-layer'),
                    style={'backgroundColor':'#a58f86'}
                ),
                fac.AntdDivider(),
                fac.AntdButton(
                    'Export Data', 
                    type='primary',
                    id='ali-button-exportData',
                    block=True,
                    icon=fac.AntdIcon(icon='antd-save'),
                    style={'backgroundColor':'#ca8269'}
                )
            ],
            direction='vertical',
            style={'width':'100%'},
        ),
        styles={'header': {'display': 'none'}},
        style={'height':'93vh', 'maxHeight': '93vh','overflowY': 'auto'},
    ), 
    id='ali-control_panel',
)

content_tabs = fac.AntdTabs(
    items=[
        {
            'key': 'ProjectInfo',
            'label': 'ProjectInfo',
            'children': fac.AntdSpace(
                [
                    dcc.Store(id='ali-store-figureScale', data={}),
                    dcc.Store(id='ali-store-leftLayout', data=None),
                    dcc.Store(id='ali-store-rightLayout', data=None),
                    fac.AntdSpace(
                        [
                            fac.AntdButton(
                                'Start Project', 
                                type='primary',
                                id='ali-start-project',
                                icon=fac.AntdIcon(icon='antd-carry-out'),
                                style={'backgroundColor':'#7d8a70'}
                            ),
                            fac.AntdPopconfirm(
                                fac.AntdButton(
                                    'Delete Project', 
                                    id='ali-delete-project',
                                    type='primary',
                                    icon=fac.AntdIcon(icon='antd-delete'),
                                    style={'backgroundColor':'#ca8269'}
                                ),
                                id='ali-delete-project-confirm',
                                locale='en-us',
                                arrow='hide',
                                okText='yes',
                                placement='bottomLeft',
                                title='Confirm Delete?'
                            ),
                            fac.AntdSpace(
                                [
                                    fac.AntdButton(
                                        'Show Bug',
                                        id='ali-button-showBug',
                                        type='primary', 
                                        icon=fac.AntdIcon(icon='antd-bug'), 
                                        style={'backgroundColor':'#a87a76', 'display':'block'}
                                    ),
                                    fac.AntdText("🚨 Error occurred during execution. Click Show Bug for details !", style={'color':'#a87a76'}, strong=True)
                                ],
                                id='ali-bug-panel',
                                size='middle',
                                style={'display':'none'}
                            ),
                            html.Div(children=None, id='refresh-ali-status', style={'display':'none'}),
                        ],
                        size='middle',
                        style={'width':'100%'},
                    ),
                    fac.AntdCard(
                        fac.AntdSteps(
                            id='ali-alignment-steps',
                            steps=[
                                {
                                    'title': 'Preprocess',
                                    'status': 'wait',
                                },
                                {
                                    'title': 'Alignment',
                                    'status': 'wait',
                                },
                                {
                                    'title': 'Postprocess',
                                    'status': 'wait',
                                }
                            ]
                        ),
                        title='Alignment Steps',
                    ),
                    fac.AntdTable(
                        columns=[
                            {'title': 'creator', 'dataIndex': 'creator', 'width':'10%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'date', 'dataIndex': 'date', 'width':'15%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'initialData', 'dataIndex': 'initialData', 'width':'45%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'xField', 'dataIndex': 'x', 'width':'10%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'yField', 'dataIndex': 'y', 'width':'10%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'zField', 'dataIndex': 'z', 'width':'10%', 'renderOptions': {'renderType': 'ellipsis'}},
                        ],
                        data=[{}],
                        style={'width': '100%', 'display':'inline-block'},
                        id='ali-table-metadata',
                        locale='en-us',
                        bordered=True,
                        pagination=False,
                    )
                ],
                direction='vertical',
                size='large',
                style={'height': 'calc(93vh - 75px)', 'width':'100%', 'overflowY': 'auto'},
            ),
        },
        {
            'key': 'Figure',
            'label': 'Figure',
            'children': fac.AntdSplitter(
                items=[
                    {
                        'children': fac.AntdCenter(
                            dcc.Graph(
                                # figure=aligned,
                                config={'displaylogo':False}, 
                                id='ali-graph-left', 
                                style={'height': '100%', 'width': '100%'},
                            ),
                            style={'height': '100%', 'width': '100%'}
                        ),
                        'size': '100%',
                        'collapsible': False,
                        'resizable': False,
                    },
                    {
                        'children': fac.AntdCenter(
                            dcc.Graph(
                                # figure=origion,
                                figure=hiddenGraph,
                                config={'displaylogo':False}, 
                                id='ali-graph-right', 
                                style={'height': '100%', 'width': '100%'},
                            ),
                            style={'height': '100%', 'width': '100%'}
                        ),
                        # 'size': '50%',
                        'collapsible': False,
                        'resizable': False,
                    }
                ],
                id='ali-splitter-figure',
                style={'height': 'calc(93vh - 65px)', 'width':'100%'}
            )
        }
    ],
    id='ali-content-tabs'
)

content_panel = html.Div(
    [
        new_project_modal,
        ManualAdjustDrawer,
        content_tabs,
        fac.AntdSpace(
            [
                fac.AntdBadge(dot=True, color='#5F9EA0'),
                fac.AntdText(
                    id='ali-live-username',
                    style={'color':'#5F9EA0'}
                ),
            ],
            id='ali-live-userspace',
            style={
                'opacity': 0,
                'position':'absolute',
                'top':'50px',
                'left':'1px',
            }
        )
    ],
    id='ali-content_panel',
    style={'width':'98%', 'height':'100%', 'position':'relative'},
)