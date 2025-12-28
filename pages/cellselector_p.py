from dash import html
import feffery_antd_components as fac
from callbacks.cellselector_c import *
from dash import dcc

new_project_modal = fac.AntdModal(
    fac.AntdSpace(
        [
            fac.AntdDivider(),
            fac.AntdSelect(
                id='sel-select-project-modal',
                placeholder='Select Segmentation Project',
                debounceWait=300,
                locale='en-us',
                allowClear=False,
                style={'width':335}
            ),
            fac.AntdSpace(
                [
                    fac.AntdSelect(
                        id='sel-select-result-modal',
                        placeholder='Select Result',
                        debounceWait=300,
                        locale='en-us',
                        allowClear=False,
                        style={'width':335}
                    ),
                    fac.AntdButton(
                        'Submit',
                        type='primary',
                        id='sel-button-submitProject-modal',
                        icon=fac.AntdIcon(icon='md-launch'),
                        style={'backgroundColor':'#d0826c', 'width':100}
                    ),
                ],
                style={'width':'100%', 'marginBottom':'30px'},
                size='middle'
            )
        ],
        size='middle',
        direction='vertical',
        style={'width':'100%'},
    ),
    id='sel-modal-newproject', 
    title='Create New Project',
    mask=False,
    width=500,
    maskClosable=False,
    visible=False
)

control_panel = html.Div(
    fac.AntdCard(
        fac.AntdSpace(
            [
                new_project_modal,
                html.Div(id='refresh-sel-status', style={'display':'none'}),
                fac.AntdButton(
                    'Create New Project', 
                    type='primary',
                    id='sel-button-newProject',
                    icon=fac.AntdIcon(icon='antd-plus'),
                    style={'backgroundColor':'#698aab', 'width': '100%'}
                ),
                fac.AntdDivider(),
                fac.AntdText('Select Project', id='sel-text-project'),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        id='sel-select-project',
                        placeholder='Select Project',
                        debounceWait=300,
                        locale='en-us',
                        allowClear=False,
                        style={'width':'100%'}
                    ),
                    id='sel-select-project-tooltip',
                    open=False,
                    title=fac.AntdText('Select Project'), 
                    color='white'
                ),
                fac.AntdText('Select Slice', id='sel-text-slice'),
                fac.AntdSelect(
                    id='sel-select-slice',
                    placeholder='Select Slice',
                    debounceWait=300,
                    locale='en-us',
                    allowClear=False,
                    style={'width':'100%'}
                ),
                fac.AntdDivider(),
                fac.AntdButton(
                    'Export Data', 
                    type='primary',
                    id='sel-button-exportData',
                    block=True,
                    icon=fac.AntdIcon(icon='antd-save'),
                    style={'backgroundColor':'#ca8269'}
                )
            ],
            size='middle',
            direction='vertical',
            style={'width':'100%'},
        ),
        styles={'header': {'display': 'none'}},
        style={'height':'93vh', 'maxHeight': '93vh','overflowY': 'auto'},
    ), 
    id='sel-control_panel',
)

content_panel = html.Div(
    fac.AntdSpace(
        [
            fac.AntdSpace(
                [
                    fac.AntdButton(
                        'Retain', 
                        type='primary',
                        id='sel-button-retain',
                        disabled=True,
                        icon=fac.AntdIcon(icon='antd-check-circle'),
                        style={'backgroundColor':'#7d8a70'}
                    ),
                    fac.AntdButton(
                        'Remove', 
                        type='primary',
                        id='sel-button-remove',
                        disabled=True,
                        icon=fac.AntdIcon(icon='antd-close-circle'),
                        style={'backgroundColor':"#b39863"}
                    ),
                    fac.AntdButton(
                        'Reset', 
                        type='primary',
                        id='sel-button-reset',
                        icon=fac.AntdIcon(icon='antd-undo'),
                        style={'backgroundColor':'#698aab'}
                    ),
                    fac.AntdPopconfirm(
                        fac.AntdButton(
                            'Delete Project', 
                            id='sel-delete-project',
                            type='primary',
                            icon=fac.AntdIcon(icon='antd-delete'),
                            style={'backgroundColor':'#a87a76'}
                        ),
                        id='sel-delete-project-confirm',
                        locale='en-us',
                        arrow='hide',
                        okText='yes',
                        placement='bottomLeft',
                        title='Confirm Delete?'
                    )
                ],
                size='middle',
                style={'marginTop':'24.5px'},
            ),
            fac.AntdCenter(
                dcc.Graph(
                    id="sel-cell-graph", 
                    config={'displaylogo':False}, 
                    style={'height':'calc(90vh - 50px)', 'width':'100%'}
                ),
            )
        ],
        direction='vertical',
        size='middle',
        style={'width':'calc(100vw - 480px)', 'height':'100%'}
    ),
    id='sel-content_panel',
    style={'width':'100%', 'height':'100%', 'display':'flex', 'justifyContent':'center'}
)