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
            fac.AntdSelect(
                id='sel-select-result-modal',
                placeholder='Select Result',
                debounceWait=300,
                locale='en-us',
                allowClear=False,
                style={'width':335}
            ),
            fac.AntdSpace(
              [
                fac.AntdText('Enable Clustering'),
                fac.AntdSwitch(checked=False, id='sel-switch-clustering-modal')
              ]
            ),
            fac.AntdTooltip(
              fac.AntdInputNumber(
                  min=0.1,
                  max=3.0,
                  step=0.1,
                  precision=1,
                  disabled=True,
                  addonBefore='Resolution',
                  placeholder='Input Resolution',
                  id='sel-input-resolution-modal',
                  style={'width': 335},
              ),
              title=fac.AntdText('Resolution controls the granularity of clustering. Higher values result in more clusters.'), 
              color='white'
            ),
            fac.AntdSpace(
                [
                  fac.AntdTooltip(
                    fac.AntdInputNumber(
                        min=0,
                        max=20,
                        step=1,
                        precision=0,
                        disabled=True,
                        addonBefore='n_iteration',
                        placeholder='Input Iteration',
                        id='sel-input-iteration-modal',
                        style={'width': 335},
                    ),
                    title=fac.AntdSpace(
                        [
                            fac.AntdText('Leiden Iterations:', strong=True),
                            fac.AntdText('🔸0: Run until convergence'),
                            fac.AntdText('🔸others: Number of iterations')
                        ],
                        direction='vertical',
                    ),
                    color='white'
                  ),
                  fac.AntdButton(
                      'Submit',
                      type='primary',
                      id='sel-button-submitProject-modal',
                      icon=fac.AntdIcon(icon='md-launch'),
                      style={'backgroundColor':'#d0826c', 'width':100}
                  ),
                ],
                size='middle',
                style={'width':'100%', 'marginBottom':'30px'},
            ),
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

clustring_status_modal = fac.AntdModal(
    fac.AntdSpace(
        [
            fac.AntdDivider(),
            fac.AntdTimeline(
                items=[
                    {
                        'content': html.Div(
                            [
                                html.Div(
                                    [
                                        fac.AntdText('Creator:', strong=True, style={'marginRight':'6px'}),
                                        fac.AntdText(type='success', id='sel-cluster-creator'),  
                                    ],
                                    style={'marginRight':'20px'}
                                ),
                                fac.AntdButton(
                                    'Show Bug',
                                    id='sel-cluster-showBug',
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
                        'content': 'Data Preparation',
                        'color':'gray',
                        'icon': fac.AntdIcon(icon='md-schedule', id='sel-cluster-step1')
                    },
                    {
                        'content': fac.AntdSpace(
                            [
                                'Generate AnnData',
                                fac.AntdProgress(percent=0, id='sel-cluster-percent', style={'width': '200px'}),
                            ],
                            size='middle'
                        ),
                        'color':'gray',
                        'icon': fac.AntdIcon(icon='md-schedule', id='sel-cluster-step2')
                    },
                    {
                        'content': 'Preprocess',
                        'color':'gray',
                        'icon': fac.AntdIcon(icon='md-schedule', id='sel-cluster-step3')
                    },
                    {
                        'content': 'Clustering',
                        'color':'gray',
                        'icon': fac.AntdIcon(icon='md-schedule', id='sel-cluster-step4')
                    },
                    {
                        'content': 'Export Results',
                        'color':'gray',
                        'icon': fac.AntdIcon(icon='md-schedule', id='sel-cluster-step5')
                    }
                ],
                id = 'sel-cluster-status'
            )   
        ],
        size='middle',
        direction='vertical',
        style={'width':'100%'},
    ),
    id='sel-modal-clusteringStatus', 
    title='Clustering Status',
    mask=False,
    width=500,
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
                fac.AntdText('Figure Configuration'),
                fac.AntdSpace(
                    [
                        fac.AntdTooltip(
                            fac.AntdSlider(id='sel-select-spotSize', min=5, max=20, defaultValue=10, style={'width': '145px'}),
                            title=fac.AntdText('select spot size'), color='white'
                        ),
                        fac.AntdTooltip(
                            fac.AntdColorPicker(
                                disabledAlpha=False,
                                locale='en-us',
                                presets=[
                                    {'colors': ["#defcff", "#ffdbfc", "#ffcac5"], 'label': 'presets'}
                                ],
                                value='#defcff',
                                id='sel-colorPicker-spotColor',
                            ),
                            title=fac.AntdText('select spot color'), color='white'
                        ),
                    ],
                    size='middle',
                    style={'width':'100%', 'marginLeft':'-4.5px'}
                ),
                fac.AntdCenter(
                    html.Div(
                        [
                            fac.AntdCheckbox(id='sel-checkbox-image', label='Image', checked=True),
                            fac.AntdCheckbox(id='sel-checkbox-spot', label='Spot', checked=True),
                        ],
                        style={
                            'display': 'flex',
                            'flex-direction': 'row',
                            'gap': '45px'
                        }
                    ),
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
            new_project_modal,
            clustring_status_modal,
            dcc.Store(id='sel-store-clusters', data=1),
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
                    ),
                    html.Div(id='refresh-sel-status', style={'display':'none'})
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