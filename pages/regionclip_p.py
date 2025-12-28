from dash import html, dcc
import feffery_antd_components as fac
from callbacks.regionclip_c import *

new_task_modal = fac.AntdModal(
    fac.AntdSpace(
        [
            fac.AntdDivider(),
            fac.AntdSpace(
                [
                    fac.AntdInput(
                        id='clip-input-taskname',
                        placeholder='Input Task Name...', 
                        autoComplete="off",
                        style={'width':'492px', 'display':'block'},
                    ),
                    fac.AntdButton(
                        'Import Task List', 
                        type='primary',
                        id='clip-button-importTaskList',
                        icon=fac.AntdIcon(icon='antd-cloud-download'),
                        style={'backgroundColor':'#5383c3'}
                    ),
                ],
                style={'width':'100%'},
            ),
            fac.AntdDraggerUpload(
                id='clip-dragger-upload',
                apiUrl='/upload/regionClip',
                fileMaxSize=100,
                fileListMaxLength=1,
                text='Upload Task List',
                locale='en-us',
                showUploadList=False,
                hint='Click or drag file here to upload',
                showErrorMessage=False,
            ),
            fac.AntdText('No file', id='clip-tasklist-filename', type='secondary'),
            fac.AntdDivider(),
            fac.AntdButton(
                'Submit',
                type='primary',
                id='clip-button-submitTaskList',
                icon=fac.AntdIcon(icon='md-launch'),
                style={'backgroundColor':'#d0826c', 'width':150, 'float':'right'}
            ),
        ],
        size='middle',
        direction='vertical',
        style={'width':'100%'},
    ),
    id='clip-modal-newtask', 
    title='Create New Task',
    mask=False,
    width=700,
    maskClosable=False,
    visible=False
)

clip_modal = fac.AntdModal(
    [
        fac.AntdDivider(),
        fac.AntdSpace(
            [
                fac.AntdInput(
                    id='clip-input-clipName',
                    placeholder='Input Clip Name...', 
                    autoComplete="off",
                    style={'width':'200px', 'display':'block'},
                ),
                fac.AntdButton(
                    'Add',
                    type='primary',
                    id='clip-button-submitClipName',
                    icon=fac.AntdIcon(icon='antd-check'),
                    style={'backgroundColor':'#d0826c', 'width':100, 'float':'right'}
                ),
            ],
            style={'width':'100%'},
        )
    ],
    id='clip-modal-inputClipName', 
    title='Set Clip Name',
    centered=True,
    mask=False,
    width=350,
    maskClosable=False,
    visible=False
)

control_panel = html.Div(
    fac.AntdCard(
        fac.AntdSpace(
            [
                new_task_modal,
                clip_modal,
                fac.AntdButton(
                    'Create New Task', 
                    type='primary',
                    id='clip-button-newtask',
                    icon=fac.AntdIcon(icon='antd-plus'),
                    style={'backgroundColor':'#698aab', 'width': '100%'}
                ),
                fac.AntdDivider(),
                fac.AntdText('Select Project'),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        id='clip-select-taskname',
                        placeholder='Select Project',
                        debounceWait=300,
                        locale='en-us',
                        allowClear=False,
                        style={'width':'100%'}
                    ),
                    id='clip-select-taskname-tooltip',
                    open=False,
                    title=fac.AntdText('Select Project'), 
                    color='white'
                ),
                dcc.Store(id='clip-store-taskname', storage_type='local'),
                html.Div(id="refresh-clip-status", style={'display':'none'}),

                fac.AntdText('Select Slice', id='clip-text-slice'),
                fac.AntdSelect(
                    id='clip-select-slice',
                    placeholder='Select Slice',
                    debounceWait=300,
                    locale='en-us',
                    allowClear=False,
                    style={'width':'100%'}
                ),
                fac.AntdSpace(
                    [
                        fac.AntdText('Select Clip', id='clip-text-clip'),
                        fac.AntdButton(
                            id='clip-add-clipName',
                            size='small',
                            shape='circle',
                            icon=fac.AntdIcon(icon='pi-plus'),
                        ),
                        fac.AntdPopconfirm(
                            fac.AntdButton(
                                id='clip-delete-clipName',
                                size='small',
                                shape='circle',
                                icon=fac.AntdIcon(icon='pi-minus')
                            ),
                            id='clip-delete-clipName-confirm',
                            locale='en-us',
                            arrow='hide',
                            okText='yes',
                            placement='bottomLeft',
                            title='Confirm Delete?'
                        )
                    ],
                    size='small'
                ),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        id='clip-select-clipname',
                        placeholder='Select Clip',
                        debounceWait=300,
                        locale='en-us',
                        allowClear=False,
                        style={'width':'100%'}
                    ),
                    id='clip-select-clipname-tooltip',
                    open=False,
                    title=fac.AntdText('Select Clip'), 
                    color='white'
                ),
                fac.AntdDivider(),
                fac.AntdText('Export Project'),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        id='clip-select-taskname-export',
                        placeholder='Export Project',
                        debounceWait=300,
                        locale='en-us',
                        allowClear=False,
                        style={'width':'100%'}
                    ),
                    id='clip-select-taskname-tooltip-export',
                    open=False,
                    title=fac.AntdText('Export Project'), 
                    color='white'
                ),
                fac.AntdText('Export Clip'),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        id='clip-select-clipname-export',
                        placeholder='Export Clip',
                        debounceWait=300,
                        locale='en-us',
                        allowClear=False,
                        style={'width':'100%'}
                    ),
                    id='clip-select-clipname-tooltip-export',
                    open=False,
                    title=fac.AntdText('Export Clip'), 
                    color='white'
                ),
                dcc.Download(id="clip-download-exportData"),
                fac.AntdButton(
                    'Export Data', 
                    type='primary',
                    id='clip-button-exportData',
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
    id='clip-control_panel',
)

content_panel = html.Div(
    fac.AntdSpace(
        [
            fac.AntdSpace(
                [
                    fac.AntdButton(
                        'Start Clip', 
                        type='primary',
                        id='clip-button-startClip',
                        icon=fac.AntdIcon(icon='antd-scissor'),
                        style={'backgroundColor':'#7d8a70'}
                    ),
                    fac.AntdPopconfirm(
                        fac.AntdButton(
                            'Delete Task', 
                            id='clip-delete-task',
                            type='primary',
                            icon=fac.AntdIcon(icon='antd-delete'),
                            style={'backgroundColor':'#ca8269'}
                        ),
                        id='clip-delete-task-confirm',
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
                                id='clip-button-showBug',
                                type='primary', 
                                icon=fac.AntdIcon(icon='antd-bug'), 
                                style={'backgroundColor':'#a87a76', 'display':'block'}
                            ),
                            fac.AntdText("🚨 Error occurred during clip. Click Show Bug for details !", style={'color':'#a87a76'}, strong=True)
                        ],
                        id='clip-bug-panel',
                        size='middle',
                        style={'display':'none'}
                    )
                ],
                size='middle',
                style={'marginTop':'24.5px'},
            ),
            fac.AntdCenter(
                fac.AntdSplitter(
                    items=[
                        {
                            'children': fac.AntdCenter(
                                dcc.Graph(
                                    id="clip-graph-original", 
                                    config={'displaylogo':False}, 
                                    style={'height':'calc(90vh - 50px)', 'width':'100%'}
                                ),
                            ),
                            'collapsible': True,
                            'defaultSize': '70%',
                        },
                        {
                            'children': fac.AntdSpace(
                                [
                                    fac.AntdRibbon(
                                        fac.AntdImage(
                                            src='',
                                            height='calc(38vh - 17px)',
                                            id='clip-stain-image',
                                            style={'visibility':'hidden'}
                                        ),
                                        color='#698aab',
                                        placement='start',
                                        text='Clipped Image',
                                        style={'width':'105px'}
                                    ),
                                    fac.AntdRibbon(
                                        fac.AntdImage(
                                            src='',
                                            height='calc(38vh - 17px)',
                                            id='clip-gem-image',
                                            style={'visibility':'hidden'}
                                        ),
                                        color='#ca8269',
                                        placement='start',
                                        text='Clipped GEM',
                                        style={'width':'105px'}
                                    ),
                                ],
                                direction='vertical',
                                size='middle',
                                style={'height':'100%', 'width':'100%','display': 'flex', 'alignItems': 'center', 'justifyContent':'center'}              
                            ),
                            'collapsible': True,
                            'defaultSize': '30%',
                        }    
                    ],
                ),
            )
        ],
        direction='vertical',
        size='middle',
        style={'width':'calc(100vw - 480px)', 'height':'100%'}
    ),
    id='clip-content_panel',
    style={'width':'100%', 'height':'100%', 'display':'flex', 'justifyContent':'center'}
)