from dash import html, dcc
import feffery_antd_components as fac
from callbacks.segmentation_c import *


new_task_modal = fac.AntdModal(
    fac.AntdSpace(
        [
            fac.AntdDivider('Data Preparation'),
            fac.AntdSpace(
                [
                    fac.AntdInput(
                        id='seg-input-taskname',
                        placeholder='Input Task Name...', 
                        style={'width':'492px', 'display':'block'},
                    ),
                    fac.AntdButton(
                        'Import Task List', 
                        type='primary',
                        id='seg-button-importTaskList',
                        icon=fac.AntdIcon(icon='antd-cloud-download'),
                        style={'backgroundColor':'#5383c3'}
                    ),
                ],
                style={'width':'100%'},
            ),
            fac.AntdDraggerUpload(
                id='seg-dragger-upload',
                apiUrl='/upload/',
                fileMaxSize=100,
                fileListMaxLength=1,
                text='Upload Task List',
                locale='en-us',
                showUploadList=False,
                hint='Click or drag file here to upload',
                showErrorMessage=False,
            ),
            fac.AntdText('No file', id='seg-tasklist-filename', type='secondary'),
            fac.AntdDivider('Model Configuration'),
            fac.AntdSpace(
                [
                    fac.AntdTooltip(
                        fac.AntdSelect(
                            id='seg-select-modelType',
                            options=['nuclei', 'cyto', 'cyto2', 'cyto3'],
                            value='cyto3',
                            prefix=fac.AntdIcon(icon='pi-cube', style={'fontSize': 10}),
                            allowClear=False,
                            style={'width':120}
                        ),
                        title=fac.AntdText(
                            fac.AntdSpace(
                                [
                                    fac.AntdText('Select Model Type:', strong=True),
                                    fac.AntdText('🔸nuclei: Nuclear model'),
                                    fac.AntdText('🔸cyto: Cytoplasm model'),
                                    fac.AntdText('🔸cyto2: Enhanced cytoplasm model'),
                                    fac.AntdText('🔸cyto3: Generalist super-model'),
                                ],
                                direction='vertical',
                            )
                        ), 
                        color='white'
                    ),
                    fac.AntdTooltip(
                        fac.AntdInputNumber(
                            id='seg-input-diameter',
                            value=0,
                            precision=0,
                            min=0,
                            addonAfter='px',
                            style={'width':120}
                        ),
                        title=fac.AntdText(
                            fac.AntdSpace(
                                [
                                    fac.AntdText('Input Estimated Diameter:', strong=True),
                                    fac.AntdText('🔸0: automatically estimated'),
                                    fac.AntdText('🔸others: cell diameter estimation')
                                ],
                                direction='vertical',
                            )
                        ), 
                        color='white'
                    ),
                    fac.AntdTooltip(
                        fac.AntdInputNumber(
                            id='seg-input-batchsize',
                            value=8,
                            precision=0,
                            min=1,
                            style={'width':120}
                        ),
                        title=fac.AntdText('Batch Size: Number of image patches processed simultaneously on the GPU'), 
                        color='white'
                    ),
                    fac.AntdCheckbox(label='UseGPU', checked=True, id='seg-checkbox-useGPU'),
                    fac.AntdButton(
                        'Submit',
                        type='primary',
                        id='seg-button-submitTaskList',
                        icon=fac.AntdIcon(icon='md-launch'),
                        style={'backgroundColor':'#d0826c', 'width':150, 'marginLeft':'-6px'}
                    ),
                ],
                size='middle',
                style={'width':'100%'},
            )
        ],
        size='middle',
        direction='vertical',
        style={'width':'100%'},
    ),
    id='seg-modal-newtask', 
    title='Create New Task',
    mask=False,
    width=700,
    maskClosable=False,
    visible=False
)

control_panel = html.Div(
    fac.AntdCard(
        fac.AntdSpace(
            [
                dcc.Interval(id="init-restore-segmentation", interval=1, max_intervals=1),
                dcc.Interval(interval=1000, disabled=True, id='segmentation-interval'),
                dcc.Interval(id="segmentation-event-loop", interval=1000),
                new_task_modal,
                fac.AntdButton(
                    'Create New Task', 
                    type='primary',
                    id='segmentation-button-newtask',
                    icon=fac.AntdIcon(icon='antd-plus'),
                    style={'backgroundColor':'#698aab', 'width': '100%'}
                ),
                fac.AntdDivider(),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        id='seg-select-taskname',
                        # options=['task1', 'task2', 'task3'],
                        # value='task1',
                        placeholder='Select Task',
                        debounceWait=300,
                        # persistence=True,
                        locale='en-us',
                        allowClear=False,
                        style={'width':'100%'}
                    ),
                    id='seg-select-taskname-tooltip',
                    open=False,
                    title=fac.AntdText('Select Task'), 
                    color='white'
                ),
                dcc.Store(id='seg-store-taskname', storage_type='local'),
                # html.Div(
                #     [
                #         fac.AntdText('Creator:', strong=True, style={'marginRight':'6px'}),
                #         fac.AntdText('zhouyb', type='success', id='seg-creator'),  
                #     ],
                #     style={'width':'100%', 'marginTop':'5px'}
                # ),
                # fac.AntdProgress(percent=100, id='annot-percent', style={'display':'block', 'width': '100%'}),

            ],
            size='middle',
            direction='vertical',
            style={'width':'100%'},
        ),
        headStyle={'display': 'none'},
        style={'height':'93vh', 'maxHeight': '93vh','overflowY': 'auto'},
    ), 
    id='segmentation-control_panel',
)

content_panel = html.Div(
    fac.AntdSpace(
        [
            fac.AntdSpace(
                [
                    fac.AntdButton(
                        'Start Task', 
                        type='primary',
                        id='seg-start-task',
                        icon=fac.AntdIcon(icon='antd-carry-out'),
                        style={'backgroundColor':'#7d8a70'}
                    ),
                    fac.AntdPopconfirm(
                        fac.AntdButton(
                            'Delete Task', 
                            id='seg-delete-task',
                            type='primary',
                            icon=fac.AntdIcon(icon='antd-delete'),
                            style={'backgroundColor':'#ca8269'}
                        ),
                        id='seg-delete-task-confirm',
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
                                id='seg-button-showBug',
                                type='primary', 
                                icon=fac.AntdIcon(icon='antd-bug'), 
                                style={'backgroundColor':'#a87a76', 'display':'block'}
                            ),
                            fac.AntdText("🚨 Error occurred during execution. Click Show Bug for details !", style={'color':'#a87a76'}, strong=True)
                        ],
                        id='seg-bug-panel',
                        size='middle',
                        style={'display':'none'}
                    )
                ],
                size='middle',
                style={'marginTop':'24.5px'},
            ),
            fac.AntdTable(
                columns=[
                    {'title': 'creator', 'dataIndex': 'creator', 'width':'21.7%', 'renderOptions': {'renderType': 'ellipsis'}},
                    {'title': 'model', 'dataIndex': 'model', 'width':'21.6%', 'renderOptions': {'renderType': 'ellipsis'}},
                    {'title': 'diameter', 'dataIndex': 'diameter', 'renderOptions': {'renderType': 'ellipsis'}},
                    {'title': 'batchsize', 'dataIndex': 'batchsize', 'renderOptions': {'renderType': 'ellipsis'}},
                    {'title': 'GPU', 'dataIndex': 'GPU', 'renderOptions': {'renderType': 'ellipsis'}},
                    {'title': 'progress', 'dataIndex': 'progress', 'width':'21.2%', 'renderOptions': {'renderType': 'mini-ring-progress'}}
                ],
                # data=[
                #     {
                #         'creator': 'zhouyb',
                #         'model': 'cyto3',
                #         'diameter': '0',
                #         'batchsize': '8',
                #         'GPU': 'True',
                #         'progress': 0.5,
                #     }
                # ],
                style={'width': '100%'},
                id='seg-table-metadata',
                locale='en-us',
                bordered=True,
                miniChartHeight=75,
                pagination=False,
            ),
            fac.AntdTable(
                columns=[
                    {'title': 'z', 'dataIndex': 'z', 'width':'8%', 'renderOptions': {'renderType': 'ellipsis'}},
                    {'title': 'image', 'dataIndex': 'image', 'width':'36%', 'renderOptions': {'renderType': 'ellipsis'}},
                    {'title': 'gem', 'dataIndex': 'gem', 'width':'36%', 'renderOptions': {'renderType': 'ellipsis'}},
                    {'title': 'registration','dataIndex': 'registration', 'width':'10%', 'renderOptions': {'renderType': 'status-badge'}},
                    {'title': 'segmentation','dataIndex': 'segmentation', 'width':'10%', 'renderOptions': {'renderType': 'status-badge'}}
                ],
                # data=[
                #     {
                #         'z': 1,
                #         'image': 'image2',
                #         'gem': 'gem1',
                #         'registration': {'status': 'success', 'text': 'success'},
                #         'segmentation': {'status': 'success', 'text': 'success'},
                #     }
                # ]*3,
                style={'width': '100%'},
                id='seg-table-tasklist',
                bordered=True,
                locale='en-us',
                pagination=False,
                maxHeight='calc(93vh - 290px)'
            )
        ],
        direction='vertical',
        size='middle',
        style={'width':'calc(100vw - 480px)', 'height':'100%'}
    ),
    style={'width':'100%', 'height':'100%', 'display':'flex', 'justifyContent':'center'}
)