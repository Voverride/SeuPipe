from dash import html
from dash import dcc
import feffery_antd_components as fac
from callbacks.expansion_c import *



new_task_modal = fac.AntdModal(
    fac.AntdSpace(
        [
            fac.AntdDivider('Data Preparation'),
            fac.AntdSpace(
                [
                    fac.AntdTooltip(
                        fac.AntdSelect(
                            id='expandtask-select-taskname',
                            placeholder='Select Project',
                            debounceWait=300,
                            locale='en-us',
                            allowClear=False,
                            style={'width':512}
                        ),
                        id='expandtask-select-taskname-tooltip',
                        title=fac.AntdText('Select Segmentation Project'),
                        color='white'
                    ),
                    fac.AntdButton(
                        'Submit',
                        type='primary',
                        id='expand-button-submitTask',
                        icon=fac.AntdIcon(icon='md-launch'),
                        style={'backgroundColor':'#d0826c', 'width':123}
                    ),
                ],
                style={'width':'100%'},
                size='middle'
            ),
            fac.AntdDivider('Spatial Parameters'),
            fac.AntdSpace(
                [
                    fac.AntdTooltip(
                        fac.AntdSelect(
                            id='expand-select-mode',
                            options=[
                                {'label':'patch by patch', 'value':'patch'},
                                {'label':'whole image', 'value':'whole'},
                            ],
                            allowClear=False,
                            value='patch',
                            style={'width':206}
                        ),
                        title=fac.AntdText(
                            fac.AntdSpace(
                                [
                                    fac.AntdText('Processing Mode:', strong=True),
                                    fac.AntdText('🔸patch by patch: Split image for memory efficiency (recommended for large images)'),
                                    fac.AntdText('🔸whole image: Process at once (requires high memory, for small images only)')
                                ],
                                direction='vertical',
                            )
                        ), 
                        color='white'
                    ),
                    fac.AntdTooltip(
                        fac.AntdInputNumber(
                            id='expand-input-patchsize',
                            value=300,
                            precision=0,
                            min=200,
                            addonBefore='patchSize',
                            style={'width':206}
                        ),
                        title=fac.AntdText('patchSize: width × height for splitting large tissue into smaller blocks. (eg. patch_size=1200 means processing the 1200×1200 spots as one block.)'), 
                        color='white'
                    ),
                    fac.AntdTooltip(
                        fac.AntdInputNumber(
                            id='expand-input-binsize',
                            value=3,
                            precision=0,
                            min=1,
                            addonBefore='binSize',
                            style={'width':206}
                        ),
                        title=fac.AntdText('binSize: number of adjacent spots (e.g., 3×3) merged into one unit to reduce noise and speed up computation'), 
                        color='white'
                    ),
                ],
                size='middle',
                style={'width':'100%'},
            ),
            
            fac.AntdDivider('Model Configuration'),
            fac.AntdSpace(
                [
                    fac.AntdTooltip(
                        fac.AntdInputNumber(
                            id='expand-input-epochs',
                            value=100,
                            precision=0,
                            min=1,
                            addonBefore='epochs',
                            style={'width':206}
                        ),
                        title=fac.AntdText('epochs: Number of times the model processes the entire dataset during training (e.g., 100 passes). More epochs improve learning but may cause overfitting if excessive'), 
                        color='white'
                    ),

                    fac.AntdTooltip(
                        fac.AntdInputNumber(
                            id='expand-input-diameter',
                            value=15,
                            precision=0,
                            min=1,
                            addonBefore='diameter',
                            style={'width':206}
                        ),
                        title=fac.AntdText('diameter: Estimated cell diameter (in spot units), used to guide the algorithm in identifying individual cell boundaries(Example: If cells span ~15 spots on average, set r_estimate=15 to help the model distinguish neighboring cells.)'), 
                        color='white'
                    ),

                    fac.AntdTooltip(
                        fac.AntdInputNumber(
                            id='expand-input-neighbor',
                            value=50,
                            precision=0,
                            min=1,
                            addonBefore='neighbors',
                            style={'width':206}
                        ),
                        title=fac.AntdText('neighbors: Number of nearby spots the model considers when predicting cell boundaries for each spot (default: 50). More neighbors capture broader context but increase compute time'), 
                        color='white'
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
    id='expand-modal-newtask', 
    title='Create New Project',
    mask=False,
    width=700,
    maskClosable=False,
    visible=False
)


control_panel = html.Div(
    fac.AntdCard(
        fac.AntdSpace(
            [
                dcc.Interval(id="init-restore-expansion", interval=1, max_intervals=1),
                new_task_modal,
                fac.AntdButton(
                    'Create New Project', 
                    type='primary',
                    id='expansion-button-newtask',
                    icon=fac.AntdIcon(icon='antd-plus'),
                    style={'backgroundColor':'#698aab', 'width': '100%'}
                ),
                fac.AntdDivider(),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        id='expand-select-taskname',
                        placeholder='Select Project',
                        debounceWait=300,
                        locale='en-us',
                        allowClear=False,
                        style={'width':'100%'}
                    ),
                    id='expand-select-taskname-tooltip',
                    open=False,
                    title=fac.AntdText('Select Project'), 
                    color='white'
                ),
                dcc.Store(id='expand-store-taskname', storage_type='local'),
            ],
            size='middle',
            direction='vertical',
            style={'width':'100%'},
        ),
        styles={'header': {'display': 'none'}},
        style={'height':'93vh', 'maxHeight': '93vh','overflowY': 'auto'},
    ), 
    id='expand-control_panel',
)

content_panel = html.Div(
    fac.AntdSpace(
        [
            fac.AntdSpace(
                [
                    fac.AntdButton(
                        'Start Project', 
                        type='primary',
                        id='exp-start-task',
                        icon=fac.AntdIcon(icon='antd-carry-out'),
                        style={'backgroundColor':'#7d8a70'}
                    ),
                    fac.AntdPopconfirm(
                        fac.AntdButton(
                            'Delete Project', 
                            id='exp-delete-task',
                            type='primary',
                            icon=fac.AntdIcon(icon='antd-delete'),
                            style={'backgroundColor':'#ca8269'}
                        ),
                        id='exp-delete-task-confirm',
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
                                id='exp-button-showBug',
                                type='primary', 
                                icon=fac.AntdIcon(icon='antd-bug'), 
                                style={'backgroundColor':'#a87a76', 'display':'block'}
                            ),
                            fac.AntdText("🚨 Error occurred during execution. Click Show Bug for details !", style={'color':'#a87a76'}, strong=True)
                        ],
                        id='exp-bug-panel',
                        size='middle',
                        style={'display':'none'}
                    )
                ],
                size='middle',
                style={'marginTop':'24.5px'},
            ),
            fac.AntdTable(
                columns=[
                    {'title': 'creator', 'dataIndex': 'creator', 'width':'15%', 'renderOptions': {'renderType': 'ellipsis'}},
                    {'title': 'mode', 'dataIndex': 'mode', 'width':'15%', 'renderOptions': {'renderType': 'ellipsis'}},
                    {'title': 'patchSize', 'dataIndex': 'patchSize', 'renderOptions': {'renderType': 'ellipsis'}},
                    {'title': 'binSize', 'dataIndex': 'binSize', 'renderOptions': {'renderType': 'ellipsis'}},
                    {'title': 'epochs', 'dataIndex': 'epochs', 'renderOptions': {'renderType': 'ellipsis'}},
                    {'title': 'diameter', 'dataIndex': 'diameter', 'renderOptions': {'renderType': 'ellipsis'}},
                    {'title': 'neighbors', 'dataIndex': 'neighbors', 'renderOptions': {'renderType': 'ellipsis'}},
                    {'title': 'progress', 'dataIndex': 'progress', 'width':'21.2%', 'renderOptions': {'renderType': 'mini-ring-progress'}}
                ],
                # data = [
                #     {
                #         'creator': 'zhouyb',
                #         'mode': 'patch',
                #         'patchSize': 300,
                #         'binSize': 3,
                #         'epochs': 100,
                #         'diameter': 16,
                #         'neighbors': 50,
                #         'progress': 5,
                #     }
                # ],
                style={'width': '100%'},
                id='exp-table-metadata',
                locale='en-us',
                bordered=True,
                miniChartHeight=75,
                pagination=False,
            ),
            fac.AntdTable(
                columns=[
                    {'title': 'z', 'dataIndex': 'z', 'width':'10%', 'renderOptions': {'renderType': 'ellipsis'}},
                    {'title': 'segmentation','dataIndex': 'segmentation', 'group': 'Segmentation', 'width':'15%', 'renderOptions': {'renderType': 'status-badge'}},
                    {'title': 'postprocess','dataIndex': 'seg_postprocess', 'group': 'Segmentation', 'width':'15%', 'renderOptions': {'renderType': 'status-badge'}},
                    {'title': 'preprocess','dataIndex': 'preprocess', 'group': 'Expansion', 'width':'15%', 'renderOptions': {'renderType': 'status-badge'}},
                    {'title': 'train','dataIndex': 'train', 'group': 'Expansion', 'width':'15%', 'renderOptions': {'renderType': 'status-badge'}},
                    {'title': 'postprocess','dataIndex': 'postprocess', 'group': 'Expansion', 'width':'15%', 'renderOptions': {'renderType': 'status-badge'}},
                    {'title': 'patchprocess','dataIndex': 'patchprocess', 'group': 'Expansion', 'width':'15%', 'renderOptions': {'renderType': 'status-badge'}},
                ],
                conditionalStyleFuncs={
                    label: """
                        (record, index) => {
                            return { style: { backgroundColor: "#FAFAFA" } }
                        }
                    """ 
                    for label in ['z', 'segmentation', 'seg_postprocess']
                },
                style={'width': '100%'},
                id='exp-table-tasklist',
                bordered=True,
                locale='en-us',
                pagination=False,
                maxHeight='calc(93vh - 335px)'
            )
        ],
        direction='vertical',
        size='middle',
        style={'width':'calc(100vw - 480px)', 'height':'100%'}
    ),
    style={'width':'100%', 'height':'100%', 'display':'flex', 'justifyContent':'center'}
)