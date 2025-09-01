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
                            placeholder='Select Task',
                            debounceWait=300,
                            locale='en-us',
                            allowClear=False,
                            style={'width':512}
                        ),
                        id='expandtask-select-taskname-tooltip',
                        title=fac.AntdText('Select Segmentation Task'),
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
    title='Create New Task',
    mask=False,
    width=700,
    maskClosable=False,
    visible=True
)


control_panel = html.Div(
    fac.AntdCard(
        fac.AntdSpace(
            [
                dcc.Interval(id="init-restore-expansion", interval=1, max_intervals=1),
                dcc.Interval(interval=1000, disabled=True, id='expansion-interval'),
                dcc.Interval(id="expansion-event-loop", interval=1000),
                new_task_modal,
                fac.AntdButton(
                    'Create New Task', 
                    type='primary',
                    id='expansion-button-newtask',
                    icon=fac.AntdIcon(icon='antd-plus'),
                    style={'backgroundColor':'#698aab', 'width': '100%'}
                ),
                fac.AntdDivider(),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        id='expand-select-taskname',
                        # options=['task1', 'task2', 'task3'],
                        # value='task1',
                        placeholder='Select Task',
                        debounceWait=300,
                        # persistence=True,
                        locale='en-us',
                        allowClear=False,
                        style={'width':'100%'}
                    ),
                    id='expand-select-taskname-tooltip',
                    open=False,
                    title=fac.AntdText('Select Task'), 
                    color='white'
                ),
                dcc.Store(id='expand-store-taskname', storage_type='local'),
            ],
            size='middle',
            direction='vertical',
            style={'width':'100%'},
        ),
        headStyle={'display': 'none'},
        style={'height':'93vh', 'maxHeight': '93vh','overflowY': 'auto'},
    ), 
    id='expand-control_panel',
)

content_panel = html.Div(
    fac.AntdEmpty(
        description=fac.AntdText('当前页面开发中', type='secondary'),
        imageStyle={'height': 250},
    )
)