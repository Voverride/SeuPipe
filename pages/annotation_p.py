from dash import html, dcc
import feffery_antd_components as fac
from callbacks.annotation_c import *

new_project_modal = fac.AntdModal(
    fac.AntdSpace(
        [
            fac.AntdDivider(),
            fac.AntdSpace(
                [
                    fac.AntdButton(
                        'Reference Data', 
                        id='ann-button-refdata',
                        type='primary', 
                        icon=fac.AntdIcon(icon='antd-cloud-download'), 
                        style={'backgroundColor':'#5383c3','width':'150px'}
                    ),
                    dcc.Store(id='ann-store-refdata', data=None),
                    fac.AntdTooltip(
                        fac.AntdSelect(
                            placeholder='label field',
                            optionFilterMode='case-insensitive',
                            allowClear=False,
                            options=[],
                            id='ann-select-label',
                            style={'width':'115px'},
                            locale='en-us',
                        ),
                        title=fac.AntdText('select label field'), color='white'
                    ),
                ],
                size='middle'
            ),
            fac.AntdSpace(
                [
                    fac.AntdButton(
                        'Query Data', 
                        id='ann-button-querydata',
                        type='primary', 
                        icon=fac.AntdIcon(icon='antd-cloud-download'), 
                        style={'backgroundColor':'#d0826c', 'width':'150px'}
                    ),
                    dcc.Store(id='ann-store-querydata', data=None),
                    fac.AntdTooltip(
                        fac.AntdSelect(
                            placeholder='x field',
                            optionFilterMode='case-insensitive',
                            allowClear=False,
                            options=[],
                            style={'width':'115px'},
                            id='ann-select-x',
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
                            style={'width':'115px'},
                            id='ann-select-y',
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
                            style={'width':'115px'},
                            id='ann-select-z',
                            locale='en-us',
                        ),
                        title=fac.AntdText('select z field'), color='white'
                    )
                ],
                size='middle'
            ),
            fac.AntdDivider(),
            fac.AntdCenter(
                [
                    fac.AntdSpace(
                        [
                            fac.AntdCheckbox(label='Remove Mitochondrial Genes (MT)', id='ann-remove-mt', checked=True, style={'marginBottom':'0.6rem'}),
                            fac.AntdCheckbox(label='Remove Ribosomal Genes (Ribo)', id='ann-remove-ribo', checked=True),
                        ],
                        direction='vertical',
                        style={'width':'50%', 'display':'inline-block'}
                    ),
                    fac.AntdSpace(
                        [
                            fac.AntdCheckbox(label='Remove Hemoglobin Genes (Hb)', id='ann-remove-hb', checked=True, style={'marginBottom':'0.6rem'}),
                            fac.AntdCheckbox(label='Use Highly Variable Genes (HVG)', id='ann-use-hvg', checked=True),
                        ],
                        direction='vertical',
                        style={'width':'50%', 'display':'inline-block'}
                    ),
                ],
                style={'width':'100%'}
            ),
            fac.AntdDivider(),
            fac.AntdSpace(
                [
                    fac.AntdTooltip(
                        fac.AntdInputNumber(
                                min=2,
                                value=2,
                                placeholder='n_layers',
                                id='ann-nlayers',
                                precision=0,
                                debounceWait=500,
                                addonBefore='n_layers',
                                style={'width':'170px'}
                        ),
                        title=fac.AntdText('Number of hidden layers. Deeper networks (2-3) model hierarchical patterns; shallow networks reduce compute cost.'), color='white'
                    ),
                    fac.AntdTooltip(
                        fac.AntdInputNumber(
                                min=1,
                                value=256,
                                placeholder='n_hiddens',
                                id='ann-nhiddens',
                                debounceWait=500,
                                step=1,
                                precision=0,
                                addonBefore='n_hiddens',
                                style={'width':'170px'}
                        ),
                        title=fac.AntdText('Number of neurons per hidden layer. Higher values capture complex gene-cell relationships but may overfit. Typical: 128-512.'), color='white'
                    ),
                    fac.AntdTooltip(
                        fac.AntdInputNumber(
                                min=1,
                                value=20,
                                placeholder='n_latent',
                                id='ann-nlatent',
                                precision=0,
                                debounceWait=500,
                                addonBefore='n_latent',
                                style={'width':'170px'}
                        ),
                        title=fac.AntdText('Dimension of latent space. Lower values (10-30) compress data better for integration; higher values preserve subtle differences.'), color='white'
                    ),  
                ],
                size='middle'
            ),
            fac.AntdSpace(
                [
                    fac.AntdTooltip(
                        fac.AntdInputNumber(
                                min=1,
                                placeholder='epochs',
                                id='ann-epochs',
                                value=100,
                                precision=0,
                                debounceWait=500,
                                addonBefore='epochs',
                                style={'width':'170px'}
                        ),
                        title=fac.AntdText('Maximum training cycles. Stop early if validation loss plateaus. Typical: 100-1000 for scRNA-seq + spatial integration'), color='white'
                    ),
                    fac.AntdTooltip(
                        fac.AntdInputNumber(
                                min=1,
                                placeholder='batch_size',
                                id='ann-batchsize',
                                value=128,
                                precision=0,
                                debounceWait=500,
                                addonBefore='batch_size',
                                style={'width':'170px'}
                        ),
                        title=fac.AntdText('Cells per training batch. Smaller sizes (64-256) help generalization; larger sizes speed up training if memory allows.'), color='white'
                    ),
                    fac.AntdTooltip(
                        fac.AntdInputNumber(
                                min=0.1,
                                value=0.2,
                                placeholder='dropout',
                                id='ann-dropout',
                                debounceWait=500,
                                step=0.1,
                                addonBefore='dropout',
                                style={'width':'170px'}
                        ),
                        title=fac.AntdText('Fraction of neurons randomly disabled during training (0.1-0.3). Prevents overfitting to reference data.'), color='white'
                    ),
                ],
                size='middle',
            ),
            fac.AntdSpace(
                [
                    fac.AntdInput(
                        placeholder='input project name',
                        id='ann-projectname',
                        autoComplete="off",
                        status='error',
                        style={'width':'357px'}
                    ),
                    fac.AntdButton(
                        'Submit', 
                        id='ann-button-submit',
                        type='primary', 
                        icon=fac.AntdIcon(icon='md-launch'), 
                        style={'backgroundColor': '#d0826c', 'width': '170px'}
                    ),
                ],
                size='middle',
                style={'marginBottom':'15px'}
            ),
        ],
        direction='vertical',
        size='middle',
    ),
    title='Create New Project',
    id='ann-modal-newproject',
    maskClosable=False,
    mask=False,
    visible=False,
    width='590px',
)

control_panel = html.Div(
    fac.AntdCard(
        fac.AntdSpace(
            [
                fac.AntdButton(
                    'Create New Project', 
                    type='primary',
                    id='ann-button-newproject',
                    icon=fac.AntdIcon(icon='antd-plus'),
                    style={'backgroundColor':'#698aab', 'width': '100%'}
                ),
                fac.AntdDivider(),
                fac.AntdText('Select Project', id='ann-text-project'),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        id='ann-select-project',
                        placeholder='Select Project',
                        debounceWait=300,
                        locale='en-us',
                        allowClear=False,
                        style={'width':'100%'}
                    ),
                    id='ann-select-project-tooltip',
                    open=False,
                    title=fac.AntdText('Select Project'), 
                    color='white'
                ),
                fac.AntdDivider(),
                fac.AntdText('Spot Configuration'),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        options=[1,2,3,4,5,6,7,8,9,10],
                        value=5,
                        id='ann-select-spotSize',
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
                                id='ann-select-borderWidth',
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
                                value='#0d0015',
                                id='ann-colorPicker-boarderColor',
                            ),
                            title=fac.AntdText('Border Color'), color='white'
                        ),
                    ],
                    size='middle',
                    style={'width':'100%'}
                ),
                fac.AntdDivider(),
                fac.AntdText('Slice Configuration'),
                fac.AntdSlider(
                    range=True,
                    min=0,
                    max=0,
                    tooltipPrefix='',
                    id='ann-slider-slicer',
                    style={'width':'100%'}
                ),
                fac.AntdButton(
                     'Slicer',
                     type='primary', 
                     id='ann-button-slicer',
                     style={'backgroundColor':'#867ba9', 'width':'100%'},
                     icon=fac.AntdIcon(icon='md-content-cut')
                 ),
                fac.AntdDivider(),
                fac.AntdButton(
                    'Export Data', 
                    type='primary',
                    id='ann-button-exportData',
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
    id='ann-control_panel',
)


content_tabs = fac.AntdTabs(
    items=[
        {
            'key': 'ProjectInfo',
            'label': 'ProjectInfo',
            'children': fac.AntdSpace(
                [
                    fac.AntdSpace(
                        [
                            fac.AntdButton(
                                'Start Project', 
                                type='primary',
                                id='ann-start-project',
                                icon=fac.AntdIcon(icon='antd-carry-out'),
                                style={'backgroundColor':'#7d8a70'}
                            ),
                            fac.AntdPopconfirm(
                                fac.AntdButton(
                                    'Delete Project', 
                                    id='ann-delete-project',
                                    type='primary',
                                    icon=fac.AntdIcon(icon='antd-delete'),
                                    style={'backgroundColor':'#ca8269'}
                                ),
                                id='ann-delete-project-confirm',
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
                                        id='ann-button-showBug',
                                        type='primary', 
                                        icon=fac.AntdIcon(icon='antd-bug'), 
                                        style={'backgroundColor':'#a87a76', 'display':'block'}
                                    ),
                                    fac.AntdText("🚨 Error occurred during execution. Click Show Bug for details !", style={'color':'#a87a76'}, strong=True)
                                ],
                                id='ann-bug-panel',
                                size='middle',
                                style={'display':'none'}
                            ),
                            html.Div(children=None, id='refresh-ann-status', style={'display':'none'}),
                        ],
                        size='middle',
                        style={'width':'100%'},
                    ),
                    fac.AntdCard(
                        fac.AntdSteps(
                            id='ann-annotation-steps',
                            steps=[
                                {
                                    'title': 'Preprocess',
                                    'status': 'wait',
                                },
                                {
                                    'title': 'Training',
                                    'status': 'wait',
                                },
                                {
                                    'title': 'Postprocess',
                                    'status': 'wait',
                                }
                            ]
                        ),
                        title='Annotation Steps',
                    ),
                    fac.AntdTable(
                        columns=[
                            {'title': 'creator', 'dataIndex': 'creator', 'width':'10%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'date', 'dataIndex': 'date', 'width':'12%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'rmMt', 'dataIndex': 'rm_mt', 'width':'8%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'rmRibo', 'dataIndex': 'rm_ribo', 'width':'8%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'rmHb', 'dataIndex': 'rm_hb', 'width':'8%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'useHVG', 'dataIndex': 'use_hvg', 'width':'8%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'nLayers', 'dataIndex': 'n_layers', 'width':'8%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'nHiddens', 'dataIndex': 'n_hiddens', 'width':'8%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'nLatent', 'dataIndex': 'n_latent', 'width':'8%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'epochs', 'dataIndex': 'epochs', 'width':'8%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'batchSize', 'dataIndex': 'batch_size', 'width':'8%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'dropout', 'dataIndex': 'dropout', 'width':'8%', 'renderOptions': {'renderType': 'ellipsis'}},
                        ],
                        data=[{}],
                        style={'width': '100%'},
                        id='ann-table-metadata',
                        locale='en-us',
                        bordered=True,
                        pagination=False,
                    ),
                    fac.AntdTable(
                        columns=[
                            {'title': 'refdata', 'dataIndex': 'refdata', 'width':'30%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'querydata', 'dataIndex': 'querydata', 'width':'30%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'labelField', 'dataIndex': 'label', 'width':'10%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'xField', 'dataIndex': 'x', 'width':'10%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'yField', 'dataIndex': 'y', 'width':'10%', 'renderOptions': {'renderType': 'ellipsis'}},
                            {'title': 'zField', 'dataIndex': 'z', 'width':'10%', 'renderOptions': {'renderType': 'ellipsis'}},
                        ],
                        data=[{}],
                        style={'width': '100%', 'display':'inline-block'},
                        id='ann-table-refquery',
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
            'key': 'Results',
            'label': 'Results',
            'children': fac.AntdSplitter(
                items=[
                    {
                        'children': fac.AntdCenter(
                            dcc.Graph(
                                config={'displaylogo':False}, 
                                id='ann-graph-result', 
                                style={'height': '100%', 'width': '100%'},
                                responsive=True
                            ),
                            style={'height': '100%', 'width': '100%'}
                        ),
                        'defaultSize': '60%',
                        'collapsible': True,
                    },
                    {
                        'children': fac.AntdCenter(
                            dcc.Graph(
                                config={'displaylogo':False}, 
                                id='ann-graph-heatmap', 
                                style={'height': '100%', 'width': '100%'},
                                responsive=True
                            ),
                            style={'height': '100%', 'width': '100%'}
                        ),
                        'defaultSize': '40%',
                        'collapsible': True,
                    },
                ],
                style={'height': 'calc(93vh - 65px)', 'width':'100%', 'overflow': 'auto'}
            )
        }
    ],
    id='ann-content-tabs'
)

content_panel = html.Div(
    [   
        new_project_modal,
        content_tabs
    ],
    id='ann-content-panel',
    style={'width':'98%', 'height':'100%', 'position':'relative'}
)