from dash import html, dcc
import feffery_antd_components as fac

row_left = html.Div(
    fac.AntdSpace(
        [
            fac.AntdDivider('Data Preparation'),
            html.Div(
                [
                    fac.AntdSpace(
                        [
                            fac.AntdButton(
                                'Reference Data', 
                                id='annotask-button-refdata',
                                type='primary', 
                                icon=fac.AntdIcon(icon='antd-cloud-download'), 
                                style={'backgroundColor':'#5383c3','width':'90%'}
                            ),
                            fac.AntdTooltip(
                                fac.AntdButton(
                                    'Query Data', 
                                    id='annotask-button-querydata',
                                    type='primary', 
                                    icon=fac.AntdIcon(icon='antd-cloud-download'), 
                                    style={'backgroundColor':'#d0826c', 'width':'90%'}
                                ), 
                                title=fac.AntdText('⚠️ Uploading a new file will overwrite the processed data. Remember to save your data before proceeding'), color='white'
                            )
                        ],
                        size='middle',
                        direction='vertical',
                        style={'width':'29%', 'folat':'left'}
                    ),
                    fac.AntdSpace(
                        [
                            fac.AntdTooltip(
                                fac.AntdSelect(
                                    placeholder='label field',
                                    optionFilterMode='case-insensitive',
                                    allowClear=False,
                                    options=[],
                                    id='annotask-select-label',
                                    style={'width':'90%'},
                                    locale='en-us',
                                ),
                                title=fac.AntdText('select label field'), color='white'
                            ),
                            fac.AntdTooltip(
                                fac.AntdSelect(
                                    placeholder='x field',
                                    optionFilterMode='case-insensitive',
                                    allowClear=False,
                                    options=[],
                                    id='annotask-select-x',
                                    style={'width':'90%'},
                                    locale='en-us',
                                ),
                                title=fac.AntdText('select x field'), color='white'
                            ),   
                        ],
                        size='middle',
                        direction='vertical',
                        style={'width':'21%', 'folat':'left'}
                    ),
                    fac.AntdSpace(
                        [
                            fac.AntdSelect(
                                style={'width':'90%', 'visibility':'hidden'},
                            ),
                            fac.AntdTooltip(
                                fac.AntdSelect(
                                    placeholder='y field',
                                    optionFilterMode='case-insensitive',
                                    allowClear=False,
                                    options=[],
                                    style={'width':'90%'},
                                    id='annotask-select-y',
                                    locale='en-us',
                                ),
                                title=fac.AntdText('select y field'), color='white'
                            ),
                        ],
                        size='middle',
                        direction='vertical',
                        style={'width':'20%', 'folat':'left'}
                    ),
                    fac.AntdSpace(
                        [
                            fac.AntdSelect(
                                style={'width':'90%', 'visibility':'hidden'},
                            ),
                            fac.AntdTooltip(
                                fac.AntdSelect(
                                    placeholder='z field',
                                    optionFilterMode='case-insensitive',
                                    allowClear=False,
                                    options=[],
                                    style={'width':'90%'},
                                    id='annotask-select-z',
                                    locale='en-us',
                                ),
                                title=fac.AntdText('select z field'), color='white'
                            ),
                        ],
                        size='middle',
                        direction='vertical',
                        style={'width':'20%', 'folat':'left'}
                    ),
                ],
                style={'width':'100%'}
            ),
            fac.AntdDivider('Preprocessing'),
            html.Div(
                [
                    fac.AntdSpace(
                        [
                            fac.AntdCheckbox(label='Remove Mitochondrial Genes (MT)', id='annotask-remove-mt', checked=True, style={'marginBottom':'0.6rem'}),
                            fac.AntdCheckbox(label='Remove Ribosomal Genes (Ribo)', id='annotask-remove-ribo', checked=True),
                        ],
                        direction='vertical',
                        style={'width':'50%', 'display':'inline-block'}
                    ),
                    fac.AntdSpace(
                        [
                            fac.AntdCheckbox(label='Remove Hemoglobin Genes (Hb)', id='annotask-remove-hb', checked=True, style={'marginBottom':'0.6rem'}),
                            fac.AntdCheckbox(label='Use Highly Variable Genes (HVG)', id='annotask-use-hvg', checked=True),
                        ],
                        direction='vertical',
                        style={'width':'50%', 'display':'inline-block'}
                    ),
                ],
                style={'width':'100%'}
            ),
            fac.AntdDivider('Model Configuration'),
            html.Div(
                [
                    fac.AntdSpace(
                        [
                            fac.AntdTooltip(
                                fac.AntdInputNumber(
                                        min=2,
                                        value=2,
                                        placeholder='n_layers',
                                        id='annotask-nlayers',
                                        precision=0,
                                        debounceWait=500,
                                        addonBefore='n_layers',
                                        style={'width':'90%'}
                                ),
                                title=fac.AntdText('Number of hidden layers. Deeper networks (2-3) model hierarchical patterns; shallow networks reduce compute cost.'), color='white'
                            ),
                            fac.AntdTooltip(
                                fac.AntdInputNumber(
                                        min=1,
                                        placeholder='epochs',
                                        id='annotask-epochs',
                                        value=100,
                                        precision=0,
                                        debounceWait=500,
                                        addonBefore='epochs',
                                        style={'width':'90%'}
                                ),
                                title=fac.AntdText('Maximum training cycles. Stop early if validation loss plateaus. Typical: 100-1000 for scRNA-seq + spatial integration'), color='white'
                            ),
                        ],
                        direction='vertical',
                        size='middle',
                        style={'width':'30%', 'float':'left'}
                    ),
                    fac.AntdSpace(
                        [
                            fac.AntdTooltip(
                                fac.AntdInputNumber(
                                        min=1,
                                        value=256,
                                        placeholder='n_hiddens',
                                        id='annotask-nhiddens',
                                        debounceWait=500,
                                        step=1,
                                        precision=0,
                                        addonBefore='n_hiddens',
                                        style={'width':'90%'}
                                ),
                                title=fac.AntdText('Number of neurons per hidden layer. Higher values capture complex gene-cell relationships but may overfit. Typical: 128-512.'), color='white'
                            ),
                            fac.AntdTooltip(
                                fac.AntdInputNumber(
                                        min=1,
                                        placeholder='batch_size',
                                        id='annotask-batchsize',
                                        value=128,
                                        precision=0,
                                        debounceWait=500,
                                        addonBefore='batch_size',
                                        style={'width':'90%'}
                                ),
                                title=fac.AntdText('Cells per training batch. Smaller sizes (64-256) help generalization; larger sizes speed up training if memory allows.'), color='white'
                            )
                        ],
                        direction='vertical',
                        size='middle',
                        style={'width':'30%', 'float':'left'}
                    ),
                    fac.AntdSpace(
                        [
                            fac.AntdTooltip(
                                fac.AntdInputNumber(
                                        min=1,
                                        value=20,
                                        placeholder='n_latent',
                                        id='annotask-nlatent',
                                        precision=0,
                                        debounceWait=500,
                                        addonBefore='n_latent',
                                        style={'width':'90%'}
                                ),
                                title=fac.AntdText('Dimension of latent space. Lower values (10-30) compress data better for integration; higher values preserve subtle differences.'), color='white'
                            ),
                            fac.AntdTooltip(
                                fac.AntdInputNumber(
                                        min=0.1,
                                        value=0.2,
                                        placeholder='dropout',
                                        id='annotask-dropout',
                                        debounceWait=500,
                                        step=0.1,
                                        addonBefore='dropout',
                                        style={'width':'90%'}
                                ),
                                title=fac.AntdText('Fraction of neurons randomly disabled during training (0.1-0.3). Prevents overfitting to reference data.'), color='white'
                            ),
                            fac.AntdButton(
                                'Start Training', 
                                id='anntask-button-train',
                                type='primary', 
                                icon=fac.AntdIcon(icon='antd-carry-out'), 
                                style={'backgroundColor':'#5383c3','width':'90%'}
                            ),
                        ],
                        direction='vertical',
                        size='middle',
                        style={'width':'30%', 'float':'left'}
                    ),
                ],
                style={'width':'100%'} 
            ),
            fac.AntdDivider('Training Progress'),
            fac.AntdTimeline(
                items=[
                    {
                        'content': html.Div(
                            [
                                html.Div(
                                    [
                                        fac.AntdText('Creator:', strong=True, style={'marginRight':'6px'}),
                                        fac.AntdText(type='success', id='annotask-timeline-creator'),  
                                    ],
                                    style={'marginRight':'20px'}
                                ),
                                fac.AntdButton(
                                    'Show Bug',
                                    id='annotask-button-showBug',
                                    size='small', 
                                    type='primary', 
                                    icon=fac.AntdIcon(icon='antd-bug'), 
                                    style={'backgroundColor':'#bb5548', 'display':'none'}
                                )
                            ],
                            style={'width':'100%', 'display':'flex'}
                        ),
                        'icon':fac.AntdAvatar(size='small'),
                    },
                    {
                        'content': 'Preprocessing',
                        'color':'gray',
                        'icon': fac.AntdIcon(icon='md-schedule', id='annot-step1')
                    },
                    {
                        'content': fac.AntdSpace(
                            [
                                html.Div(
                                    [
                                        fac.AntdText('Training', style={'marginRight':'10px'}),  
                                    ],
                                ),
                                fac.AntdProgress(percent=0, id='annot-percent', style={'display':'block', 'width': '90%'}),
                            ],
                            style={'width':'100%', 'display':'block'},
                            size='middle'
                        ),
                        'color':'gray',
                        'icon': fac.AntdIcon(icon='md-schedule', id='annot-step2')
                    },
                    {
                        'content': 'Prediction',
                        'color':'gray',
                        'icon': fac.AntdIcon(icon='md-schedule', id='annot-step3')
                    }
                ],
                id = 'annotation-timeline-status',
                style={'width': '100%', 'clear':'both'}
            )
        ],
        direction='vertical',
        size='middle',
    ),
    style={'float':'left'}
)

row_right = dcc.Graph(
    config={'displaylogo':False}, 
    id='annotation-graph-queryumap',
    style={'display':'none'},
    responsive=True
)

AnnotationTask = fac.AntdDrawer(
    fac.AntdSplitter(
        items=[
            {
                'children': html.Div(
                    row_left,
                    style={'height': '100%', 'paddingLeft':'10px', 'display': 'grid', 'placeItems': 'center' }
                ),
                'collapsible': True,
            },
            {
                'children': fac.AntdCenter(
                    row_right,
                    style={'height': '100%', 'width': '100%'}
                ),
                'collapsible': True,
            },
        ],
        style={'height': '100%', 'width': '100%'},
    ),
    title='Create New Task', 
    id='annotask-drawer',
    visible=False, 
    mask=False,
    maskClosable=False,
    height='100%',
    placement='bottom',
    containerId='annotation-content-panel',
)