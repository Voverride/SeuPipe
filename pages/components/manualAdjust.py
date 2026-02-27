from dash import html
import feffery_antd_components as fac
from dash import dcc
from utils.colors import redScaleColor, greenScaleColor, blueScaleColor

ManualAdjustDrawer = fac.AntdDrawer(
    [
        fac.AntdPopupCard(
            fac.AntdSpace(
                [
                    html.Div(
                        id='ma-div-legend',
                    ),
                    html.Div(
                        fac.AntdButton(
                            'load more', type='link', ghost=True, 
                            icon=fac.AntdIcon(icon='antd-down'),
                            iconPosition='end',
                            id='ma-button-loadMore'
                        ),
                        id='ma-div-loadMore',
                        style={'textAlign':'center', 'display':'block', 'legendNum':0}
                    )
                ],
                style={
                    'maxHeight': '50vh',
                    'overflow': 'auto',  
                    'width':'210px',
                },
                direction='vertical'
            ),
            id='ma-card-legend',
            title='',
            visible=False,
            draggable=True,
            width='250px',
            transitionType='none',
            destroyOnClose=False,
            forceRender=True,
            closable=False,
            style={
                'top': '7vh',
                'left': '450px',
                'position': 'fixed',
            },
        ),
        fac.AntdSpace(
            [
                fac.AntdTitle('Spot Properties', level=5, style={'margin-top':1}),
                fac.AntdSpace(
                    [
                        fac.AntdTooltip(
                            fac.AntdSelect(
                                placeholder='Spot Size=3',
                                options=[1,2,3,4,5,6,7,8,9,10],
                                value=5,
                                allowClear=False,
                                id='ma-selecter-spotSize',
                                locale='en-us',
                                style={'width': '140px'},
                            ),
                            title=fac.AntdText('Spot Size'), color='white'
                        ),
                        fac.AntdTooltip(
                            fac.AntdSelect(
                                placeholder='Border Width=0',
                                options=[0,1],
                                allowClear=False,
                                value=1,
                                id='ma-selecter-borderWidth',
                                locale='en-us',
                                style={'width': '140px'},
                            ),  
                            title=fac.AntdText('Border Width'), color='white'
                        ),
                        fac.AntdTooltip(
                            fac.AntdColorPicker(
                                id='ma-ColorPicker-boarderColor',
                                showText=True,
                                disabledAlpha=False,
                                locale='en-us',
                                presets=[
                                    {'colors': ['#0d0015', '#474a4d', '#1c305c', '#640125'], 'label': 'presets'}
                                ],
                                value='#0d0015'
                            ),
                            title=fac.AntdText('Border Color'), color='white'
                        ),
                    ],
                    style={'width':'100%'}
                ),

                fac.AntdTitle('Slice Filter', level=5),
                fac.AntdSpace(
                    [
                        fac.AntdSlider(
                            id='ma-slider-z',
                            range=True,
                            min=0,
                            max=0,
                            tooltipPrefix='',
                            style={'width':279}
                        ),
                        dcc.Store(id='ma-slider-znum'),
                        dcc.Store(id='ma-slider-znums'),
                        fac.AntdButton(
                            'Slicer',
                            id='ma-button-slicer',
                            type='primary', 
                            style={'backgroundColor':'#ca8269'},
                            icon=fac.AntdIcon(icon='md-content-cut')
                        ),
                    ],
                    size='middle'
                ),
                fac.AntdTitle('Adjust Alignment', level=5),
                
                fac.AntdSpace(
                    [
                        fac.AntdSelect(
                            placeholder='Active Slice',
                            id='ma-selecter-activeSlice',
                            locale='en-us',
                            style={'width': '117px'},
                        ),
                        fac.AntdTooltip(
                            fac.AntdColorPicker(
                                id='ma-ColorPicker-actSlice',
                                disabledAlpha=False,
                                disabled=False,
                                locale='en-us',
                                presets=[
                                    {'colors': ['#ff6570', '#fa8c35', '#f1f907', '#06b16a', '#0bffd2', '#257cff', '#8c00ff', '#474a4d', '#afafb0', '#ffffff'], 'label': 'presets'}
                                ],
                                value='#f1f907'
                            ),
                            title=fac.AntdText('Active Slice Color'), color='white'
                        ),
                        fac.AntdIcon(
                            icon='md-trending-flat'
                        ),
                        fac.AntdSelect(
                            placeholder='Reference Slice',
                            id='ma-selecter-referenceSlice',
                            locale='en-us',
                            style={'width': '135px'},
                        ),
                        fac.AntdTooltip(
                            fac.AntdColorPicker(
                                id='ma-ColorPicker-refSlice',
                                disabledAlpha=False,
                                disabled=False,
                                locale='en-us',
                                presets=[
                                    {'colors': ['#ff6570', '#fa8c35', '#f1f907', '#06b16a', '#0bffd2', '#257cff', '#8c00ff', '#474a4d', '#afafb0', '#ffffff'], 'label': 'presets'}
                                ],
                                value='#afafb0'
                            ),
                            title=fac.AntdText('Reference Slice Color'), color='white'
                        ),
                    ],
                    size='middle'
                ),

                fac.AntdSpace(
                    [
                        fac.AntdTooltip(
                            fac.AntdSelect(
                                options=['Custom', 'Field', 'Gene'],
                                value = 'Custom',
                                id='ma-selecter-colorMode',
                                allowClear=False,
                                disabled=False,
                                locale='en-us',
                                style={'width': '117px', 'marginRight':'0px'},
                            ),
                            title=fac.AntdText('Select Color Mode'), color='white'
                        ),  
                        fac.AntdTooltip(
                            fac.AntdSelect(
                                placeholder='Pick Fields',
                                id='ma-selecter-colorField',
                                optionFilterMode='regex',
                                allowClear=False,
                                debounceWait=300,
                                locale='en-us',
                                disabled=True,
                                style={'width': '122px'},
                            ),
                            title=fac.AntdText('Select Colored Field'), color='white'
                        ), 
                        fac.AntdTooltip(
                            fac.AntdSelect(
                                placeholder='Pick Genes',
                                id='ma-selecter-colorGene',
                                allowClear=False,
                                debounceWait=300,
                                optionFilterMode='regex',
                                locale='en-us',
                                disabled=True,
                                style={'width': '122px'},
                            ),
                            title=fac.AntdText('Select Colored Gene'), color='white'
                        ),     
                    ],
                    size='middle'
                ),
                fac.AntdSelect(
                    defaultValue='C1',
                    allowClear=False,
                    disabled=True,
                    locale='en-us',
                    id='ma-selecter-colorRange',
                    options=[
                        {
                            'label': 'C1',
                            'value': 'C1',
                            'colors': redScaleColor,
                        },
                        {
                            'label': 'C2',
                            'value': 'C2',
                            'colors': greenScaleColor,
                        },
                        {
                            'label': 'C3',
                            'value': 'C3',
                            'colors': blueScaleColor,
                        }
                    ],
                    colorsMode='sequential',
                    style={'width': '392px'},
                ),
                fac.AntdDivider(),
                fac.AntdSpace(
                    [
                        fac.AntdSpace(
                        [
                            fac.AntdTooltip(
                                fac.AntdInputNumber(
                                        min=0,
                                        value=5,
                                        placeholder='Step Size',
                                        id = 'ma-inputNum-stepSize',
                                        precision=5,
                                        debounceWait=300,
                                        style={'width': 126, 'margin-right':21},
                                ),
                                title=fac.AntdText('Step Size'), color='white'
                            ),
                            fac.AntdButton(
                                id = 'ma-button-unclockwise',
                                type='primary', 
                                style={'width':50, 'backgroundColor':'#5F9EA0'},
                                icon=fac.AntdIcon(icon='pi-arrow-counter-clockwise')
                            ),
                            fac.AntdButton(
                                id='ma-button-up',
                                type='primary', 
                                style={'width':50, 'backgroundColor':'#5F9EA0'},
                                icon=fac.AntdIcon(icon='antd-arrow-up')
                            ),
                            fac.AntdButton(
                                id='ma-button-clockwise',
                                type='primary',
                                style={'width':50, 'backgroundColor':'#5F9EA0'},
                                icon=fac.AntdIcon(icon='pi-arrow-clockwise')
                            ),
                            fac.AntdTooltip(
                                fac.AntdButton(
                                    id='ma-button-reFocus',
                                    type='primary',
                                    # shape = 'circle',
                                    style={'backgroundColor':'#ca8269'},
                                    icon=fac.AntdIcon(icon='antd-aim')
                                ),
                                title = fac.AntdText('Reset View'), color='white'
                            ),
                        ],
                        size='middle',
                        ),
                        fac.AntdSpace(
                        [
                            fac.AntdTooltip(
                                fac.AntdInputNumber(
                                    min=0,
                                    id = 'ma-inputNum-rotationAngle',
                                    placeholder='Rotation Angle',
                                    value = 2,
                                    precision=5,
                                    debounceWait=300,
                                    style={'width': 126, 'margin-right':21},
                                ),
                                title=fac.AntdText('Rotation Angle'), color='white'
                            ),

                            fac.AntdButton(
                                id = 'ma-button-left',
                                type='primary',
                                style={'width':50, 'backgroundColor':'#5F9EA0'},
                                icon=fac.AntdIcon(icon='antd-arrow-left')
                            ),
                            fac.AntdButton(
                                id = 'ma-button-down',
                                type='primary', 
                                style={'width':50, 'backgroundColor':'#5F9EA0'},
                                icon=fac.AntdIcon(icon='antd-arrow-down')
                            ),
                            fac.AntdButton(
                                id = 'ma-button-right',
                                type='primary', 
                                style={'width':50, 'backgroundColor':'#5F9EA0'},
                                icon=fac.AntdIcon(icon='antd-arrow-right'),
                            ),
                            fac.AntdTooltip(
                                fac.AntdButton(
                                    id='ma-button-showGrid',
                                    type='primary',
                                    style={'backgroundColor':'#949495'},
                                    icon=fac.AntdIcon(icon='antd-number')
                                ),
                                title=fac.AntdText('Switch Grid'), color='white'
                            ), 
                        ],
                        size='middle',
                        )
                    ],
                    size='middle',
                    direction='vertical',
                ),
                fac.AntdTitle('Sync Slices', level=5),
                fac.AntdSpace(
                    [
                        fac.AntdSpace(
                            [
                                fac.AntdButton(
                                    'Pick Above',
                                    id='ma-button-pickAbove',
                                    type='primary', 
                                    style={'backgroundColor':'#507ea4'},
                                    icon=fac.AntdIcon(icon='antd-caret-up')
                                ),
                                fac.AntdTooltip(
                                    fac.AntdButton(
                                        id='ma-button-autoPickAbove',
                                        icon=fac.AntdIcon(id='ma-icon-autoPickAbove', icon='antd-check'),
                                        style={'marginLeft': '-10px'}
                                    ),
                                    title=fac.AntdText('Auto Pick Above'), color='white'
                                ),
                                fac.AntdButton(
                                    'Pick Below',
                                    id='ma-button-pickBelow',
                                    type='primary', 
                                    style={'backgroundColor':'#895b8a', 'marginLeft': '50px'},
                                    icon=fac.AntdIcon(icon='antd-caret-down')
                                ),
                                fac.AntdTooltip(
                                    fac.AntdButton(
                                        id='ma-button-autoPickBelow',
                                        icon=fac.AntdIcon(id='ma-icon-autoPickBelow', icon='antd-close'),
                                        style={'marginLeft': '-10px'}
                                    ),
                                    title=fac.AntdText('Auto Pick Below'), color='white'
                                ),
                                # fac.AntdButton(
                                #     'Apply',
                                #     id='ma-button-apply',
                                #     type='primary', 
                                #     style={'backgroundColor':'#a16d5d', 'width':110},
                                #     icon=fac.AntdIcon(icon='pi-broom')
                                # ),
                            ],
                            size='middle',
                        ),  
                        fac.AntdTable(
                            columns=[
                                {'title': 'slices', 'dataIndex': 'slices', 'align':'center'}
                            ],
                            # data=[
                            #     {'slices': 'slice'+str(i)} for i in range(1, 15)
                            # ],
                            style = {'width':389},
                            # filterOptions={'slices': {'filterMode': 'keyword'}},
                            id = 'ma-table-SyncSlices',
                            rowSelectionType='checkbox',
                            bordered=True,
                            pagination=False,
                            locale='en-us',
                            size='small',
                            maxHeight='25vh'
                        )
                    ],
                    size='middle',
                    direction='vertical', 
                )
                
            ],
            style={'width':'100%'},
            direction='vertical',
            size='middle'
        ), 
    ],
    title='Manual Adjust', 
    id='gs-drawer-manual',
    forceRender=True,
    visible=False,
    mask=False,
    width='450px',
    placement='left'
)