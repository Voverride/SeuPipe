from dash import html, dcc
import feffery_antd_components as fac
from .components.annotationTask import AnnotationTask
from callbacks.annotation_c import *

control_panel = html.Div(
    fac.AntdCard(
        fac.AntdSpace(
            [
                dcc.Interval(id="init-restore-annotation", interval=1, max_intervals=1),
                dcc.Interval(interval=1000, disabled=True, id='annotation-interval'),
                dcc.Interval(id="annotation-event-loop", interval=1000),
                fac.AntdButton(
                    'Create New Task', 
                    type='primary',
                    id='annotation-button-newtask',
                    icon=fac.AntdIcon(icon='antd-plus'),
                    style={'backgroundColor':'#698aab', 'width': '100%'}
                ),
                fac.AntdDivider('Spot Properties', innerTextOrientation='left'),
                # fac.AntdTitle('Spot Properties', level=5, style={'margin-top':1}),
                fac.AntdTooltip(
                    fac.AntdSelect(
                        options=[1,2,3,4,5,6,7,8,9,10],
                        value=5,
                        id='annotation-select-spotSize',
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
                                id='annotation-select-borderWidth',
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
                                id='annotation-colorPicker-boarderColor',
                            ),
                            title=fac.AntdText('Border Color'), color='white'
                        ),
                    ],
                    size='middle',
                    style={'width':'100%'}
                ),
                # fac.AntdTitle('Slice Filter', level=5),
                fac.AntdDivider('Slice Filter', innerTextOrientation='left'),
                fac.AntdSlider(
                    range=True,
                    min=0,
                    max=0,
                    tooltipPrefix='',
                    id='annotation-slider-slicer',
                    style={'width':'100%'}
                ),
                fac.AntdButton(
                     'Slicer',
                     type='primary', 
                     id='annotation-button-slicer',
                     style={'backgroundColor':'#867ba9', 'width':'100%'},
                     icon=fac.AntdIcon(icon='md-content-cut')
                 ),
                fac.AntdDivider(),
                fac.AntdButton(
                    'Export Data', 
                    type='primary',
                    id='annotation-button-exportData',
                    block=True,
                    icon=fac.AntdIcon(icon='antd-save'),
                    style={'backgroundColor':'#ca8269'}
                )
            ],
            size='middle',
            direction='vertical',
            style={'width':'100%'},
        ),
        headStyle={'display': 'none'},
        style={'height':'93vh', 'maxHeight': '93vh','overflowY': 'auto'},
    ), 
    id='annotation-control_panel',
)

content_panel = html.Div(
    [
        AnnotationTask,
        fac.AntdSplitter(
            items=[
                {
                    'children': fac.AntdCenter(
                        dcc.Graph(
                            config={'displaylogo':False}, 
                            id='annotation-graph-result', 
                            style={'display':'none'},
                            responsive=True
                        ),
                        style={'height': '100%', 'width': '100%'}
                    ),
                    'defaultSize': '60%',
                    'collapsible': True,
                },
                {
                    # 'children': fac.AntdSplitter(
                    #     items=[
                    #         {
                    #             'children': fac.AntdCenter(
                    #                 dcc.Graph(
                    #                     config={'displaylogo':False}, 
                    #                     # id='annotation-graph-refumap', 
                    #                     style={'width': '100%', 'height':'90%', 'margin':'auto', 'maxWidth':'50vw'},
                    #                     responsive=True
                    #                 ),
                    #                 style={'height': '100%', 'width': '100%'}
                    #             ),
                    #             'collapsible': True,
                    #         },
                    #         {
                    #             'children': fac.AntdCenter(
                    #                 dcc.Graph(
                    #                     config={'displaylogo':False}, 
                    #                     # id='annotation-graph-queryumap', 
                    #                     style={'width': '100%', 'height':'90%', 'margin':'auto', 'maxWidth':'50vw'},
                    #                     responsive=True
                    #                 ),
                    #                 style={'height': '100%', 'width': '100%'}
                    #             ),
                    #             'collapsible': True,
                    #         },
                    #     ],
                    #     layout='vertical',
                    #     style={'height': '100%', 'width': '100%'},
                    # ),
                    'children': fac.AntdCenter(
                        dcc.Graph(
                            config={'displaylogo':False}, 
                            id='annotation-graph-heatmap', 
                            style={'display':'none'},
                            responsive=True
                        ),
                        style={'height': '100%', 'width': '100%'}
                    ),
                    'defaultSize': '40%',
                    'collapsible': True,
                },
            ],
            style={'width':'100%', 'height':'100%'}
        )
    ],
    id='annotation-content-panel',
    style={'width':'95%', 'height':'100%', 'position':'relative'}
)