from dash import html
import feffery_antd_components as fac
from callbacks.passcode_c import *

control_panel = html.Div(
    fac.AntdCard(
        fac.AntdEmpty(image='simple', description='Empty Toolbox'),
        headStyle={'display': 'none'},
        style={'height':'93vh','display': 'flex', 'alignItems': 'center', 'justifyContent': 'center'},
    ), 
    id='passcode-control-panel'
)

content_panel = html.Div(
    fac.AntdSpace(
        [
            fac.AntdModal(
                [
                    fac.AntdSpace(
                        [
                            fac.AntdInput(
                                id='input-passcode', 
                                placeholder='Please Input Passcode'
                            ),
                            fac.AntdSelect(
                                options=[
                                    {'label': 'Editable', 'value': 'Editable'},
                                    {'label': 'ReadOnly', 'value': 'ReadOnly'},
                                ],
                                id='select-permission',
                                value='Editable',
                                allowClear=False,
                                style={'width':'100%'}
                            ),
                            fac.AntdSelect(
                                options=[
                                    {'label': 'Active', 'value': 'Active'},
                                    {'label': 'Disabled', 'value': 'Disabled'},
                                ],
                                id='select-status',
                                value='Active',
                                allowClear=False,
                                style={'width':'100%'}
                            ),
                        ],
                        direction='vertical',
                        size='middle',
                        style={'width':'80%', 'margin': '0 auto', 'display': 'flex', 'marginTop':'30px', 'marginBottom':'30px'}
                    ),
                    html.Div(
                        fac.AntdSpace(
                            [
                                fac.AntdButton('Cancel', id='create-cancel', type='primary', style={'backgroundColor':'#ab6953'}),
                                fac.AntdButton('Submit', id='create-submit', type='primary', style={'backgroundColor':'#507ea4'})
                            ],
                            size='middle'
                        ),
                        style={'text-align':'right'}
                    )
                ], 
                id='create-passcode', 
                title='Create New Passcode',
                maskClosable=False,
                visible=False
            ),
            fac.AntdSpace(
                [
                    fac.AntdButton(
                        'New Passcode', 
                        type='primary',
                        id='new-passcode',
                        icon=fac.AntdIcon(icon='antd-plus'),
                        style={'backgroundColor':'#698aab'}
                    ),
                    fac.AntdButton(
                        'Delete', 
                        type='primary',
                        id='delete-passcode',
                        icon=fac.AntdIcon(icon='antd-delete'),
                        style={'backgroundColor':'#b88884'}
                    ),
                ],
                size='middle',
                style={'marginBottom':'2vh'},
            ),
            fac.AntdTable(
                id='passcode-table',
                maxHeight='70vh',
                pagination=False,
                rowSelectionType='checkbox',
                filterOptions={'Passcode': {'filterMode': 'keyword'}},
                columns=[
                    {'title': 'id','dataIndex': 'id', 'width':'6%', 'renderOptions': {'renderType': 'ellipsis'}},
                    {'title': 'Passcode','dataIndex': 'Passcode', 'editable': True,'renderOptions': {'renderType': 'ellipsis-copyable'}},
                    {'title': 'Permission','dataIndex': 'Permission', 'renderOptions': {'renderType': 'select'}},
                    {'title': 'Status','dataIndex': 'Status', 'renderOptions': {'renderType': 'select'}},
                ],
                locale='en-us',
                bordered=True,
            ),
        ],
        direction='vertical',
        style={'width':'100%', 'marginTop':'5vh'},
    ),
    style={'width':'95%', 'height':'100%'},
)