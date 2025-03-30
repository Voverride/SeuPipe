from dash import html
import feffery_antd_components as fac
from dash.dependencies import Input, Output, State
from dash import callback, set_props, no_update
from dash.exceptions import PreventUpdate
from dataManager.workspace import set_workpase
from controller.auth import verify_modify_permission
from controller.notice import set_aside_notice
from controller.alignment_ctl import read_alignment_file, export_alignment_file
from controller.annotation_ctl import read_annotask_refdata, read_annotask_querydata, export_annotation_file
import os

class FileSelecter:
    def __init__(self, id:str):
        self.boxid = id+'1'
        self.fpathid = id+'2'
        self.submitid = id+'3'
        self.tableid = id+'4'
        self.confirmid = id+'5'
        self.replaceid = id+'6'
        self.cancleid = id+'7'
        self.filenameid = id+'8'
        self.annotask = None # ref or query
        self.box = html.Div(
            [   
                fac.AntdModal(
                    fac.AntdSpace(
                        [
                            fac.AntdSpace(
                                [
                                    fac.AntdInput(
                                        id=self.fpathid,
                                        value='/',
                                        persistence=True,
                                        placeholder='please input file path',
                                        prefix=fac.AntdIcon(icon='di-linux'),
                                        style={'width': '660px'},
                                    ),
                                    fac.AntdButton(
                                        'submit',
                                        type='primary',
                                        id=self.submitid,
                                        disabled=True,
                                        style={'width': '90px'},
                                        icon=fac.AntdIcon(icon='md-launch')
                                    ),
                                ],
                            ),
                            html.Div(
                                [
                                    fac.AntdModal(
                                        [ 
                                            html.P(
                                               [
                                                   fac.AntdIcon(icon='fc-high-priority', style={'fontSize': 20}),
                                                   fac.AntdText(' The file '),
                                                   fac.AntdText('', id=self.filenameid, code=True),
                                                   fac.AntdText(' already exists in this location. Do you want to replace it?'),
                                               ] 
                                            ),
                                            html.Div(
                                                fac.AntdSpace(
                                                    [
                                                        fac.AntdButton('Cancle', id=self.cancleid),
                                                        fac.AntdButton('Replace', id=self.replaceid, type='primary')
                                                    ],
                                                    size='middle'
                                                ),
                                                style={'text-align':'right'}
                                            )
                                        ], 
                                        id=self.confirmid, 
                                        title='File Already Exists',
                                        locale = 'en-us',
                                        maskClosable=False,
                                        closable=False,
                                    ),
                                    fac.AntdTable(
                                        id=self.tableid,
                                        columns=[
                                            {
                                                'title': '',
                                                'dataIndex': 'fileInfo',
                                                'renderOptions': {'renderType': 'button'},
                                                'align': 'left'
                                            }
                                        ],
                                        bordered=True,
                                        pagination=False,
                                        maxHeight='45vh'
                                    )
                                ], style={'height':'50vh'}
                            )
                        ],
                        direction='vertical'
                    ),
                    id=self.boxid, 
                    title='Import Data', 
                    width='800px',
                    visible=False
                )
            ],
        )
        
        empty_fileList = [
            {
                'fileInfo': {
                    'content': 'there was nothing found',
                    'type': 'link',
                    'icon':'fc-high-priority',
                    'disabled':True
                }
            }
        ]

        @callback(
            Output(self.confirmid, 'visible', allow_duplicate=True),
            Output(self.boxid, 'visible', allow_duplicate=True),
            Input(self.replaceid, 'nClicks'),
            State(self.fpathid, 'value'),
            State('main-title-header', 'children'), 
            running=[
                (Output(self.replaceid, 'loading'), True, False),
                (Output(self.cancleid, 'disabled'), True, False),
                (Output(self.replaceid, 'style'), {'backgroundColor':'#a86965'}, None),
            ],
            prevent_initial_call=True,
        )
        def confirm_overwrite(nc, path, page):
            if nc:
                status = False
                if page=='Alignment':
                    status = export_alignment_file(path)
                elif page=='Annotation':
                    status = export_annotation_file(path)
                if status:
                    return False, False
            raise PreventUpdate

        @callback(
            Output(self.confirmid, 'visible', allow_duplicate=True),
            Input(self.cancleid, 'nClicks'),
            prevent_initial_call=True,
        )
        def cancel_confirm(nc):
            if nc:
                return False
            return no_update

        @callback(
            Output('main-loading-area', 'children'),
            Output(self.confirmid, 'visible'),
            Output(self.filenameid, 'children'),
            Input(self.submitid, 'nClicks'),
            State(self.boxid, 'title'),
            State(self.fpathid, 'value'),
            State('main-title-header', 'children'), 
            running=[
                (Output(self.submitid, 'loading'), True, False),
                (Output(self.submitid, 'style'), {'backgroundColor':'#a86965'}, None),
            ]
        )
        def submit_filePath(nClicks, title, path, page):
            visible = no_update
            filename = no_update
            if nClicks:
                if title=='Import Data':
                    if page=='Alignment':
                        status = read_alignment_file(path)
                        if status:
                            set_props(self.boxid, {'visible':False})
                    if page=='Annotation':
                        type = self.annotask
                        status = False
                        if type=='ref':
                            status = read_annotask_refdata(path)
                        elif type=='query':
                            status = read_annotask_querydata(path)
                        if status:
                            set_props(self.boxid, {'visible':False})

                elif title=='Export Data':
                    if page=='Alignment':
                        if os.path.exists(path):
                            visible = True
                            filename = path
                        else:
                            status = export_alignment_file(path)
                            if status:
                                set_props(self.boxid, {'visible':False})
                    elif page=='Annotation':
                        if os.path.exists(path):
                            visible = True
                            filename = path
                        else:
                            status = export_annotation_file(path)
                            if status:
                                set_props(self.boxid, {'visible':False})
                elif title=='Select Workspace':
                    set_workpase(path)
                    set_props(self.boxid, {'visible':False})
                    set_aside_notice('Spatpy workspace has been set to '+path+'SpatpyWorkspace', 'success', duration=20)
                    
            return no_update, visible, filename

        @callback(
            Output(self.submitid, 'disabled'), 
            Input(self.fpathid, 'value'),
            State(self.boxid, 'title'),
        )
        def update_submitButton_by_filePath(path:str, title:str):
            if title=='Export Data':
                if path.endswith('.h5ad'):
                    return False
                return True
            elif title=='Import Data':
                if os.path.exists(path) and path.endswith('.h5ad'):
                    return False
                return True
            elif title=='Select Workspace':
                if os.access(path, os.W_OK):
                    return False
                return True

        @callback(
            Output(self.fpathid, 'value'),
            Input(self.tableid, 'nClicksButton'),
            State(self.tableid, 'clickedContent'),
            State(self.fpathid, 'value'),
            prevent_initial_call=True,
        )
        def update_filePath_by_button_click(n_clicks, content, filePath):
            folderPath = os.path.dirname(filePath)
            newPath = os.path.join(folderPath, content)
            if os.path.isdir(newPath):
                return newPath+'/'
            return newPath

        @callback(
            Output(self.tableid, 'data'), 
            Output(self.tableid, 'columns'), 
            Input(self.fpathid, 'value'),
            Input(self.boxid, 'title'),
            State(self.tableid, 'columns'),
        )
        def update_fileList_by_inputPath(path, title, columns):
            filename = os.path.basename(path).lower()
            folderpath = os.path.dirname(path)
            if os.access(folderpath, os.R_OK):
                columns[0]['title']=folderpath
                fileList = [file for file in os.listdir(folderpath) if filename in file.lower() and os.access(os.path.join(folderpath, file), os.R_OK)]
                if title=='Select Workspace':
                    fileList = [file for file in fileList if os.path.isdir(os.path.join(folderpath, file))]
                fileList.sort()
                data = [
                    {
                        'fileInfo': {
                            'content': file,
                            'type': 'link',
                            'icon':'fc-folder' if os.path.isdir(os.path.join(folderpath, file)) else 'fc-file'
                        }
                    } for file in fileList
                ]
                if len(data)==0:
                    data = empty_fileList
                return data, columns
            else:
                columns[0]['title']=''
                return empty_fileList, columns

    def set_annotask(self, type:str)->None:
        self.annotask = type

    def get_annotask(self)->str:
        return self.annotask

    def get_boxid(self) -> str:
        return self.boxid
    
    def get_fpathid(self)->str:
        return self.fpathid
    
    def get_submitid(self)->str:
        return self.submitid
    
    def get_tableid(self)->str:
        return self.tableid
    
    def open_import_box(self):
        permission = verify_modify_permission()
        if permission:
            set_props(self.boxid, {'visible':True})
            set_props(self.boxid, {'closable':True})
            set_props(self.boxid, {'maskClosable':True})
            set_props(self.boxid, {'keyboard':True})
            set_props(self.boxid, {'title':'Import Data'})

    def open_export_box(self):
        permission = verify_modify_permission()
        if permission:
            set_props(self.boxid, {'visible':True})
            set_props(self.boxid, {'closable':True})
            set_props(self.boxid, {'maskClosable':True})
            set_props(self.boxid, {'keyboard':True})
            set_props(self.boxid, {'title':'Export Data'})
    
    def open_workspace_box(self):
        set_props(self.boxid, {'visible':True})
        set_props(self.boxid, {'closable':False})
        set_props(self.boxid, {'maskClosable':False})
        set_props(self.boxid, {'keyboard':False})
        set_props(self.boxid, {'title':'Select Workspace'})

    def get_box(self) -> html.Div:
        return self.box
    
fileSelecter = FileSelecter('fileSelecter')