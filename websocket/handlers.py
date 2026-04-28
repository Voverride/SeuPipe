from controller.auth import *
from setting import setting
from dash import set_props, Patch
import json
from dataManager.coordinate_d import get_coordinate, get_initialfig

async def serverSendUpdateManualAdjustStatus(key, websocket, ws, ms, **kwargs):
    """
    服务端发送更新手动调整任务状态信息
    """
    project = kwargs.get('project', None)
    operations = kwargs.get('operations', None)
    username = kwargs.get('username', None)
    data = dict(
        key=key,
        operations=operations,
        username=username,
    )
    for connected in ws.clients:
        if connected.curMenuItem==setting.alignment:
            ws_project = connected.params.get('project', None)
            if project==ws_project:
                await connected.send(json.dumps(data))

def clientParseUpdateManualAdjustStatus(message, websocket, ws, ms, **kwargs):
    """
    客户端解析更新手动调整任务状态信息
    """
    operations = message.get('operations', None)
    username = message.get('username', None)
    if operations is not None:
        patch = Patch()
        project = operations.get('project', None)
        x = operations.get('xfield', None)
        y = operations.get('yfield', None)
        sliceidx = operations.get('sliceidx', [])
        initFig = get_initialfig(project)
        coordinate = get_coordinate(project)
        for idx in sliceidx:
            obsIndex = initFig['data'][idx]['customdata'].flatten()
            coords = coordinate.loc[obsIndex, [x, y]]
            x_coords = coords[x].tolist()
            y_coords = coords[y].tolist()
            patch['data'][int(idx)]['x'] = x_coords
            patch['data'][int(idx)]['y'] = y_coords
        set_props('ali-graph-left', dict(figure=patch))
        if username is not None:
            set_props('ali-live-username', dict(children=username))

async def serverSendUpdateAlignmentStatus(key, websocket, ws, ms, **kwargs):
    """
    服务端发送更新对齐任务状态信息
    """
    project = kwargs.get('project', None)
    data = dict(
        key=key,
        update=True,
    )
    for connected in ws.clients:
        if connected.curMenuItem==setting.alignment:
            ws_project = connected.params.get('project', None)
            if project==ws_project:
                await connected.send(json.dumps(data))

def clientParseUpdateAlignmentStatus(message, websocket, ws, ms, **kwargs):
    """
    客户端解析更新对齐任务状态信息
    """
    update = message.get('update', False)
    if update:
        set_props('refresh-ali-status', dict(children=None))

async def serverSendUpdateAnnotationStatus(key, websocket, ws, ms, **kwargs):
    """
    服务端发送更新细胞注释任务状态信息
    """
    project = kwargs.get('project', None)
    data = dict(
        key=key,
        update=True,
    )
    for connected in ws.clients:
        if connected.curMenuItem==setting.annotation:
            ws_project = connected.params.get('project', None)
            if project==ws_project:
                await connected.send(json.dumps(data))

def clientParseUpdateAnnotationStatus(message, websocket, ws, ms, **kwargs):
    """
    客户端解析更新细胞注释任务状态信息
    """
    update = message.get('update', False)
    if update:
        set_props('refresh-ann-status', dict(children=None))

async def serverSendUpdateClusteringStatus(key, websocket, ws, ms, **kwargs):
    """
    服务端发送更新聚类任务状态信息
    """
    project = kwargs.get('project', None)
    data = dict(
        key=key,
        update=True,
    )
    for connected in ws.clients:
        if connected.curMenuItem==setting.cellSelector:
            ws_project = connected.params.get('project', None)
            if project==ws_project:
                await connected.send(json.dumps(data))

def clientParseUpdateClusteringStatus(message, websocket, ws, ms, **kwargs):
    """
    客户端解析更新聚类任务状态信息
    """
    update = message.get('update', False)
    if update:
        set_props('refresh-sel-status', dict(children=None))

async def serverSendUpdateSegmentationStatus(key, websocket, ws, ms, **kwargs):
    """
    服务端发送更新分割任务状态信息
    """
    project = kwargs.get('project', None)
    data = dict(
        key=key,
        update=True,
    )
    for connected in ws.clients:
        if connected.curMenuItem==setting.segmentation:
            ws_project = connected.params.get('project', None)
            if project==ws_project:
                await connected.send(json.dumps(data))

def clientParseUpdateSegmentationStatus(message, websocket, ws, ms, **kwargs):
    """
    客户端解析更新分割任务状态信息
    """
    update = message.get('update', False)
    if update:
        set_props('seg-refresh-status', dict(children=None))

async def serverSendUpdateCellSelectorStatus(key, websocket, ws, ms, **kwargs):
    """
    服务端发送更新细胞选择任务状态信息
    """
    project = kwargs.get('project', None)
    slice = kwargs.get('slice', None)
    data = dict(
        key=key,
        update=True,
    )
    for connected in ws.clients:
        if connected.curMenuItem==setting.cellSelector:
            ws_project = connected.params.get('project', None)
            ws_slice = connected.params.get('slice', None)
            if project==ws_project and slice==ws_slice:
                await connected.send(json.dumps(data))

def clientParseUpdateCellSelectorStatus(message, websocket, ws, ms, **kwargs):
    """
    客户端解析更新细胞选择任务状态信息
    """
    update = message.get('update', False)
    if update:
        set_props('refresh-sel-status', dict(children=None))

async def serverSendUpdateRegionClipStatus(key, websocket, ws, ms, **kwargs):
    """
    服务端发送更新区域剪切任务状态信息
    """
    taskName = kwargs.get('taskName', None)
    slice = kwargs.get('slice', None)
    clipName = kwargs.get('clipName', None)
    data = dict(
        key=key,
        update=True,
    )
    for connected in ws.clients:
        if connected.curMenuItem=='RegionClip':
            ws_taskName = connected.params.get('taskName', None)
            ws_sliceName = connected.params.get('sliceName', None)
            ws_clipName = connected.params.get('clipName', None)
            if taskName==ws_taskName and slice==ws_sliceName and clipName==ws_clipName:
                await connected.send(json.dumps(data))

def clientParseUpdateRegionClipStatus(message, websocket, ws, ms, **kwargs):
    """
    客户端解析更新区域剪切任务状态信息
    """
    update = message.get('update', False)
    if update:
        set_props('refresh-clip-status', dict(children=None))

async def serverSendUpdateExpansionStatus(key, websocket, ws, ms, **kwargs):
    """
    服务端发送更新胞域扩增任务状态信息
    """
    data = dict(
        key=key,
        update=True
    )
    for connected in ws.clients:
        if connected.curMenuItem=='Expansion':
            await connected.send(json.dumps(data))


def clientParseUpdateExpansionStatus(message, websocket, ws, ms, **kwargs):
    """
    客户端解析更新胞域扩增任务状态信息
    """
    update = message.get('update', False)
    if update:
        set_props('init-restore-expansion', dict(disabled=False))

def clinetSendUpdateParams(key, websocket, ws, ms, **kwargs):
    """
    客户端发送更新参数信息
    """
    params = kwargs.get('params', None)
    if params is None:
        return
    data = dict(
        key=key,
        params=params,
    )
    set_props('ws', dict(send=json.dumps(data)))

async def serverParseUpdateParams(message, websocket, ws, ms, **kwargs):
    """
    服务端解析更新参数信息
    """
    params = message.get('params', None)
    if params is not None:
        websocket.params = params

def clinetSendUpdateMenuItem(key, websocket, ws, ms, **kwargs):
    """
    客户端发送更新菜单项信息
    """
    menuItem = kwargs.get('menuItem', None)
    if menuItem is None:
        return
    data = dict(
        key=key,
        menuItem=menuItem,
    )
    set_props('ws', dict(send=json.dumps(data)))

async def serverParseUpdateMenuItem(message, websocket, ws, ms, **kwargs):
    """
    服务端解析更新菜单项信息
    """
    menuItem = message.get('menuItem', None)
    if menuItem is not None:
        websocket.curMenuItem = menuItem

async def serverSendSyncLoginState(key, websocket, ws, ms, **kwargs):
    """
    服务端发送同步登录状态信息
    """
    host, port = websocket.remote_address
    user = search_user(usrhost=host)
    username = user[0]['usrname'] if user else ''
    data = dict(
        key=key,
        syncLogin=True,
    )
    for connected in ws.clients:
        curHost, curPort = connected.remote_address
        if curPort==port:
            continue
        curUser = search_user(usrhost=curHost)
        if not curUser or curUser[0]['usrname']==username:
            await connected.send(json.dumps(data))
            
def clientParseSyncLoginState(message, websocket, ws, ms, **kwargs):
    """
    客户端解析同步登录状态信息
    """
    syncLogin = message.get('syncLogin', False)
    if syncLogin:
        set_props('main-refresh', dict(href='/'))

def clinetSendLoginState(key, websocket, ws, ms, **kwargs):
    """
    客户端发送用户登录状态信息
    """
    data = dict(
        key=key,
        login=True,
    )
    set_props('ws', dict(send=json.dumps(data)))

async def serverParseLoginState(message, websocket, ws, ms, **kwargs):
    """
    服务端解析用户登录状态信息
    """
    login = message.get('login', False)
    if login:
        await ms.serverSend(ms._syncLoginState_s2c, websocket, ws)

async def serverSendLogoutClient(key, websocket, ws, ms, **kwargs):
    """
    服务端发送用户登出通知
    """
    host = websocket.remote_address[0]
    data = dict(
        key=key,
        logout=True,
    )
    for connected in ws.clients:
        if connected.remote_address[0] == host:
            await connected.send(json.dumps(data))

def clientParseLogoutClient(message, websocket, ws, ms, **kwargs):
    """
    客户端解析用户登出通知
    """
    logout = message.get('logout', False)
    if logout:
        set_props('main-refresh', dict(href='/'))

def clinetSendLogoutState(key, websocket, ws, ms, **kwargs):
    """
    客户端发送用户登出状态信息
    """
    data = dict(
        key=key,
        logout=True,
    )
    set_props('ws', dict(send=json.dumps(data)))

async def serverParseLogoutState(message, websocket, ws, ms, **kwargs):
    """
    解析用户登出状态信息
    """
    logout = message.get('logout', False)
    if logout:
        await ms.serverSend(ms._logoutClient_s2c, websocket, ws)
async def serverSendLoginState(key, websocket, ws, ms, **kwargs):
    """
    服务端发送用户登录状态信息
    """
    host = websocket.remote_address[0]
    isLogin = verify_host(host)
    data = dict(
        key=key,
        isLogin=isLogin,
    )
    await websocket.send(json.dumps(data))

def clientParseLoginState(message, websocket, ws, ms, **kwargs):
    """
    客户端解析用户登录状态信息
    """
    isLogin = message.get('isLogin', False)
    if isLogin:
        restore_usrinfo()
    else:
        set_props('login-box', dict(visible=True))
    

