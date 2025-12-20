from controller.auth import *
import json

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

def clinetSendUpdateMenuItem(key, websocket, ws, ms, **kwargs):
    """
    客户端发送更新菜单项信息
    """
    menuItem = kwargs.get('menuItem', None)
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
    

