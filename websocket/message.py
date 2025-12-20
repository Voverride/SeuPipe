import json
from websocket.handlers import *

class Message:
    def __init__(self):
        self._loginState_s2c = 'login_state_s2c'
        self._logout_c2s = 'logout_c2s'
        self._logoutClient_s2c = 'logout_client_s2c'
        self._login_c2s = 'login_c2s'
        self._syncLoginState_s2c = 'sync_login_state_s2c'
        self._updateMenuItem_c2s = 'update_menuitem_c2s'
        self._updateExpansionStatus_s2c = 'update_expansion_status_s2c'
        self._handlers = {
            # 建立ws连接时检查登录状态
            self._loginState_s2c: {
                'send': serverSendLoginState,
                'parse': clientParseLoginState
            },
            # 登出时客户端发送登出状态
            self._logout_c2s: {
                'send': clinetSendLogoutState,
                'parse': serverParseLogoutState
            },
            # 服务端通知客户端登出
            self._logoutClient_s2c: {
                'send': serverSendLogoutClient,
                'parse': clientParseLogoutClient
            },
            # 登录时客户端发送登录状态
            self._login_c2s: {
                'send': clinetSendLoginState,
                'parse': serverParseLoginState
            },
            # 同步登录状态, 注销相同账号异地登录状态
            self._syncLoginState_s2c: {
                'send': serverSendSyncLoginState,
                'parse': clientParseSyncLoginState
            },
            # 更新每个websocket的菜单项
            self._updateMenuItem_c2s: {
                'send': clinetSendUpdateMenuItem,
                'parse': serverParseUpdateMenuItem
            },
            # 更新胞域扩增任务状态
            self._updateExpansionStatus_s2c: {
                'send': serverSendUpdateExpansionStatus,
                'parse': clientParseUpdateExpansionStatus
            }
        }
    
    async def serverSend(self, key, websocket=None, ws=None, **kwargs):
        """
        服务端发送消息
        """
        if key not in self._handlers:
            return
        await self._handlers[key]['send'](key, websocket, ws, self, **kwargs)
    
    def clientParse(self, message, websocket=None, ws=None, **kwargs):
        """
        客户端解析消息
        """
        msg_data = json.loads(message['data'])
        key = msg_data.get('key', None)
        if key is None or key not in self._handlers:
            return
        self._handlers[key]['parse'](msg_data, websocket, ws, self, **kwargs)

    def clinetSend(self, key, websocket=None, ws=None, **kwargs):
        """
        客户端发送消息
        """
        if key not in self._handlers:
            return
        self._handlers[key]['send'](key, websocket, ws, self, **kwargs)

    async def serverParse(self, message, websocket=None, ws=None, **kwargs):
        """
        服务端解析消息
        """
        msg_data = json.loads(message)
        key = msg_data.get('key', None)
        if key is None or key not in self._handlers:
            return
        await self._handlers[key]['parse'](msg_data, websocket, ws, self, **kwargs)


ms = Message()