import asyncio
import websockets
from utils.commonfuc import get_local_ip
from websocket.message import ms

class WebSocket:
    def __init__(self):
        self.clients = set()  # 存储所有ws连接

    def run(self):
        """
        启动WebSocket服务
        """
        async def startWebSocketServer():
            host, port = get_local_ip(), 8765
            print(f"✅ Starting WebSocket server on ws://{host}:{port}")
            async with websockets.serve(self.handleConnect, host, port):
                await asyncio.Future()
        asyncio.run(startWebSocketServer())
    
    async def handleConnect(self, websocket):
        """
        处理每个WebSocket连接
        """
        websocket.curMenuItem = None
        self.clients.add(websocket)
        await ms.serverSend(ms._loginState_s2c, websocket, self)
        try:
            async for message in websocket:
                await ms.serverParse(message, websocket, self)
        finally:
            self.clients.remove(websocket)
    
    async def notifyUpdateExpansionStatus(self, operation, key, value):
        """
        通知所有处于扩增任务页面的用户更新状态
        """
        await ms.serverSend(ms._updateExpansionStatus_s2c, ws=self)

ws = WebSocket()