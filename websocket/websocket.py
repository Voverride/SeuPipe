import asyncio
import websockets
from utils.commonfuc import get_local_ip

class WebSocket:
    def __init__(self):
        self.connected_clients = {}  # 存储所有连接的客户端, key:username, value:websocket
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
        print(f"新连接: {websocket.remote_address}")
        try:
            async for message in websocket:
                # 广播给所有人
                for client in self.connected_clients.values():
                    await client.send(f"Broadcast: {message}")  # 👈 发给所有人
        finally:
            # 清理断开的连接
            # connected_clients.remove(websocket)
            pass

ws = WebSocket()