import asyncio
import websockets
from utils.commonfuc import get_local_ip
from websocket.message import ms
from setting import setting

class WebSocket:
    def __init__(self):
        self.clients = set()  # 存储所有ws连接

    def run(self):
        """
        启动WebSocket服务
        """
        async def startWebSocketServer():
            host, port = get_local_ip(), setting.ws_port
            async with websockets.serve(
                self.handleConnect, host, port,
                compression=None,
                open_timeout=30,
                ping_interval=30,
                ping_timeout=10,
                close_timeout=5,
                max_size=5 * 1024 * 1024,
                max_queue=32
            ):
                await asyncio.Future()
        asyncio.run(startWebSocketServer())
    
    async def handleConnect(self, websocket):
        """
        处理每个WebSocket连接
        """
        try:
            websocket.curMenuItem = None
            websocket.params = {}
            self.clients.add(websocket)
            await ms.serverSend(ms._loginState_s2c, websocket, self)
            async for message in websocket:
                await ms.serverParse(message, websocket, self)
        finally:
            self.clients.remove(websocket)

    async def notifyUpdateManualAdjustStatus(self, operation, key, value, **kwargs):
        """
        通知处于对齐任务页面的用户更新手动调整状态
        """
        await ms.serverSend(ms._updateManualAdjustStatus_s2c, ws=self, **kwargs)

    async def notifyUpdateAlignmentStatus(self, operation, key, value, **kwargs):
        """
        通知处于对齐任务页面的用户更新状态
        """
        await ms.serverSend(ms._updateAlignmentStatus_s2c, ws=self, **kwargs)

    async def notifyUpdateAnnotationStatus(self, operation, key, value, **kwargs):
        """
        通知处于细胞注释任务页面的用户更新状态
        """
        await ms.serverSend(ms._updateAnnotationStatus_s2c, ws=self, **kwargs)
    
    async def notifyUpdateCellSelectorStatus(self, operation, key, value, **kwargs):
        """
        通知处于细胞选择页面的用户更新状态
        """
        await ms.serverSend(ms._updateCellSelectorStatus_s2c, ws=self, **kwargs)
    
    async def notifyUpdateExpansionStatus(self, operation, key, value, **kwargs):
        """
        通知处于扩增任务页面的用户更新状态
        """
        await ms.serverSend(ms._updateExpansionStatus_s2c, ws=self)

    async def notifyUpdateSegmentationStatus(self, operation, key, value, **kwargs):
        """
        通知处于分割任务页面的用户更新状态
        """
        await ms.serverSend(ms._updateSegmentationStatus_s2c, ws=self, **kwargs)
    
    async def notifyUpdateRegionClipStatus(self, operation, key, value, **kwargs):
        """
        通知处于区域裁剪任务页面的用户更新状态
        """
        await ms.serverSend(ms._updateRegionClipStatus_s2c, ws=self, **kwargs)

    async def notifyUpdateClusteringStatus(self, operation, key, value, **kwargs):
        """
        通知处于聚类任务页面的用户更新状态
        """
        await ms.serverSend(ms._updateClusteringStatus_s2c, ws=self, **kwargs)

ws = WebSocket()