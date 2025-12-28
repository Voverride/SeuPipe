import dash
from pages.main_p import main_layout
import argparse
import warnings
from numba import errors
import tensorflow as tf
from tensorflow.python.ops.numpy_ops import np_config
gpus = tf.config.list_physical_devices('GPU')
if gpus:
    tf.config.experimental.set_memory_growth(gpus[0], True)
np_config.enable_numpy_behavior()
warnings.filterwarnings("ignore", category=errors.NumbaWarning)
from controller.segmentation_ctl import parse_tasklist
from controller.regionclip_ctl import parse_regionclip_tasklist
from utils.commonfuc import get_local_ip
from flask import request
from io import TextIOWrapper, BytesIO
from controller.notice import global_error_handler
from dash.dependencies import Input, Output, State
from websocket.websocket import ws
import os
import threading

app = dash.Dash(
    __name__, 
    title='SeuPipe',
    update_title=None,
    suppress_callback_exceptions=True,
    on_error=global_error_handler,
    prevent_initial_callbacks='initial_duplicate'
)

app.layout = main_layout

app.clientside_callback(
    """(nClicks, collapsed) => {
        return [!collapsed, collapsed ? 'antd-arrow-left' : 'antd-arrow-right'];
    }""",
    [
        Output('main-sider-collapse', 'collapsed'),
        Output('main-icon-menuItem', 'icon'),
    ],
    Input('main-button-trigger', 'nClicks'),
    State('main-sider-collapse', 'collapsed'),
    prevent_initial_call=True,
)

app.clientside_callback(
    """
    function () {
        document.addEventListener('keydown', function(e) {
            let iniCode = 0
            if(e.ctrlKey) {
                iniCode = 1000
            }
            keyCode = JSON.stringify(iniCode+e.keyCode)
            dash_clientside.set_props("key-pressed-events", {data: keyCode})
        });
        return dash_clientside.no_update;
    }
    """,
    Output('SeuPipe', 'id'),
    Input('SeuPipe', 'id'),
)

@app.server.route('/upload/', methods=['POST'])
def upload():
    '''
    构建细胞分割清单上传服务
    :return:
    '''
    try:
        file = request.files['file']
        text_stream = TextIOWrapper(BytesIO(file.read()), encoding='utf-8')
        lines = text_stream.readlines()
        parse_tasklist(lines)
    except Exception as e:
        raise e
    return {'filename': ''}

@app.server.route('/upload/regionClip', methods=['POST'])
def uploadRegionClip():
    '''
    构建图像裁剪清单上传服务
    :return:
    '''
    try:
        file = request.files['file']
        text_stream = TextIOWrapper(BytesIO(file.read()), encoding='utf-8')
        lines = text_stream.readlines()
        parse_regionclip_tasklist(lines)
    except Exception as e:
        raise e
    return {'filename': ''}

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Run a Dash app with custom parameters.")
    parser.add_argument('--host', type=str, default=get_local_ip(), help="Host address to run the app on (default: auto-detected network IP).")
    parser.add_argument('--port', type=int, default=8088, help="Port to run the app on (default: 8088).")
    parser.add_argument('--debug', action='store_true', default=False, help="Enable or disable debug mode (default: False).")
    
    args = parser.parse_args()
    reload = True if args.debug else False

    if os.environ.get("WERKZEUG_RUN_MAIN") is None:
        ws_thread = threading.Thread(target=ws.run, daemon=True)
        ws_thread.start()

    app.run(
        host=args.host, 
        port=args.port,
        threaded=True,
        debug=args.debug, 
        use_reloader=reload
    )