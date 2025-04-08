import dash
from pages.main_p import main_layout
import argparse
from controller.notice import global_error_handler
from dash.dependencies import Input, Output, State

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

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Run a Dash app with custom parameters.")
    parser.add_argument('--host', type=str, default="127.0.0.1", help="Host address to run the app on (default: 127.0.0.1).")
    parser.add_argument('--port', type=int, default=8088, help="Port to run the app on (default: 8088).")
    parser.add_argument('--debug', type=bool, default=False, help="Enable or disable debug mode (default: False).")
    
    args = parser.parse_args()
    reload = True if args.debug else False

    app.run(
        host=args.host, 
        port=args.port,
        threaded=True,
        debug=args.debug, 
        use_reloader=reload
    )