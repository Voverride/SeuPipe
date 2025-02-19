import dash
import traceback
import json
import feffery_antd_components as fac

def set_head_notice(content:str, type:str)->None:
    """
    Set the header notice with specified content and type.

    This function updates the 'main-notice-area' component by setting its children
    to an Ant Design Message component with the provided content and type.

    Parameters:
    - content (str): The text message to be displayed in the notice.
    - type (str): The type of the notice, which determines its appearance (e.g., success, error, warning).

    Returns:
    - None
    """
    dash.set_props('main-notice-area', {'children':fac.AntdMessage(content=content, type=type)})


def global_error_handler(err:Exception)->None:
    """
    This function is designed to be used as a global error handler in a Dash
    application. When an error occurs, it creates an Ant Design notification
    with the error message, traceback, and input information. The notification
    will be displayed in the main error area for 30 seconds.

    Args:
        err (Exception): The error that occurred during the Dash callback.

    Returns:
        None
    """
    callback_context = dash.ctx
    style={
        'position': 'fixed',
        'top': '1vh',
        'right': '0',
        'width': '20vw',
        'zIndex': 9999,
        'display': 'flex',
    }
    info = str(err)
    input = json.dumps(callback_context.triggered)
    trace = traceback.format_exc()
    dash.set_props(
        'main-error-area', dict(style=style)
    )
    dash.set_props(
        'main-error-area-info', dict(children=info)
    )
    dash.set_props(
        'main-error-area-input', dict(children=input)
    )
    dash.set_props(
        'main-error-area-traceback', dict(children=trace)
    )