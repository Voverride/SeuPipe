from dash import html
import feffery_antd_components as fac

control_panel = html.Div(
    ''
)

content_panel = html.Div(
    fac.AntdEmpty(
        description=fac.AntdText('当前页面开发中', type='secondary'),
        imageStyle={'height': 250},
    )
)