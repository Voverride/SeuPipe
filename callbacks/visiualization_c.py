from controller.visiualization_ctl import *
from dash import callback, Input, Output, State, no_update, Patch
from controller.auth import *

datasets = visData.get_dataset_list()
initial_dataset = datasets[0] if len(datasets)>0 else None

@callback(
    Output('visiualization-graph-result', 'figure', allow_duplicate=True),
    Output('visiualization-store-relayoutData', 'data'),
    Input('visiualization-graph-result', 'relayoutData'),
    prevent_initial_call=True
)
def backup_relayoutData(relayoutData):
    """
    备份 relayoutData 数据
    """
    if relayoutData is None:
        return no_update, no_update
    fig_patch = Patch()
    update_relayout(fig_patch, relayoutData)
    return fig_patch, relayoutData


@callback(
    Output('visiualization-graph-result', 'figure', allow_duplicate=True),
    Input('visiualization-select-spotSize', 'value'),
    Input('visiualization-select-borderWidth', 'value'),
    Input('visiualization-colorPicker-boarderColor', 'value'),
    State('visiualization-select-dataset', 'value'),
    prevent_initial_call=True
)
def update_figure_bySpotConfig(spotSize, borderWidth, borderColor, dataset_name):
    """
    基于点配置集更新图
    """
    fig = visData.get_3d_scatter(
        dataset_name=dataset_name,
        marker_size=spotSize,
        boarder_width=borderWidth,
        boarder_color=borderColor,
    )
    fig_patch = no_update
    def update_patch(fig, patch, use_border=False):
        for i in range(len(fig['data'])):
            patch['data'][i]['marker']['size'] = spotSize
            if use_border:
                patch['data'][i]['marker']['line']['width'] = borderWidth
                patch['data'][i]['marker']['line']['color'] = borderColor
            
    if fig is not None:
        fig_patch = Patch()
        update_patch(fig, fig_patch, use_border=True)
    return fig_patch

@callback(
    Output('visiualization-graph-result', 'figure'),
    Output('visiualization-graph-result', 'style'),
    Input('visiualization-select-dataset', 'value'),
    State('visiualization-select-spotSize', 'value'),
    State('visiualization-select-borderWidth', 'value'),
    State('visiualization-colorPicker-boarderColor', 'value'),
    State('visiualization-graph-result', 'relayoutData'),
    running=[
        (Output('visiualization-control-dataset-divider', 'children'), 'Loading...', 'Dataset Selection'),
    ]
)
def update_figure_byDataset(dataset_name, spotSize, borderWidth, borderColor, relayoutData):
    """
    基于数据集更新图
    """
    if not dataset_name:
        return no_update, no_update
    fig = visData.get_3d_scatter(
        dataset_name=dataset_name,
        marker_size=spotSize,
        boarder_width=borderWidth,
        boarder_color=borderColor,
    )
    if relayoutData:
        layout_update = {}

        if 'scene.camera' in relayoutData:
            layout_update.setdefault('scene', {})['camera'] = relayoutData['scene.camera']

        if 'scene.aspectratio' in relayoutData:
            layout_update.setdefault('scene', {})['aspectmode'] = 'manual'
            layout_update['scene']['aspectratio'] = relayoutData['scene.aspectratio']

        if layout_update:
            fig.update_layout(layout_update)
    return fig, {'height':'100%', 'width': '100%', 'display':'block'}

@callback(
    Output('visiualization-select-dataset', 'options'),
    Input('visiualization-select-dataset-tooltip', 'open'),
)
def update_dataset_list(open):
    """
    更新数据集列表
    """
    if open:
        return visData.get_dataset_list()
    return no_update