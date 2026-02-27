from controller.lassion_ctl import *
from dash import callback, Input, Output, State, no_update, Patch
from controller.auth import *

datasets = lasData.get_dataset_list()
initial_dataset = datasets[0] if len(datasets)>0 else None

@callback(
    Input('lasso-button-download', 'nClicks'),
    State('lasso-graph-result', 'selectedData'),
    State('lasso-select-dataset', 'value'),
    running=[
        (Output('lasso-button-download', 'loading'), True, False),
    ]
)
def download_csv(nClicks, selectedData, dataset_name):
    """
    下载选中的点
    """
    if nClicks:
        if not selectedData:
            pointIndex = None
        if selectedData:
            points = selectedData.get('points', None)
            if points is None or len(points) == 0:
                pointIndex = None
            else:
                pointIndex = [point['customdata'] for point in points]
        select_cell(pointIndex, dataset_name)


@callback(
    Output('lasso-graph-result', 'figure', allow_duplicate=True),
    Input('lasso-select-spotSize', 'value'),
    Input('lasso-select-borderWidth', 'value'),
    Input('lasso-colorPicker-boarderColor', 'value'),
    State('lasso-select-dataset', 'value'),
    prevent_initial_call=True
)
def update_figure_bySpotConfig(spotSize, borderWidth, borderColor, dataset_name):
    """
    基于点配置集更新图
    """
    fig = lasData.get_3d_scatter(
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
    Output('lasso-graph-result', 'figure'),
    Output('lasso-graph-result', 'style'),
    Input('lasso-select-dataset', 'value'),
    State('lasso-select-spotSize', 'value'),
    State('lasso-select-borderWidth', 'value'),
    State('lasso-colorPicker-boarderColor', 'value'),
    running=[
        (Output('lasso-control-dataset-divider', 'children'), 'Loading...', 'Dataset Selection'),
    ]
)
def update_figure_byDataset(dataset_name, spotSize, borderWidth, borderColor):
    """
    基于数据集更新图
    """
    if not dataset_name:
        return no_update, no_update
    fig = lasData.get_3d_scatter(
        dataset_name=dataset_name,
        marker_size=spotSize,
        boarder_width=borderWidth,
        boarder_color=borderColor,
    )
    return fig, {'height':'100%', 'width': '100%', 'display':'block'}

@callback(
    Output('lasso-select-dataset', 'options'),
    Input('lasso-select-dataset-tooltip', 'open'),
)
def update_dataset_list(open):
    """
    更新数据集列表
    """
    if open:
        return lasData.get_dataset_list()
    return no_update