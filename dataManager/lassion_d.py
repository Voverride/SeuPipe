from anndata import AnnData
from utils.colors import get_color_map
import plotly.express as px
import pickle
import os

class LassionData:
    def __init__(self):
        self.figures = {}
    
    def get_dataset_list(self):
        dataset_names = os.listdir('dataManager/vispkl')
        dataset_names.sort()
        return dataset_names
    
    def get_dataset_data(self, dataset_name:str):
        with open(os.path.join('dataManager/vispkl', dataset_name), 'rb') as f:
            df = pickle.load(f)
        return df
    
    def set_3d_scatter(self, dataset_name:str, fig):
        self.figures[dataset_name] = fig
    
    def get_3d_scatter(self, dataset_name:str, marker_size, boarder_width, boarder_color):
        if not dataset_name:
            return None
        if dataset_name in self.figures:
            fig = self.figures[dataset_name]
        else:
            fig = self.plot_3d_scatter(dataset_name)
            self.set_3d_scatter(dataset_name, fig)
        fig.update_traces(
            marker=dict(
                size=marker_size,
                line=dict(color=boarder_color, width=boarder_width)
            )
        )
        return fig

    def plot_3d_scatter(self, dataset_name, marker_size:int=5, boarder_width:int=1, boarder_color:str='#0d0015', grid_color:str = '#5F9EA0'):
        """
        绘制三维散点图
        """
        df = self.get_dataset_data(dataset_name)
        label_field = 'celltype'
        fields = set(df[label_field].unique())
        cmap = get_color_map(fields, type='COLORS_60')
        df = df.sort_values(by=label_field)
        fig = px.scatter(
            df, x='x', y='y', 
            color=label_field, 
            color_discrete_map=cmap,
            custom_data=[df.index.tolist()]
        )

        fig.update_traces(
            marker=dict(
                size=marker_size,
                line=dict(color=boarder_color, width=boarder_width)
            )
        )
        fig.update_layout(
            autosize=True,
            legend=dict(
                itemsizing='constant',
                title=None
            ),
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
        )
        return fig

lasData = LassionData()
