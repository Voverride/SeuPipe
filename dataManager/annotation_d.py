from anndata import AnnData
from utils.colors import get_color_map

class AnnotationData:
    def __init__(self):
        self._refdata = None
        self._refname = None
        self._querydata = None
        self._queryname = None
        self._labelfield = None
        self._xfield = None
        self._yfield = None
        self._zfield = None
        self._resultfig = None
        self._refumap = None
        self._queryumap = None
        self._heatmap = None
        self._anntstatus = {
            'creator':None,
            'thread':None,
            'rmMt':True,
            'rmHb':True,
            'rmRibo':True,
            'useHvg':True,
            'epochs':100,
            'batch_size':128,
            'n_layers':2,
            'dropout':0.2,
            'n_hiddens': 256,
            'n_latent': 20, 
            'exception':None,
            'steps':{
                1:{'running':False, 'complete':False, 'failed':False},
                2:{'running':False, 'complete':False, 'failed':False, 'percent':0, 'epoch': '0/0'},
                3:{'running':False, 'complete':False, 'failed':False},
            }
        }
    def set_annstatus_props(self, creator:str, thread, rmMt:bool, rmHb:bool, rmRibo:bool, useHvg:bool, epochs:int, batch_size:int, n_layers:int, dropout:float, n_hiddens:int, n_latent:int)->None:
        self._anntstatus['creator'] = creator
        self._anntstatus['thread'] = thread
        self._anntstatus['rmMt'] = rmMt
        self._anntstatus['rmHb'] = rmHb
        self._anntstatus['rmRibo'] = rmRibo
        self._anntstatus['useHvg'] = useHvg
        self._anntstatus['epochs'] = epochs
        self._anntstatus['batch_size'] = batch_size
        self._anntstatus['n_layers'] = n_layers
        self._anntstatus['dropout'] = dropout
        self._anntstatus['n_hiddens'] = n_hiddens
        self._anntstatus['n_latent'] = n_latent
    def get_anntstatus(self)->dict:
        return self._anntstatus

    def set_refumap(self, fig):
        self._refumap = fig

    def get_refumap(self):
        return self._refumap
    
    def get_queryumap(self):
        return self._queryumap
    
    def set_heatmap(self, fig):
        self._heatmap = fig
    def get_heatmap(self):
        return self._heatmap
    def set_queryumap(self, fig):
        self._queryumap = fig
    def set_resultfig(self, fig)->None:
        self._resultfig = fig
    def get_resultfig(self):
        return self._resultfig
    def reset_anntstatus(self)->None:
        self._anntstatus = {
            'creator':None,
            'thread':None,
            'rmMt':True,
            'rmHb':True,
            'rmRibo':True,
            'useHvg':True,
            'epochs':100,
            'batch_size':128,
            'n_layers':2,
            'dropout':0.2,
            'n_hiddens': 256,
            'n_latent': 20, 
            'exception':None,
            'steps':{
                1:{'running':False, 'complete':False, 'failed':False},
                2:{'running':False, 'complete':False, 'failed':False, 'percent':0, 'epoch': '0/0'},
                3:{'running':False, 'complete':False, 'failed':False},
            }
        }

    def get_cmap(self)->dict:
        ref = self.get_refdata()
        label_field = self.get_labelfield()
        if label_field is None:
            return None
        fields = set(ref.obs[label_field].unique())
        return get_color_map(fields, type='COLORS_60')

    def get_zfield_min(self)->float:
        """
            Retrieves the minimum value of the z-field in the AnnData object.

            Returns:
                float: Minimum value of the z-field, or None if the AnnData or z-field is not set.
        """
        if self._querydata is None or self._zfield is None:
            return None
        col = self._querydata.obs[self._zfield]
        return col.min()
    def get_zfield_max(self)->float:
        """
            Retrieves the maximum value of the z-field in the AnnData object.

            Returns:
                float: Maximum value of the z-field, or None if the AnnData or z-field is not set.
        """
        if self._querydata is None or self._zfield is None:
            return None
        col = self._querydata.obs[self._zfield]
        return col.max()
    def reset_allprops(self)->None:
        """
            Resets all internal properties of the object to their default values.
        """
        self.reset_refprops()
        self.reset_queryprops()
    def reset_refprops(self)->None:
        self._refdata = None
        self._refname = None
        self._labelfield = None
    def reset_queryprops(self)->None:
        self._querydata = None
        self._queryname = None
        self._xfield = None
        self._yfield = None
        self._zfield = None
    def set_labelfield(self, field:str)->None:
        self._labelfield = field
    
    def get_project_name(self)->str:
        refname = self.get_refname()
        queryname = self.get_queryname()
        if refname and queryname:
            return refname+'_to_'+queryname
        return None

    def set_xfield(self, field:str)->None:
        self._xfield = field

    def set_yfield(self, field:str)->None:
        self._yfield = field

    def set_zfield(self, field:str)->None:
        self._zfield = field

    def get_labelfield(self)->str:
        return self._labelfield
    
    def get_xfield(self)->str:
        return self._xfield
    
    def set_refname(self, name)->None:
        self._refname = name
    
    def set_queryname(self, name)->None:
        self._queryname = name

    def get_refname(self)->str:
        return self._refname
    
    def get_queryname(self)->str:
        return self._queryname

    def get_yfield(self)->str:
        return self._yfield

    def get_zfield(self)->str:
        return self._zfield
    def set_refdata(self, adata:AnnData)->None:
        self._refdata = adata
    
    def get_refdata(self)->AnnData:
        return self._refdata
    
    def set_querydata(self, adata:AnnData)->None:
        self._querydata = adata

    def get_querydata(self)->AnnData:
        return self._querydata
    
    def get_refdata_fields(self)->list:
        ad = self.get_refdata()
        if ad is None:
            return []
        obs_fields = ad.obs.columns.tolist()
        obs_fields.sort()
        return obs_fields
    
    def get_querydata_fields(self)->list:
        ad = self.get_querydata()
        if ad is None:
            return []
        obs_fields = ad.obs.columns.tolist()
        obs_fields.sort()
        return obs_fields

    

annData = AnnotationData()