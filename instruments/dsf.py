import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class DSF(dict):
    """
    Written to parse .xls outputs from
        Applied Biosystems(TM) StepOnePlus(TM) Real-Time PCR System
    
    Should also be relatively easy to parse their proprietary .eds run files
    Its just a gzipped archieve with xml files inside
    
    I only tested this for doing differential scanning fluorimetry (thermal shift assay) on proteins
    My own stuff used SyPRO Orange Dye (just tell the machine its ROX)
    
    Curve-fitting assumes a single global unfolding event
    """
    
    def __init__(self, eds_file):
        eds_name     = lambda x: eds_file.split('/')[-1]
        eds_sheet    = lambda x: pd.read_excel(eds_file, x, skiprows=7)
        loc_readings = lambda x: [i for i in x.columns if i.startswith('Reading ')]
        readings     = lambda x: x[['Well Location']+loc_readings(x)].set_index('Well Location').T

        let2ind      = lambda x: 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'.index(x)
        num2ind      = lambda x: int(x) - 1
        grid2ind     = lambda x: (let2ind(x[0]), num2ind(x[1]))

        well_locs    = lambda x: np.array([grid2ind(i) for i in x.columns.to_list()])
        well_topleft = lambda x: np.array((x[:,0].min(), x[:,1].min()))
        well_botrght = lambda x: np.array((x[:,0].max(), x[:,1].max()))

        transpose    = lambda *x: list(map(list, list(zip(*x))))
        checkwells   = lambda *x: all(len(set(i))==1 for i in transpose(*[i.columns.values for i in x]))
        checkshape   = lambda *x: all(len(set(i))==1 for i in transpose(*[i.shape for i in x]))

        reads_data   = readings(eds_sheet('Melt Region Normalized Data'))
        reads_temp   = readings(eds_sheet('Melt Region Temperature Data'))

        assert checkwells(reads_data, reads_temp)
        assert checkshape(reads_data, reads_temp)

        self['Name'         ]  = eds_name(eds_file)
        self['Data'         ]  = np.array([i for i in zip(reads_temp.values.T, reads_data.values.T)])
        self['Well Location']  = reads_data.columns.values
        self['Well Indices' ]  = well_locs(reads_temp)
        self['Well Indices' ] -= well_topleft(self['Well Indices'])
        self['Grid Size'    ]  = well_botrght(self['Well Indices'])+1

    def label_axes(self, x=[], y=[]):
        name   = lambda *x: ' & '.join(i for i in x if i!='')
        xs, ys = self['Grid Size'][::-1]
        x , y  = ['']*xs if len(x)==0 else x, ['']*ys if len(y)==0 else y
        assert len(x)==xs and len(y)==ys
        self['Labels'] = [name(y[i[0]],x[i[1]]) for i in self['Well Indices']]

    def samples(self):
        get_ind = (lambda x: x[0]) if (self['Grid Size'][1]==1) else (lambda x: tuple(x))
        for n in range(self['Data'].shape[0]):
            yield {
                'Sample' : self['Well Location'][n] if 'Labels' not in self else (
                         '%s : %s' % (self['Well Location'][n], self['Labels'][n])),
                'Data'   : self['Data'][n],
                'ind'    : get_ind(self['Well Indices'][n]),
            }
        
    def plot_raw(self, scale=[5,2]):
        print(self['Name'])
        figsize = self['Grid Size'][::-1]*scale
        fig, ax = plt.subplots(*self['Grid Size'], figsize=figsize, sharex=True)
        
        for _ in self.samples():
            x, y  = _['Data']
            ax[_['ind']].set_title(_['Sample'])
            ax[_['ind']].plot(x, y)
            
            ax[_['ind']].set_yticklabels([])
            ax[_['ind']].set_yticks([])
            
        ax[_['ind']].set_xlim(x.min(), x.max())
        return plt.show()

