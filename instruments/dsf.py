from itertools import groupby
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
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

    def plot_fit(self, scale=[5,2]):
        """
        melt curves are fit using the boltzman sigmoid by the following algorithm
            (1) melt region indices are detected by the first instance of a positive first derivative
                ! Currently trying other solutions
            (2) flourescense is normalized based on the melt region as defined by in (1)
            (3) indices of the melt region are expanded outwards by 2 datapoints
            (4) the expanded melt region is fit to a boltzman sigmoid equation
                ! Since the data was normalized, baseline and plateau are assumed 0 and 1
        """
        enends  = lambda x: (x[0][0], x[-1][0])
        mranges = lambda x: [enends(list(i)) for k, i in groupby(enumerate(inc), lambda x: x[1]) if k]
        longest = lambda x: x[np.argmax([i[1]-i[0] for i in x])] # for detection of tm shift, reconsider
        clmask  = lambda x: x if sum(x)!=0 else np.array([True]*x.shape[0])
        
        refzero = lambda x, y: x - y.min()
        refstnd = lambda x, y: refzero(x, y)/refzero(x, y).max()
        
        pushb   = lambda x: 0 if x-2<0 else  x-2 
        pusht   = lambda x, y: y-1 if x+2>y else  x+2
        push    = lambda x, y: (pushb(x[0]), pusht(x[1], y.shape[0]))

        print(self['Name'])
        figsize = self['Grid Size'][::-1]*scale
        fig, ax = plt.subplots(*self['Grid Size'], figsize=figsize, sharex=True)

        for _ in self.samples():
            x, y   = _['Data']
            
            d1     = np.gradient(y)
            inc    = clmask(d1>0)
            ystnd  = refstnd(y, y[inc])
            
            # t, b   = mranges(inc)[0]   
            t, b   = longest(mranges(inc))
            et, eb = push((t, b), inc) # extended rise for curve fitting
            
            ax[_['ind']].set_title(_['Sample'])
            ax[_['ind']].plot(x, ystnd, color='red', lw=0.35)
            # ax[_['ind']].plot(x[et:eb], ystnd[et:eb])
            
            ###### 
            
            model      = lambda x, tm, slope: (1/(1+np.exp((tm-x)/slope)))
            popt, pcov = curve_fit(model, 
                                   x[et:eb],
                                   ystnd[et:eb], 
                                   p0 = (x[t:b].mean(), 0.5),
                                   method='lm', 
                                   maxfev=9999)
            fx = np.arange(x.min(), x.max(), 0.1)
            fy = np.array([model(i, *popt) for i in fx])
            ax[_['ind']].plot(fx, fy, color='darkgrey', lw=2)
            
            if popt[0] < x.max() and popt[0] > x.min() :
                ax[_['ind']].plot([popt[0]]*2, (-2,2), color='black', ls=':')
                ax[_['ind']].text(popt[0], 0.9, '%.2f  ' % popt[0], ha='right', va='top')
            
            ###### 

            ax[_['ind']].set_ylim(-0.1, 1.1)
            ax[_['ind']].set_yticklabels([])
            ax[_['ind']].set_yticks([])
        ax[_['ind']].set_xlim(x.min(), x.max())
        return plt.show()
        
