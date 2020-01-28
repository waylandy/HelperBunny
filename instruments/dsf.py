from itertools import groupby
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

from IPython.display import display
np.set_printoptions(precision=6, suppress=True)

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
        def gradient(x, y):
            frame  = lambda x: [(x[n-1 if n>1 else 0:n+2]) for n, xi in enumerate(x)] 
            avg    = lambda x: sum(x)/len(x)
            slope  = lambda x, y: -avg([(xj-xi)/(yi-yj) for (xi,xj),(yi, yj) in zip(zip(x,x[1:]),zip(y,y[1:]))])
            turn   = lambda x, y: (y[1]>y[0] and y[1]>y[2]) or (y[1]<y[0] and y[1]<y[2])
            der    = lambda x, y: slope(x,y) if len(x)!=3 else 0 if turn(x,y) else slope(x,y)
            return np.array([der(xf, yf) for xf, yf in zip(frame(x),frame(y))])
            
        def detector(x, y):
            ends   = lambda x: (x[0][0], x[-1][0])
            rise   = lambda x: np.array([(0,0)]+[ends(list(i)) for k, i in groupby(enumerate(x),lambda x:x[1]>=0) if k])
            melt   = lambda x: x[np.subtract(*x.T[::-1])>1]
            melts  = melt(rise(gradient(x, y)))
            try:
                ml, mr = melts[np.argmax([y[j]-y[i] for i, j in melts])]
            except:
                ml, mr = 0, y.shape[0]
            return x[ml:mr+1], y[ml:mr+1]
        
        def rescale(x, y, xm, ym):
            b, t   = ym.min(), ym.max()
            return (y-b)/(t-b)
        
        report  = []
        figsize = self['Grid Size'][::-1]*scale
        fig, ax = plt.subplots(*self['Grid Size'], figsize=figsize, sharex=True)
        fig.suptitle(self['Name'], fontsize="xx-large")
        for _ in self.samples():
            ax[_['ind']].set_title(_['Sample'])
            x, y = _['Data']
            
            inc_   = lambda x: [j>i for i,j in zip(x, x[1:])]+[True]
            inc    = inc_(x)
            xm, ym = detector(x[inc], y[inc])
            yr     = rescale(x, y, xm, ym)
            
            ################################################################## 
            model      = lambda x, tm, slope: (1/(1+np.exp((tm-x)/slope)))
            norm       = lambda x: (x-x.min())/(x.max()-x.min())
            popt, pcov = curve_fit(model, 
                                   xm,
                                   norm(ym), 
                                   p0 = ((xm[0]+xm[-1])/2, 0.5),
                                   method='lm', 
                                   maxfev=9999,
                                  )
            fx = np.arange(x.min(), x.max(), 0.1)
            fy = np.array([model(i, *popt) for i in fx])
            ax[_['ind']].plot(fx, fy, color='darkgrey', lw=2.4)
            
            if popt[0] < x.max() and popt[0] > x.min() :
                ax[_['ind']].plot([popt[0]]*2, (yr.min()-0.1,yr.max()+0.1), color='black', ls=':')
                ax[_['ind']].text(popt[0], yr.min()+(yr.max()-yr.min()+0.2)*0.75,
                                  '%.2f  ' % popt[0], ha='right', va='top')
                report.append((_['Sample'], round(popt[0], 2)))
            ##################################################################

            ax[_['ind']].plot(x, yr, color='blue', lw=0.35)
            ax[_['ind']].set_xlim(x.min(), x.max())
            # ax[_['ind']].set_ylim(-0.1,1.1)
            ax[_['ind']].set_ylim(yr.min()-0.1,yr.max()+0.1)
            ax[_['ind']].set_yticklabels([])
            ax[_['ind']].set_yticks([])
            ax[_['ind']].tick_params(axis='x', length=0)
            
        display(pd.DataFrame(report, columns=[self['Name'],'Tm']).set_index(self['Name']))
        return plt.show()

