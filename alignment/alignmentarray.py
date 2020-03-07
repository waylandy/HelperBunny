import sys
from itertools import groupby
import numpy as np

from scipy.stats import entropy
from scipy.spatial.distance import jensenshannon

from  .. import fig
from .exporter import AlignmentExporter

class AlignmentArray(np.ndarray):
    """
        uses references as array element dtype
    """

    def __new__(cls, ndarray, names=None):
        obj = np.asarray(ndarray, dtype=object).view(cls)
        return obj

    def __array_finalize__(self, obj):
        if obj is None: 
            return

############################################ low level routines

    def ispos(self):
        # boolean array of whether each column is a position in the alignment
        # makes judgement using ONLY the first sequence
        ispos = lambda x: False if len(x)!=1 else not x.islower()
        return np.array(list(map(ispos, self[0])))

    def cols(self):
        # gives the position number for each column in the array
        # "positions" will be int, "insertions" will be float
        pos = self.ispos().astype(object)
        pos[pos==True] = np.arange(1,1+pos[pos==True].shape[0])
        i = 0
        for n, v in enumerate(pos):
            if type(v)==bool:
                pos[n]=i+0.5
            else:
                i=v
        return pos

    def positions(self, f=0, fxn=None):
        # remove insertion columns from the alignment array
        return self[:,self.ispos()]

############################################ editting the alignment

    def make_insertion(self, start, end):
        # defines the start to end indices will be redefined as inserts
        # index locations start and end are INCLUSIVE
        # if you define the start of an insertion next to an existing insertion, this will automatically extend the boundary
        # this makes sure 2 insertions will never be next to each other
        assert len(self.shape)==2
        reform    = lambda x: ''.join(i.lower() for i in x if i!='-')
        compress  = lambda x: np.array([reform(i) for i in x], dtype=object)
        hasi      = frozenset(n for n, i in enumerate(self.ispos()) if not i)
        while start-1 in hasi:
            start-= 1
        while end+1 in hasi:
            end  += 1
        col       = np.arange(start, end+1)
        r, *d     = col
        self[:,r] = compress(self[:,col])
        if start==end:
            return self    
        return np.delete(self, d, axis=1)

    def define_insertions(self, gap=0.5):
        # automatically defines inserts based on the proportion of gaps at each alignment position
        # ADD OPTION TO KEEP COLUMNS BASED ON LOW ENTROPY POSITIONS
        assert len(self.shape) == 2
        insert = lambda x: (x=='-').mean(0)
        mask   = insert(self)>gap
        insi   = [list(ind) for b, ind in groupby(enumerate(mask), lambda x: x[1])]
        insi   = [(i[0][1], (i[0][0], i[-1][0])) for i in insi if i[0][1]]
        for b, i in insi[::-1]:
            self = self.make_insertion(*i)
        return self

############################################ creating a position array for statistics

    def PositionArray(self, v=1):
        # remove insertions from the alignment array and return a PositionArray
        mask = self.ispos()
        tran = (mask==True).sum(), self.shape[1]
        if v>0:
            sys.stderr.write('PositionArray           : Found %s POSITIONS from %s COLUMNS\n'% tran)
            sys.stderr.write('PositionArray           : Removing %s INSERTIONS from alignment\n' % (tran[1]-tran[0]))
        data = np.array(self[:,self.ispos()], dtype='|U1')
        return PositionArray(data)

############################################ writing files

    def export(self, name, output=None, format=None):
        ext = lambda x: x.split('.')[-1]
        if format==None:
            if output==None:
                format = 'cfa'
            else:
                format = ext(output)
        if output==None:
            output=sys.stdout
            cl = False
        else:
            output=open(output, 'w')
            cl = True

        exp = AlignmentExporter()
        if   format in ('cfa'):
            exp.to_cfa(self, name, output)
        elif format in ('fasta', 'fa'):
            exp.to_fasta(self, name, output)
        elif format in ('cma', 'mma'):
            exp.to_xma(self, name, output)

        if cl:
            output.close()

############################################ visualizing the alignment

    def show_logo(self, **kwargs):
        return self.PositionArray(v=0).show_logo(**kwargs)

    def show_aln(self, **kwargs):
        return self.PositionArray(v=0).show_aln(**kwargs)

    def show(self, **kwargs):
        return self.PositionArray(v=0).show(**kwargs)

############################################    ############################################

alphabet = 'ARNDCQEGHILKMFPSTWYV-'

background = ('ARNDCQEGHILKMFPSTWYV', np.array( # BLOSUM62 DISTRIBUTION
            [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062,
             0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]))

class PositionArray(np.ndarray):
    """
        uses single char as array element dtype
    """

    def __new__(cls, ndarray, names=None):
        obj = np.asarray(ndarray).view(cls)
        return obj

    def __array_finalize__(self, obj):
        if obj is None: 
            return

    def gaps(self, gap='-'):
        # returns %gaps at each alignment position
        return (self==gap).mean(axis=0)

    def frequency(self, alphabet=alphabet, pseudocount=0):
        # returns position frequency matrix as int AND the alphabet used
        return pseudocount+np.array([(i==self).sum(axis=0) for i in alphabet]), alphabet

    def probability(self, alphabet=alphabet, pseudocount=0):
        # returns position probability matrix as float AND the alphabet used
        ans, _ = self.frequency(alphabet, pseudocount=pseudocount)
        return ans/ans.sum(0), alphabet

    def weight(self, background=background, pseudocount=0):
        # returns position weight matrix normalized against background as float AND the alphabet used
        ans, _ = self.probability(alphabet=background[0], pseudocount=pseudocount)
        return (ans.T/background[1]).T, background[0]

    def kldivergence(self, *arg, background=background):
        # calculates kullback leibler divergence of each column to the background; does not yield bidirectional equality
        # if argument is given, compares each column of the 2 alignments, number of columns must be same in both
        l = len(arg)
        if l>1:
            raise Exception('Provide only 1 argument to be compared against or none to compare against predetermined background (blosum62)')
        elif l==1:
            assert arg[0].shape[1], self.shape[1] 
            b = self.probability(alphabet=background[0])[0].T
            q = arg[0].probability(alphabet=background[0])[0].T # idea: allow direct passing of a probability matrix to speed up calculation
            return np.array([entropy(*i) for i in zip(b, q)]) # self is considered the background
        elif l==0:
            p = self.probability(alphabet=background[0])[0].T
            return np.array([entropy(i,background[1]) for i in p])

    def jsdivergence(self, *arg, background=background):
        # calculates jensen shannon divergence of each column to the background; will yield bidirectional equality
        # if argument is given, compares each column of the 2 alignments, number of columns must be same in both
        l = len(arg)
        if l>1:
            raise Exception('Provide only 1 argument to be compared against or none to compare against predetermined background (blosum62)')
        elif l==1:
            assert arg[0].shape[1], self.shape[1] 
            b = self.probability(alphabet=background[0])[0].T
            q = arg[0].probability(alphabet=background[0])[0].T
            return np.array([jensenshannon(*i) for i in zip(b, q)])
        elif l==0:
            p = self.probability(alphabet=background[0])[0].T
            return np.array([jensenshannon(i,background[1]) for i in p])

    def consensus(self, norm=False, background=background):
        cons  = lambda x: ''.join(background[0][i] for i in x.argmax(0))
        if norm:
            return cons(self.weight(background=background)[0])
        else:
            return cons(self.probability(alphabet=background[0])[0])

############################################ visualizing the alignment

    def show_logo(self, **kwargs):
        # takes in a position array, arguments are passed to SequenceLogoViewer
        plot = fig.SequenceLogoViewer(self, **kwargs)
        fig.show(plot)

    def show_aln(self, **kwargs):
        # takes in a position array, arguments are passed to AlignmentViewer
        plot = fig.AlignmentViewer(self, **kwargs)
        fig.show(plot)

    def show(self, plot_width=1000):
        # takes in a position array
        plot1   = fig.SequenceLogoViewer(self, plot_width=plot_width, plot_height=230)
        plot2   = fig.AlignmentViewer(self, plot_width=plot_width, scale=20)
        plot1.x_range = plot2.x_range
        plot1.xaxis.visible = False
        plot = fig.gridplot([[plot1], [plot2]], toolbar_location=None)
        fig.show(plot)






















