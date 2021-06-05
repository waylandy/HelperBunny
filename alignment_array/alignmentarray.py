import sys
from itertools import groupby

import numpy as np
from scipy.stats import entropy
from scipy.spatial.distance import jensenshannon

from .sequencelogo import SequenceLogo

AMINO       = 'ARNDCQEGHILKMFPSTWYV-'

BLOSUM62_BG = 'ARNDCQEGHILKMFPSTWYV', np.array(
            #  A      R      N      D      C      Q      E      G      H      I
            [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062,
            #  L      K      M      F      P      S      T      W      Y      V
             0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072])

class AlignmentArray(np.ndarray):
    
    #################################################################################
    #####  inherit np.ndarray                                                   #####
    #################################################################################
    
    def __new__(cls, ndarray, names=None):
        obj = np.asarray(ndarray, dtype=object).view(cls)
        return obj

    def __array_finalize__(self, obj):
        # indexing and ufuncs pass thru this
        if obj is None: 
            return
        self.info = getattr(obj, 'info', None)
        
    def __array_wrap__(self, out_arr, context=None):
        # only ufuncs pass thru this
        # TO DO: test behavior on all possible ufuncs applicable to sequence alignments
 
        # documentation claims i should return:
        # super(AlignmentArray, self).__array_wrap__(self, out_arr, context)
        # this breaks code after second run (need more testing)

        if (out_arr.dtype == object):
            return out_arr
        else:
            return out_arr.view(np.ndarray)

    #################################################################################
    #####  basic functions & properties                                         #####
    #################################################################################

    def is_position(self):
        # returns boolean array of whether each column is a position in the alignment
        ispos = lambda x: False if len(x)!=1 else not x.islower()
        return np.array(list(map(ispos, self[0])))

    def remove_inserts(self):
        # returns alignment array without inserts; needed for doing hard math
        return self[:,self.is_position()]
    
    def alignment_positions(self): # cols
        # returns the alignment position number for each column in the array
        is_pos    = self.is_position()
        positions = np.cumsum(is_pos).astype(float)
        return positions + (~is_pos * 0.5)

    def get_positions(self, *args):
        # returns columns based on the alignment positions <int>
        if len(args)==0:
            raise Exception('didnt say which alignment positions you wanted')
        else:
            ndx = {int(j):i for i,j in enumerate(self.alignment_positions()) if j%1!=0.5}
            return self[:,[ndx[i] for i in args]]
        
    def get_range(self, start, end):
        # returns columns (with inserts) based on the range of alignment positions <int>
        assert start <= end
        assert {int}==set(map(type,(start,end)))
        aln_pos = self.alignment_positions()
        get_ind = lambda x: np.where(x==aln_pos)[0][0]
        return self[:,get_ind(start):1+get_ind(end)]

    #################################################################################
    #####  alignment editing                                                    #####
    #################################################################################

    def define_inserts(self, gap=0.5):
        # redefine inserts based on the proportion of gaps at each alignment position
        # lower gap cutoff = more positions removed
        inserts = ((self=='-').mean(0) > gap) + ~self.is_position()
        grouper = groupby(enumerate(inserts), lambda x: x[1])
        ranges  = [(k, list(i[0] for i in g)) for k, g in grouper]
        unalign = lambda x: ''.join(filter(lambda _: _!='-', x)).lower()
        squash  = lambda x: np.array([[unalign(i)] for i in x], dtype=object)
        proc    = lambda x: squash(self[:,x[1]]) if x[0] else np.array(self[:,x[1]], dtype=object)
        new_aln = AlignmentArray(np.concatenate(tuple(map(proc,ranges)), axis=1))
        assert all(unalign(i)==unalign(j) for i, j in zip(new_aln,self))
        return new_aln
    
    #################################################################################
    #####  position arrays for statistics                                       #####
    #################################################################################
    
    def frequency(self, alphabet=AMINO):
        # returns position frequency matrix
        assert all(self.is_position())
        F = np.array([(i==self).sum(axis=0) for i in alphabet])
        return F, np.array(list(alphabet))

    def probability(self, alphabet=AMINO, pseudocount=0):
        # returns position probability matrix
        assert all(self.is_position())
        F = np.array([(i==self).sum(axis=0) for i in alphabet]) + pseudocount
        P = F / F.sum(0)
        return P, np.array(list(alphabet))

    def entropy(self, alphabet=AMINO, pseudocount=0):
        # returns the shannon entropy of each position
        assert all(self.is_position())
        F = np.array([(i==self).sum(axis=0) for i in alphabet]) + pseudocount
        P = F / F.sum(0)
        S = -(P * np.log(P, where=P!=0)).sum(0)
        return S, np.array(list(alphabet))
    
    def nats(self, alphabet=AMINO, pseudocount=0):
        # returns the information in each position in nats
        assert all(self.is_position())
        F = np.array([(i==self).sum(axis=0) for i in alphabet]) + pseudocount
        P = F / F.sum(0)
        I = np.log(len(alphabet))+(P * np.log(P, where=P!=0)).sum(0)
        return I

    def bits(self, alphabet=AMINO, pseudocount=0):
        # returns the information in each position in bits
        assert all(self.is_position())
        F = np.array([(i==self).sum(axis=0) for i in alphabet]) + pseudocount
        P = F / F.sum(0)
        I = np.log2(len(alphabet))+(P * np.log2(P, where=P!=0)).sum(0)
        return I

    def kldivergence(self, *arg, background=BLOSUM62_BG):
        # calculates kullback leibler divergence of each column to the background; does not yield bidirectional equality
        # if argument is given, compares each column of the 2 alignments, number of columns must be same in both
        assert all(self.is_position())
        if len(arg)>1:
            raise Exception('Provide only 1 argument to be compared against or none to compare against the BLOSUM62 background')
        elif len(arg)==1:
            assert arg[0].shape[1], self.shape[1] 
            b = self.probability(alphabet=background[0])[0].T
            q = arg[0].probability(alphabet=background[0])[0].T
            return np.array([entropy(*i) for i in zip(b, q)])
        elif len(arg)==0:
            p = self.probability(alphabet=background[0])[0].T
            return np.array([entropy(i,background[1]) for i in p])

    def jsdivergence(self, *arg, background=BLOSUM62_BG):
        # calculates jensen shannon divergence of each column to the background; will yield bidirectional equality
        # if argument is given, compares each column of the 2 alignments, number of columns must be same in both
        assert all(self.is_position())
        if len(arg)>1:
            raise Exception('Provide only 1 argument to be compared against or none to compare against the BLOSUM62 background')
        elif len(arg)==1:
            assert arg[0].shape[1], self.shape[1] 
            b = self.probability(alphabet=background[0])[0].T
            q = arg[0].probability(alphabet=background[0])[0].T
            return np.array([jensenshannon(*i) for i in zip(b, q)])
        elif len(arg)==0:
            p = self.probability(alphabet=background[0])[0].T
            return np.array([jensenshannon(i,background[1]) for i in p])

    ####################################################################################
    #####  visualizations                                                          #####
    ####################################################################################

    def logo(self):
        return SequenceLogo(self)

    def show_logo(self, **kwargs):
        return SequenceLogo(self).plot(**kwargs)
    







