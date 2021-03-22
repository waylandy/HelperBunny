import sys
from itertools import groupby

import numpy as np
from scipy.stats import entropy
from scipy.spatial.distance import jensenshannon

from ..widgets import *

alphabet = 'ARNDCQEGHILKMFPSTWYV-'

background = ('ARNDCQEGHILKMFPSTWYV', np.array( # BLOSUM62 DISTRIBUTION
            [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062,
             0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]))

class AlignmentArray(np.ndarray):
    """
        uses references as array element dtype
    """

    ####################################################################################
    #####  inherit np.ndarray                                                      #####
    ####################################################################################
    
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
 
        if (out_arr.dtype == object):
            return out_arr
        else:
            return out_arr.view(np.ndarray)

        # documentation claims i should return:
        # super(AlignmentArray, self).__array_wrap__(self, out_arr, context)
        # this breaks code after second run (need more testing)
    
    ####################################################################################
    #####  alignment properties                                                    #####
    ####################################################################################

    def is_position(self): # ispos
        # boolean array of whether each column is a position in the alignment
        # makes judgement using ONLY the first sequence
        ispos = lambda x: False if len(x)!=1 else not x.islower()
        return np.array(list(map(ispos, self[0])))

    def position_number(self): # cols
        # gives the position number for each column in the array
        # "positions" will be int, "insertions" will be float
        pos = self.is_position().astype(object)
        pos[pos==True] = np.arange(1,1+pos[pos==True].shape[0])
        i = 0
        for n, v in enumerate(pos):
            if type(v)==bool:
                pos[n]=i+0.5
            else:
                i=v
        return pos

    def positions(self, *args): # CHECK IF THIS IS BEING USED ANYWHERE ELSE
        # returns alignment positions (NOT index) if given
        # if no positions given, returns all positions
        if len(args)==0:
            return self[:,self.is_position()]
        else:
            obj = set(map(type,args))
            assert 1==len(obj) and obj.pop()==int
            mask = [i in args for i in self.position_number()]
            return self[:,mask]

    ####################################################################################
    #####  alignment editting                                                      #####
    ####################################################################################

    def make_insertion(self, start, end):
        # defines the start to end indices will be redefined as inserts
        # index locations start & end are INCLUSIVE
        # if you define the start of an insertion next to an existing insertion, this will automatically extend the boundary
        # this makes sure 2 insertions will never be next to each other
        try:
            del self.position_array
        except:
            pass
        assert len(self.shape)==2
        reform    = lambda x: ''.join(i.lower() for i in x if i!='-')
        compress  = lambda x: np.array([reform(i) for i in x], dtype=object)
        hasi      = frozenset(n for n, i in enumerate(self.is_position()) if not i)
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
    
    ####################################################################################
    #####  position arrays for statistics                                          #####
    ####################################################################################
    #####
    #####  always update position_array if __getitem__ or __setitem__ is invoked
    #####  currently no strict enforcing to rebuild position_array if invoked
    #####
    
    def get_position_array(self, v=1):
        # remove insertions from the alignment array and return a PositionArray
        mask = self.is_position()
        tran = (mask==True).sum(), self.shape[1]
        if v>0:
            sys.stderr.write('PositionArray           : Found %s POSITIONS from %s COLUMNS; removing %s INSERTIONS\n'% (*tran, tran[1]-tran[0]))
        self.position_array = np.array(self[:,self.is_position()], dtype='|U1')
        return self.position_array
    
    def check_position_array(self):
        # check the position array
        # can still fail if you use __setitem__ at any point
        try:
            if shape.shape!=self.position_array.shape:
                self.get_position_array()
        except:
            self.get_position_array()
    
    def gaps(self, gap='-'):
        # CORE STAT FUNCTION
        # returns %gaps at each alignment position
        self.check_position_array()
        return (self.position_array==gap).mean(axis=0)

    def frequency(self, alphabet=alphabet, pseudocount=0):
        # CORE STAT FUNCTION
        # returns position frequency matrix as int AND the alphabet used
        self.check_position_array()
        return pseudocount+np.array([(i==self.position_array).sum(axis=0) for i in alphabet]), alphabet

    def probability(self, alphabet=alphabet, pseudocount=0):
        # DERIVED CALCULATION: frequency
        # returns position probability matrix as float AND the alphabet used
        ans, _ = self.frequency(alphabet, pseudocount=pseudocount)
        return ans/ans.sum(0), alphabet

    def weight(self, background=background, pseudocount=0):
        # DERIVED CALCULATION: probability
        # returns position weight matrix normalized against background as float & the alphabet used
        ans, _ = self.probability(alphabet=background[0], pseudocount=pseudocount)
        return (ans.T/background[1]).T, background[0]

    def kldivergence(self, *arg, background=background):
        # DERIVED CALCULATION: probability
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
        # DERIVED CALCULATION: probability
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
        # DERIVED CALCULATION: weight & probability
        cons  = lambda x: ''.join(background[0][i] for i in x.argmax(0))
        if norm:
            return cons(self.weight(background=background)[0])
        else:
            return cons(self.probability(alphabet=background[0])[0])
        
    ####################################################################################
    #####  writing files and exporting                                             #####
    ####################################################################################

    def sequences(self, *arg):
        # used to write files
        n = list(map(str,range(self.shape[0]))) if len(arg)==0 else arg[0]
        assert len(n)==self.shape[0]
        for n, s in zip(n, self):
            yield n, s
   
    def generate_fasta(self, *args, aligned=True):
        for i, j in self.sequences(*args):
            j = ''.join(k for k in ''.join(j) if not k.islower()) if aligned else ''.join(j).upper().replace('-','')
            k = '>%s\n%s\n\n' % (i, j)
            yield k

    def generate_cfa(self, *args):
        for i, j in self.sequences(*args):
            k = '>%s\n%s\n\n' % (i, ''.join(j))
            yield k

    ####################################################################################
    #####  visualizations                                                          #####
    ####################################################################################

    def show_logo(self, **kwargs):
        return SequenceLogoScroll(self, **kwargs).show()







