import sys
from itertools import groupby

import numpy as np
from scipy.stats import entropy
from scipy.spatial.distance import jensenshannon

from .utils import (
    iter_fasta,
    partition_a2m_string,
    )
from .logo import (
    SequenceLogo,
    plot_interactive,
    )

AMINO       = 'ARNDCQEGHILKMFPSTWYV-'

BLOSUM62_BG = 'ARNDCQEGHILKMFPSTWYV', np.array([
    #  A      R      N      D      C      Q      E      G      H      I
     0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062,
    #  L      K      M      F      P      S      T      W      Y      V
     0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072,
    ])

class AlignmentArray(np.ndarray):
    
    def __new__(cls, ndarray, names=None):
        obj = np.asarray(ndarray, dtype=object).view(cls)
        return obj

    def __array_finalize__(self, obj):
        if obj is None: 
            return
        self.info = getattr(obj, 'info', None)
        
    def __array_wrap__(self, out_arr, context=None):
        return out_arr if (out_arr.dtype == object) else out_arr.view(np.ndarray)
    
    #################################################################################
    
    @staticmethod
    def from_a2m_file(a2m_file, flanking=True, verbose=True):
        n_seq = sum(1 for i in iter_fasta(a2m_file))
        n_pos = len(partition_a2m_string(next(iter_fasta(a2m_file))[1]))
        left, right = (0, n_pos) if flanking else (1, n_pos - 1)
        aligned = np.full((n_seq, n_pos if flanking else n_pos - 2), '', dtype=object)
        headers = np.full(n_seq, '', dtype=object)
    
        for n, (h, s) in enumerate(iter_fasta(a2m_file)):
            aligned[n] = partition_a2m_string(s)[left:right]
            headers[n] = h
    
            if verbose:
                n += 1
                if n % 1000 == 0:
                    sys.stderr.write(f'{n} / {n_seq} ({n_pos})\r')

        aligned = AlignmentArray(aligned)
        positions = aligned.is_position()
        if verbose:
            sys.stderr.write(f'{n} / {n_seq} ({positions.sum()} + {(~positions).sum()})\r')
        
        aligned = aligned[:, ~(aligned == '').all(0)]
        positions = aligned.is_position()
        if verbose:
            sys.stderr.write(f'{n} / {n_seq} ({positions.sum()} + {(~positions).sum()})\n')
            
        return headers, aligned
    
    @staticmethod
    def from_a2m_list(a2m_list, flanking=True, verbose=True):
        n_seq = len(a2m_list)
        n_pos = len(partition_a2m_string(a2m_list[0]))
        left, right = (0, n_pos) if flanking else (1, n_pos - 1)
        aligned = np.full((n_seq, n_pos if flanking else n_pos - 2), '', dtype=object)
    
        for n, s in enumerate(a2m_list):
            aligned[n] = partition_a2m_string(s)[left:right]
    
            if verbose:
                n += 1
                if n % 1000 == 0:
                    sys.stderr.write(f'{n} / {n_seq} ({n_pos})\r')
        
        aligned = AlignmentArray(aligned)
        positions = aligned.is_position()
        if verbose:
            sys.stderr.write(f'{n} / {n_seq} ({positions.sum()} + {(~positions).sum()})\r')
        
        aligned = aligned[:, ~(aligned == '').all(0)]
        positions = aligned.is_position()
        if verbose:
            sys.stderr.write(f'{n} / {n_seq} ({positions.sum()} + {(~positions).sum()})\n')
        
        return aligned
    
    #################################################################################

    def is_position(self):
        return np.array([i.isupper() or i == '-' for i in self[0]])

    def remove_inserts(self):
        return self[:,self.is_position()]
    
    def get_positions(self):
        is_position = self.is_position()
        positions = np.cumsum(is_position) + (~is_position * 0.5)
        return np.array([i if i % 1 == 0.5 else int(i) for i in positions], dtype=object)

    def get_residue_positions(self, *args):
        mapping = {j: i  for i, j in enumerate(self.get_positions()) if j % 1 != 0.5}
        return self[:, [mapping[i] for i in args]]

    def get_residue_range(self, start, end):
        mapping = {j: i  for i, j in enumerate(self.get_positions()) if j % 1 != 0.5}
        if start == end == None:
            raise Exception()
        elif start == None:
            return self[:, :1 + mapping[end]]
        elif end == None:
            return self[:, mapping[start]:] 
        else:
            return self[:, mapping[start]: 1 + mapping[end]]
    
    #################################################################################

    def trim_alignment(self, occupancy=0.5, **kwargs):
        blueprint = np.array([i.islower() or i == '' for i in self[0]])
        blueprint += (self == '-').mean(0) > occupancy
        return self.squeeze_alignment(blueprint, **kwargs)
    
    def squeeze_alignment(self, blueprint, **kwargs):
        squeeze = lambda x: ''.join(i.lower().replace('-','') if j else i for i, j in zip(x, blueprint))
        return AlignmentArray.from_a2m_list(list(map(squeeze, self)), **kwargs)

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
    
    #################################################################################

    def show(self, **kwargs):
        return plot_interactive(self, **kwargs)
