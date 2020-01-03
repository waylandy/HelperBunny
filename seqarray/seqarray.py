import sys
from itertools import groupby
import numpy as np
import pandas as pd

"""
    CFA format notes:
     1) sequences provided are assumed to be the FULL SEQUENCE
     2) duplicate identifiers are allowed in case of multiple profile hits
     3) sequences must be provided as a single unbroken line after the identifier
          a) feature sequences can be provided below and must be the same length
          b) each sequence must have the same number of feature sequences
"""

class CFA_Iterator:
    """
        Iterator object for CFA sequences
        Should also be capable of reading feature-linked CFA (not tested)
        Upon iteration, returns [name, seq]
    """
    def __init__(self, cfa):
        self._cfafile = cfa
        self.isname = lambda l: l[0]=='>'

    def __iter__(self):
        self._r = open(self._cfafile)
        self._i = groupby(self._r, self.isname)
        self._s = False
        return self

    def __next__(self):
        while True:
            head, group = next(self._i)
            if not head and self._s:
                self._s = False
            elif head:
                line = next(group)
            else:
                group = [i.rstrip() for i in group if not i.isspace()]
                assert len(set([len(i) for i in group])) == 1
                return line.rstrip()[1:], group

class CFA(np.ndarray):
    """
        Representation of CFA sequences as arrays
        Should also be capable of reading feature-linked CFA (not tested)
        Attributes:
            ar    : sequence array (try subclassing np.ndarray object)
            names : names array
        Example:
            ca = CFA_Array(cfa_file)
            ca.ar
    """
    def __new__(cls, ndarray, names=None):
        obj = np.asarray(ndarray, dtype=object).view(cls)
        obj.names = names
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.names = getattr(obj, 'names', None)

    # def fmt(st, maxsize=12, fillchar='.', align='left'):
    #     stl     = len(st)
    #     if stl>maxsize:
    #         maxsize -=2
    #         digit = lambda x : len(str(x))
    #         ceil  = lambda x : int(np.ceil(x))
    #         floor = lambda x : int(np.floor(x))
    #
    #         mid = digit(stl - maxsize + digit(stl))
    #         topbot = maxsize-mid
    #         top, bot = ceil(topbot/2), floor(topbot/2)
    #         mids = st[top:-bot]
    #
    #         return '%s(%s)%s' % (st[:top], len(mids), st[-bot:])
    #     #, sum((len(st[:top]), digit(len(mids)), len(st[-bot:])))
    #     else:
    #         if align=='left':
    #             return st.ljust(maxsize, fillchar)
    #         if align=='right':
    #             return st.rjust(maxsize, fillchar)

    # def _update(self):
    #     ispos = lambda x: False if len(x)!=1 else not x.islower()

    def _frame(self, x):
        return self if len(self.shape)!=3 else self[:,:,x]

    def aln(self, f=0, fxn=None):
        ispos = lambda x: False if len(x)!=1 else not x.islower()
        fxn   = fxn if fxn!=None else ispos
        test  = self[0,:] if len(self.shape)!=3 else self[0,:,f]
        mask  = [fxn(x) for x in test]
        return self._frame(f)[:,mask]

    def out(self, f=0):
        self._frame(f)
        return 0
    # def out(self, f=0, aln=True):
    #     aln  = None if aln else lambda x: True
    #     bind = lambda x: [''.join(i) for i in x]
    #     return bind(self.aln(f=f, fxn=aln))

def from_cfa(cfa, v=1):
    ispos = lambda x: x.isupper() or x=='-'
    isins = lambda x: x.islower()

    def aln2ar(line):
        x, pre = [''], 'a'
        for now in line:
            if   ispos(pre) and ispos(now):
                x += ['',now]
            elif isins(pre) and isins(now):
                x[-1] += now
            elif isins(pre) and ispos(now):
                x += [now]
            elif ispos(pre) and isins(now):
                x += [now]
            else:
                raise Exception('Unrecognized case: %s%s'%(pre,now))
            pre=now
        return x + [''] if ispos(x[-1][0]) else x

    def align(ref, target):
        assert len(''.join(ref)) == len(target)
        it = iter(target)
        return [''.join([next(it) for i in range(len(b))]) for b in ref]

    def data2ar(data):
        ref = aln2ar(data[0])
        if len(data)==1:
            return ref
        else:
            return [ref] + [align(ref, i) for i in data[1:]]

    name, data = [], []

    for i, (n, d) in enumerate(CFA_Iterator(cfa)):
        name  += [n]
        data  += [data2ar(d)]
        if i%200==0 and v>0:
            sys.stderr.write('[ CFA_Array ] Importing Row           : %s\r' % (1+i))
    if v>0:
        sys.stderr.write(    '[ CFA_Array ] Importing Row           : %s\n' % (1+i))
    ar = CFA(np.array(data, dtype=object), names=np.array(name, dtype=object))
    ar = ar if len(ar.shape)!=3 else ar.transpose(0,2,1)
    return ar

    
