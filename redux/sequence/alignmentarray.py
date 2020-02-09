import sys
import numpy as np

from ..formatIter.cfa import cfaReader

class AlignmentArray(np.ndarray):
    """
        Representation of CFA sequences as arrays
        Should also be capable of reading feature-linked CFA (not heavily tested)
        Should also be capable of reading FA (bad optimization)
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


"""

    Optimized constructor for the AlignmentArray object
       The AlignmentArray object itself was not optimized for the code below
       I just copied and pasted it from the previous version

"""

isins = lambda x: x.islower()
def aln2array(inseq):
    x, preins = [''], True
    for now in inseq:
        nowins = isins(now)
        if nowins:
            if preins:
                x[-1] += now
            else:
                x += [now]
        else:
            if preins:
                x += [now]
            else:
                x += ['',now]
        preins = nowins
    return x if isins(x[-1][0]) else x+['']

align = lambda x, y: [''.join([next(y) for _ in range(len(z))]) for z in x]
def all2array(inseqs):
    out, *fts = inseqs
    out = [aln2array(out)]
    for _ in fts:
        out += [align(out[0], iter(_))]
    return out
    
def cfa2array(infile, v=1):
    name, data = [], []
    for i, (n, d) in enumerate(cfaReader(infile)):
        name  += [n]
        data  += [all2array(d)]
        if i%200==0 and v>0:
            sys.stderr.write('[ CFA_Array ] Importing Row           : %s\r' % (1+i))
    if v>0:
        sys.stderr.write(    '[ CFA_Array ] Importing Row           : %s\n' % (1+i))
    ar = AlignmentArray(np.array(data, dtype=object), names=np.array(name, dtype=object))
    return ar if len(ar.shape)!=3 else ar.transpose(0,2,1)




