import sys
import numpy as np

from .alignmentarray     import AlignmentArray     as alnarray2d
from .alignmentarrayplus import AlignmentArrayPlus as alnarray3d

from ..read.cfa import cfaReader
from ..read.fasta import fastaReader

def from_cfa(infile, v=1):
    isins = lambda x: x.islower()
    align = lambda x, y: [''.join([next(y) for _ in range(len(z))]) for z in x]
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
    def all2array(inseqs):
        out, *fts = inseqs
        out = [aln2array(out)]
        for _ in fts:
            out += [align(out[0], iter(_))]
        return out

    name, data = [], []
    for i, (n, d) in enumerate(cfaReader(infile)):
        name  += [n]
        data  += [all2array(d)]
        if i%200==0 and v>0:
            sys.stderr.write('Importing Row           : %s\r' % (1+i))
    if v>0:
        sys.stderr.write(    'Importing Row           : %s\n' % (1+i))
    return from_list(name, data)

def from_fasta(infile, v=1):
    name, data = [], []
    for i, (n, d) in enumerate(fastaReader(infile)):
        name  += [n]
        data  += [list(d)]
        if i%200==0 and v>0:
            sys.stderr.write('Importing Row           : %s\r' % (1+i))
    if v>0:
        sys.stderr.write(    'Importing Row           : %s\n' % (1+i))
    return from_list(name, data)

def from_list(name, data):
    name = np.array(name, dtype=object)
    data = np.array(data, dtype=object)
    return from_ndarray(name, data)

def from_ndarray(name, data):
    if len(data.shape)==3:
        data = data if len(data.shape)!=3 else data.transpose(0,2,1)
        if data.shape[2]==1:
            data = data[:,:,0]
        else:
            return finalize(name, data)
    if len(data.shape)==2:
        return finalize(name, data)
    else:
        raise AssertionError('Could not create alignment array, unexpected format received')
        
def finalize(name, data, optimize=True):
    dim = len(data.shape)
    if optimize and 2==dim:
        # optimization: delete dead insertion columns
        ispos  = lambda x: False if len(x)!=1 else not x.islower()
        isdead = lambda x: {''}==set(x)
        ins    = np.array([n for n, b in enumerate(map(ispos, data[0])) if not b])
        if len(ins) !=0:
            dead   = [n for n, c in zip(ins, data.T[ins]) if isdead(c)]
            data   = np.delete(data, dead, axis=1)
    
    if dim==2:
        return name, alnarray2d(data)
    if dim==3:
        return name, alnarray3d(data)

def AlignmentArray(input, format=None):
    if type(input)==str:
        fmt    = lambda x: x.split('.')[-1].lower()
        format = fmt(input) if format == None else format
    if format in ['cfa']:
        return from_cfa(input)
    if format in ['fasta','fa']:
        return from_fasta(input)
    else:
        raise Exception('Format not recognized')

