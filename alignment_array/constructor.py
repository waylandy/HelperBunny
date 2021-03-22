import sys
import numpy as np

from .sequenceheaders import SequenceHeaders

from .alignmentarray import AlignmentArray
from .alignmentarrayfeaturized import AlignmentArrayFeaturized

from ..sequence_parser.cfa import read_cfa
from ..sequence_parser.fasta import read_fasta


def from_cfa(file, v=1):
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
    for i, (n, d) in enumerate(read_cfa(file)):
        name  += [n]
        data  += [all2array(d)]
        if i%200==0 and v>0:
            sys.stderr.write('Importing Row           : %s\r' % (1+i))
    if v>0:
        sys.stderr.write(    'Importing Row           : %s\n' % (1+i))
    return from_list(name, data)

def from_fasta(file, v=1):
    name, data = [], []
    for i, (n, d) in enumerate(read_fasta(file)):
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
            return name, data
    if len(data.shape)==2:
        return name, data
    else:
        raise AssertionError('Could not create alignment array, unexpected format received')
        
def optimize(name, data, remove_unused_cols=True):
    dim = len(data.shape)
    if remove_unused_cols and 2==dim:
        # optimization: delete dead insertion columns
        ispos  = lambda x: False if len(x)!=1 else not x.islower()
        isdead = lambda x: {''}==set(x)
        ins    = np.array([n for n, b in enumerate(map(ispos, data[0])) if not b])
        if len(ins) !=0:
            dead   = [n for n, c in zip(ins, data.T[ins]) if isdead(c)]
            data   = np.delete(data, dead, axis=1)
    
    if dim==2:
        return SequenceHeaders(name), AlignmentArray(data)
    if dim==3:
        return SequenceHeaders(name), AlignmentArrayFeaturized(data)

def constructor(*args, format=None, remove_unused_cols=True, v=1):
    """
    USAGE 1: FROM FILES
        name, aln = AlignmentArray('file.fasta', remove_unused_cols=False)

    USAGE 2: FROM LISTS 
        name   = ['a','b','c']
        aln    = [
                    ['A','','A','A'],
                    ['A','','A','A'],
                    ['A','','A','A']
                 ]
        name, aln = AlignmentArray(name, aln, remove_unused_cols=False)
    """
    # first check if the input is a filename
    if len(args)==1 and type(args[0])==str:
        input   = args[0]
        fmt     = lambda x: x.split('.')[-1].lower()
        format  = fmt(input) if format == None else format
    if format in ['cfa','a2m']:
        return optimize(*from_cfa(input, v=v), remove_unused_cols=remove_unused_cols)
    if format in ['fasta','fa']:
        return optimize(*from_fasta(input, v=v), remove_unused_cols=remove_unused_cols)

    # now check if the input is 2 lists: name and alignment list
    # new plan : just try to force it into alignment array
    try:
        assert len(args)==2
        args = tuple(arg.tolist() if 'tolist' in dir(arg) else arg for arg in args)
        assert all(type(i) in (list,tuple) for i in args)
        return optimize(*from_list(*args), remove_unused_cols=remove_unused_cols)
    except:
        pass

    # at this point, give up on trying to parse it
    else:
        raise Exception('Failed to create array representation.')


