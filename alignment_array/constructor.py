import sys
import numpy as np

from .sequenceheaders import SequenceHeaders

from .alignmentarray import AlignmentArray
from .alignmentarrayfeaturized import AlignmentArrayFeaturized

#from ..sequence_parser.cfa import read_cfa
from ..sequence_parser.a2m import read_a2m
from ..sequence_parser.fasta import read_fasta


def vectorize_aligned_sequence(seq, flanking=True):
    is_insert  = lambda x: x.islower()
    now_insert = True
    was_insert = True
    A          = ['']
    for x in seq:
        now_insert = is_insert(x)
        if now_insert:
            if was_insert:
                A[-1] += x
            else:
                A += [x]
        else:
            if was_insert:
                A += [x]
            else:
                A += ['',x]
        was_insert = now_insert
    if not is_insert(A[-1][0]):
        A += ['']
    return A if flanking else A[1:-1]

def a2m_to_lists(file, flanking=True, quiet=False):
    report_row = lambda x: sys.stderr.write('Importing Row : %s\r' % (1+n))
    names, aln = [], []
    for n, (name, seq) in enumerate(read_a2m(file)):
        names += [name]
        aln   += [vectorize_aligned_sequence(seq, flanking=flanking)]
        if not quiet and n%200==0: report_row(n)
    report_row(n)
    return names, aln

def fasta_to_lists(file, quiet=False):
    report_row = lambda x: sys.stderr.write('Importing Row : %s\r' % (1+n))
    name, data = [], []
    for i, (n, d) in enumerate(read_fasta(file)):
        name  += [n]
        data  += [list(d)]
        if not quiet and n%200==0: report_row(n)
    report_row(n)
    return names, aln

def constructor(*args, format=None, remove_unused_cols=True, v=1):
    is_file  = lambda x: len(x)==1 and type(x[0])==str
    is_list  = lambda x: len(x) in [1,2] and all(type(i)==list for i in x)
    get_ext  = lambda x: x.split('.')[-1].lower()
    
    if is_file(args):
        file = next(iter(args))
        ext  = get_ext(file)
        if   ext in ['a2m']:
            name, data = a2m_to_lists(file)
        elif ext in ['fasta','fa']:
            name, data = fasta_to_lists(file)
        else:
            raise Exception('unrecognized file extension: use a2m or fasta')
    elif is_list(args):
        if len(args)==1:
            name, data = (list(range(len(args[0])+1))[1:], args[0])
        elif len(args)==2:
            assert 1==len(set(map(len,args)))
            name, data = args
        else:
            raise Exception('this shouldnt be possible')
    else:
        raise Exception('unrecognized case')
    
    name = np.array(name, dtype=object)
    data = np.array(data, dtype=object)
    
    
    if data.ndim == 1:
        raise Exception('unexpected shape. is input aligned?')
    elif data.ndim == 2:
        data = data[:,~(data=='').all(0)]
        return SequenceHeaders(name), AlignmentArray(data)
    elif data.ndim == 3:
        raise Exception('featurized alignments in development')
    else:
        raise Exception('unexpected shape.')


