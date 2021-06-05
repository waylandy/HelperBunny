import sys
import numpy as np

from ..parsers.sequence_fasta import read_fasta

from .headerarray import SequenceHeaders
from .alignmentarray import AlignmentArray

def gen_array(infile, **kwargs):
    
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
    
    def fast_check(infile, max_iter=420):
        get_vlen = lambda x: len(vectorize_aligned_sequence(x))
        reader   = read_fasta(infile)
        vlen     = get_vlen(next(reader)[1])
        for n, (h, s) in enumerate(read_fasta(infile)):
            if n > max_iter:
                break
            if vlen != get_vlen(s):
                raise Exception(f'Input file "{infile}" does not seem properly aligned')
    
    def parse_file(infile, flanking=True, remove_unused_cols=True, max_seqs=-1, quiet=False):
        printrow   = lambda x: sys.stderr.write('Importing Row : %s\r' % (1+n))
        reporter   = lambda x: None if quiet else printrow(x)
        exceeded   = lambda x: x == max_seqs - 1
        reader     = read_fasta(infile)
        names, aln = [], []
        
        for n, (name, seq) in enumerate(reader):
            names += [name]
            aln   += [vectorize_aligned_sequence(seq, flanking=flanking)]
            if n % 1000 == 999:
                reporter(n)
            if exceeded(n):
                break
        reporter(n)
        
        names = np.array(names, dtype=object)
        
        if remove_unused_cols:
            is_empty = lambda x: all(i=='' for i in x)
            index    = {i for i, j in enumerate(map(is_empty,zip(*aln))) if not j}
            aln_T    = [j for i, j in enumerate(zip(*aln)) if i in index]
            aln      = np.array(aln_T, dtype=object).T
        else:
            aln      = np.array(aln, dtype=object)
        
        return names, aln
    
    fast_check(infile)
    names, aln = parse_file(infile, **kwargs)
    return SequenceHeaders(names), AlignmentArray(aln)
