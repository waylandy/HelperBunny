import sys
import random
from itertools import groupby

import numpy as np

def iter_fasta(fasta_file):
    is_header = lambda x: x.startswith('>')
    compress  = lambda x: ''.join(_.strip() for _ in x)
    reader    = iter(groupby(open(fasta_file), is_header))
    reader    = iter(groupby(open(fasta_file), is_header)) if next(reader)[0] else reader
    for key, group in reader:
        if key:
            for header in group:
                header = header[1:].strip()
        else:
            sequence = compress(group)
            if sequence != '':
                yield header, sequence

def read_fasta(fasta_file):
    headers, sequences = zip(*iter_fasta(fasta_file))
    return np.array(headers, dtype=object), np.array(sequences, dtype=object)

#################################################################################

class A2M_Map:
    def __init__(self, a2m_string):
        self.a2m_string = a2m_string
        pos, gap, ins = np.array([[1,0,0] if i.isupper() else [0,1,0] if i == '-' else [0,0,1] for i in self.a2m_string], dtype=bool).T
        self.aligned_pos = np.cumsum(gap + pos) + ins * 0.5
        self.residue_pos = np.cumsum(ins + pos) + gap * 0.5
        self.aligned_pos_to_index = {int(j): i for i, j in enumerate(self.aligned_pos) if j % 1 != 0.5}

    def get_aligned_residue_positions(self, *args):
        return [self.a2m_string[self.aligned_pos_to_index[i]] for i in args]
        
    def get_aligned_residue_range(self, start, end):
        return self.a2m_string[self.aligned_pos_to_index[start]: 1 + self.aligned_pos_to_index[end]]

def partition_a2m_string(a2m_string):
    current, previous = True, True
    partitions = ['']
    
    for i in a2m_string:
        current = i.islower()
        
        if current:
            if previous:
                partitions[-1] += i
            else:
                partitions += [i]
        else:
            if previous:
                partitions += [i]
            else:
                partitions += ['', i]
        
        previous = current
    if not previous:
        partitions += ['']
    return partitions

