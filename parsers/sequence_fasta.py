import io
import os
from itertools import groupby

def _read_fasta(file):
    """ 
    iterable  : works for reading fasta, fa, a3m, a2m
    yields    : (header, sequence)
    behaviors : ignores headers without sequences, will read blank headers
    """
    is_header = lambda x: x.startswith('>')
    compress  = lambda x: ''.join(_.strip() for _ in x)
    reader    = iter(groupby(open(file), is_header))
    reader    = iter(groupby(open(file), is_header)) if next(reader)[0] else reader
    for key, group in reader:
        if key:
            for header in group:
                header = header[1:].strip()
        else:
            sequence = compress(group)
            if sequence != '':
                yield header, sequence

def read_fasta(arg):
    if os.path.exists(arg): 
        # if arg is an alignment file
        handle = open(arg)
    elif type(arg)==str: 
        # assume that arg the contents of a file
        handle = io.StringIO()
        handle.write(arg)
        handle.seek(0)
    else:
        raise Exception('Invalid input')

    is_header = lambda x: x.startswith('>')
    compress  = lambda x: ''.join(_.strip() for _ in x)
    for n, (key, group) in enumerate(groupby(handle, is_header)):
        if key:
            for header in group:
                header = header[1:].strip()
        elif n==0:
            continue
        else:
            sequence = compress(group)
            if sequence != '':
                yield header, sequence
