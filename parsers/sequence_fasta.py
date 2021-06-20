from itertools import groupby

def read_fasta(file):
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

