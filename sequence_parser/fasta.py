from itertools import groupby

def read_fasta(file):
    """
    generator: yields entries in a fasta file, (identifier, sequence)
    """
    f1 = lambda x: x.startswith('>')
    f2 = lambda x: ''.join(i.strip() for i in x)
    it = iter(groupby(open(file), f1))
    it = iter(groupby(open(file), f1)) if next(it)[0] else it
    for k, g in it:
        if not k:
            s = f2(g)
            if s!='':
                yield n, s
        if k:
            for i in g:
                n = i.strip()[1:]
