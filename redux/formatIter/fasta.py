from itertools import groupby

class fastaReader:
    """
    An iterator to read FA files
    returns (name, seq)
    """
    def __init__(self, infile):
        self.infile = infile
    
    def __iter__(self):
        return self.walk(fx_name=lambda x: x[1:])
        
    def walk(self, fx_name=None, fx_seq=None):
        fx_name = fx_name if fx_name != None else lambda x: x
        fx_seq  = fx_seq  if fx_seq  != None else lambda x: x
        isseq   = lambda x: x[1].startswith('>')
        names   = lambda x: [(i[0], i[1].strip()) for i in x]
        flip    = lambda x: list(zip(*list(x)))
        aggseq  = lambda x: fx_seq(''.join(i.strip() for i in flip(x)[1]))
        first   = True
        for k, g in groupby(enumerate(open(self.infile)), isseq):
            if k:
                first=False
                g = names(g)
                for n, i in names(g)[:-1]:
                    yield fx_name(i), fx_seq('')
                now = names(g)[-1][1]
            else:
                if first:
                    continue
                yield fx_name(now), aggseq(g)
