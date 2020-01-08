import sys
from itertools import groupby

"""
    XMA format notes:
     An attempt at a reading a generalized format of CMA, MMA, SMA, etc.
"""

class XMA_Block:
    """
        An iterator to alignment blocks from XMA
        Upon iteration, returns (num, name, aln)
            num   : string; sequence numbers, probably useless
            name  : size 2 tuple; sequence name and random XMA-related tags
            aln   : size 3 tuple; fragments of the alignment
    """
    def __init__(self, name, it):
        self.name, self.it = name, iter(it)
        return

    def __iter__(self):
        return self

    def _excise(self, s, border=')('):
        b1, b2 = border
        b1s, b2s = True, True
        l, p = ['','',''], 0
        for e in s:
            if e==b1 and b1s:
                b1s = False
                p+=1
                continue
            elif e==b2 and b2s and not b1s:
                b2s = False
                p+=1
                continue
            l[p]+=e
        l = [''.join([y for y in x if y not in '{\n})(*']) for x in l]
        l[0], l[2] = l[0].lower(), l[2].lower()
        return l

    def _restore(self, name):
        n = name.find(')|}')
        if n==-1:
            return name, ''
        adi = name[name[:n+3].find('{|'):n+3]
        return name.replace(adi, '').strip(), adi

    def __next__(self):
        while True:
            l = next(self.it)
            if l[0]=='$':
                num = l[:-2]
            elif l[0]=='>':
                name = self._restore(l[1:].strip())
            elif '{(' in l:
                aln = self._excise(l)
                return num, name, aln

class XMA:
    """
        An iterator object to read CMA, MMA, etc...
        Upon iteration, returns XMA_Block
        Example:
            XMA('file.mma').to_cfa('file.cfa')
    """
    def __init__(self, xmafile):
        self.xmafile = xmafile

    def __iter__(self):
        self.it = groupby(open(self.xmafile), lambda x: x.startswith('[0_'))
        return self

    def __next__(self):
        name = lambda x: x.split('=')[1].split('(')[0]
        while True:
            k, g = next(self.it)
            k, g = (name(next(g)), next(self.it)[1]) if k else (k, g)
            return XMA_Block(k, g)

    def to_cfa(self, *args):
        h = sys.stdout if len(args)==0 else open(args[0], 'w')
        alnl = lambda x: len(x)-sum(1 for i in x if i.islower())
        lens = set(alnl(next(iter(x))[2][1]) for x in self)
        assert len(lens)==1
        for i in self:
            for j in i:
                h.write('>%s XMA:%s\n%s\n\n' % (j[1][0], i.name, ''.join(j[2])))
