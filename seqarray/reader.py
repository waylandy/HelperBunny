import sys
from itertools import groupby
import pandas as pd

"""
    XMA format notes:
     An attempt at a generalized format of CMA, MMA, SMA, etc.
"""

class FA_Iterator:
    """
    An iterator to read FA files
    returns (name, seq)
    """
    def __init__(self, fafile, parsetags=False):
        self._fafile    = fafile
        self._tags      = parsetags
        self.isname     = lambda l: l[0]=='>'

    def __iter__(self):
        self._r = open(self._fafile)
        self._i = groupby(self._r, self.isname)
        self._s = False
        return self

    def _parse_tags(self, seq):
        seq = seq.split()
        chunk, k = [seq[0]], True
        for n, (k, group) in enumerate(groupby(seq[1:], lambda x: '=' in x or ':' in x)):
            if k == True:
                for g in group:
                    chunk += [g]
            else:
                if n != 0:
                    chunk[-1] += ' '+' '.join(list(group))
                else:
                    chunk += [' '.join(list(group))]
        return chunk

    def __next__(self):
        while True:
            head, group = next(self._i)
            if not head and self._s:
                self._s = False
            elif head:
                line = next(group)
            else:
                name = self._parse_tags(line.rstrip()[1:]) if self._tags else line.rstrip()[1:]
                seq  = ''.join(l.rstrip() for l in group)
                return name, seq

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

"""
    Other useful stuff
"""

class CDD:
    def __init__(self, cdd_file):
        self.cdd_file = cdd_file
        self._build()
        return

    def __getitem__(self, key):
        return self.df.xs(key, level='Query')

    def _build(self):
        it, i = iter(open(self.cdd_file)), -1
        while not next(it).startswith('Query'): i+=1
        self.df = pd.read_csv(self.cdd_file, sep='\t', skiprows=i)
        self.df['mid']   =  (self.df['From'].values + self.df['To'].values)/2
        index, query     = list(zip(*[i.split('#')[1].split(' - ') for i in self.df['Query']]))
        self.df['Query'] = query
        self.df['Index'] = list(map(int, index))
        self.df.sort_values(by=['Index', 'mid'], ascending=[1, 1], inplace=True)
        self.df.set_index(['Query', 'Accession'], inplace=True)
        self.df = self.df[['From', 'To', 'Short name', 'Superfamily', 'Hit type', 'PSSM-ID', 'E-Value', 'Incomplete']]
