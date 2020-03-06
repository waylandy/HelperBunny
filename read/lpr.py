import sys
import re
from itertools import groupby
import numpy as np
import pandas as pd

"""
    BPPS output parsers:
     these objects help parse BPPS outputs
     these are not guaranteed to work as I cant find documented standards
"""

class LPR:
    """
        Reads LPR outputs from the BPPS programs
        Attributes:
            data  : dataframe representation of bpps
        Examples:
            lpr        = LPR(lpr_file)
            cfa        = CFA_Array(cfa_file)
            score, cat = lpr.fit(cfa.ar[:,1::2])
            lpr[1]
    """
    def __init__(self, lprfile, verbose=1):
        self.v       = verbose
        self.lprfile = lprfile
        self._populate()

    def __getitem__(self, key):
        if type(key)==int:
            return self.data.xs(key, level='category')
        elif len(key)==2:
            return self.data.xs(key, level=('category','rank'))

    def _populate(self):
        def iter_pats(block):
            block = ''.join(list(block))
            ispat = lambda x: bool(re.search(".*\(.*\)\:.*\%\).*\%\)", x))
            nopar = lambda x: re.sub(r'\([^)]*\)', '', x)
            parse = lambda x: [''.join(i[1]) for i in groupby(x, lambda x: x.isalpha())]
            for pat in filter(ispat, block.split('\n')):
                pat, info = [x.strip() for x in nopar(pat).split(':')]
                info, (patres, patpos) = info.split(), parse(pat)
                yield [patres, patpos] + info

        def iter_lpr():
            breakpad = '::::::::::'
            catbreak = lambda x: x.startswith(breakpad) and x.rstrip().endswith(breakpad)
            catnum   = lambda x: int(x.split(' ')[3][:-1])
            with open(self.lprfile) as r:
                for key, group in groupby(r, catbreak):
                    if key:
                        cat = catnum(next(group))
                    else:
                        yield cat, group
        multidex, data = [], []
        for key, group in iter_lpr():
            for rank, pat in enumerate(iter_pats(group)):
                multidex.append((key, rank+1))
                data.append(pat)
        columns = ['residue','position','fg-match','fg-diverge','bg-match','bg-diverge','nats','binaryurn','%info','wtnumseq']
        self.data = pd.DataFrame(data, index=list(map(list, list(zip(*multidex)))), columns=columns)
        for c in self.data.columns:
            self.data[c] = pd.to_numeric(self.data[c], errors='ignore')
        self.data.index.names = 'category', 'rank'
        td = dict((cat, self.data.loc[cat]['nats'].values.sum()) for cat in set(i[0] for i in self.data.index.values))
        t  = np.array([td[i[0]] for i in self.data.index.values])
        self.data['%info'] = (self.data['nats'].values/t).astype(np.float16)
        self.data = self.data.sort_values(['category','rank'])

    def fit(self, ar):
        cats = sorted(list(set(i[0] for i in self.data.index.values)))
        try:
            self._scord
        except:
            self._scord = dict((cat, self.data.loc[cat][['residue', 'position', '%info']].values) for cat in cats)
        if len(ar.shape)==1:
            ar = np.array([ar])
        # assume indices = position - 1
        scoreone = lambda c, x: sum(n for r, i, n in self._scord[c] if x[i-1] in r)
        scoreall = lambda x: [scoreone(i, x) for i in cats]
        x = np.zeros((ar.shape[0], len(cats)), dtype=np.float16)
        for n, s in enumerate(ar):
            x[n] = scoreall(s)
            if n%200==0 and self.v>0:
                sys.stderr.write('[ LPR       ] Calculating Likelihoods : %s\r' % n)
        if self.v>0:
            sys.stderr.write(    '[ LPR       ] Calculating Likelihoods : %s\n' % n)
        return x.T, np.array(cats)

    def vectorize(self, group, f='%info', aln=None):
        aln     = aln if aln!=None else self.data['position'].values.max()
        x , y   = 1 + np.arange(aln), np.zeros(aln)
        fx, fy  = self[group][['position',f]].values.T
        fx      = fx.astype(int)
        y[fx-1] = fy
        return x, y
