import sys
import re
from itertools import groupby

import numpy as np
import pandas as pd

def read_lpr(file):
    """
    generator: yields each category in an lpr file
    """
    col = ['residue','position','fg-match','fg-diverge','bg-match','bg-diverge','nats','binaryurn','%info','wtnumseq']
    f = lambda x: x.endswith('::::::::::\n')
    q = lambda x: re.sub(r'\([^)]*\)', '', x)
    s = lambda x: re.split('(\d+)',x[0])[:2]+x[1:]
    p = re.compile(r".*\(.*\)\:.*\%\).*\%\)")
    i = iter(groupby(open(file), f))
    i = iter(groupby(open(file), f)) if next(i)[0] else i
    for k, g in i:
        n = next(g).strip('\n :') if k else n
        if not k:
            g = map(q,filter(lambda x: p.match(x), g))
            g = [s(j.replace(':','').split()) for j in g]
            g = pd.DataFrame(g, columns=col)
            yield int(n.split()[-1]), g

def read_con(file):
    """
    generator: yield the constraints and weights
    """
    b, col = '', ['residue','position','weight']
    for l in open(file):
        if '}' in l and b!='':
            n, b = (b+l).split('=')
            b = b[1+b.index('{'):b.index('}')]
            b = [i.strip().split('\t') for i in b.strip().split('\n')]
            yield n.strip(), pd.DataFrame([i for i in b if len(i)==3], columns=col)
            b=''
        elif all(i in l for i in '={') or b!='':
            b+=l

def read_hpt(file):
    """
    parser: returns the names and topology of a hpt file
    """
    r = [i.split(' ') for i in filter(lambda x: x.startswith('+'), open(file))]
    r = [(i, int(j), k[:-2]) for i, j, k in [[i]+j.split('.') for i, j in r]]
    names, hpt = dict(i[-2:] for i in r), next(zip(*r))
    hpt = [[names[n+1] for n, _ in enumerate(i) if _=='+'] for i in hpt]
    return names, hpt

class Constraints:
    """
    represents a set of constraints defining one cluster
    """
    def __init__(self, name, residues, positions, weights):
        self.name      = name
        self.residues  = residues
        self.positions = positions
        self.weights   = weights
        self.norm      = weights/weights.sum()
        
    def __repr__(self):
        n = ('Set %s' % self.name) if type(self.name)==int else self.name
        return '<Constraints: "%s" containing %s constraints>' % (n, len(self.positions))
    
    def fit(self, x):
        # assume indices = position - 1
        x = x.positions()
        s = np.zeros(x.shape[0])
        for r, n, w in zip(self.residues, self.positions-1, self.norm):
            s += w*np.array([a==x[:,n] for a in r]).any(0)
        return s

class SequenceConstraints:
    """
    represents constraint sets defining many clusters
    """
    def __init__(self, file):
        extension        = file.split('.')[-1].lower()
        self.constraints = []
        
        if extension=='lpr':
            self._parse_lpr(file)
            
        elif extension=='con':
            self._parse_con(file)
        
        else:
            raise Exception('File extension not recognized. Use "lpr" or "con" format.')
    
    def __iter__ (self):
        return iter(self.constraints)
    
    def __len__(self):
        return sum(1 for i in self)
    
    def __repr__(self):
        r = [j for i in self for j in i.positions]
        return '<SequenceConstraints: container of %s constraint sets (%s,%s)>' % (len(self), min(r), max(r))
    
    def _parse_lpr(self, file):
        for g, d in read_lpr(file):
            self.constraints += [Constraints(g,
                d['residue' ].values.astype(object),
                d['position'].values.astype(int   ),
                d['nats'    ].values.astype(float ))]

    def _parse_con(self, file):
        for g, d in read_con(file):
            self.constraints += [Constraints(g,
                d['residue' ].values.astype(object),
                d['position'].values.astype(int   ),
                d['weight'  ].values.astype(float ))]
    
    def label(self, label):
        if type(label)==str:
            if label.endswith('.hpt'):
                names, hpt = read_hpt(label)
        elif type(label)==dict:
            names = label
        else:
            raise Exception('Provide labels through a dictionary or an hpt file.')
        for c in self.constraints:
            c.name = names[c.name]
        return self
    
    def fit(self, x):
        t, s, a = len(self), [], []
        for n, c in enumerate(self):
            sys.stderr.write('SequenceConstraints : fitting %s out of %s\r' % (1+n, t))
            s += [c.fit(x)]
            a += [c.name]
        return a, np.array(s).astype(np.float16)

    def export(self, file):
        w = open(file, 'w')
        for c in self:
            b = ['\t'+'\t'.join(map(str,i)) for i in zip(c.residues, c.positions, c.weights)]
            w.write('%s = {\n%s\n}\n\n' % (c.name, '\n'.join(b)))
        w.close()


