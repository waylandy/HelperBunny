import sys
import numpy as np

"""
    BPPS output parsers:
     these objects help parse BPPS outputs
     these are not guaranteed to work as I cant find documentated standards
"""

class HPT:
    """
        Reads HPT outputs from the BPPS programs
        Attributes:
            n     : total number of groups, not including randoms
            names : dictionary mapping group numbers to names
            group : dictionary mapping group hierarchies
    """
    def __init__(self, hptfile):
        self.hptfile = hptfile
        self._populate()

    def _populate(self):
        r  = open(self.hptfile)
        it = iter(r)
        l  = next(it)
        while l[0]!='!':
            l = next(it).rstrip()
        self.n     = len(l)
        self.names = {}
        self.group = {}
        for i in range(self.n):
            tree, nem = next(it).rstrip()[:-1].split(' ')
            num , nam = nem.split('.')
            hir       = [1+n for n, _ in enumerate(tree) if _=='+']
            assert hir[-1]==int(num)
            self.names[hir[-1]] = nam
            self.group[hir[-1]] = hir
        r.close()

    def get(self, q):
        a = []
        for i in self.group:
            if q in self.group[i]:
                a += self.group[i]
        a = sorted(i for i in list(set(a)) if i>=q)
        return np.array([(i, self.names[i]) for i in a], dtype=object)

    def show(self):
        nal = max(len(str(i)) for i in self.names)
        nul = 1+max(len(self.names[i]) for i in self.names)
        for i in self.group:
            a, b = str(i).ljust(nal,' '), self.names[i].ljust(nul, ' ')
            c    = ' > '.join(self.names[_] for _ in self.group[i])
            sys.stdout.write('[ %s ] %s: %s\n' % (a,b,c))
