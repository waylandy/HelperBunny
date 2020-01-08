import sys
# from .seqarray import *




"""


MOST LIKELY BROKEN AT THE MOMENT



"""



class AlignmentMap:
    """
        Single alignment map for a CFA sequence
        To be used with CFA maps
    """
    def __init__(self, seq):
        self._build(seq)
        self._seq = seq

    def _build(self, seq):
        isres = lambda x: x != '-'
        isaln = lambda x: not x.islower()
        self._apos, self._rpos = [], []
        caln, cres = 0, 0
        for r in seq:
            if isres(r):
                cres += 1
                if isaln(r):
                    caln += 1
                    self._apos += [caln]
                    self._rpos += [cres]
                else:
                    self._apos += [caln+0.5]
                    self._rpos += [cres]
            else:
                caln += 1
                self._apos += [caln]
                self._rpos += [cres+0.5]

    def res(self, x):
        i = self._rpos.index(x)
        return self._apos[i], self._seq[i]

    def aln(self, x):
        i = self._apos.index(x)
        return self._rpos[i], self._seq[i]

    def frag2pos(self, frag):
        try:
            self._raw
        except:
            self._raw = filter(lambda x: x!='-', self._seq).upper()
        s  = len(frag)
        v  = []
        rv =  [i for i in range(len(self._raw)-s+1) if self._raw[i:i+s] == frag]
        for r in rv:
            d = self.res(r)
            v += [(r, d[0], d[1])]
        return v

class CFA_Map:
    """
        Maps residue number to alignment number and vice versa
        For efficiency/speed, CFA_Map uses lazy loading; maps are generated in a need-based manner
        Attributes:
            names : ndarray of sequence names
            seqs  : sequences in same order as names
            maps  : dictionary of AlignmentMap objects
        Example:
            cm = CFA_Map(cfa_file)
            cm.aln(137, lambda x: 'VEGFR' in x)
    """

    def __init__(self, cfa):
        self._from_file(cfa)
        self.maps  = {}

    def _from_file(self, cfa):
        self.names, self.seqs = list(zip(*[(name, seq[0]) for name, seq in CFA_Iterator(cfa)]))
        self.names, self.seqs = np.array(self.names, dtype=object), np.array(self.seqs, dtype=object)

    def search(self, fxn):
        return [n for n, e in enumerate(self.names) if fxn(e)]

    def _parse_query(self, query):
        if type(query)==int:
            return [query]
        if callable(query):
            return self.search(query)
        try:
            if all(type(q)==int for q in query):
                return query
        except:
            pass
        return None

    def res(self, x, query):
        query = self._parse_query(query)
        ans   = []
        for q in query:
            if q not in self.maps:
                self.maps[q] = AlignmentMap(self.seqs[q])

            ans += [self.maps[q].res(x)]
        return ans

    def aln(self, x, query):
        query = self._parse_query(query)
        ans   = []
        for q in query:
            if q not in self.maps:
                self.maps[q] = AlignmentMap(self.seqs[q])

            ans += [self.maps[q].aln(x)]
        return ans

    def frag2pos(self, frag, query):
        query = self._parse_query(query)
        ans   = []
        for q in query:
            if q not in self.maps:
                self.maps[q] = AlignmentMap(self.seqs[q])

            ans += [self.maps[q].frag2pos(frag)]
        return ans
