import sys
from .seqarray import *

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
        For efficiency/speed, CFA_Map implements lazy loading; maps are generated in a need-based manner
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

class CFA_Overlay:
    """
        Overlays features found in ffa onto cfa
    """
    def __init__(self, ffa_file, cfa_file):
        self.cfa_file = cfa_file
        self._build_ref(ffa_file)
        self._reg = lambda x: x.replace('-','').upper()

    def _build_ref(self, ffa_file):
        self.ref = {}
        for i in CFA_Iterator(ffa_file):
            if i[0] in self.ref:
                sys.stderr.write('[ CFA_Overlay ] DUPLICATE ENTRY    : %s\n' % i[0])
            else:
                self.ref[i[0]] = i[1]

    def _overlay(self, ref, aln):
        seq, ft = ref[0], [iter(i) for i in ref[1:]]
        if len(self._reg(aln))!=len(seq):
            sys.stderr.write(    '[ CFA_Overlay ] MISMATCH           : \n' )
        return [aln]+[''.join(next(it) if i!='-' else '-' for i in aln) for it in ft]

    def to_cfa(self, *args):
        h = sys.stdout if len(args)==0 else open(args[0], 'w')
        for i in CFA_Iterator(self.cfa_file):
            n = i[0][::-1]
            n = n[n.index(':eliforp ')+9:][::-1]
            if n not in self.ref:
                sys.stderr.write('[ CFA_Overlay ] REFERENCE NOT FOUND: %s\n' % i[0])
            else:
                h.write('>%s\n%s\n\n' % (i[0], '\n'.join(self._overlay(self.ref[n], i[1][0]))))

class CFA_Names(dict):
    """
        Parses a name vector using a given standard format
    """
    def __init__(self, names, format=[]):
        self.names = names
        if 'up' in format:
            self.uniprot()
        if 'mf' in format:
            self.mapfaster()
        self._build()

    def _build(self):
        for k in self:
            self[k] = np.array(self[k], dtype=object)

    def uniprot(self):
        self['id'       ] = []
        self['name'     ] = []
        self['UniProtKB'] = []
        self['OS'       ] = []
        self['OX'       ] = []
        self['GN'       ] = []
        for name in self.names:
            try:
                fr  = []
                it  = iter(name.split())
                fr += next(it).split('|')+['']
                i   = next(it)
                while '=' not in i:
                    fr[3]+=i+' '
                    i     = next(it)
            except:
                fr  = ['','','','']
            self['id'       ] += [fr[1]]
            self['name'     ] += [fr[3]]
            self['UniProtKB'] += ['|'.join(fr[:3])]
            try:
                fr = []
                if i.startswith('OS='):
                    i     = i[3:]
                    fr   += ['']
                    while '=' not in i:
                        fr[-1]+= i+' '
                        i     = next(it)
                    fr[-1] = fr[-1][:-1]
                else:
                    fr += ''
                fr += [i[3:] if i.startswith('OX=') else '']
                i   = next(it)
                fr += [i[3:] if i.startswith('GN=') else '']
            except:
                fr  = ['','','','']
            self['OS'] += [fr[0]]
            self['OX'] += [fr[1]]
            self['GN'] += [fr[2]]

    def mapfaster(self):
        self['profile'] = []
        for name in self.names:
            assert ' profile:' in name
            self['profile'] += [name.split(' profile:')[-1].split()[0]]
