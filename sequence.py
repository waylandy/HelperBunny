import sys
import re
from itertools import groupby
import numpy as np
import pandas as pd

"""
    XMA format notes:
     An attempt at a generalized format of CMA, MMA, SMA, etc.
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

"""
    CFA format notes:
     1) sequences provided are assumed to be the FULL SEQUENCE
     2) duplicate identifiers are allowed in case of multiple profile hits
     3) sequences must be provided as a single unbroken line after the identifier
          a) feature sequences can be provided below and must be the same length
          b) each sequence must have the same number of feature sequences
"""

class CFA_Iterator:
    """
        Iterator object for CFA sequences
        Should also be capable of reading feature-linked CFA (not tested)
        Upon iteration, returns [name, seq]
    """
    def __init__(self, cfa):
        self._cfafile = cfa
        self.isname = lambda l: l[0]=='>'

    def __iter__(self):
        self._r = open(self._cfafile)
        self._i = groupby(self._r, self.isname)
        self._s = False
        return self

    def __next__(self):
        while True:
            head, group = next(self._i)
            if not head and self._s:
                self._s = False
            elif head:
                line = next(group)
            else:
                group = [i.rstrip() for i in group if not i.isspace()]
                assert len(set([len(i) for i in group])) == 1
                return line.rstrip()[1:], group

class CFA(np.ndarray):
    """
        Representation of CFA sequences as arrays
        Should also be capable of reading feature-linked CFA (not tested)
        Attributes:
            ar    : sequence array (try subclassing np.ndarray object)
            names : names array
        Example:
            ca = CFA_Array(cfa_file)
            ca.ar
    """
    def __new__(cls, ndarray, names=None):
        obj = np.asarray(ndarray, dtype=object).view(cls)
        obj.names = names
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.names = getattr(obj, 'names', None)

    # def fmt(st, maxsize=12, fillchar='.', align='left'):
    #     stl     = len(st)
    #     if stl>maxsize:
    #         maxsize -=2
    #         digit = lambda x : len(str(x))
    #         ceil  = lambda x : int(np.ceil(x))
    #         floor = lambda x : int(np.floor(x))
    #
    #         mid = digit(stl - maxsize + digit(stl))
    #         topbot = maxsize-mid
    #         top, bot = ceil(topbot/2), floor(topbot/2)
    #         mids = st[top:-bot]
    #
    #         return '%s(%s)%s' % (st[:top], len(mids), st[-bot:])
    #     #, sum((len(st[:top]), digit(len(mids)), len(st[-bot:])))
    #     else:
    #         if align=='left':
    #             return st.ljust(maxsize, fillchar)
    #         if align=='right':
    #             return st.rjust(maxsize, fillchar)

    # def _update(self):
    #     ispos = lambda x: False if len(x)!=1 else not x.islower()

    def _frame(self, x):
        return self if len(self.shape)!=3 else self[:,:,x]

    def aln(self, f=0, fxn=None):
        ispos = lambda x: False if len(x)!=1 else not x.islower()
        fxn   = fxn if fxn!=None else ispos
        test  = self[0,:] if len(self.shape)!=3 else self[0,:,f]
        mask  = [fxn(x) for x in test]
        return self._frame(f)[:,mask]

    def out(self, f=0):
        self._frame(f)
        return 0
    # def out(self, f=0, aln=True):
    #     aln  = None if aln else lambda x: True
    #     bind = lambda x: [''.join(i) for i in x]
    #     return bind(self.aln(f=f, fxn=aln))

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
                self.maps[q] = self.AlignmentMap(self.seqs[q])

            ans += [self.maps[q].res(x)]
        return ans

    def aln(self, x, query):
        query = self._parse_query(query)
        ans   = []
        for q in query:
            if q not in self.maps:
                self.maps[q] = self.AlignmentMap(self.seqs[q])

            ans += [self.maps[q].aln(x)]
        return ans

    def frag2pos(self, frag, query):
        query = self._parse_query(query)
        ans   = []
        for q in query:
            if q not in self.maps:
                self.maps[q] = self.AlignmentMap(self.seqs[q])

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

def from_cfa(cfa, v=1):
    ispos = lambda x: x.isupper() or x=='-'
    isins = lambda x: x.islower()

    def aln2ar(line):
        x, pre = [''], 'a'
        for now in line:
            if   ispos(pre) and ispos(now):
                x += ['',now]
            elif isins(pre) and isins(now):
                x[-1] += now
            elif isins(pre) and ispos(now):
                x += [now]
            elif ispos(pre) and isins(now):
                x += [now]
            else:
                raise Exception('Unrecognized case: %s%s'%(pre,now))
            pre=now
        return x + [''] if ispos(x[-1][0]) else x

    def align(ref, target):
        assert len(''.join(ref)) == len(target)
        it = iter(target)
        return [''.join([next(it) for i in range(len(b))]) for b in ref]

    def data2ar(data):
        ref = aln2ar(data[0])
        if len(data)==1:
            return ref
        else:
            return [ref] + [align(ref, i) for i in data[1:]]

    name, data = [], []

    for i, (n, d) in enumerate(CFA_Iterator(cfa)):
        name  += [n]
        data  += [data2ar(d)]
        if i%200==0 and v>0:
            sys.stderr.write('[ CFA_Array ] Importing Row           : %s\r' % (1+i))
    if v>0:
        sys.stderr.write(    '[ CFA_Array ] Importing Row           : %s\n' % (1+i))
    ar = CFA(np.array(data, dtype=object), names=np.array(name, dtype=object))
    ar = ar if len(ar.shape)!=3 else ar.transpose(0,2,1)
    return ar

"""
    BPPS output parsers:
     these objects help parse BPPS outputs
     these are not guaranteed to work as I cant find documentated standards
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
