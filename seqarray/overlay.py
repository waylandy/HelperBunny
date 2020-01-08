import sys


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
