import sys
import os
import subprocess
from itertools import groupby

class DSSP:
    """
        Read DSSP output files; Dumps to featured-linked fa files (ffa)
        Provide a file handle to output, not a filename
        Example:
            DSSP(dssp_file).tocfa()

    """
    def __init__(self, dsspfile, ss=7):
        assert ss in [3, 7]
        self.ss       = ss
        self.tr       = dict(['GH', 'HH', 'IH', 'EE', 'BE', 'SL', 'TL', 'CL'])
        self.name     = os.path.splitext(os.path.basename(dsspfile))[0]
        self._build(open(dsspfile))

    def _build(self, h):
        # mostly taken from PDB.DSSP
        self.data, start = [], 0
        for l in h.readlines():
            sl = l.split()
            if len(sl) < 2:
                continue
            if sl[1] == "RESIDUE":
                start = 1
                continue
            if not start:
                continue
            if l[9] == " ":
                self.data += [l[13] if l[14]!='*' else l[14]]
                continue
            dssp_index = int(l[:5])
            resseq = int(l[5:10])
            icode = l[10]
            chainid = l[11]
            aa = l[13]
            ss = l[16]
            if ss == " ":
                ss = "C"
            try:
                NH_O_1_relidx = int(l[38:45])
                NH_O_1_energy = float(l[46:50])
                O_NH_1_relidx = int(l[50:56])
                O_NH_1_energy = float(l[57:61])
                NH_O_2_relidx = int(l[61:67])
                NH_O_2_energy = float(l[68:72])
                O_NH_2_relidx = int(l[72:78])
                O_NH_2_energy = float(l[79:83])
                acc = int(l[34:38])
                phi = float(l[103:109])
                psi = float(l[109:115])
            except ValueError as exc:
                if l[34] != " ":
                    shift = l[34:].find(" ")
                    NH_O_1_relidx = int(l[38 + shift : 45 + shift])
                    NH_O_1_energy = float(l[46 + shift : 50 + shift])
                    O_NH_1_relidx = int(l[50 + shift : 56 + shift])
                    O_NH_1_energy = float(l[57 + shift : 61 + shift])
                    NH_O_2_relidx = int(l[61 + shift : 67 + shift])
                    NH_O_2_energy = float(l[68 + shift : 72 + shift])
                    O_NH_2_relidx = int(l[72 + shift : 78 + shift])
                    O_NH_2_energy = float(l[79 + shift : 83 + shift])
                    acc = int((l[34 + shift : 38 + shift]))
                    phi = float(l[103 + shift : 109 + shift])
                    psi = float(l[109 + shift : 115 + shift])
                else:
                    raise ValueError(exc)
            self.data+=[(chainid, aa, ss, phi, psi)]

    def to_cfa(self, *args):
        h = sys.stdout if len(args)==0 else args[0]
        for i, g in groupby(self.data, lambda x: x!='*'):
            if not i:
                continue
            g   = list(g)
            ch  = g[0][0]

            gap = ['.' if i!='!' else '!' for i in g]
            gpi = [n for n, i in enumerate(gap) if i=='!']
            gpm = [i+1 for i in gpi] + [i-1 for i in gpi]
            gap = ''.join('!' if i in gpm else '.' for i in range(len(gap)) if i not in gpi)
            aa  = ''.join(i[1] for i in g if i!='!')
            ss  = ''.join(i[2] for i in g if i!='!')
            if len(gap)!=len(aa):
                sys.stderr.write('ERROR PARSING: %s' % self.name)
                continue
            h.write('>%s_%s\n%s\n%s\n%s\n\n'%(self.name, ch, aa, ss, gap))

def tmalign(pdb_a, pdb_b, bin='TMalign'):
    """
    Captures the output of TM-align; tested on TM-align (Version 20170708)
    Link: https://zhanglab.ccmb.med.umich.edu/TM-align/
    Use the fortran version!
    """

    result={}
    aln=False
    for l in subprocess.check_output([bin, pdb_a, pdb_b]).split('\n'):
        if   'Length of Chain_1' in l:
            result['Length of Chain_1'] = int(l.split()[3])
        elif 'Length of Chain_2' in l:
            result['Length of Chain_2'] = int(l.split()[3])
        elif 'Aligned length='   in l:
            for a in l.split(', '):
                a = a.split('=')
                if   'Alig'   in a[0]:
                    result['Aligned length'] = a[1]
                elif 'RMSD'   in a[0]:
                    result['RMSD'] = a[1]
                elif 'Seq_ID' in a[0]:
                    result['Seq_ID'] = a[2]
        elif 'normalized by length of Chain_1)' in l:
            result['TM-score1'] = float(l.split()[1])
        elif 'normalized by length of Chain_2)' in l:
            result['TM-score2'] = float(l.split()[1])
        elif '(":" denotes' in l:
            aln=1
        elif   aln==1:
            result['alignment']=[l]
            aln+=1
        elif aln==2:
            result['alignment'].append(l)
            aln+=1
        elif aln==3:
            result['alignment'].append(l)
            aln+=1
    return result
