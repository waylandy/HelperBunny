import os
import sys
from subprocess import call

import numpy as np

class CD_HIT:
    def __init__(self):
        return

    def get_word_size(self, c):
        if c<0.4:
            raise Exception('Identity parameter must be between 1.0 - 0.5')
        return 5 if c>= 0.7 else 4 if c>= 0.6 else 3 if c>= 0.5 else 2

    def run_cdhit(self, AlignmentArray, bin='cdhit', c=0.98, T=20, M=0):
        char    = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
        rand    = ''.join(np.random.choice(list(char), 6))
        A       = AlignmentArray

        queries = c if type(c)==list else [c]
        queries = [(c, self.get_word_size(c)) for c in queries]

        nseqs   = A.shape[0]
        digits  = int(np.ceil(np.log10(1+nseqs)))
        unalign = lambda x: ''.join(x).replace('-','').upper()
        pad_int = lambda x: 'NDX'+str(x).zfill(digits)
        infile  = f'cdhit_{rand}.fasta'

        try:
            with open(infile, 'w') as w:
                for n, s in enumerate(map(unalign,A)):
                    w.write(f'>{pad_int(n)}\n{s}\n')

            indices = {}
            is_head = lambda x: x.startswith('>')
            get_ind = lambda x: int(x[4:])
            reader  = lambda x: map(get_ind,filter(is_head,open(x)))
            for c, n in queries:
                output = f'cdhit_{rand}_{c}.fasta'
                cmd    = f'{bin} -c {c} -n {n} -T {T} -M {M} -i {infile} -o {output}'
                sys.stderr.write(cmd+'\n')
                call(cmd.split())

                os.remove(f'{output}.clstr')

                ind = np.zeros(nseqs, dtype=bool)
                ind[list(reader(output))] = True
                indices[c] = ind

                os.remove(output)

            os.remove(infile)
            return indices
        
        except KeyboardInterrupt:
            try:
                os.remove(infile)
                os.remove(output)
                os.remove(f'{output}.clstr')
            except:
                pass
        raise Exception('Failed!')
