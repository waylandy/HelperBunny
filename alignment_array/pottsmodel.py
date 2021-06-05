import os
import sys
from subprocess import call
from itertools import groupby

import numpy as np

class Potts:
    def __init__(self):
        self.H     = None
    
    def __repr__(self):
        if type(self.H)!=np.ndarray:
            return f'< Potts: no data >'
        else:
            return f'< Potts: size={self.H.shape}, nseq={self.nseqs} >'
    
    def get_fn_apc(self, H):
        # calculate frobenius norm
        fn        = np.linalg.norm(H[:,:,:20,:20], axis=(2,3)) # index assumes that gap is last
        
        # do average product correction
        avg_sites = fn.mean(0)
        avg_total = fn.mean()
        fn_apc    = np.zeros(fn.shape)
        for i in np.arange(fn.shape[0]):
            for j in np.arange(i+1,fn.shape[0]):
                fn_apc[i, j] = fn[i, j] - avg_sites[i] * (avg_sites[j] / avg_total)
        fn_apc   += fn_apc.T
        
        return fn, fn_apc
    
    def from_ccmpred(self, AlignmentArray, bin='ccmpred'):
        if type(self.H)==np.ndarray:
            raise Exception('refusing to overwrite existing model')

        char        = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
        rand        = ''.join(np.random.choice(list(char), 6))
        A           = AlignmentArray.remove_inserts()
        
        file_psicov = f'ccmpred_{rand}.psicov'
        file_fn_apc = f'ccmpred_{rand}.mat'
        file_params = f'ccmpred_{rand}.potts'
        
        try:
            with open(file_psicov, 'w') as w:
                for sequence in A:
                    w.write(''.join(sequence)+'\n')
            cmd = f'{bin} -n 5000 -e 0.01 -r {file_params} {file_psicov} {file_fn_apc}'
            sys.stderr.write(f'{cmd}\n')
            call(cmd.split())
            
            get_line  = lambda x: list(map(float, x.strip().split('\t')))
            get_mat   = lambda x: np.array(list(map(get_line,x)))
            get_pair  = lambda x: tuple(map(int,x[1:].split())) if len(x.split())==3 else None
            reader    = groupby(open(file_params), lambda x: x.startswith('#'))
            fields    = get_mat(next(reader)[1])
            pairs     = []
            for i in reader:
                try:
                    pairs.append([get_pair(next(i[1])), get_mat(next(reader)[1])])
                except StopIteration:
                    break
            sites     = 1 + max(j for i in pairs for j in i[0])
            states    = pairs[0][1].shape[0]
            couplings = np.zeros((sites,sites,states,states))
            for (i, j), c in pairs:
                    couplings[i, j] = c
                    couplings[j, i] = c
                    
            for file in (file_psicov,file_fn_apc,file_params):
                if os.path.isfile(file):
                    os.remove(file)
            
            # verify ccmpred alphabet in source code
            self.states = np.array(list('ARNDCQEGHILKMFPSTWYV-'))
            
            self.nsites = A.shape[1]
            self.nseqs  = A.shape[0]
            self.H      = couplings
            self.fields = fields
            
            self.fn, self.fn_apc = self.get_fn_apc(self.H)
        
        except KeyboardInterrupt:
            for file in (file_psicov,file_fn_apc,file_params):
                if os.path.exists(file):
                    os.remove(file_psicov)
            raise Exception('Model was not generated')
        
        return self

    def dump(self, filename):
        if type(self.H)!=np.ndarray:
            raise Exception('no data to save')
        data = {'H':self.H, 'fields':self.fields, 'nseqs':self.nseqs, 'states':self.states}
        np.savez_compressed(filename, **data)
        
    def load(self, filename):
        if type(self.H)==np.ndarray:
            raise Exception('refusing to overwrite existing model')
        data = np.load(filename)
        self.H      = data['H']
        self.fields = data['fields']
        self.nseqs  = data['nseqs']
        self.states = data['states']
        self.nsites = self.H.shape[0]
        self.fn, self.fn_apc = self.get_fn_apc(self.H)
        return self

