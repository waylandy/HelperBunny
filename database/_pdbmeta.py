import os
import sys
import gzip
import tarfile
import requests
from ftplib import FTP

import pandas as pd

LIBHB    = '%s/.libhb/' % os.path.expanduser('~')

class pdbmeta:
    def __init__(self, path=LIBHB):
        self.path = path
        if not os.path.exists(self.path):
            os.makedirs(self.path)
    
        self.retrieve_entries()
        self.retrieve_pdb_chain_uniprot()
        
    def download(self, url, path='', unzip=True):
        if url.startswith('https://'):
            output  = path+url.split('/')[-1]
            request = requests.get(url, stream=True)
            with open(output, 'wb') as f:
                for chunk in request.iter_content(chunk_size=1024):
                    f.write(chunk)
                    f.flush()
        if url.startswith('ftp://'):
            ftp, d = url[6:].split('/', 1)
            output = path+d.split('/')[-1]
            ftp    = FTP(ftp)
            ftp.login()
            with open(output, 'wb') as fp:
                ftp.retrbinary('RETR %s' % d, fp.write)
            ftp.quit()
        if unzip and url.endswith('.gz'):
            with gzip.open(output, 'rb') as r:
                output = output[:-3]
                with open(output, 'wb') as w:
                    for l in r:
                        w.write(l)
            os.remove(output+'.gz')
        return output

    def retrieve_entries(self):
        url = 'ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx'
        self.entries = self.download(url, path=self.path)
        
        with open(self.entries, 'r') as r:
            it     = iter(r)
            header = next(it).strip().split(', ')
            dat    = [i.strip().split('\t') for i in it][1:]
        def num(i):
            try:    return float(i)
            except: return float('nan')
        self.entries = pd.DataFrame(dat, columns=header)
        self.entries['IDCODE'] = [i.lower() for i in self.entries['IDCODE']]
        self.entries['RESOLUTION'] = [num(i)  for i in self.entries['RESOLUTION']]
        return self.entries
    
    def retrieve_pdb_chain_uniprot(self):
        url = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz'
        self.pdb_chain_uniprot = self.download(url, path=self.path)

        with open(self.pdb_chain_uniprot, 'r') as r:
            it     = iter(r)
            next(it)
            header = next(it).strip().split('\t')
            dat    = [i.strip().split('\t') for i in it]
        self.pdb_chain_uniprot = pd.DataFrame(dat, columns=header)
        return self.pdb_chain_uniprot
