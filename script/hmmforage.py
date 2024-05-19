import re
import os
import sys
import shutil
import hashlib
import sqlite3
from itertools import groupby

import numpy as np
import pandas as pd

import pyhmmer

"""

db = HMMforage.run(
    'database.fa',
    'profiles.hmm',
    'hits.db',
    e_value = 0.01,
    db_size = None,
    threads = 6,
    buffer_size = 10000000,
    verbose = True,
)

db.get_profiles()
db.export_a2m('hits.a2m', n_pad=0, profile=None)
db.export_parquet('hits.parquet', block_size=500)

"""

class HMMforage:
    @staticmethod
    def get_size_fasta(fasta_file, verbose=True):
        it = filter(lambda x: x.startswith('>'), open(fasta_file, 'r'))
        for i, _ in enumerate(it, 1):
            if verbose and (i % 200 == 0):
                sys.stderr.write(f'{i}\r')
        if verbose:
            sys.stderr.write(f'{i}\n')
        return i        
    
    @staticmethod
    def iter_chunks_fasta(fasta_file, n_chars=0):
        n_chars = n_chars if n_chars > 0 else 1
        with pyhmmer.easel.SequenceFile(fasta_file, digital=True, alphabet=pyhmmer.easel.Alphabet.amino()) as seqs:
            buffer, buffer_size = [], 0
            for seq in seqs:
                buffer += [seq]
                buffer_size += seq.sequence.shape[0]
                if buffer_size > n_chars:                
                    yield buffer
                    buffer, buffer_size = [], 0
            if buffer_size > 0:
                yield buffer
    
    @staticmethod
    def hashify_sequences(seqs):
        hashmap = {}
        for seq in seqs:
            header = f'{seq.name.decode()} {seq.description.decode()}'
            sequence = seq.textize().sequence
            hash = hashlib.sha256((header+sequence).encode()).hexdigest()
            hashmap[hash] = (header, sequence)
            seq.name = hash.encode()
            seq.description = b''
        return hashmap, seqs
    
    @staticmethod
    def iter_hmmsearch_hits(hmms, seqs, **params):
        for hits in pyhmmer.hmmer.hmmsearch(hmms, seqs, **params):
            for hit in hits:
                for domain in hit.domains:
                    a2m = [
                        (domain.alignment.hmm_from - 1) * '-',
                        domain.alignment.target_sequence,
                        (hits.searched_nodes - domain.alignment.hmm_to) * '-',
                    ]
                    yield {
                        'hash'          : hit.name.decode(),
                        'hmm_name'      : domain.alignment.hmm_name.decode(),
                        'hmm_accession' : domain.alignment.hmm_accession.decode(),
                        'env_from'      : domain.env_from,
                        'env_to'        : domain.env_to,
                        'evalue'        : domain.i_evalue,
                        'score'         : domain.score,
                        'a2m'           : ''.join(a2m),
                    }
    
    @staticmethod
    def run(fasta_input, hmm_input, db_output, e_value=0.01, db_size=None, threads=6, buffer_size=1000, verbose=True):
        db = HMMhits(db_output, overwrite=True)
        db_size = HMMforage.get_size_fasta(fasta_input, verbose=verbose) if db_size == None else db_size
        n_seqs, n_hits, n_doms = 0, 0, 0

        for seqs in HMMforage.iter_chunks_fasta(fasta_input, n_chars=buffer_size):
            hashmap, seqs = HMMforage.hashify_sequences(seqs)
            n_seqs += len(seqs)
        
            with pyhmmer.plan7.HMMFile(hmm_input) as hmms:
                it = HMMforage.iter_hmmsearch_hits(hmms, seqs, cpus=threads, E=e_value, Z=db_size)
                hits = list(filter(lambda x: e_value >= x['evalue'], it))
                db.insert_hits(hits)
                n_doms += len(hits)
                
                hashes = [(hash, *hashmap[hash]) for hash in set(hit['hash'] for hit in hits)]
                sequences = [{'hash': a, 'header': h, 'sequence': s} for a, h, s in hashes]
                db.insert_sequences(sequences)
                n_hits += len(sequences)
        
            if verbose:
                sys.stderr.write(f'{n_hits} ({n_doms}) / {n_seqs}\r')
        
        return db

class HMMhits:
    def __init__(self, db_file, overwrite=False):
        if os.path.exists(db_file) and overwrite:
            os.remove(db_file)
        
        self.conn = sqlite3.connect(db_file)
        
        if overwrite:
            self._build()
        
    def _build(self):
        cur = self.conn.cursor()
        cur.execute("""
            CREATE TABLE sequences (
                hash      TEXT,
                header    TEXT,
                sequence  TEXT
            )""")
        cur.execute("""
            CREATE TABLE hits (
                hash           TEXT,
                hmm_name       TEXT,
                hmm_accession  TEXT,
                env_from       INTEGER,
                env_to         INTEGER,
                evalue         REAL,
                score          REAL,
                a2m            TEXT
            )""")
        self.conn.commit()
    
    def insert_sequences(self, data):
        if len(data) != 0:
            cur = self.conn.cursor()
            for row in data:
                cur.execute('INSERT INTO sequences VALUES (?,?,?)', (
                    row['hash'],
                    row['header'],
                    row['sequence'],
                ))
            self.conn.commit()
            
    def insert_hits(self, data):
        if len(data) != 0:
            cur = self.conn.cursor()
            for row in data:
                cur.execute('INSERT INTO hits VALUES (?,?,?,?,?,?,?,?)', (
                    row['hash'],
                    row['hmm_name'],
                    row['hmm_accession'],
                    row['env_from'],
                    row['env_to'],
                    row['evalue'],
                    row['score'],
                    row['a2m'],
                ))
            self.conn.commit()

    def is_valid(self):
        cur = self.conn.cursor()
        cur.execute("""
            SELECT 
                hits.hash
            FROM 
                hits
            WHERE
                hits.hash NOT IN (SELECT sequences.hash FROM sequences)
            GROUP BY
                hits.hash
            """)
        return cur.fetchone() == None 

    def get_profiles(self):
        cur = self.conn.cursor()
        cur.execute("""
            SELECT 
                hits.hmm_name,
                COUNT(hits.a2m)
            FROM 
                hits
            GROUP BY
                hits.hmm_name
            """)
        return dict(cur.fetchall())

    def export_a2m(self, a2m_file, n_pad=10, profile=None):
        assert self.is_valid()
        
        profiles = self.get_profiles()
        profile = list(profiles)[0] if (profile is None) and (len(profiles) == 1) else profile
        assert profile in profiles
        
        cur = self.conn.cursor()
        cur.execute(f"""
            SELECT 
                sequences.header,
                sequences.sequence,
                hits.a2m,
                hits.env_from,
                hits.env_to
            FROM 
                hits
            JOIN
                sequences ON sequences.hash = hits.hash
            WHERE
                hits.hmm_name = '{profile}'
            ORDER BY
                hits.hash      ASC,
                hits.env_from  ASC
            """)
        
        def _get_match(domain, sequence, env):
            rank = lambda x: abs(x[1] - env[1]) + abs(x[0] - env[0])
            matches = sorted([(1 + i.start(), i.end()) for i in re.finditer(domain, sequence)], key=rank)
            assert len(matches) > 0
            return matches[0]
        
        with open(a2m_file, 'w') as w:
            n_pad = 0 if n_pad < 0 else n_pad
            for header, sequence, a2m, *env in cur.fetchall():
                match = _get_match(a2m.replace('-','').upper(), sequence, env)
                left, right = sequence[:match[0]-1], sequence[match[1]:]
                padded = left[len(left)-n_pad:].lower() + a2m + right[:n_pad].lower()
                assert padded.replace('-','').upper() in sequence
                w.write(f'>{header}\n{padded}\n\n')
    
    def export_parquet(self, parquet_file, block_size=10000):
        assert self.is_valid()
        
        cur = self.conn.cursor()
        cur.execute("""
            SELECT 
                sequences.hash,
                sequences.header,
                sequences.sequence,
                hits.hmm_name,
                hits.hmm_accession,
                hits.env_from,
                hits.env_to,
                hits.evalue,
                hits.score,
                hits.a2m
            FROM 
                hits
            JOIN
                sequences ON sequences.hash = hits.hash
            ORDER BY
                hits.hash ASC
            """)

        def _get_match(domain, sequence, env):
            rank = lambda x: abs(x[1] - env[1]) + abs(x[0] - env[0])
            matches = sorted([(1 + i.start(), i.end()) for i in re.finditer(domain, sequence)], key=rank)
            assert len(matches) > 0
            return matches[0]

        def _iter_query(cursor):
            for k, group in groupby(cursor, lambda x: x[0]):
                group = sorted(set(map(tuple,group)), key=lambda x: x[5])
                for hash, header, sequence, hmm_name, hmm_accession, env_from, env_to, evalue, score, a2m in group:
                    hit_start, hit_end = _get_match(a2m.replace('-','').upper(), sequence, (env_from, env_to))
                    yield header, hmm_name, hmm_accession, hit_start, hit_end, evalue, score, a2m
        
        if os.path.isdir(parquet_file):
            shutil.rmtree(parquet_file)
        if os.path.isfile(parquet_file):
            os.remove(parquet_file)
        os.makedirs(parquet_file)
        
        buffer = []
        cols = ['header', 'hmm_name', 'hmm_accession', 'hit_start', 'hit_end', 'evalue', 'score', 'a2m']
        for n, row in enumerate(_iter_query(cur), 1):
            buffer += [row]
            if n % block_size == 0:
                pd.DataFrame(buffer, columns=cols).to_parquet(f'{parquet_file}/{n:09}')
                buffer = []
        pd.DataFrame(buffer, columns=cols).to_parquet(f'{parquet_file}/{n:09}')

def cluster_ranges(starts, ends, evalues):
    ind = np.argsort(starts)
    assert (starts <= ends).all()
    starts, ends, evalues = starts[ind], ends[ind], evalues[ind]
    
    maxacc = np.maximum.accumulate([0, *ends])[:-1]
    cluster = np.cumsum(starts > maxacc)

    # overlap = lambda a1, a2, b1, b2: max(0, min(a2, b2) - max(a1, b1))
    # for c in np.unique(cluster):
    #     mask = cluster == c
    #     ranges = np.array([starts[mask], ends[mask]]).T
    #     sizes = ranges[:,1] - ranges[:,0]
    #     dist = np.zeros([mask.sum(), mask.sum()], dtype=int)
    #     for ai, (a1, a2) in enumerate(ranges):
    #         for bi, (b1, b2) in enumerate(ranges[:ai+1]):
    #             dist[ai, bi] = dist[bi, ai] = overlap(a1, a2, b1, b2)
    #     distfrac = (dist / sizes)
    
    rank = evalues.argsort()
    iter_masks = (cluster == c for c in np.unique(cluster))
    keep = np.array([(evalues == evalues[c].min()) * c for c in iter_masks]).any(0)

    rev = np.argsort(ind)
    return cluster[rev], keep[rev]
