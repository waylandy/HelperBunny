import os
import random
import datetime
from subprocess import call, DEVNULL

import numpy as np

"""
# https://github.com/weizhongli/cdhit/releases

import pandas as pd
pd.set_option('display.max_rows', 50)

sequences = sequences
blacklist = sequences[:2]

data = {'sequences': sequences}

seqids, clusters, centroids = run_cdhit_series(sequences, [0.8, 0.9, 0.7])
data.update(dict(zip([f'cdhit${int(i)}' for i in seqids * 100], clusters.T)))

seqids, whitelist = run_cdhit2d_series(sequences, blacklist, [0.85, 0.95, 0.75, 0.65])
data.update(dict(zip([f'cdhit2d${int(i)}' for i in seqids * 100], whitelist.T)))

pd.DataFrame(data).sort_values(['cdhit2d$65', 'cdhit2d$75', 'cdhit2d$85', 'cdhit2d$95'])

"""

get_random_string = lambda x: ''.join(random.choice('abcdefghijklmnopqrstuvwxyz0123456789') for i in range(x))
get_datetime_string = lambda: datetime.datetime.now().strftime("%y%m%d-%H%M%S")
get_temp_name = lambda: f'{get_datetime_string()}-{get_random_string(6)}'

def _seqid_to_wordsize(seqid):
    if seqid < 0.5 or seqid > 1.0:
        raise Exception()
    else:
        return 5 if seqid >= 0.7 else 4 if seqid >= 0.6 else 3

def _write_temp_fasta(sequences, temp_fasta, prefix=''):
    with open(temp_fasta, 'w') as w:
        for n, s in enumerate(sequences):
            w.write(f'>{prefix}{n:09d}\n{s}\n\n')

def _run_cdhit(fasta_file, cdhit_file, seqid, threads=0, memory=0, binary='cd-hit', verbose=False):
    clstr_file = f'{cdhit_file}.clstr'
    cmd = f'{binary} -i {fasta_file} -o {cdhit_file} -n {_seqid_to_wordsize(seqid)} -c {seqid} -M {memory} -T {threads} -d 0'
    call(cmd.split(), stdout=None if verbose else DEVNULL)
    assert os.path.exists(cdhit_file) and os.path.exists(clstr_file)
    return clstr_file

def _run_cdhit2d(fasta_file, exclude_file, cdhit_file, seqid, threads=0, memory=0, binary='cd-hit-2d', verbose=False):
    clstr_file = f'{cdhit_file}.clstr'
    cmd = f'{binary} -i {exclude_file} -i2 {fasta_file} -o {cdhit_file} -n {_seqid_to_wordsize(seqid)} -c {seqid} -M {memory} -T {threads} -d 0'
    call(cmd.split(), stdout=None if verbose else DEVNULL)
    assert os.path.exists(cdhit_file) and os.path.exists(clstr_file)
    return clstr_file

def run_cdhit_series(sequences, seqids, **kwargs):
    seqids = sorted(set(seqids))[::-1]
    temp_name = get_temp_name()
    temp_files = [f'{temp_name}-{n}.fasta' for n in range(len(seqids) + 1)]
    _write_temp_fasta(sequences, temp_files[0], prefix='')

    n_seqs, n_series = len(sequences), len(seqids)
    clusters = np.full((n_seqs, n_series), -1, dtype=int)
    centroids = np.full((n_seqs, n_series), False, dtype=bool)
    
    for n, (temp_fasta, temp_cdhit, seqid) in enumerate(zip(temp_files[:-1], temp_files[1:], seqids)):

        cluster = {}
        clstr_file = _run_cdhit(temp_fasta, temp_cdhit, seqid, **kwargs)
        for line in open(clstr_file) :
            if line.startswith('>'):
                label = int(line.split('>Cluster ')[1])
            else:
                cluster[int(line.strip().split('aa, >')[1].split('... ')[0])] = label
        clusters[:, n] = [cluster[i] if i in cluster else -1 for i in range(n_seqs)]
        os.remove(clstr_file)

        stream = filter(lambda x: x.startswith('>'), open(temp_cdhit))
        centroids[[int(i[1:]) for i in stream], n] = True
        os.remove(temp_fasta)
    
    os.remove(temp_cdhit)
        
    for n in range(clusters.shape[1] - 1):
        succ = dict(clusters[(clusters[:, n: n + 2] != -1).all(1), n: n + 2])
        clusters[:, n + 1] = [succ[i] if i != -1 else -1 for i in clusters[:, n]]

    return np.array(seqids[::-1]), clusters[:,::-1], centroids[:,::-1]

def run_cdhit2d_series(sequences, blacklist, seqids, **kwargs):
    seqids = sorted(set(seqids))[::-1]
    temp_name = get_temp_name()
    temp_files = [f'{temp_name}-{n}.fasta' for n in range(len(seqids) + 1)]
    _write_temp_fasta(sequences, temp_files[0], prefix='i2_')

    temp_exclu = f'{temp_name}.fasta'
    _write_temp_fasta(blacklist, temp_exclu, prefix='i1_')

    n_seqs, n_series = len(sequences), len(seqids)
    whitelist = np.full((n_seqs, n_series), True, dtype=bool)
    
    for n, (temp_fasta, temp_cdhit, seqid) in enumerate(zip(temp_files[:-1], temp_files[1:], seqids)):
        clstr_file = _run_cdhit2d(temp_fasta, temp_exclu, temp_cdhit, seqid)
        stream = filter(lambda x: 'aa, >i2_' in x, open(clstr_file))
        whitelist[[int(i.split('aa, >i2_')[1].split('... ')[0]) for i in stream], n] = False
        os.remove(clstr_file)
        os.remove(temp_fasta)

    os.remove(temp_exclu)
    os.remove(temp_cdhit)

    for n in range(whitelist.shape[1] - 1):
        whitelist[whitelist[:, n] == False, n + 1] = False

    return np.array(seqids[::-1]), whitelist[:,::-1]

