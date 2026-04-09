import re
import os
import sys
import shutil
import random
import datetime
from itertools import groupby

from subprocess import Popen, call, PIPE, STDOUT

import numpy as np
import pandas as pd
import polars as pl
# import matplotlib.pyplot as plt
# import seaborn as sns

"""
# https://github.com/soedinglab/CCMpred

import os
import numpy as np
import pandas as pd
from lib.sequence import ccmpred

os.makedirs("data/test-ccmpred", exist_ok=True)
df = pd.read_parquet("data/up200402_tk261.parquet")
a3m_strings = df["a2m"].tolist()[:100]

couplings_npz = "data/test-ccmpred/up200402_tk261.npz"
ccmpred.get_couplings_npz(a3m_strings, couplings_npz)
couplings = np.load(couplings_npz)

alignment = a3m_strings[0]
coupling_scores = ccmpred.get_coupling_scores_from_a3m_string(alignment, couplings['couplings'], couplings['states'])
display(coupling_scores.shape)

position1, position2 = 126, 146
ccmpred.plot_pair_from_couplings(position1, position2, couplings['couplings'], couplings['states'])

"""

###

is_a2m = lambda x: bool(re.search("^[A-Za-z-.]+$", x))
count_postions = lambda x: sum(1 for i in re.finditer("[A-Z-]", x))
a2m_to_aln = lambda x: x.translate(str.maketrans("", "", "abcdefghijklmnopqrstuvwxyz."))

###

def iter_fasta(file):
    is_header = lambda x: x.startswith('>')
    compress  = lambda x: ''.join(_.strip() for _ in x)
    reader    = iter(groupby(open(file), is_header))
    reader    = iter(groupby(open(file), is_header)) if next(reader)[0] else reader
    for key, group in reader:
        if key:
            for header in group:
                header = header[1:].strip()
        else:
            sequence = compress(group)
            if sequence != '':
                yield header, sequence

def get_polars_from_fasta(fasta_file):
    schema = pl.Schema((("header", pl.String), ("sequence", pl.String)))
    table = pl.DataFrame(iter_fasta(fasta_file), orient="row", schema=schema)
    return table

def validate_a2m(a2m_file):
    header, sequence = next(iter_fasta(a2m_file))
    n_positions = count_postions(sequence)
    for n_sequences, (header, sequence) in enumerate(iter_fasta(a2m_file), 1):
        assert n_positions == count_postions(sequence)
        assert is_a2m(sequence)
    return n_positions, n_sequences

###

def get_fn_apc(couplings):
    # calculate frobenius norm
    fn        = np.linalg.norm(couplings[:,:,:20,:20], axis=(2,3)) # index assumes that gap is last
    
    # do average product correction
    avg_sites = fn.mean(0)
    avg_total = fn.mean()
    fn_apc    = np.zeros(fn.shape)
    for i in np.arange(fn.shape[0]):
        for j in np.arange(i+1,fn.shape[0]):
            fn_apc[i, j] = fn[i, j] - avg_sites[i] * (avg_sites[j] / avg_total)
    fn_apc   += fn_apc.T
    
    return fn, fn_apc

def get_couplings_npz(
    a3m_strings,
    couplings_npz,
    device  = 0,    #  -d DEVICE   Calculate on CUDA device number DEVICE (set to -1 to use CPU) [default: 0]
    threads = 1,    #  -t THREADS  Calculate using THREADS threads on the CPU (automatically disables CUDA if available) [default: 1]
    numiter = 50,   #  -n NUMITER  Compute a maximum of NUMITER operations [default: 50]
    epsilon = 0.01, #  -e EPSILON  Set convergence criterion for minimum decrease in the last K iterations to EPSILON [default: 0.01]
    lastk   = 5,    #  -k LASTK    Set K parameter for convergence criterion to LASTK [default: 5]
    idthres = 0.8,  #  -w IDTHRES  Set sequence reweighting identity threshold to IDTHRES [default: 0.8]
    lfactor = 0.2,  #  -l LFACTOR  Set pairwise regularization coefficients to LFACTOR * (L-1) [default: 0.2]
    ccmpred = "ccmpred",
    ):

    ### staging
    # get_random_string = lambda x: ''.join(random.choice("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") for _ in range(x))
    # get_datetime_string = lambda: datetime.datetime.now().strftime("%y%m%d-%H%M%S")
    # temp_name = f'{get_datetime_string()}-{get_random_string(8)}'
    temp_name = f'/tmp/temp_files'

    os.makedirs(temp_name, exist_ok=True)
    psicov_file = f"{temp_name}/alignment.psicov"
    mat_file = f"{temp_name}/output.mat"
    npy_file = f"{temp_name}/output.npy"
    # raw_file = f"{temp_name}/output.raw"

    a2m_to_aln = lambda x: x.translate(str.maketrans("", "", "abcdefghijklmnopqrstuvwxyz."))
    aln_strings = [a2m_to_aln(i) for i in a3m_strings]
    assert 1 == len(set(map(len,aln_strings)))
    with open(psicov_file, "w") as w:
        w.write('\n'.join(aln_strings))
    
    ### subprocess
    params = {
        "d": device,
        "t": threads,
        "n": numiter,
        "e": epsilon,
        "k": lastk,
        "w": idthres,
        "l": lfactor,
        "y": npy_file,
        # "r": raw_file,
        }
    if device == -1:
        del params["threads"]
    args = ' '.join([f"-{k} {v}" for k, v in params.items()])
    cmd = f"{ccmpred} {args} {psicov_file} {mat_file}"

    proc = Popen(f"stdbuf -o0 {cmd}".split(), stdout=PIPE, stderr=STDOUT, text=True)
    for buffer in proc.stdout:
        sys.stdout.write(buffer)
        sys.stdout.flush()
    proc.wait()

    ### parse output
    n_sites = len(aln_strings[0])
    couplings = np.load(npy_file).reshape(n_sites, n_sites, 21, 21)
    fn, fn_apc = get_fn_apc(couplings)
    shutil.rmtree(temp_name)
    
    ### write output : mixed format
    np.savez_compressed(couplings_npz, **{
        'n_seqs'    : len(aln_strings),
        'n_sites'   : n_sites,
        'states'    : list('ARNDCQEGHILKMFPSTWYV-'),
        'couplings' : couplings.astype(np.float32),
        'fn_apc'    : fn_apc.astype(np.float32),
        })
    return temp_name

###

def get_coupling_scores_from_a3m_string(a3m_string, couplings, states):

    a3m_to_aln = lambda x: ''.join(i for i in x if not i.islower())
    aln_string = a3m_to_aln(a3m_string)
    n_sites = len(aln_string)
    assert couplings.shape[:2] == (n_sites, n_sites)
    assert couplings.shape[2:] == (len(states), len(states))
    
    state_to_index = {i: n for n, i in enumerate(states)}
    sites = np.arange(n_sites)
    indices = np.array([state_to_index.get(i, -1) for i in aln_string])

    scores = couplings[sites[:,None], sites[None,:], indices[:,None], indices[None,:]]
    scores[:, indices == -1] = scores[indices == -1, :] = np.nan
    return scores

def plot_pair_from_couplings(position1, position2, couplings, states):
    pairs = pd.DataFrame(couplings[position1 - 1, position2 - 1], index=states, columns=states)
    abs_max = pairs.abs().values.max()

    plt.figure(figsize=(12, 9))
    sns.heatmap(pairs, vmin=-abs_max, vmax=abs_max, cmap="seismic", annot_kws={"size": 8},
        linewidths=1, linecolor='lightgrey', annot=True, fmt="0.2f", cbar=False)
    plt.xticks(np.arange(21) + 0.5, states, rotation=0, ha='center', size=14, family="monospace")
    plt.yticks(np.arange(21) + 0.5, states, rotation=0, va='center', size=14, family="monospace")
    plt.tick_params(axis="both", width=0)
    plt.ylabel(f"position {position1}", loc="top")
    plt.xlabel(f"position {position2}", loc="left")
    plt.show()