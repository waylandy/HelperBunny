import os
import re
import shutil
import random
import datetime
from subprocess import call, PIPE, STDOUT

import numpy as np

from Bio import SeqIO

"""
# https://www.bioinformatics.org/cd-hit/cd-hit-user-guide

# apt-get update && apt-get install -y cd-hit

import cdhit

sequences = [
    "EVAKDADLVIEAIPEIFDLKKKVFSEIEQYCP",
    "QEAARIGLVNEVVPQERFWDRVMEVANRLAGPP",
    "EEAKNLGLVAEVFPQERFWDEVMKLAREVAELPPA",
    "NEEAKEIGLVLDYVPDDVFMDEVMKIAKQIAKNAP",
    "KAKEIGLVLEVFPEDVFWDKALELARRVAAMPEC",
    "EAEMGFVLHVYPKDRFWEEVDEYAREIAKMPAP",
    "DEVADADLVIEAIPEIFDLRVFSEIEQYAP",
    "ETAKNLGLVAEVFPEEDFMEKVIEFAKNLTELPPG",
    "ENEAK",
    ]
sequences_blacklist = [
    "EVAKDADLVIEAIPEIFDLKKKVFSEIEQYCP",
    "DEVADADLVIEAIPEIFDLRVFSEIEQYAP",
    "ETAKNLGLVAEVFPEEDFMEKVIEFAKNLTELPPG",
    ]

output = cdhit.get_cdhit_clusters(sequences, percent_id=90)
display(output)

output = cdhit.get_cdhit2d_clusters(sequences, sequences_blacklist, percent_id=90)
display(output)

"""

def get_cdhit_clusters(
    sequences,          # a3m or unaligned proteins
    percent_id = 90,    # percent identity threshold, default 90%
    threads    = 4,     # number of threads, default 1; with 0, all CPUs will be used
    memory     = 0,     # memory limit (in MB) for the program, default 800; 0 for unlimitted;
    bandwidth  = 20,    # bandwidth of alignment, default 20
    min_length = 10,    # length of throwaway sequences, default 10
    tolerence  = 2,     # tolerance for redundance, default 2
    slow_mode  = False, # accurate but slow mode
    cdhit      = "cd-hit",
    ):

    ### parameters
    assert min_length > 3
    sequence_id = percent_id / 100
    assert sequence_id >= 0.5 and sequence_id <= 1.0
    word_size = 5 if sequence_id >= 0.7 else 4 if sequence_id >= 0.6 else 3

    ### staging
    get_random_string = lambda x: ''.join(random.choice("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") for _ in range(x))
    get_datetime_string = lambda: datetime.datetime.now().strftime("%y%m%d-%H%M%S")
    temp_name = f'{get_datetime_string()}-{get_random_string(8)}'

    os.makedirs(temp_name, exist_ok=True)
    input_fasta = f"{temp_name}/input.fasta"
    output_fasta = f"{temp_name}/output.fasta"
    output_clstr = f'{temp_name}/output.fasta.clstr'

    reformat = lambda x: ''.join(i if i != '*' else "X" for i in x.upper() if i.isalpha() or i=="*")
    with open(input_fasta, 'w') as w:
        for n, sequence in enumerate(sequences):
            w.write(f">{n}\n{reformat(sequence)}\n")
    
    ### subprocess
    cmd = f'{cdhit} -i {input_fasta} -o {output_fasta} -c {sequence_id} -b {bandwidth} -M {memory} -T {threads} -n {word_size} -l {min_length} -t {tolerence} -d 0 -g {int(slow_mode)}'
    # print(cmd)
    call(cmd.split(), stdout=PIPE, stderr=STDOUT)

    ### parse output : clusters
    index_to_cluster = {}
    pattern1 = re.compile(r"(?<=>Cluster )\d+$")
    pattern2 = re.compile(r"(?<=aa, >)\d+(?=... )")
    for line in open(output_clstr):
        if line.startswith('>'):
            cluster = int(pattern1.search(line).group())
        else:
            index_to_cluster[int(pattern2.search(line).group())] = cluster
    labels = np.array([index_to_cluster.get(n, -1) for n, i in enumerate(sequences)])
    
    ### parse output : centroids
    parser = SeqIO.FastaIO.SimpleFastaParser(open(output_fasta))
    indices = {int(h) for h, s in parser}
    centroids = np.array([n in indices for n, i in enumerate(sequences)])
    shutil.rmtree(temp_name)

    ### format output : column oriented
    output = {
        # "sequence" : sequences,
        "label"    : labels,
        "centroid" : centroids,
        }
    return output

def get_cdhit2d_clusters(
    sequences,           # a3m or unaligned proteins
    sequences_blacklist, # a3m or unaligned proteins
    percent_id = 90,     # percent identity threshold, default 90%
    threads    = 4,      # number of threads, default 1; with 0, all CPUs will be used
    memory     = 0,      # memory limit (in MB) for the program, default 800; 0 for unlimitted;
    bandwidth  = 20,     # bandwidth of alignment, default 20
    min_length = 10,     # length of throwaway sequences, default 10
    tolerence  = 2,      # tolerance for redundance, default 2
    slow_mode  = False,  # accurate but slow mode
    cdhit2d    = "cd-hit-2d",
    ):

    ### parameters
    assert min_length > 3
    sequence_id = percent_id / 100
    assert sequence_id >= 0.5 and sequence_id <= 1.0
    word_size = 5 if sequence_id >= 0.7 else 4 if sequence_id >= 0.6 else 3

    ### staging
    get_random_string = lambda x: ''.join(random.choice("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") for _ in range(x))
    get_datetime_string = lambda: datetime.datetime.now().strftime("%y%m%d-%H%M%S")
    temp_name = f'{get_datetime_string()}-{get_random_string(8)}'

    os.makedirs(temp_name, exist_ok=True)
    input_fasta = f"{temp_name}/input.fasta"
    blacklist_fasta = f"{temp_name}/blacklist.fasta"
    output_fasta = f"{temp_name}/output.fasta"
    output_clstr = f'{temp_name}/output.fasta.clstr'

    reformat = lambda x: ''.join(i if i != '*' else "X" for i in x.upper() if i.isalpha() or i=="*")
    with open(input_fasta, 'w') as w:
        for n, sequence in enumerate(sequences):
            w.write(f">{n}\n{reformat(sequence)}\n")
    with open(blacklist_fasta, 'w') as w:
        for n, sequence in enumerate(sequences_blacklist):
            w.write(f">x{n}\n{reformat(sequence)}\n")

    ### subprocess
    cmd = f'{cdhit2d} -i {blacklist_fasta} -i2 {input_fasta} -o {output_fasta} -c {sequence_id} -b {bandwidth} -M {memory} -T {threads} -n {word_size} -l {min_length} -t {tolerence} -d 0 -g {int(slow_mode)}'
    # print(cmd)
    call(cmd.split(), stdout=PIPE, stderr=STDOUT)

    ### parse output : clusters
    index_to_cluster, clusters_blacklist = {}, set()
    pattern1 = re.compile(r"(?<=>Cluster )\d+$")
    pattern2 = re.compile(r"(?<=aa, >)x?\w+(?=... )")
    for line in open(output_clstr):
        if line.startswith('>'):
            cluster = int(pattern1.search(line).group())
        else:
            match = pattern2.search(line).group()
            if match.startswith("x"):
                clusters_blacklist.add(cluster)
            else:
                index_to_cluster[int(pattern2.search(line).group())] = cluster
    labels = np.array([index_to_cluster.get(n, -1) for n, i in enumerate(sequences)])

    ### parse output : centroids
    parser = SeqIO.FastaIO.SimpleFastaParser(open(output_fasta))
    indices = {int(h) for h, s in parser}
    whitelist = np.array([n in indices for n, i in enumerate(sequences)])
    shutil.rmtree(temp_name)

    ### format output : column oriented
    output = {
        # "sequence"    : sequences,
        "label"       : labels,
        "whitelisted" : whitelist,
        }
    return output
