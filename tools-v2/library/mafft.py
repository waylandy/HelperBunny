import os
import re
import shutil
import random
import datetime
from subprocess import call, PIPE, STDOUT, DEVNULL

from Bio import SeqIO

"""
# https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html

# apt-get update && apt-get install -y mafft

import mafft

sequences = [
    "IGLVLDDDDYVPDD*",
    "IGLVLYVPDD",
    "AAAAEVAKDADLVIEAIPAAAAAAAEIFDLK*",
    "QEAARIGLVNEVVPQERFWDRVMEVANRLAGPP-.",
    "AAAAAKELGLVAEVFPQERFWGEVMKLAAAAA",
    "NEEAKEIGLVLDYVPDDVFMDEVMKIAKQIAKNAP",
    "DEVADADLVIEAIPEIFDLRVFSEIEQYAP",
    "ETAKNLGLVAEVFPEEDFMEKVIEFAKNLTELPPG",
    ]
aligned_sequences = [
    'MEEAKELGLVAEVFPQERFWGEVMKLARWMAHV',
    'MEEAKEIGLVLDYVPDDVFMDGVMKIARWMAHV',
    'MEKAKEIGLVLEVFPEDVFWADALKLARWMAHV',
    ]

output = mafft.get_mafft_alignment(sequences, method="einsi")
display(output)

output = mafft.get_mafft_alignment_from_aligned(sequences, aligned_sequences, method="multipair")
display(output)

"""

def get_mafft_alignment(
    sequences,
    threads = 4,
    method = "fft-ns-2",
    mafft = "mafft",
    ):

    ### staging
    get_random_string = lambda x: ''.join(random.choice("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") for _ in range(x))
    get_datetime_string = lambda: datetime.datetime.now().strftime("%y%m%d-%H%M%S")
    temp_name = f'{get_datetime_string()}-{get_random_string(8)}'

    os.makedirs(temp_name, exist_ok=True)
    input_fasta = f"{temp_name}/input.fasta"
    output_aln = f"{temp_name}/output.aln"

    is_sequence  = lambda x: bool(re.search("^([A-Za-z*]+|[A-Za-z-.]+)$", x))
    sequence_to_alpha = lambda x: x.translate(str.maketrans("abcdefghijklmnopqrstuvwxyz*", "ABCDEFGHIJKLMNOPQRSTUVWXYZX", "-."))
    with open(input_fasta, 'w') as w:
        for n, sequence in enumerate(sequences):
            if is_sequence(sequence):
                w.write(f">{n}\n{sequence_to_alpha(sequence)}\n")
    
    ### subprocess
    match method:
        case 'fft-ns-2': # guide tree, progressive alignment, re-tree 2x
            cmd = f"{mafft} --reorder --retree 2 --thread {threads} {input_fasta}"
        case 'ginsi': # sequences are fully alignable with similar lengths
            cmd = f"{mafft} --reorder --maxiterate 1000 --globalpair --thread {threads} {input_fasta}"
        case 'linsi': # sequences have single alignable region and unalignable flanks
            cmd = f"{mafft} --reorder --maxiterate 1000 --localpair --thread {threads} {input_fasta}"
        case 'einsi': # sequences have multiple alignable regions and unalignable regions
            cmd = f"{mafft} --reorder --maxiterate 1000 --genafpair --thread {threads} {input_fasta}"
        case _:
            raise Exception('')
    with open(output_aln, 'wb') as wb:
        call(cmd.split(), stdout=wb, stderr=DEVNULL)
    
    ### parse output
    index_to_aligned = {int(h): s for h, s in SeqIO.FastaIO.SimpleFastaParser(open(output_aln))}
    aligned_strings = [index_to_aligned.get(n, None) for n, i in enumerate(sequences)]
    shutil.rmtree(temp_name)

    ### output : column oriented
    output = {
        # "sequence"   : sequences,
        "aln_string" : aligned_strings,
        }
    return output

def get_mafft_alignment_from_aligned(
    sequences,
    aligned_sequences,
    threads = 4,
    method = "6merpair",
    mafft = "mafft",
    ):

    ### staging
    get_random_string = lambda x: ''.join(random.choice("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") for _ in range(x))
    get_datetime_string = lambda: datetime.datetime.now().strftime("%y%m%d-%H%M%S")
    temp_name = f'{get_datetime_string()}-{get_random_string(8)}'

    os.makedirs(temp_name, exist_ok=True)
    input_fasta = f"{temp_name}/input.fasta"
    aligned_fasta = f"{temp_name}/aligned.fasta"
    output_aln = f"{temp_name}/output.aln"
    output_map = f"{input_fasta}.map"

    is_sequence = lambda x: bool(re.search("^([A-Za-z*]+|[A-Za-z-.]+)$", x))
    sequence_to_alpha = lambda x: x.translate(str.maketrans("abcdefghijklmnopqrstuvwxyz*", "ABCDEFGHIJKLMNOPQRSTUVWXYZX", "-."))
    alphas = [sequence_to_alpha(i) for i in sequences]
    with open(input_fasta, 'w') as w:
        for n, alpha in enumerate(alphas):
            if is_sequence(alpha):
                w.write(f">{n}\n{sequence_to_alpha(alpha)}\n")
    
    a2m_to_aln = lambda x: x.translate(str.maketrans("", "", "abcdefghijklmnopqrstuvwxyz."))
    assert all(map(is_sequence, aligned_sequences))
    aln_strings = [a2m_to_aln(i) for i in aligned_sequences]
    n_positions = len(aln_strings[0])
    assert all(n_positions == len(i) for i in aln_strings)
    with open(aligned_fasta, 'w') as w:
        for n, aln_string in enumerate(aln_strings):
            w.write(f">a{n}\n{aln_string}\n")

    ### subprocess
    match method:
        case 'auto':
            cmd = f"{mafft} --reorder --thread {threads} --{method} --addfull {input_fasta} --mapout {aligned_fasta}"
        case 'multipair':
            cmd = f"{mafft} --reorder --thread {threads} --{method} --addfull {input_fasta} --mapout {aligned_fasta}"
        case '6merpair':
            cmd = f"{mafft} --reorder --thread {threads} --{method} --addfull {input_fasta} --mapout {aligned_fasta}"
        case _:
            raise Exception('')
    with open(output_aln, 'wb') as wb:
        call(cmd.split(), stdout=wb, stderr=DEVNULL)

    ### parse output : map
    index_to_a2m_strings = {}
    for line in open(output_map):
        if line.startswith('#'):
            continue
        elif line.startswith('>'):
            index = int(line[1:])
            index_to_a2m_strings[index] = []
        else:
            char, residue, position = line.strip().split(", ")
            index_to_a2m_strings[index] += [[char, int(residue), int(position) if position != "-" else "-"]]
    for key, value in index_to_a2m_strings.items():
        for n, item in enumerate(value):
            if item[2] != "-":
                index_nterm = n
                break
        for n, item in enumerate(value[::-1]):
            if item[2] != "-":
                index_cterm = len(value) - n
                break
        index_to_a2m_strings[key] = ''.join(i[0] for i in value[:index_nterm]).lower() # nterm
        position_previous = 0
        for char, residue, position in value[index_nterm:index_cterm]:
            if position == '-':
                index_to_a2m_strings[key] += char.lower()
            else:
                index_to_a2m_strings[key] += '-' * (position - position_previous -1) + char.upper()
                position_previous = position
        index_to_a2m_strings[key] += '-' * (n_positions - position)
        index_to_a2m_strings[key] += ''.join(i[0] for i in value[index_cterm:]).lower()  # cterm
    a2m_strings = [index_to_a2m_strings.get(n, None) for n, i in enumerate(sequences)]
    aln_strings = [a2m_to_aln(i) for i in a2m_strings]
    assert all(len(i) == n_positions for i in aln_strings)
    shutil.rmtree(temp_name)

    ### output : column oriented
    output = {
        # "sequence"    : sequences,
        "a2m_string"  : a2m_strings,
        "a2m_trimmed" : aln_strings,
        }
    return output

