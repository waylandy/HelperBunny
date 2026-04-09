import inspect
from multiprocessing import Pool
import numpy as np

import parasail

"""
# https://github.com/jeffdaily/parasail-python
# https://github.com/jeffdaily/parasail

# pip3 install parasail==1.3.4

import _parasail

sequence1 = "ADLIEVAKDADLVIEAIPEIFDLKKKVFSEIEQYCPDHTIFATNTSS"
sequence2 = "MLKRNIKPEEAKNLGLVAEVFPQERFWDEVMKLAREVAELPPLAVKY"
result = _parasail.get_pairwise_alignment(sequence1, sequence2, method="sg")

print(result["alignment1"])
print(result["comparison"])
print(result["alignment2"])


sequence_pairs = [
    ("EVAKDADLVIEAIPEIFDLKKKVFSEIEQYCP", "EEAKNLGLVAEVFPQERFWDEVMKLAREVAELPP"),
    ("QEAARIGLVNEVVPQERFWDRVMEVANRLAGPP", "EEAKNLGLVAEVFPQERFWDEVMKLAREVAELPPA"),
    ("NEEAKEIGLVLDYVPDDVFMDEVMKIAKQIAKNAP", "KAKEIGLVLEVFPEDVFWDKALELARRVAAMPEC"),
    ("EAEMGFVLHVYPKDRFWEEVDEYAREIAKMPAP", "DEVADADLVIEAIPEIFDLRVFSEIEQYAP"),
    ("ETAKNLGLVAEVFPEEDFMEKVIEFAKNLTELPPG", "ENEAKELGLVAKIFPAEELWEETMKVAKSLAQ"),
    ]
output = _parasail.get_pairwise_alignment(sequence_pairs, threads=5)
display(output)

for result in output:
    print(result["alignment1"])
    print(result["comparison"])
    print(result["alignment2"])
    print()

sequences = [
    "EVAKDADLVIEAIPEIFDLKKKVFSEIEQYCP",
    "QEAARIGLVNEVVPQERFWDRVMEVANRLAGPP",
    "NEEAKEIGLVLDYVPDDVFMDEVMKIAKQIAKNAP",
    "EAEMGFVLHVYPKDRFWEEVDEYAREIAKMPAP",
    "ETAKNLGLVAEVFPEEDFMEKVIEFAKNLTELPPG",
    ]
output = _parasail.get_distance_matrix(sequences)
display(output)

"""

def _get_pairwise_alignment(sequence1, sequence2, gap_open=11, gap_extend=1, method="sg", alphabet="auto"):

    ### parameters
    sequence_to_unaligned = lambda x: x.translate(str.maketrans("abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", "-."))
    unaligned1 = sequence_to_unaligned(sequence1)
    unaligned2 = sequence_to_unaligned(sequence2)

    match alphabet:
        case "auto":
            is_nucl = set(unaligned1.upper() + unaligned2.upper()).issubset({"A", "C", "G", "N", "T"})
            matrix = parasail.dnafull if is_nucl else parasail.blosum62
        case "prot":
            matrix = parasail.blosum62
        case "nucl":
            matrix = parasail.dnafull
        case _:
            raise Exception('alphabet = {"auto", "prot", "nucl"}')

    ### process
    match method:
        case "nw": # needleman-wunch
            result = parasail.nw_trace_striped_64(unaligned1, unaligned2, gap_open, gap_extend, matrix)
        case "sg": # semi-global
            result = parasail.sg_trace_striped_64(unaligned1, unaligned2, gap_open, gap_extend, matrix)
        case "sw": # smith-waterman
            result = parasail.sw_trace_striped_64(unaligned1, unaligned2, gap_open, gap_extend, matrix)
        case _:
            raise Exception('method = {"nw", "sg", "sw"}')
    
    n_positions = len(result.traceback.comp)
    n_identical = result.traceback.comp.count("|")
    n_similar = result.traceback.comp.count(":")
    n_gaps = result.traceback.comp.count(" ")

    segment1 = result.traceback.query.replace('-','')
    segment2 = result.traceback.ref.replace('-','')
    start1 = 1 + unaligned1.index(segment1)
    start2 = 1 + unaligned2.index(segment2)
    end1 = start1 + len(segment1) - 1
    end2 = start2 + len(segment2) - 1
    # start1, end1 = 1 + result.cigar.beg_query, 1 + result.end_query
    # start2, end2 = 1 + result.cigar.beg_ref, 1 + result.end_ref
    assert segment1 == unaligned1[start1-1:end1]
    assert segment2 == unaligned2[start2-1:end2]

    ### output : row oriented
    output = {
        # "sequence1": unaligned1,
        # "sequence2": unaligned2,
        "alignment1": result.traceback.query, 
        "alignment2": result.traceback.ref,
        "start1": start1,
        "end1": end1,
        "start2": start2,
        "end2": end2,
        "comparison": result.traceback.comp,
        "cigar": result.cigar.decode.decode("utf-8"),
        "identity": n_identical / n_positions,
        "similarity": (n_identical + n_similar) / n_positions,
        "gap": n_gaps / n_positions,
        "score": result.score, # alignment score
        }
    return output

def get_pairwise_alignment(*args, threads=1, **kwargs):

    match len(args):

        case 1:
            sequence_pairs, = args # list of sequence pairs

            signature = inspect.signature(_get_pairwise_alignment)
            default_keys = [k for k, v in signature.parameters.items() if v.default != inspect._empty]
            default_kwargs = {k: v.default for k, v in signature.parameters.items() if v.default != inspect._empty}

            kwargs = [kwargs.get(k, default_kwargs[k]) for k in default_keys]
            inputs = ([*i, *kwargs] for i in sequence_pairs)
            with Pool(threads) as pool:
                outputs = pool.starmap(_get_pairwise_alignment, inputs)
            
            return outputs
        
        case 2:
            sequence1, sequence2 = args  # two sequences

            return _get_pairwise_alignment(sequence1, sequence2, **kwargs)
        
        case _:
            raise Exception()

def get_distance_matrix(sequences, **kwargs):

    n_sequences = len(sequences)
    indices = [i for i in range(n_sequences)]
    index_pairs = [(i, j) for i in indices for j in indices[:i]]
    sequences_pairs = ((sequences[i], sequences[j]) for i, j in index_pairs)

    outputs = get_pairwise_alignment(sequences_pairs, **kwargs)
    distance_matrix = np.zeros((n_sequences, n_sequences))
    for (i, j), k in zip(index_pairs, outputs):
        distance_matrix[i,j] = distance_matrix[j,i] = k["distance"]

    ### format output : mixed
    output = {
        # "sequence"        : sequences,
        "distance_matrix" : distance_matrix,
        }
    return output

def get_distance_matrix(sequences, **kwargs):

    n_sequences = len(sequences)
    indices = [i for i in range(n_sequences)]
    index_pairs = [(i, j) for i in indices for j in indices[:i]]
    sequences_pairs = ((sequences[i], sequences[j]) for i, j in index_pairs)

    outputs = get_pairwise_alignment(sequences_pairs, **kwargs)
    distance_matrix1 = np.zeros((n_sequences, n_sequences))
    distance_matrix2 = np.zeros((n_sequences, n_sequences))
    for (i, j), k in zip(index_pairs, outputs):
        distance_matrix1[i,j] = distance_matrix1[j,i] = k["identity"]
        distance_matrix2[i,j] = distance_matrix2[j,i] = k["similarity"]

    ### format output : mixed
    output = {
        # "sequence"   : sequences,
        "identity"   : distance_matrix1,
        "similarity" : distance_matrix2,
        }
    return output
