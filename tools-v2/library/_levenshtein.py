import inspect
from multiprocessing import Pool
import numpy as np

import Levenshtein

"""
# https://rawgit.com/ztane/python-Levenshtein/master/docs/Levenshtein.html
# pip3 install python-Levenshtein==0.27.3

import _levenshtein

sequence1 = "ADLIEVAKDADLVIEAIPEIFDLKKKVFSEIEQYCPDHTIFATNTSS"
sequence2 = "MLKRNIKPEEAKNLGLVAEVFPQERFWDEVMKLAREVAELPPLAVKY"
output = _levenshtein.get_pairwise_distance(sequence1, sequence2, metric="levenshtein")
display(output)

sequence_pairs = [
    ("EVAKDADLVIEAIPEIFDLKKKVFSEIEQYCP", "EEAKNLGLVAEVFPQERFWDEVMKLAREVAELPP"),
    ("QEAARIGLVNEVVPQERFWDRVMEVANRLAGPP", "EEAKNLGLVAEVFPQERFWDEVMKLAREVAELPPA"),
    ("NEEAKEIGLVLDYVPDDVFMDEVMKIAKQIAKNAP", "KAKEIGLVLEVFPEDVFWDKALELARRVAAMPEC"),
    ("EAEMGFVLHVYPKDRFWEEVDEYAREIAKMPAP", "DEVADADLVIEAIPEIFDLRVFSEIEQYAP"),
    ("ETAKNLGLVAEVFPEEDFMEKVIEFAKNLTELPPG", "ENEAKELGLVAKIFPAEELWEETMKVAKSLAQ"),
    ]
output = _levenshtein.get_pairwise_distance(sequence_pairs, metric="levenshtein", threads=5)
display(output)

sequences = [
    "EVAKDADLVIEAIPEIFDLKKKVFSEIEQYCP",
    "QEAARIGLVNEVVPQERFWDRVMEVANRLAGPP",
    "NEEAKEIGLVLDYVPDDVFMDEVMKIAKQIAKNAP",
    "EAEMGFVLHVYPKDRFWEEVDEYAREIAKMPAP",
    "ETAKNLGLVAEVFPEEDFMEKVIEFAKNLTELPPG",
    ]
output = _levenshtein.get_distance_matrix(sequences, threads=4, metric="levenshtein")
display(output)

"""

def _get_pairwise_distance(sequence1, sequence2, metric="levenshtein", func_kwargs={}):

    ### parameters
    sequence_to_unaligned = lambda x: x.translate(str.maketrans("abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", "-."))
    unaligned1 = sequence_to_unaligned(sequence1)
    unaligned2 = sequence_to_unaligned(sequence2)

    match metric:
        case "levenshtein":
            func = Levenshtein.distance
        case "ratio":
            func = Levenshtein.ratio
        case "hamming":
            func = Levenshtein.hamming
        case "jaro":
            func = Levenshtein.jaro
        case "jaro_winkler":
            func = Levenshtein.jaro_winkler
        case _:
            raise Exception('alphabet = {"levenshtein", "ratio", "hamming", "jaro", "jaro_winkler"}')
    
    ### format output : row oriented
    output = {
        # "sequence1": unaligned1,
        # "sequence2": unaligned2,
        "distance": func(unaligned1, unaligned2, **func_kwargs),
        }
    return output

def get_pairwise_distance(*args, threads=1, **kwargs):

    match len(args):

        case 1:
            sequence_pairs, = args # list of sequence pairs

            signature = inspect.signature(_get_pairwise_distance)
            default_keys = [k for k, v in signature.parameters.items() if v.default != inspect._empty]
            default_kwargs = {k: v.default for k, v in signature.parameters.items() if v.default != inspect._empty}

            kwargs = [kwargs.get(k, default_kwargs[k]) for k in default_keys]
            inputs = ([*i, *kwargs] for i in sequence_pairs)
            with Pool(threads) as pool:
                outputs = pool.starmap(_get_pairwise_distance, inputs)
            
            return outputs
        
        case 2:
            sequence1, sequence2 = args  # two sequences

            return _get_pairwise_distance(sequence1, sequence2, **kwargs)
        
        case _:
            raise Exception()

def get_distance_matrix(sequences, **kwargs):

    n_sequences = len(sequences)
    indices = [i for i in range(n_sequences)]
    index_pairs = [(i, j) for i in indices for j in indices[:i]]
    sequences_pairs = ((sequences[i], sequences[j]) for i, j in index_pairs)

    outputs = get_pairwise_distance(sequences_pairs, **kwargs)
    distance_matrix = np.zeros((n_sequences, n_sequences))
    for (i, j), k in zip(index_pairs, outputs):
        distance_matrix[i,j] = distance_matrix[j,i] = k["distance"]

    ### format output : mixed
    output = {
        # "sequence" : sequences,
        "distance" : distance_matrix,
        }
    return output
