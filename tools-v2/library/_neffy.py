import os
import math
import string
import random
import datetime

import neffy

"""
# https://maryam-haghani.github.io/NEFFy/usage_guide.html#python

# pip3 install neffy==0.1.1

import _neffy

aln_strings = [
    "-EVAKDADLVIEAIPE--IFDLKKKVFSEIEQYCP-",
    "--EVADADLVIEAIPE--IFDL--RVFSEIEQYAP-",
    "--EAKNLGLVAEVFPQERFWDEVMKLAREVAELP--",
    "-EEAKNLGLVAEVFPQERFWDEVMKLAREVAELPPA",
    "-ETAKNLGLVAEVFPEEDFMEKVIEFAKNLTELPPG",
    "ENEAKELGLVAKIFPA--LWEETMKVAKSLAQ----",
    "---AKEIGLVLEVFP---FWDKALELARRVAAM---",
    "-QEAARIGLVNEVVPQ--FWDRVMEVANRLAGPP--",
    "NEEAKEIGLVLDYVPDDVFMDEVMKIAKQIAKNAP-",
    ]
output = _neffy.get_sequence_weights(aln_strings)
display(output)

a2m_strings = [
    'eVAKDADLVIEAIPEIFDLKKKVFSEIEQYCp',
    'EVADADLVIEAIPEIFDL--RVFSEIEQYAp',
    'EAKNLGLVAEVFPQerFWDEVMKLAREVAELP',
    'eEAKNLGLVAEVFPQerFWDEVMKLAREVAELPpa',
    'eTAKNLGLVAEVFPEedFMEKVIEFAKNLTELPpg',
    'enEAKELGLVAKIFPALWEETMKVAKSLAQ--',
    '-AKEIGLVLEVFP-FWDKALELARRVAAM-',
    'qEAARIGLVNEVVPQFWDRVMEVANRLAGPP',
    'neEAKEIGLVLDYVPDdvFMDEVMKIAKQIAKNAp',
    ]
output = _neffy.get_sequence_weights(a2m_strings)
display(output)

"""

def get_sequence_weights(
    a3m_strings,
    alphabet = "protein",               # Enum to specify the type of sequences in the MSA ("protein", "rna", or "dna") 
    check_validation = True,            # Validate the input MSA file based on alphabet or not 
    threshold = 0.8,                    # Similarity threshold for sequence weighting, must be between 0 and 1
    omit_query_gaps = True,             # Omit gap positions of query sequence from all sequences for NEFF computation
    is_symmetric = True,                # Consider gaps in number of differences when computing sequence similarity cutoff (asymmetric) or not (symmetric)
    non_standard_option = "AsStandard", # Enum to handle non-standard residues of the specified alphabet ("AsStandard", "ConsiderGap", "ConsiderGapInCutoff") 
    depth = "inf",                      # Depth of MSA to be used in NEFF computation (starting from the first sequence)
    gap_cutoff = 1,                     # Threshold for considering a position as gappy and removing that (between 0 and 1; 1 = no gappy position)
    ):

    ### parameters
    match alphabet:
        case "protein":
            alphabet = neffy.Alphabet.Protein
        case "rna":
            alphabet = neffy.Alphabet.RNA
        case "dna":
            alphabet = neffy.Alphabet.DNA
        case _:
            raise Exception('alphabet = {"protein", "rna", "dna"}')

    match non_standard_option:
        case "AsStandard":
            non_standard_option = neffy.NonStandardOption.AsStandard
        case "ConsiderGap":
            non_standard_option = neffy.NonStandardOption.ConsiderGap
        case "ConsiderGapInCutoff":
            non_standard_option = neffy.NonStandardOption.ConsiderGapInCutoff
        case _:
            raise Exception('non_standard_option = {"AsStandard", "ConsiderGap", "ConsiderGapInCutoff"}')

    ### staging
    get_random_string = lambda x: ''.join(random.choice("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") for _ in range(x))
    get_datetime_string = lambda: datetime.datetime.now().strftime("%y%m%d-%H%M%S")
    temp_name = f'{get_datetime_string()}-{get_random_string(8)}'
    aln_file = f"{temp_name}.a3m"

    with open(aln_file, "w") as w:
        for n, sequence in enumerate(a3m_strings):
            aln_string = ''.join(i for i in sequence if not i.islower())
            w.write(f">{n}\n{aln_string}\n")

    ### process
    msa_length, msa_depth, weights = neffy.compute_neff(file=aln_file, only_weights=True, alphabet=alphabet, 
        check_validation=check_validation, threshold=threshold, norm=neffy.Normalization.No_Normalization, 
        omit_query_gaps=omit_query_gaps, is_symmetric=is_symmetric, non_standard_option=non_standard_option, 
        depth=depth, gap_cutoff=gap_cutoff)
    os.remove(aln_file)
    assert len(weights) == len(a3m_strings)

    ### output : mixed
    output = {
        "msa_length": msa_length,
        "msa_depth": msa_depth,
        "neff": sum(weights) * (1 / math.sqrt(msa_length)),  # normalize by sqrt of sequence length
        "weights": weights,
        }
    return output
