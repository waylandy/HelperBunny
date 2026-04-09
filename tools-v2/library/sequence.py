import re
import sys
import random
import string

assert re._MAXCACHE > 100

######### alphabets

iupac_dna              = "GATC"
iupac_dna_extended     = "GATCRYWSMKHBVDN"

iupac_rna              = "GAUC"
iupac_rna_extended     = "GAUCRYWSMKHBVDN"

iupac_protein          = "ACDEFGHIKLMNPQRSTVWY"
iupac_protein_extended = "ACDEFGHIKLMNPQRSTVWYBXZJUO" # same as ascii upper

ascii_letters          = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
ascii_lowercase        = "abcdefghijklmnopqrstuvwxyz"
ascii_uppercase        = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

######### validators : format

is_sequence  = lambda x: bool(re.search("^([A-Za-z*]+|[A-Za-z-.]+)$", x))  #  sequence     : unaligned or a2m         #  
is_unaligned = lambda x: bool(re.search("^[A-Za-z*]+$", x))                #  |  unaligned : upper, lower, asterisk   #  
is_alpha     = lambda x: bool(re.search("^[A-Za-z]+$", x))                 #  |  |  alpha  : upper, lower             #  
is_a2m       = lambda x: bool(re.search("^[A-Za-z-.]+$", x))               #  |  a2m       : upper, lower, dash, dot  #  or "aligned"
is_aln       = lambda x: bool(re.search("^[A-Z-]+$", x))                   #  |  |  aln    : upper, dash              #  or "trimmed"

######### conversions

def sequence_to_alpha(sequence):
    assert is_sequence(sequence)
    sequence_to_alpha = lambda x: x.translate(str.maketrans("abcdefghijklmnopqrstuvwxyz*", "ABCDEFGHIJKLMNOPQRSTUVWXYZX", "-."))
    return sequence_to_alpha(sequence)

def sequence_to_unaligned(sequence):
    assert is_sequence(sequence)
    sequence_to_unaligned = lambda x: x.translate(str.maketrans("abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", "-."))
    return sequence_to_unaligned(sequence)

def a2m_to_aln(sequence):
    assert is_a2m(sequence)
    a2m_to_aln = lambda x: x.translate(str.maketrans("", "", "abcdefghijklmnopqrstuvwxyz."))
    return a2m_to_aln(sequence)

def a2m_to_a2m_unpadded(sequence):
    assert is_a2m(sequence)
    a2m_to_a2m_unpadded = lambda x: x.replace('.', '')
    return a2m_to_a2m_unpadded(sequence)

######### validators : composition

def is_protein(sequence, mode=1):
    alpha = sequence_to_alpha(sequence)
    match mode:
        case 0:
            return bool(re.search("^[ACDEFGHIKLMNPQRSTVWY]+$", alpha.upper()))
        case 1:
            return bool(re.search("^[ACDEFGHIKLMNPQRSTVWYX]+$", alpha.upper()))
        case 2:
            return bool(re.search("^[ACDEFGHIKLMNPQRSTVWYBXZJUO]+$", alpha.upper()))
        case _:
            raise Exception("")

def is_dna(sequence, mode=1):
    alpha = sequence_to_alpha(sequence)
    match mode:
        case 0:
            return bool(re.search("^[GATC]+$", alpha.upper()))
        case 1:
            return bool(re.search("^[GATCN]+$", alpha.upper()))
        case 2:
            return bool(re.search("^[GATCRYWSMKHBVDN]+$", alpha.upper()))
        case _:
            raise Exception("")

def is_rna(sequence, mode=1):
    alpha = sequence_to_alpha(sequence)
    match mode:
        case 0:
            return bool(re.search("^[GAUC]+$", alpha.upper()))
        case 1:
            return bool(re.search("^[GAUCN]+$", alpha.upper()))
        case 2:
            return bool(re.search("^[GAUCRYWSMKHBVDN]+$", alpha.upper()))
        case _:
            raise Exception("")

######### statistics

def count_postions(sequence):
    assert is_a2m(sequence)
    count_postions = lambda x: sum(1 for i in re.finditer("[A-Z-]", x))
    return count_postions(sequence)

def count_occupied(sequence):
    assert is_a2m(sequence)
    count_occupied = lambda x: sum(1 for i in re.finditer("[A-Z]", x))
    return count_occupied(sequence)

def count_gaps(sequence):
    assert is_a2m(sequence)
    count_gaps = lambda x: x.count("-")
    return count_gaps(sequence)

######### parsing

def get_a2m_partitions(sequence, regex=False):
    assert is_a2m(sequence)
    
    if regex:
        pattern = re.compile(
            "^[a-z]*(?=[A-Z-])"             # insert n-term
            "|(?<=[A-Z-])[a-z]*?(?=[A-Z-])" # insert
            "|[A-Z-]"                       # position or gap
            "|(?<=[A-Z-])[a-z]*$"           # insert c-term
            "|^[a-z]*$"                     # insert only
            )
        matches = list(pattern.finditer(sequence))
        return [i.group() for i in matches]
    
    else:
        assert all(i.isalpha() for i in sequence if not i in {'-', '.'})
        current, previous = True, True
        partitions = ['']
        
        for i in sequence:
            current = i.islower()
            if current:
                if previous:
                    partitions[-1] += i
                else:
                    partitions += [i]
            else:
                if previous:
                    partitions += [i]
                else:
                    partitions += ['', i]
            previous = current
        
        if not previous:
            partitions += ['']
        
        return partitions


# split a2m into 3 parts : (nterm, aligned, cterm)
get_a2m_sections = lambda x: re.search("^([a-z]*)(?=[A-Z-])([A-Za-z-.]*?)([a-z]*)$", x).groups()

