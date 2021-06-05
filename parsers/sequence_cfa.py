from itertools import groupby

# CURRENTLY NOT BEING USED

"""
    CFA format standards:
     1) sequences provided are ASSUMED to be the FULL SEQUENCE
     2) duplicate identifiers are allowed in case of multiple profile hits
     3) sequences must be provided as a single unbroken line after the identifier
          a) the first sequence MUST be amino acids with
             uppercase positions + lowercase insertions
             dashes are treated as positions
          b) feature sequences must have the same length
          c) each header must have the same number of feature sequences
"""

def read_cfa(file):
    """
    generator: yields entries in a cfa file, (identifier, sequence-feature list)

    Should also be capable of reading feature-linked CFA more testing
    This object offers no direct error checking for formatting
    """
    f1 = lambda x: x.startswith('>')
    f2 = lambda x: tuple(i.strip() for i in x if not i.isspace())
    it = iter(groupby(open(file), f1))
    it = iter(groupby(open(file), f1)) if next(it)[0] else it
    for k, g in it:
        if not k:
            s = f2(g)
            yield n, s
        if k:
            n = list(g)[-1].strip()[1:]
