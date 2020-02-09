from itertools import groupby

"""
    CFA format notes:
     1) sequences provided are ASSUMED to be the FULL SEQUENCE
     2) duplicate identifiers are allowed in case of multiple profile hits
     3) sequences must be provided as a single unbroken line after the identifier
          a) the first sequence MUST be amino acids with
             uppercase positions + lowercase insertions
             dashes are treated as positions
          b) feature sequences must have the same length
          c) each header must have the same number of feature sequences
"""


class cfaReader:
    """
        Iterator object for CFA sequences
        Should also be capable of reading feature-linked CFA more testing
        This object offers no direct error checking for formatting
        Upon iteration, returns ( name, tuple(seqs) )
    """
    def __init__(self, infile):
        self.infile = infile
    
    def __iter__(self):
        return self.walk(fx_name=lambda x: x[1:])
        
    def walk(self, fx_name=None, fx_seq=None):
        fx_name = fx_name if fx_name != None else lambda x: x
        fx_seq  = fx_seq  if fx_seq  != None else lambda x: x
        
        isseq   = lambda x: x.startswith('>')
        fmtnam  = lambda x: next(x).strip()
        preseq  = lambda x: tuple(i.strip() for i in x if not i.isspace())
        first   = True
        for k, g in groupby(open(self.infile), isseq):
            if k:
                first = False
                name  = fmtnam(g)
            else:
                if first:
                    continue
                yield fx_name(name), fx_seq(preseq(g))
