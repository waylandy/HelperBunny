
import sys
from itertools import groupby

"""
    CFA format notes:
     1) sequences provided are assumed to be the FULL SEQUENCE
     2) duplicate identifiers are allowed in case of multiple profile hits
     3) sequences must be provided as a single unbroken line after the identifier
          a) feature sequences can be provided below and must be the same length
          b) each sequence must have the same number of feature sequences
"""

class CFA_Iterator:
    """
        Iterator object for CFA sequences
        Should also be capable of reading feature-linked CFA (not tested)
        Upon iteration, returns [name, seq]
    """
    def __init__(self, cfa):
        self._cfafile = cfa
        self.isname = lambda l: l[0]=='>'

    def __iter__(self):
        self._r = open(self._cfafile)
        self._i = groupby(self._r, self.isname)
        self._s = False
        return self

    def __next__(self):
        while True:
            head, group = next(self._i)
            if not head and self._s:
                self._s = False
            elif head:
                line = next(group)
            else:
                group = [i.rstrip() for i in group if not i.isspace()]
                assert len(set([len(i) for i in group])) == 1
                return line.rstrip()[1:], group
