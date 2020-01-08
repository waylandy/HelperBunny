from itertools import groupby

class FA_Iterator:
    """
    An iterator to read FA files
    returns (name, seq)
    """
    def __init__(self, fafile, parsetags=False):
        self._fafile    = fafile
        self._tags      = parsetags
        self.isname     = lambda l: l[0]=='>'

    def __iter__(self):
        self._r = open(self._fafile)
        self._i = groupby(self._r, self.isname)
        self._s = False
        return self

    def _parse_tags(self, seq):
        seq = seq.split()
        chunk, k = [seq[0]], True
        for n, (k, group) in enumerate(groupby(seq[1:], lambda x: '=' in x or ':' in x)):
            if k == True:
                for g in group:
                    chunk += [g]
            else:
                if n != 0:
                    chunk[-1] += ' '+' '.join(list(group))
                else:
                    chunk += [' '.join(list(group))]
        return chunk

    def __next__(self):
        while True:
            head, group = next(self._i)
            if not head and self._s:
                self._s = False
            elif head:
                line = next(group)
            else:
                name = self._parse_tags(line.rstrip()[1:]) if self._tags else line.rstrip()[1:]
                seq  = ''.join(l.rstrip() for l in group)
                return name, seq
