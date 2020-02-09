import sys
import re
from itertools import groupby


"""
    XMA format notes:
     My best attempt abstracting the format of CMA, MMA, SMA, etc.
"""


# NOTE 200207 : iteration of xmaBlock can only occur once, try and fix?
#               possible workaround is to re-iterate xmaBlock to rebuild xmaReader


class xmaBlock:
    """
        An iterator to alignment blocks from XMA
        Ignored information:
          The dollar sign line above each sequence name line
          The numbers in brackets inside the sequence name line
        Upon iteration, return (name, aln)
            name  : string; sequence name
            aln   : size 3 tuple; alignment with flanking sequences
    """
    def __init__(self, name, stream):
        head = re.compile(r'^\[0_\((.+?)\)\=(.+?)\((.+?)\)(.+)$')
        div  = lambda x: head.search(x).groups()
        _, self.name, seqs, info = div(name)
        self._stream = stream
        self.npos    = next(self._stream).count('*')
        assert self.npos>0
    
    def __iter__(self):
        return self.walk(fx_name=None, fx_seq=None)
    
    def walk(self, fx_name=None, fx_seq=None):
        fx_name  = fx_name if fx_name != None else lambda x: x
        fx_seq   = fx_seq  if fx_seq  != None else lambda x: x
        fmnam    = re.compile(r'^(.+?) +{\|(.+?)\((.+?)\)\|}(.+) +$')
        frstlast = lambda x: ' '.join((x[0],x[-1]))
        def namer(x):
            r    = fmnam.search(x)
            return x.strip() if r == None else frstlast(r.groups())
        fmseq    = re.compile(r'^{\((.*?)\)(.+?)\((.*?)\)}\*$')
        seqer    = lambda x: fmseq.search(x).groups()
        walker   = groupby(self._stream, lambda x: x.startswith('$'))
        while not next(walker):
            continue
        for k, g in walker:
            if k:
                continue
            yield fx_name(namer(next(g))), fx_seq(seqer(next(g)))

class xmaReader:
    """
        An iterator object to read XMA formatted files
        Upon iteration, returns xmaBlock

        Converts xma-formatted alignments to cfa using xmaReader.to_cfa()
          appends tag to sequence name to denote which xma block the sequence came from
          includes some checks to ensure formatting
          requires 2 iterations of xmaReader: one for checking, one for outputting
        Example:
            xmaReader('file.mma').to_cfa('file.cfa')
    """

    def __init__(self, infile):
        self.infile = infile
        
    def __iter__(self):
        return self.walk()

    def walk(self):
        isheader = lambda x: x.startswith('[0_')
        for k, g in groupby(open(self.infile), isheader):
            if k:
                name = next(g).strip()
            else:
                yield xmaBlock(name, g)

    def to_cfa(self, *args):
        h = sys.stdout if len(args)==0 else open(args[0], 'w')
        assert 1 == len(set(i.npos for i in self))
        ispos  = lambda x: not x.islower()
        poslen = lambda x: sum(1 for i in x[1][1] if ispos(i))
        names  = []
        alen   = poslen(next(iter(next(iter(self)))))
        for block in self:
            assert block.name not in names
            names += [block.name]
            for s in block:
                assert poslen(s)==alen
        outfmt = lambda x, y: '%s profile=%s\n%s%s%s\n\n' % (x[0],y,x[1][0].lower(),x[1][1],x[1][2].lower())
        for block in self:
            for s in block:
                h.write(outfmt(s, block.name))

