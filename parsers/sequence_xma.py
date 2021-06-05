import sys
import re
from itertools import groupby


"""
    XMA format notes:
    this is my best attempt at abstracting the cma, mma, and sma formats
    all formatting standards are inferred as i cant find the relevant documentation
    highly recommend converting to cfa format for easier parsing
        xma to cfa:
            removes top header, bracket accessions, becomes semi-compliant with fasta format
        cfa to xma:
            might require removal of flanking regions (unsure of standards)
"""


# NOTE 200207 : iteration of xmaBlock can only occur once, try and fix?
#               possible workaround is to re-iterate xmaBlock to rebuild xmaReader


        # mapgap is closed-source, but it seems to set a max sequence header size so headers may get cut (needs testing)

class XMAAlignmentBlock:
    """
        UNTESTED
        check if the addition of brackets in the names is properly isolated

        all programs that work with these formats are closed source, so many aspects are inferred
        known sequence alterations when aligning (rungaps/mapgaps):
            random brackets are sometimes added after the first space in the sequence header
            there is length for max sequence headers, long headers may be truncated in the output
            some programs can not parse the parenthetical flanking sequences, and their generation is not default

        Upon iteration, return (name, aln)
            name  : string; sequence name
            aln   : size 3 tuple; alignment with flanking sequences
    """
    def __init__(self, stream, name, size, params):
        self._stream   = stream
        self.name      = name
        self.size      = size
        self.params    = params
        self.positions = next(self._stream).count('*')
    
    def __repr__(self):
        x = self.name, self.size, self.positions
        return '<XMAAlignmentBlock: "%s" containing %s sequences and %s positions>' % x
    
    def __iter__(self):
        return self.walk()
    
    def walk(self, fx_name=None, fx_seq=None):
        # removes the bracket notation if it got added to the sequence header, else leave it alone
        # in most cases this should allow revert header to the original
        fmnam    = re.compile(r'^(.+?) +{\|(.+?)\((.+?)\)\|}(.+) +$')
        frstlast = lambda x: ' '.join((x[0],x[-1]))
        def namer(x):
            r    = fmnam.search(x)
            return x.strip() if r == None else frstlast(r.groups())
        # separates the sequence into flanking segments and the middle alignment part
        fmseq    = re.compile(r'^{\((.*?)\)(.+?)\((.*?)\)}\*$')
        seqer    = lambda x: fmseq.search(x).groups()
        walker   = groupby(self._stream, lambda x: x.startswith('$'))
        while not next(walker):
            continue
        for k, g in walker:
            if k:
                continue
            yield namer(next(g))[1:], seqer(next(g))

def read_xma(file):
    """
        An iterator object to read XMA formatted files
    """
    fr = re.compile(r'^\[\d+_\((.+?)\)\=(.+?)\((.+?)\)(.+)$')
    f1 = lambda x: bool(re.match('^\[\d+_\(', x))
    f2 = lambda x: fr.search(x).groups()
    it = iter(groupby(open(file), f1))
    it = iter(groupby(open(file), f1)) if next(it)[0] else it
    for k, g in it:
        if k:
            n = f2(next(g))
        else:
            yield XMAAlignmentBlock(g, n[1], n[2], n[3])
