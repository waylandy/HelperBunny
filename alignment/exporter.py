import sys

class AlignmentExporter:
    def to_cfa(self, alnarray, name, output=sys.stdout):
        # DOES NOT WORK ON 3 dimensional arrays right now
        assert len(alnarray.shape)==2
        fname  = lambda x: '>%s\n' % x
        fseq   = lambda x: '%s\n\n' % ''.join(x)

        for s, n in zip(alnarray, name):
            output.write(fname(n))
            output.write(fseq(s))
    
    def to_fasta(self, alnarray, name, output):
        assert len(alnarray.shape)==2
        fname  = lambda x: '>%s\n' % x
        fseq   = lambda x: '%s\n\n' % ''.join(x).upper()
        
        for s, n in zip(alnarray, name):
            output.write(fname(n))
            output.write(fseq(s))
    
    def to_xma(self, alnarray, name, head='xmaoutput', output=sys.stdout):
        nseq, npos = alnarray.positions().shape
        head = ('[0_(1)=%s(%s){go=0,gx=0,pn=0.0,lf=0,rf=0}:\n(%s)%s\n\n' % (head, nseq, npos, '*'*npos))
        output.write(head)
        ispos    = lambda x: False if len(x)!=1 else not x.islower()
        ispos    = [i[0] for i in enumerate(map(ispos, alnarray[0])) if i[1]==True]
        alnarray = alnarray[:,min(ispos):max(ispos)+1]
        nres     = lambda x: sum(1 for i in x for j in i if j!='-')
        fseq     = lambda x: '{()%s()}*' % ''.join(x)
        for i, (n, s) in enumerate(zip(name, alnarray)):
            header = '$%s=%s(%s)\n>%s\n%s\n\n' % (i+1,nres(s),npos,n,fseq(s))
            output.write(header)
        output.write('_0].')

