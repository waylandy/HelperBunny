import os
import subprocess
from datetime import datetime
from IPython.display import Image, display

class WebLogo:
    """
        Plots weblogos from CFA arrays; requires weblogo3 and ghostscript
        Example:
            wl = WebLogo()
            wl.plot(CFA_Array(cfa_file.ar)

        BUG: ghostrscript aparently cant find too wide a vector ~52 aln positions
             use a raster is its too long or output in fragments
    """
    def __init__(self, bin='weblogo'):
        self._temp  = lambda: datetime.now().strftime("%y%m%d-%H%M%S-%f")
        self.bin    = bin

    def plot(self, aln, name=None, save=False, mode=None):
        ar2fa = lambda x : '\n'.join([''.join(i) for i in x])
        aln   = aln if isinstance(aln, str) else ar2fa(aln)
        if name==None:
            name = self._temp()
        tempfile = name+'.WEBLOGO_TEMP.fa'
        outfile  = name
        mode     ='ss' if mode==None and set(aln).issubset(set('BEST-HIGC\n')) else 'aa'
        ssmode   = '' if mode!='ss' else " -a GHICSTBE -C red GHI helix -C black CST loop -C blue BE sheet"
        with open(tempfile, 'w') as w:
            w.write(aln)
        cmd = "%s -f %s -o %s -n 500 --scale-width NO --errorbars NO --fineprint . -F jpeg --resolution 150 -Y NO%s" % (
            self.bin, tempfile, outfile+'.jpeg', ssmode)
        subprocess.call(cmd.split())
        display(Image(filename=outfile+'.jpeg'))
        os.remove(outfile+'.jpeg')
        if save:
            cmd = '%s -f %s -o %s -n 500 --scale-width NO --errorbars NO --fineprint . -F eps -Y NO%s' % (
                self.bin, tempfile, outfile+'.temp.eps', ssmode)
            subprocess.call(cmd.split())
            cmd = 'gs -o %s -dNoOutputFonts -sDEVICE=eps2write %s' % (
                outfile+'.eps', outfile+'.temp.eps')
            subprocess.call(cmd.split())
            os.remove(outfile+'.temp.eps')
        os.remove(tempfile)
