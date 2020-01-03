import os
import subprocess
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from IPython.display import Image, display

class WebLogo:
    """
        Plots weblogos from CFA arrays; requires weblogo3 and ghostscript
        Example:
            wl = WebLogo()
            wl.plot(CFA_Array(cfa_file.ar)
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

def draw_ss(ss, x, ypos=20, height=2, ax=None):
    ypos -= height/2
    for s, p in zip(ss, x):
        if s == '-':
            continue
        elif s in 'EB':  # blue
            pat=patches.Rectangle(width=1,height=height,fc='#0000ff',lw=0,xy=(p-0.5, ypos))
        elif s in 'GHI': # red
            pat=patches.Rectangle(width=1,height=height,fc='#ff0000',lw=0,xy=(p-0.5, ypos))
        elif s in 'STC': # black
            pat=patches.Rectangle(width=1,height=height/4,fc='#000000',lw=0, xy=(p-0.5,ypos+(height*(3/8))))
        if ax == None:
            plt.gca().add_patch(pat)
        else:
            ax.add_patch(pat)
