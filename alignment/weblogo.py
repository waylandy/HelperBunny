import os 
import io
from datetime import datetime

import weblogo
from IPython.display import Image

class WebLogo:
    def __init__(self, alnarray, tmp='/tmp/', **kwargs):
        time          = lambda: datetime.now().strftime("%y%m%d-%H%M%S-%f")
        self.tempfile = '%shelperbunny.weblogo.%s' % (tmp,time())
        self.alnarray = alnarray.positions()
        self._build(**kwargs)
        
    def _build(self, title='', scale_width=True, show_yaxis=False, show_xaxis=True, stacks_per_line=100):
        fmt    = lambda x: '>_\n%s\n' % ''.join(x)
        plain  = self.tempfile+'.plain'
        wl_eps = self.tempfile+'.eps'
        wl_eps = 'eps.eps'
        with open(plain, 'w') as w:
            for i in self.alnarray:
                w.write(fmt(i))
        
        seqs     = weblogo.read_seq_data(open(plain, 'r'))
        logodata = weblogo.LogoData.from_seqs(seqs)
        
        logooptions                 = weblogo.LogoOptions()
        logooptions.title           = title
        logooptions.show_yaxis      = show_yaxis
        logooptions.show_xaxis      = show_xaxis
        logooptions.scale_width     = scale_width
        logooptions.stacks_per_line = stacks_per_line
        logooptions.fineprint       = ''
        logooptions.show_errorbars  = False
        
        logoformat  = weblogo.LogoFormat(logodata, logooptions) 
        
        self.logodata   = logodata
        self.logoformat = logoformat
        
    def show(self):
        png = weblogo.logo_formatter.png_print_formatter(self.logodata, self.logoformat)
        return Image(png)
    
        # eps = weblogo.eps_formatter(logodata, logoformat)
    
    
        
            
aar, name = hb.AlignmentArray('input/mj.fa')
input = aar.positions()[:,:30]




wl = WebLogo(aar.positions()[:,:50])
wl.show()




import io
from datetime import datetime

from PIL import Image

from helperbunny.ext.weblogo import weblogo 



name, aln    = hb.AlignmentArray(child[0])
alnarray     = aln.define_insertions(gap=.5).PositionArray()[:,:30]

wfasta       = lambda x:'\n'.join(map(lambda x:'>%s\n%s'%x,enumerate(map(lambda x:''.join(x),x))))
tempfasta    = io.StringIO(wfasta(alnarray))
seqs         = weblogo.read_seq_data(tempfasta)
logodata     = weblogo.LogoData.from_seqs(seqs)
logooptions  = weblogo.LogoOptions(
                    title           = '',
                    fineprint       = '',
                    show_yaxis      = False,
                    show_xaxis      = False,
                    scale_width     = True,
                    stacks_per_line = 50,
                    show_errorbars  = False
                                )
logoformat   = weblogo.LogoFormat(logodata, logooptions)
output       = weblogo.png_print_formatter(logodata, logoformat)

outstream    = io.BytesIO(output)

# img = Image.open(outstream).convert('RGBA')
# arr = np.array(img.convert('RGBA'))[::-1]
# arr[(arr==255).sum(2)==4,3] = 0
# col = np.where(arr[:,:,3].sum(axis=0)!=0)[0]
# arr = arr[:,col.min():1+col.max()]

from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource
from bokeh.models.glyphs import ImageURL
from bokeh.io import show

plot   = figure()
plot.image_rgba(image=[arr], x=0, y=0, dw=3, dh=0.5)
show(plot)

