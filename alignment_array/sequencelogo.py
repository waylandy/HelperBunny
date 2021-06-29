from subprocess import PIPE, Popen
from io import BytesIO
import numpy as np
from PIL import Image

from ..ext.weblogo import weblogo

from bokeh.plotting import figure, show
from bokeh.models import Range1d, Label, Title, SingleIntervalTicker, LinearAxis

import panel as pn
import panel.widgets as pnw

pn.extension()


params = {'title'           : '',
          'fineprint'       : '',
          'show_yaxis'      : False,
          'show_xaxis'      : True,
          'number_interval' : 999999,
          'tic_length'      : 1,
          'scale_width'     : True,
          'stacks_per_line' : 999999,
          'show_errorbars'  : False,
          'logo_font'       : "ArialMT"}

class SequenceLogo:
    def __init__(self, AlignmentArray):
        A             = AlignmentArray.remove_inserts()
        alphabet      = weblogo.seq.unambiguous_protein_alphabet
        counts        = np.array([(A==i).sum(0) for i in alphabet]).T
        self.logodata = weblogo.LogoData.from_counts(alphabet, counts)
        self.nseq, self.npos = A.shape
    
    def writer(self, stream, *args):
        nargs = len(args)
        if nargs == 0:
            return stream
        if nargs == 1:
            with open(args[0], 'wb') as w:
                w.write(stream)
            return None
        raise Exception('Failed to write stream.')
    
    def print_png(self, *args, **kwargs):
        logooptions = weblogo.LogoOptions(**kwargs)
        logoformat  = weblogo.LogoFormat(self.logodata, logooptions)
        stream      = weblogo.png_print_formatter(self.logodata, logoformat)
        return self.writer(stream, *args)
    
    def print_png_cropped(self, *args):
        params = {'title'           : '',
                  'fineprint'       : '',
                  'show_yaxis'      : False,
                  'show_xaxis'      : True,
                  'number_interval' : 999999,
                  'tic_length'      : 1,
                  'scale_width'     : True,
                  'stacks_per_line' : 999999,
                  'show_errorbars'  : False,
                  'logo_font'       : "ArialMT"}
        buffer = BytesIO()
        stream = BytesIO(self.print_png(**params))
        img    = Image.open(stream).convert('RGBA')
        img    = img.resize((np.array(img.size)/2).astype(int))        
        img    = np.array(img.convert('RGBA'))[::-1]
        img[(img[:,:,:3]==255).sum(2)==3,3]=0 # transparency
        a0     = np.where(img[:,:,3].sum(0)!=0)[0]
        a1     = np.where(img[:,:,3].sum(1)!=0)[0]
        img    = img[a1.min():1+a1.max(),a0.min():1+a0.max()]
        img    = img[:,img[:,:,3][0]!=0] # optional fix to remove superstretched elements?
        Image.fromarray(img[::-1]).save(buffer, format='png')
        return self.writer(buffer.getvalue(), *args)
    
    def print_eps(self, *args, **kwargs):
        logooptions = weblogo.LogoOptions(**kwargs)
        logoformat  = weblogo.LogoFormat(self.logodata, logooptions)
        stream      = weblogo.eps_formatter(self.logodata, logoformat)
        return self.writer(stream, *args)
    
    def print_eps_traced(self, *args, **kwargs):
        stream      = self.print_eps(**kwargs)
        gs          = 'gs'
        cmd         = f'{gs} -dNoOutputFonts -dEPSCrop -sDEVICE=eps2write -sstdout=%stderr -sOutputFile=- -'
        p = Popen(cmd.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate(stream)
        return self.writer(stdout, *args)

    def plot(self, plot_width=1200, plot_height=200, scale=36, title='', stats=True):
        dim  = (self.nseq, self.npos)
        
        plot = figure(plot_width=plot_width, 
              plot_height=plot_height, 
              x_axis_type= None,
              x_range=Range1d(0,plot_width/scale, bounds=(-4,5+dim[1])),
              y_range=Range1d(0,1.09, ),
              tools='xpan,reset',
              toolbar_location=None)

        title = Label(text     = title,
                      x        = 0,
                      x_offset = 20,
                      x_units  = 'screen', 
                      y        = plot.plot_height, 
                      y_offset = -55,
                      y_units  = 'screen',
                      background_fill_color = 'white',
                      text_font_size        = '12pt',
                      text_font             = 'monospace')
        plot.add_layout(title)

        stats = Label(text      = 'seq=%s pos=%s' % dim,
                      x        = plot.plot_width,
                      x_offset = -14,
                      x_units  = 'screen', 
                      y        = plot.plot_height, 
                      y_offset = -46,
                      y_units  = 'screen',
                      background_fill_color = 'white',
                      text_font_size        = '10pt',
                      text_align            = 'right',
                      text_font             = 'monospace')
        if stats:
            plot.add_layout(stats)
        
        ticker = SingleIntervalTicker(interval=5, num_minor_ticks=5)
        xaxis = LinearAxis(ticker=ticker)
        plot.add_layout(xaxis, 'below')
        plot.xaxis.major_label_text_font_style = "bold"
        plot.xaxis.major_label_text_font_size  = "14pt"
        plot.xaxis.major_label_text_font       = "monospace"
        plot.grid.visible  = False
        plot.yaxis.visible = False
        plot.xaxis.visible = True
        img = Image.open(BytesIO(self.print_png_cropped())).convert('RGB')
        img = np.array(img.convert('RGBA'))[::-1]
        plot.image_rgba(image=[img], x=.5, y=0, dw=dim[1], dh=1)
        return show(plot)
