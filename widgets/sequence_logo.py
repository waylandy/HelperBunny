import numpy as np

from io import StringIO, BytesIO
from bokeh.models import Range1d, Label, Title
from PIL import Image

from bokeh.plotting import figure, show
from bokeh.models import SingleIntervalTicker, LinearAxis

from ..ext.weblogo import weblogo

def SequenceLogo(AlignmentArray, output='png', **kwargs):
    A           = AlignmentArray.remove_inserts()
    alphabet    = weblogo.seq.unambiguous_protein_alphabet
    counts      = np.array([(A==i).sum(0) for i in alphabet]).T
    logodata    = weblogo.LogoData.from_counts(alphabet, counts)
    logooptions = weblogo.LogoOptions(**kwargs)
    logoformat  = weblogo.LogoFormat(logodata, logooptions)
    if output == 'eps':
        return BytesIO(weblogo.eps_formatter(logodata, logoformat))
        # Adobe vector formats can also be re-traced with ghostscript
        # gs -o output.eps -dNoOutputFonts -dEPSCrop -sDEVICE=eps2write input.eps    
    elif output == 'png':
        return BytesIO(weblogo.png_print_formatter(logodata, logoformat))
    else:
        raise Exception('Choose either "png" (raster) or "eps" (vector) format')

class SequenceLogoScroll:
    def __init__(self, AlignmentArray, plot_width=1200, plot_height=200, scale=36, title='', stat=True):
        params = {'title'           : '',
                  'fineprint'       : '',
                  'show_yaxis'      : False,
                  'show_xaxis'      : True,
                  'number_interval' : 999999,
                  'tic_length'      : 1,
                  'scale_width'     : True,
                  'stacks_per_line' : 999999,
                  'show_errorbars'  : False}

        b    = SequenceLogo(AlignmentArray, output='png', **params)
        img  = Image.open(b).convert('RGBA')
        img  = img.resize((np.array(img.size)/2).astype(int))
        img  = np.array(img.convert('RGBA'))[::-1]
        img[(img[:,:,:3]==255).sum(2)==3,3]=0 # transparency
        a0   = np.where(img[:,:,3].sum(0)!=0)[0]
        a1   = np.where(img[:,:,3].sum(1)!=0)[0]
        img  = img[a1.min():1+a1.max(),a0.min():1+a0.max()]
        img  = img[:,img[:,:,3][0]!=0] # optional fix to remove superstretched elements?

        dim  = (AlignmentArray.shape[0], sum(AlignmentArray.is_position()))

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

        stat = Label(text      = 'seq=%s pos=%s' % dim,
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
        if stat: plot.add_layout(stat)

        ticker = SingleIntervalTicker(interval=5, num_minor_ticks=5)
        xaxis = LinearAxis(ticker=ticker)
        plot.add_layout(xaxis, 'below')
        plot.xaxis.major_label_text_font_style = "bold"
        plot.xaxis.major_label_text_font_size  = "14pt"
        plot.xaxis.major_label_text_font       = "monospace"
        plot.grid.visible  = False
        plot.yaxis.visible = False
        plot.xaxis.visible = True

        plot.image_rgba(image=[img], x=.5, y=0, dw=dim[1], dh=1)
        
        self.plot  = plot
        self.title = title
        self.shape = np.array(dim)

    def show(self):
        show(self.plot)

