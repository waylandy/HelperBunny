import numpy as np

from io import StringIO, BytesIO
from bokeh.models import Range1d, Label, Title
from PIL import Image

from bokeh.plotting import figure, show
from bokeh.models import SingleIntervalTicker, LinearAxis

from ..ext.weblogo import weblogo

def SequenceLogoPlot(data, **kwargs):
    wfasta      = lambda x:'\n'.join(map(lambda x:'>%s\n%s'%x,enumerate(map(lambda x:''.join(x),x))))
    fasta       = StringIO(wfasta(data))
    seqs        = weblogo.read_seq_data(fasta)
    logodata    = weblogo.LogoData.from_seqs(seqs)
    logooptions = weblogo.LogoOptions(**kwargs)
    logoformat  = weblogo.LogoFormat(logodata, logooptions)
    return logodata, logoformat

def SequenceLogoViewer(data, plot_width=1000, plot_height=230, scale=60, title='', stat=True):
    params = {'title'           : '',
              'fineprint'       : '',
              'show_yaxis'      : False,
              'show_xaxis'      : True,
              'number_interval' : 999999,
              'resolution'      : 600,
              'tic_length'      : 1,
              'scale_width'     : True,
              'stacks_per_line' : 999999,
              'show_errorbars'  : False}
    stream = BytesIO(weblogo.png_print_formatter(*SequenceLogoPlot(data, **params)))

    img = Image.open(stream).convert('RGBA')
    img = np.array(img.convert('RGBA'))[::-1]
    img[(img[:,:,:3]==255).sum(2)==3,3]=0 # transparency
    a0  = np.where(img[:,:,3].sum(0)!=0)[0]
    a1  = np.where(img[:,:,3].sum(1)!=0)[0]
    img = img[a1.min():1+a1.max(),a0.min():1+a0.max()]

    plot               = figure(plot_width=plot_width, 
                                plot_height=plot_height, 
                                x_axis_type= None,
                                x_range=Range1d(0,plot_width/scale),
                                y_range=Range1d(0,1.09),
                                tools='xpan,reset',
                                toolbar_location=None)

    title = Label(text = title,
                  x        = 0,
                  x_offset = 35,
                  x_units  = 'screen', 
                  y        = plot.plot_height, 
                  y_offset = -64,
                  y_units  = 'screen',
                  background_fill_color = 'white',
                  background_fill_alpha = 50,
                  text_font_size        = '15pt',
                  text_font             = 'monospace')
    plot.add_layout(title)

    stat = Label(text      = 'seq=%s pos=%s'%data.shape,
                  x        = plot.plot_width,
                  x_offset = -14,
                  x_units  = 'screen', 
                  y        = plot.plot_height, 
                  y_offset = -46,
                  y_units  = 'screen',
                  background_fill_color = 'white',
                  background_fill_alpha = 50,
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

    plot.image_rgba(image=[img], x=.5, y=0, dw=data.shape[1], dh=1)
    return plot

