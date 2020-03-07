import numpy as np

from io import StringIO, BytesIO
from bokeh.models import Range1d
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

def SequenceLogoViewer(data, plot_width=1000, plot_height=230, scale=60):
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
                                tools="xpan,reset")

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

