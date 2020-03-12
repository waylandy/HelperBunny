import sys
import numpy as np

from bokeh.plotting import figure, show
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d, SingleIntervalTicker, LinearAxis, Label
from bokeh.models.glyphs import Text, Rect

aa_color = { # based on default colors from weblogo3
    'R': 'blue',
    'K': 'blue',
    'D': 'blue',
    'E': 'blue',
    'N': 'blue',
    'Q': 'blue',
    'S': 'green',
    'G': 'green',
    'H': 'green',
    'T': 'green',
    'A': 'green',
    'P': 'green',
    'Y': 'black',
    'V': 'black',
    'M': 'black',
    'C': 'black',
    'L': 'black',
    'F': 'black',
    'I': 'black',
    'W': 'black',
    '-': 'white'
    }

def AlignmentViewer(data, scale=17, maxrows=20, color=aa_color, plot_width=1000, highlight=False, labels=[]):
    # data is a position array
    osh    = data.shape[0]
    if len(labels)!=0:
        assert len(labels)==osh
    labels = labels[:maxrows] if 0!=len(labels) else []
    data   = data[:maxrows] if data.shape[0] > maxrows else data
    x, y   = data.shape[::-1]
    x, y   = np.meshgrid(np.arange(1, 1+x), np.arange(0,y,1))
    x, y   = x.ravel(), -1-y.flatten()
    recty  = y + 0.5
    h      = 1/data.shape[0]
    text   = data.flatten()   
    dic    = dict(x      = x, 
                  y      = y, 
                  recty  = recty, 
                  text   = text,
                  color  = [color[i] if i in color else 'white' for i in text]
                  )

    source      = ColumnDataSource(dic)

    plot_height = data.shape[0]*scale
    fontsize    = scale/1.7
    plot   = figure(title            = None, 
                    plot_width       = plot_width, 
                    plot_height      = plot_height, 
                    x_axis_type      = None,
                    x_range          = Range1d(0, plot_width/scale), 
                    y_range          = Range1d(-plot_height/scale,0), 
                    tools            = "xpan,reset", 
                    min_border       = 0, 
                    toolbar_location = 'below',
                    background_fill_color = 'black',
                    background_fill_alpha = 0.12,
                   )

    ############################ VERY BAD IMPLEMENTATION, REWORK TO USE COLUMNDATASOURCE
    if 0!=len(labels):
        for n, lab in enumerate(labels):
            seqids = Label(text    = lab,
                          x        = 0,
                          x_offset = 16,
                          x_units  = 'screen', 
                          y        = -1-n, 
                          y_units  = 'data',
                          # background_fill_color = 'white',
                          # background_fill_alpha = 50,
                          text_font_size        = '10pt',
                          text_font             = 'monospace'
                          )
            plot.add_layout(seqids)
    ############################


    ticker = SingleIntervalTicker(interval=5, num_minor_ticks=5)
    xaxis = LinearAxis(ticker=ticker)
    plot.add_layout(xaxis, 'above')

    if highlight:
        rect  = Rect(x="x", y="recty",  
                     width=1, height=1,
                     fill_color="color",
                     fill_alpha=0.2, line_color=None)
        plot.add_glyph(source, rect)

    text  = Text(x="x", y="y", 
                 text="text", text_align='center',
                 text_color="color", text_font="monospace",
                 text_font_size="%spt" % fontsize) #, text_font_style='bold')
    plot.add_glyph(source, text)

    plot.grid.visible  = False
    plot.yaxis.visible = False
    plot.xaxis.major_label_text_font_style = "bold"
    plot.xaxis.major_label_text_font_size  = "12pt"
    plot.xaxis.major_label_text_font       = "monospace"
    sys.stderr.write('AlignmentViewer         : Showing %s of %s total sequences' % (data.shape[0], osh))
    return plot

