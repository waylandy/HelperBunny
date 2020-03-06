import numpy as np

from bokeh.plotting import figure, show
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot

import panel as pn
import panel.widgets as pnw

pn.extension()

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

def AlignmentViewer(data, plot_height=None, plot_width=1000, maxrows=30, color=aa_color):
    data   = data[:maxrows] if data.shape[0] > maxrows else data
    x, y   = data.shape[::-1]
    x, y   = np.meshgrid(np.arange(1, 1+x), np.arange(0,y,1))
    x, y   = x.ravel(), -1-y.flatten()
    recty  = y + 0.5
    h      = 1/data.shape[0]
    text   = data.flatten()
    source = ColumnDataSource(dict(x     = x.tolist(), 
                                   y     = y.tolist(), 
                                   recty = recty.tolist(), 
                                   text  = text.tolist(),
                                   color = [color[i] if i in color else 'white' for i in text]
                                  )
                             )
    
    # plot_height = data.shape[0]*15+50 if plot_height==None else plot_height
    plot_height = data.shape[0]*15+50

    plot   = figure(title            = None, 
                    plot_width       = plot_width, 
                    plot_height      = plot_height, 
                    x_range          = Range1d(0, plot_width/17), 
                    y_range          = Range1d(-plot_height/17,0), 
                    tools            = "xpan,reset", 
                    min_border       = 0, 
                    toolbar_location = 'below',
                    background_fill_color = 'black',
                    background_fill_alpha = 0.25,
                   )
    
    # rect  = Rect(x="x", y="recty",  
    #              width=1, height=1,
    #              fill_alpha=0.4, line_color=None)
    # plot.add_glyph(source, rect)
    
    text  = Text(x="x", y="y", 
                 text="text", text_align='center',
                 text_color="color", text_font="monospace",
                 text_font_size="10pt") #, text_font_style='bold')
    plot.add_glyph(source, text)
    
    # plot.axis.visible = False
    plot.grid.visible = False
    plot.yaxis.visible = False
    plot.xaxis.major_label_text_font_style = "bold"
    # plot.yaxis.minor_tick_line_width = 0
    # plot.yaxis.major_tick_line_width = 0
    return show(plot)

