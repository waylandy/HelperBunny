from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource, Range1d
from bokeh.models.glyphs import Text, Rect

def BarGraph(x, y, names=None, plot_width=1000, plot_height=230):
    n      = ['' for i in x] if names==None else names
    source = ColumnDataSource(dict(x=x,y=y,m=y/2,n=[' '+_ for _ in n]))
    plot   = figure(title            = None, 
                    plot_width       = plot_width, 
                    plot_height      = plot_height, 
                    x_axis_type      = None,
                    y_range          = Range1d(0, 1.09*max(source.data['y'])), 
                    tools            = "xpan,reset", 
                    min_border       = 0, 
                    toolbar_location = 'below',
                    background_fill_color = 'black',
                    background_fill_alpha = 0.12,
                   )

    rect  = Rect(x="x", y="m", height="y", width=0.9, fill_color="red",
                 fill_alpha=0.2, line_color=None)
    plot.add_glyph(source, rect)

    text  = Text(x="x", y=0, angle=90, angle_units='deg', text="n", text_align='left',
                 text_color="black", text_font="monospace", text_font_size="12pt")
    plot.add_glyph(source, text)
    plot.grid.visible  = False
    plot.yaxis.visible = False
    return plot
