from bokeh.plotting import show
from bokeh.layouts import gridplot
from bokeh.models import Label, Title

from .sequencelogo import SequenceLogoViewer

def CompareLogo(*args, sync=True, **kwargs):
    # can also pass plot_width, scale, plot_height 
    plot = []
    for n, (name, seq) in enumerate(args):
        try:
            del logo
        except:
            pass
        logo = SequenceLogoViewer(seq, **kwargs)
        plot.append([logo])
        title = Label(text = name,
                      x        = 0,
                      x_offset = 35,
                      x_units = 'screen', 
                      y        = plot[n][0].plot_height, 
                      y_offset = -64,
                      y_units = 'screen',
                      background_fill_color = 'white',
                      background_fill_alpha = 50,
                      text_font_size        = '15pt',
                      text_font             = 'monospace',
                     )
        #title = Title(text=name,align='center',text_font='monospace',text_font_size='15pt')
        plot[n][0].add_layout(title, 'left')
    if sync:
        ip   = iter(plot)
        base = next(ip)[0].x_range
        for i in ip:
            i[0].x_range = base
    mplot = gridplot(plot, toolbar_location=None)
    show(mplot)
