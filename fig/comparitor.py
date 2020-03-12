from bokeh.plotting import show
from bokeh.layouts import gridplot
from bokeh.models import Label, Title

from .sequencelogo import SequenceLogoViewer

def CompareLogo(*args, sync=True, **kwargs):
    # can also pass plot_width, scale, plot_height 
    plot = []
    for n, (name, seq) in enumerate(args):
        try:
            seq = seq.PositionArray(v=0)
        except:
            pass
        logo = SequenceLogoViewer(seq, title=name, **kwargs)
        plot.append([logo])
    if sync:
        ip   = iter(plot)
        base = next(ip)[0].x_range
        for i in ip:
            i[0].x_range = base
    mplot = gridplot(plot, toolbar_location=None)
    show(mplot)
