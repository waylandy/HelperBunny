from bokeh.plotting import figure, output_file, show
from bokeh.layouts import gridplot
from bokeh.models import Label, Title

from .sequencelogo import SequenceLogoViewer
from .bargraph import BarGraph

def CompareLogo(*args, sync=True, **kwargs):
    # can also pass plot_width, scale, plot_height 
    plot = []
    base = 0
    for n, (name, seq) in enumerate(args):
        try:
            if 'PositionArray' in dir(seq):
                seq  = seq.PositionArray(v=0)
                logo = SequenceLogoViewer(seq, title=name, **kwargs)
                plot.append([logo])
                base = logo.x_range if base==0 else base
            else:

                if len(seq)==3: # 3 args doesn't seem to work
                    p = BarGraph(*seq[:2], names=seq[2], **kwargs)
                if len(seq)==2:
                    p = BarGraph(*seq[:2], **kwargs)
                plot.append([p])

        except:
            pass

    if sync:
        for i in plot:
            i[0].x_range = base
    mplot = gridplot(plot, toolbar_location=None)
    show(mplot)
