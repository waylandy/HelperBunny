from io import BytesIO
from subprocess import PIPE, Popen

import numpy as np
from PIL import Image

from bokeh.plotting import (
    figure,
    show,
    )
from bokeh.models import (
    Range1d, 
    Label, 
    Title, 
    SingleIntervalTicker, 
    LinearAxis
    )

import panel as pn
import panel.widgets as pnw

import weblogo

pn.extension()

WEBLOGO_KWARGS = {
    'title'           : '',
    'fineprint'       : '',
    'show_yaxis'      : False,
    'show_xaxis'      : True,
    'number_interval' : 999999,
    'tic_length'      : 1,
    'scale_width'     : True,
    'stacks_per_line' : 999999,
    'show_errorbars'  : False,
    'logo_font'       : "ArialMT"
    }

class SequenceLogo:
    def __init__(self, logodata, n_seqs=None):
        # logo.logodata.length
        self.logodata = logodata
        self.n_seqs = n_seqs
    
    def _write(self, stream, *args):
        nargs = len(args)
        if nargs == 0:
            return stream
        elif nargs == 1:
            with open(args[0], 'wb') as w:
                w.write(stream)
        else:
            raise Exception('Failed to write stream.')
    
    @staticmethod
    def from_alignment_array(alignment_array):
        A = alignment_array.remove_inserts()
        alphabet = weblogo.seq.unambiguous_protein_alphabet
        logodata = weblogo.LogoData.from_counts(alphabet, np.array([(A == i).sum(0) for i in alphabet]).T)
        return SequenceLogo(logodata, n_seqs = A.shape[0])
    
    def write_eps(self, *args, binary='gs', trace=False, **kwargs):
        kwargs = WEBLOGO_KWARGS if len(kwargs) == 0 else kwargs
        logooptions = weblogo.LogoOptions(**kwargs)
        logoformat = weblogo.LogoFormat(self.logodata, logooptions)
        buffer = weblogo.eps_formatter(self.logodata, logoformat)

        if trace:
            return self._write(buffer, *args)
        
        else:
            # # https://ghostscript.readthedocs.io/en/latest/Use.html
            # cmd = f'{binary} -dNoOutputFonts -dEPSFitPage -sDEVICE=eps2write -sstdout=%stderr -sOutputFile=- -'
            # proc = Popen(cmd.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
            # stdout, stderr = proc.communicate(buffer)
            return self._write(buffer, *args)

    def write_png(self, *args, crop=True, **kwargs):
        kwargs = WEBLOGO_KWARGS if len(kwargs) == 0 else kwargs
        image = Image.open(BytesIO(self.write_eps(**kwargs)))
        image.load(scale=3.0)

        buffer = BytesIO()
        image.save(buffer, format='png')
        img = np.array(Image.open(buffer).convert('RGBA'))
        base = np.where(img.mean(axis=(1,2)) != 255)[0].max()
        img = img[:base + 1,255 != img[base,:,:].mean(1),:]
        img[(img[:, :, :3] == 255).all(2), 3] = 0

        buffer = BytesIO()
        Image.fromarray(img).save(buffer, format='png')
        return self._write(buffer.getvalue(), *args)

    def to_display(self, **kwargs):
        kwargs = WEBLOGO_KWARGS if len(kwargs) == 0 else kwargs
        return Image.open(BytesIO(self.write_png(**kwargs)))

def plot_interactive(alignment_array, figsize=(1400, 160), scale=32, chunk_size=100, title='', stats=False):
    alignment_array = alignment_array.remove_inserts()
    n_seq, n_pos = alignment_array.shape
    width, height = figsize
    
    plot = figure(
        width                 = width, 
        height                = height, 
        x_axis_type           = None,
        x_range               = Range1d(0, width / scale, bounds=(0.5, max(n_pos, width / scale))),
        y_range               = Range1d(0, 1.08, ),
        tools                 = 'xpan,reset',
        toolbar_location      = None,
        )
    plot.add_layout(Label(
        text                  = title,
        x                     = 0,
        x_offset              = 20,
        x_units               = 'screen', 
        y                     = height, 
        y_offset              = -55,
        y_units               = 'screen',
        background_fill_color = 'white',
        text_font_size        = '12pt',
        # text_font             = 'monospace',
        ))
    plot.add_layout(Label(
        text                  = f'seq={n_seq}; pos={n_pos}' if stats else '',
        x                     = width,
        x_offset              = -14,
        x_units               = 'screen', 
        y                     = height, 
        y_offset              = -52,
        y_units               = 'screen',
        background_fill_color = 'white',
        text_font_size        = '12pt',
        text_align            = 'right',
        # text_font             = 'monospace',
        ))
    
    ticker = SingleIntervalTicker(interval=5, num_minor_ticks=5)
    xaxis = LinearAxis(ticker=ticker)
    plot.add_layout(xaxis, 'below')
    plot.xaxis.major_label_text_font_style = "bold"
    plot.xaxis.major_label_text_font_size  = "12pt"
    plot.grid.visible  = False
    plot.yaxis.visible = False
    plot.xaxis.visible = True

    chunks = np.arange(n_pos) // chunk_size
    for index in set(chunks):
        chunk = np.where(chunks == index)[0]
        buffer = Image.open(BytesIO(SequenceLogo.from_alignment_array(alignment_array[:,chunk]).write_png()))
        uint8 = np.array(buffer.convert('RGBA'))[::-1]
        uint32 = np.zeros(uint8.shape[:2], dtype=np.uint32)        
        view = uint32.view(dtype=np.uint8).reshape((*uint8.shape[:2], 4))
        view[:] = uint8
        plot.image_rgba(image=[uint32], x=chunk.min() + 0.5, y=0, dw=chunk.shape[0], dh=1)
    
    show(plot)

