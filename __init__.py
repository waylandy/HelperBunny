
"""

HelperBunny for Python3
Author : Wayland Yeung

Release numbers are basically the date updated
Edits from other contributers on the same day should append their initials to __release__

"""

import sys

__release__ = 'X'
__version__ = __release__

sys.stderr.write("""

Helper Bunny pre-Alpha [ Release %s ]
import sys 
sys.dont_write_bytecode = True

import numpy as np
np.set_printoptions(suppress=True)

Dependencies (update as i go along?):
sudo python3 -m pip install numpy pandas scipy bokeh panel pillow

""".strip() % __version__)






from .sequence_parser.fasta import read_fasta
from .sequence_parser.cfa import read_cfa
from .sequence_parser.xma import read_xma

from .alignment_array.sequenceheaders import SequenceHeaders

from .alignment_array.constraints import SequenceConstraints
from .alignment_array.constructor import constructor as AlignmentArray

from .database.taxonomy import TaxonomyDatabase
















