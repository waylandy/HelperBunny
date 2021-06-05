
"""

HelperBunny for Python3
Author : Wayland Yeung

Release numbers are basically the date updated
Edits from other contributers on the same day should append their initials to __release__

"""

import sys

__release__ = '??????'
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

from .parsers.sequence_fasta import read_fasta
from .parsers.sequence_xma   import read_xma

from .alignment_array.constructor     import gen_array
from .alignment_array.alignmentarray  import AlignmentArray
from .alignment_array.headerarray     import SequenceHeaders

from .alignment_array.pottsmodel      import Potts

from .alignment_array.constraints     import SequenceConstraints

from .database.taxonomy               import TaxonomyDatabase
















