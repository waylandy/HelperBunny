
"""


    HelperBunny for Python3
    Author : Wayland Yeung



"""

import sys


__release__ = '200116'
__version__ = __release__


"""
Release numbers are basically the date updated
Edits from other contributers on the same day should append their initials to __release__

Keep release and version number synced
"""


sys.stderr.write("""Helper Bunny pre-Alpha [ Release %s ]
For development, turn off module compiling: sys.dont_write_bytecode = True
Other helful functions: "np.set_printoptions(suppress=True)"
""" % __version__)


"""
from importlib import reload
import sys
sys.dont_write_bytecode = True
import HelperBunny as hb
reload(hb)
"""

from .reader import *
from .seqarray import *
from .wrapper import *

from .lib.math import *
from .lib.plot import *
from .lib.potpourri import *


"""
git add .
git commit -m 'message'
git push
"""

aa = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
      'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
      'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
      'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
