
"""
    HelperBunny for Python3
    Author : Wayland Yeung
"""

import sys

__release__ = '200102'

sys.stderr.write("""Helper Bunny pre-Alpha [ Release %s ]
For development, turn off module compiling: sys.dont_write_bytecode = True
Other helful functions: "np.set_printoptions(suppress=True)"
""" % __release__)

"""
from importlib import reload
import sys
sys.dont_write_bytecode = True
import HelperBunny as hb
reload(hb)
"""

from .potpourri import *
from .math import *
from .reader import *
from .sequence import *
from .wrapper import *
from .plot import *

"""
git add .
git commit -m 'message'
git push
"""
