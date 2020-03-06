
"""


    HelperBunny for Python3
    Author : Wayland Yeung



"""

import sys

__release__ = '200305'
__version__ = __release__





"""
Release numbers are basically the date updated
Edits from other contributers on the same day should append their initials to __release__

Keep release and version number synced
"""

sys.stderr.write("""Helper Bunny pre-Alpha [ Release %s ]
For development, turn off module compiling: sys.dont_write_bytecode = True
Other helful functions: "np.set_printoptions(suppress=True)"

Dependencies (update as i go along?):
sudo python3 -m pip install numpy pandas scipy bokeh panel weblogo



""" % __version__)





from .read import *


from .alignment.constructor import AlignmentArray as AlignmentArray


