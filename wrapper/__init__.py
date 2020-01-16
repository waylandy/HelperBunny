
"""
Description
"""

from .tmalign import tmalign
from .weblogo3 import WebLogo

# subprocess.check_output
# subprocess.call

def _run(bin, *args, opt='-', **kwargs):
    args   = [str(i) for i in args]
    kwargs = ['%s%s %s'%(opt, k, kwargs[k]) for k in kwargs]
    cmd    = [str(bin)]+args+kwargs

