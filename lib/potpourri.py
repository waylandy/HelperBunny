import re
import os
from datetime import datetime

def nsort(list, inplace=False):
    """
    natural sorting to read like a human
    """
    import re

    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    if not inplace:
        return sorted(list, key = alphanum_key)
    if inplace:
        return list.sort(key = alphanum_key)

def rwalk(rootdir):
    """
    simple generator for recursive walk
    """
    for subdir, dirs, files in os.walk(rootdir):
        for f in files:
            yield os.path.join(subdir, f)

def transpose(list_of_lists):
    """
    transposes 90 degrees, works on lists
    """
    return list(map(list, list(zip(*list_of_lists))))

now = lambda : datetime.now().strftime("%y%m%d-%H%M%S-%f")
day = lambda : datetime.now().strftime("%y%m%d")

# def moving_avg(x, w=3):
#     s=(w*2)+1
#     return [sum(x[n-w if n>w else 0:n+w+1])/s for n, xi in enumerate(x)]
