import re
import os
from itertools import islice

def chunker(iterable, n=10):
    it = iter(iterable)
    return iter(lambda: tuple(islice(it, n)), ())

def recursive_walk(rootdir):
    for subdir, dirs, files in os.walk(rootdir):
        for f in files:
            yield os.path.join(subdir, f)

def natural_sort(items):
    convert = lambda x: int(x) if x.isdigit() else x.lower()
    segment = lambda x: [convert(c) for c in re.split('([0-9]+)', x)]
    return sorted(items, key=segment)

