import re
import os
from itertools import islice

def recursive_walk(root_dir):
    for subdir, dirs, files in os.walk(root_dir):
        for file in files:
            yield f'{subdir}/{file}'

def natural_sort(items):
    convert = lambda x: int(x) if x.isdigit() else x.lower()
    return sorted(items, key=lambda x: [convert(i) for i in re.split('([0-9]+)', x)])

def chunker(iterable, n=10):
    it = iter(iterable)
    return iter(lambda: tuple(islice(it, n)), ())

