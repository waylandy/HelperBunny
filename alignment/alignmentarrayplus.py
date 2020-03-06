import sys
from itertools import groupby
import numpy as np

# from .alignmentviewer import AlignmentViewer
# from .alignmentexporter import AlignmentExporter

class AlignmentArrayPlus(np.ndarray):

    def __new__(cls, ndarray, names=None):
        obj = np.asarray(ndarray, dtype=object).view(cls)
        return obj

    def __array_finalize__(self, obj):
        if obj is None: 
            return


