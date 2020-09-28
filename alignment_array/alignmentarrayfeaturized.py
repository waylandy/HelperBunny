import numpy as np

class AlignmentArrayFeaturized(np.ndarray):
    """
        uses references as array element dtype
    """

    def __new__(cls, ndarray, names=None):
        obj = np.asarray(ndarray, dtype=object).view(cls)
        return obj

    def __array_finalize__(self, obj):
        if obj is None: 
            return
