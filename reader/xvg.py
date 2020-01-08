import numpy as np

def xvg(xvg):
    x,y=[],[]
    r = open(xvg)
    for l in r:
        if l[0] not in '#@':
            l=l.split()
            x.append(float(l[0]))
            y.append(float(l[1]))
    r.close()
    return np.array(x), np.array(y)
