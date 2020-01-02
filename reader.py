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

def xpm(xpm, dtype=str, trans=False):
    atr   = ['title','legend','x-label','y-label','type','x-axis','y-axis']
    atr   = dict(list(zip(atr, ['']*len(atr))))
    comment   = lambda x: x.startswith('/*') and x.endswith('*/\n')
    unquote   = lambda x: x[1+x.find('"'):x.rfind('"')]
    uncomment = lambda x: x[2+x.find('/*'):x.rfind('*/')]
    colorchk  = lambda x: ' c #' in x
    colorkey  = lambda x: x[1+x.find('"'):x.find('c #')].strip()
    colorval  = lambda x: dtype(unquote(uncomment(x)))
    keypair   = lambda x: (colorkey(x), colorval(x))
    togcom, matread, togkey, passes, mat = True, False, 0, 0, []
    ixpm = gzip.open(xpm) if xpm.endswith('.gz') else open(xpm)
    for l in ixpm:
        if comment(l):
            if togcom:
                for k in atr:
                    _k = k+':'
                    if _k in l:
                        val = l[l.find(_k)+len(_k):].lstrip().rstrip('*/\n')
                        if atr[k] == val:
                            togcom = False
                            break
                        atr[k] += val
                continue
            continue
        elif l.startswith('static char'):
            # assume next line is the settings
            x, y, c, z = list(map(int, unquote(next(ixpm)).split()))
            assert z==1
            # assume next lines are the keys
            tr = dict(keypair(next(ixpm)) for n in range(c))
            matread = True
        elif matread:
            # assume uninterupted matrix lines
            d = ''.join(unquote(l if n == 0 else next(ixpm)) for n in range(y))
            d = [tr[i] for i in  d] if trans else d
            mat += [np.reshape(list(d), (y,x))]
            matread = False     
    ixpm.close()
    for k in atr:
        atr[k] = unquote(atr[k]).strip()
    return mat[0] if len(mat)==1 else np.array(mat), atr, tr
