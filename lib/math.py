import numpy as np
from scipy import stats

def mode(x):
    i, c = np.unique(x, return_counts=True)
    return i[c.argmax()]

def seq2pmf(m, items='ACDEFGHIKLMNPQRSTVWY-X'):
    # probability mass function
    return np.array([[(i==t).sum() for t in items] for i in m.T])

def kld(p, q):
    # kullback leibler divergence
    return stats.entropy(p, q)

def jsd(p, q):
    # jensen shannon divergence
    p, q = p.flatten(), q.flatten()
    p, q = p/p.sum(), q/q.sum()
    m = (p + q) / 2
    return (kld(p, m) + kld(q, m)) / 2

def contrast_jsd(cfa_array, mask, cutoff=0.1): # confirm dependencies on the CFA array class?
    p = pmf(cfa_array[ mask], items='ACDEFGHIKLMNPQRSTVWY-')
    q = pmf(cfa_array[~mask], items='ACDEFGHIKLMNPQRSTVWY-')
    y = np.array([jsd(a,b) for a,b in zip(p,q)])
    y[((cfa_array[mask]!='-').sum(0)/mask.sum())<cutoff] = 0
    return np.arange(y.shape[0])+1, y

def contrast_kld(cfa_array, mask, cutoff=0.1): # confirm dependencies on the CFA array class?
    p = pmf(cfa_array[ mask], items='ACDEFGHIKLMNPQRSTVWY-')
    q = pmf(cfa_array[~mask], items='ACDEFGHIKLMNPQRSTVWY-')
    y = np.array([kld(a,b) for a,b in zip(p,q)])
    y[((cfa_array[mask]!='-').sum(0)/mask.sum())<cutoff] = 0
    return np.arange(y.shape[0])+1, y

def dihedral(*args):
    """
    Calculate dihedral from 4 coordinates
    """
    if len(args) == 4:
        p = np.array(args)
    elif len(args) == 1:
        p = np.array(args[0])
    else:
        raise ArithmeticError("That don't look like 4 coordinates to me...")
    b = p[:-1] - p[1:]
    b[0] *= -1
    v = np.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
    v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
    b1 = b[1] / np.linalg.norm(b[1])
    x = np.dot(v[0], v[1])
    m = np.cross(v[0], b1)
    y = np.dot(m, v[1])
    return np.degrees(np.arctan2( y, x ))

def gibbs(h, kB=0.0019872041):
    """
    NOT SURE IF THIS IS ACCURATE
    Calculate Gibbs entropy from vector (probability mass function).

    Boltzmann constant: 0.0019872041 kcal/(mol*K)
    https://en.wikipedia.org/wiki/Entropy_(statistical_thermodynamics)
    """
    if not isinstance(h, np.ndarray):
        h = np.array(h)

    v = h.flatten() / h.sum()
    return -kB * sum([i*np.log(i) for i in v[v!=0]])

def boltzmann(a, kB=0.0019872041, temperature=310):
    """
    CONFIRM THIS IS CORRECT??? SHOULD PSEUDOCOUNTS BE USED???
    Transforms a probability mass function to the Boltzmann distribution.

    Boltzmann constant: 0.0019872041 kcal/(mol*K)
    Temperature should be in Kelvin
    https://en.wikipedia.org/wiki/Boltzmann_distribution
    """
    if not isinstance(a, np.ndarray):
        a = np.array(a)
    a = a / a.sum()

    with np.errstate(divide='ignore'):
        B = -kB*temperature*np.log(a)
    m = B[B!=np.inf].max()
    B[B==np.inf] = m
    return B - m

def essential_dynamics(xyz):
    """
    NEEDS MORE TESTING & ABSTRACTION
        build a class later
        current implementation wastes potential intermediate data that is never used
    diagonalizes the cartesian coordinates, then cross correlate the eigenvectors
    ASSUMES SYSTEM IS TRANSPOSED TO THE CENTER OF GRAVITY????
    try using rasterized pcolormesh????
    """
    from sklearn import decomposition
    xyz  = np.array([sele.positions.flatten() for i in u.trajectory])
    pca = decomposition.PCA()
    pca.fit(xyz)
    modes      = 20
    n_atoms    = int(pca.components_.shape[0]/3)
    components = pca.components_[:modes,:].T
    var        = pca.explained_variance_
    shape      = (modes, n_atoms, 3)
    a1         = (components*var[:modes]).T.reshape(shape)
    a2         = components.T.reshape(shape)
    covariance = np.tensordot(a1.transpose(2, 0, 1),
                              a2.transpose(0, 2, 1),
                              axes=([0, 1], [1, 0]))
    diag       = np.power(covariance.diagonal(), 0.5)
    D          = np.outer(diag, diag)
    xcorr      = covariance/D
    return xcorr
