from itertools import cycle
import numpy as np
import networkx as nx
from Bio import PDB

"""

VERY POOR INTEGRATION AT THIS POINT

"""


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

class PDBHeader(dict):
    def __init__(self, pdbfile):
        # check out tags for other cysteines
        self._build(pdbfile)

    def _build(self, pdbfile):
        self['title' ] = ''
        self['modres'] = {}
        self['hetnam'] = {}
        self['hetsyn'] = {}
        
        with open(pdbfile, 'r') as r:
            for l in r:
                tag = l[:6]
                if   tag == 'TITLE ':
                    self['title' ] += l[6:]
                elif tag == 'MODRES':
                    name, orig = (l[12:15], l[24:27])
                    self['modres'][name] = orig
                elif tag == 'HETNAM':
                    name, desc = l[11:14], l[15:]
                    if name not in self['hetnam']:
                        self['hetnam'][name] = desc.rstrip()
                    else:
                        self['hetnam'][name] += desc.rstrip()
                elif tag == 'HETSYN':
                    name, desc = l[11:14], l[15:]
                    if name not in self['hetsyn']:
                        self['hetsyn'][name] = desc.rstrip()
                    else:
                        self['hetsyn'][name] += desc.rstrip()
        self['title' ] = ' '.join(self['title'].split())

class AtomMask(dict):
    def __init__(self, structure):
        self.structure = structure
        self.atoms     = np.array(list(structure.get_atoms()), dtype=object)

    def map_atom(self, name, fxn):
        self[name] = np.array([fxn(i) for i in self.atoms], dtype=np.bool)
        
    def map_residue(self, name, fxn):
        flat       = lambda x: [z for y in x for z in y]
        self[name] = [[fxn(i)]*len(i) for i in self.structure.get_residues()]
        self[name] = np.array(flat(self[name]), dtype=np.bool)
    
class Polypeptide(dict):
    def __init__(self, structure, isres):
        self._build(structure, isres)
        
    def _build(self, struct, fxn):
        stagger = lambda x: zip(x, x[1:])
        peptide = lambda x, y: (x['C']-y['N'])<1.7
        connect = lambda x, y: peptide(x, y) if all(j in i for i in (x,y) for j in ('C','N')) else False
        for chain in struct:
            self[chain] = [[]]
            for r1, r2 in stagger([i for i in chain if fxn(i)]):
                self[chain][-1] += [r1]
                if not connect(r1, r2):
                    self[chain] += [[]]
            self[chain][-1] += [r2]
    
    def apply(self, *args):
        la = len(args)
        if la > 2 or la < 1:
            return None
        elif la == 1:
            res, chn = args[0], lambda x: x
        else:
            res, chn = args
        return dict((chn(k), [[res(r) for r in i] for i in self[k]]) for k in self)

class NeighborList(nx.Graph):
    def __init__(self, atoms, data=None, val=None, cutoff=1.7, **attr):
        super(NeighborList, self).__init__()
        self.val, self.is_built = val, True
        
        self.kdtree  = PDB.NeighborSearch(atoms)
        self._build(cutoff)
        
    def _build(self, cutoff):
        for a, b in self.kdtree.search_all(cutoff, level='A'):
            if not self.has_node(a):
                self.add_node(a)
            if not self.has_node(b):
                self.add_node(b)
            # add blacklisting/whitelisting in the future?
            self.add_edge(a, b)

    def walk(self, here, path=()):
        peek = lambda x, y: tuple(i for i in x if i.name==y)
        path = iter(path)
        that = next(path)
        while True:
            this  = that
            try:
                that  = next(path)
            except: 
                pass
            where = peek(self.neighbors(here), that)
            nch   = len(where)
            
            yield here
            
            if   nch>1:
                print('ERROR')
                raise
            elif nch<1:
                return
            here = where[0]

class Structure:
    def __init__(self, pdbfile, model=0):
        self.name     = pdbfile
        self.struct   = PDB.PDBParser(QUIET=True).get_structure('', pdbfile)[model]
        
        self.header   = PDBHeader(pdbfile)
        self.resDIC   = {**RES3to1,**dict((k,RES3to1[self.header['modres'][k]]) for k in self.header['modres'])}
        isres         = lambda x: x.resname in self.resDIC # should work in all PDB3.0 formats
        self.pp       = Polypeptide(self.struct, isres)
        
        self.mask     = AtomMask(self.struct)
        self.mask.map_atom('H', lambda x: x.name=='H')
        
        self.covalent = NeighborList(self.mask.atoms[~self.mask['H']])
        
    def seq(self): # lambdas here can probably be abstracted into the polypeptide class, work on it later
        flat     = lambda x: [j for i in x for j in i]
        flatpp   = lambda x: dict((k, flat(x[k])) for k in x)
        seq      = flatpp(self.pp.apply(lambda x: self.resDIC[x.resname], lambda x: x.id))
        seq      = dict((k, ''.join(seq[k])) for k in seq)
        return seq
    
    def dihedrals(self):
        _dihedral = lambda x: dihedral([i.coord for i in x]) if len(x)==4 else np.nan
        _serial   = lambda x: [i.get_serial_number() for i in x] if len(x)==4 else [0]*4
        _info     = lambda x: (*_serial(x), _dihedral(x))
        phi_atoms = lambda x: list(self.covalent.walk(x['C'], path=('C','CA','N','C')))[::-1]
        psi_atoms = lambda x: list(self.covalent.walk(x['N'], path=('N','CA','C','N')))
        fxn       = lambda x, y: lambda z: y(x(z))
        
#         display(
#             self.pp.apply(fxn(psi_atoms, _info), lambda x: x.id)
#         )
        
        display(
            self.pp.apply(fxn(phi_atoms, _info), lambda x: x.id)
        )
        return 



