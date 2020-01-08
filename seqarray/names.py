# from .seqarray import *

class CFA_Names(dict):
    """
        Parses a name vector using a given standard format
    """
    def __init__(self, names, format=[]):
        self.names = names
        if 'up' in format:
            self.uniprot()
        if 'mf' in format:
            self.mapfaster()
        self._build()

    def _build(self):
        for k in self:
            self[k] = np.array(self[k], dtype=object)

    def uniprot(self):
        self['id'       ] = []
        self['name'     ] = []
        self['UniProtKB'] = []
        self['OS'       ] = []
        self['OX'       ] = []
        self['GN'       ] = []
        for name in self.names:
            try:
                fr  = []
                it  = iter(name.split())
                fr += next(it).split('|')+['']
                i   = next(it)
                while '=' not in i:
                    fr[3]+=i+' '
                    i     = next(it)
            except:
                fr  = ['','','','']
            self['id'       ] += [fr[1]]
            self['name'     ] += [fr[3]]
            self['UniProtKB'] += ['|'.join(fr[:3])]
            try:
                fr = []
                if i.startswith('OS='):
                    i     = i[3:]
                    fr   += ['']
                    while '=' not in i:
                        fr[-1]+= i+' '
                        i     = next(it)
                    fr[-1] = fr[-1][:-1]
                else:
                    fr += ''
                fr += [i[3:] if i.startswith('OX=') else '']
                i   = next(it)
                fr += [i[3:] if i.startswith('GN=') else '']
            except:
                fr  = ['','','','']
            self['OS'] += [fr[0]]
            self['OX'] += [fr[1]]
            self['GN'] += [fr[2]]

    def mapfaster(self):
        self['profile'] = []
        for name in self.names:
            assert ' profile:' in name
            self['profile'] += [name.split(' profile:')[-1].split()[0]]
