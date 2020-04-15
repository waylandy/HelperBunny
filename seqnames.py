import numpy as np

class SeqNames(dict):
    def __init__(self, names, tags=['up', 'rf'], query=[]):
        self.names = names
        if 'up' in tags:
            self.parse_uniprot()
        if 'rf' in tags:
            self.parse_profiles()
        for q in query:
            self.parse_query(q)
    
    def map_names(self, fxn=None):
        def septag(name, div='=:'):
            div  = frozenset(i for i in div)
            out  = [[]]
            for i in name.split(' '):
                if not any(s in i for s in div):
                    out[-1] += [i]
                else:
                    out     += [[i]]
            return [' '.join(i) for i in out]
        fxn = septag if fxn==None else fxn
        for name in self.names:
            yield fxn(name)
    
    def parse_uniprot(self):
        for tag in ['id','name','OS','OX','GN','PE','SV']:
            self[tag] = []
        split = lambda x: x.split(' ',1) if ' ' in x else [x,'']
        for name in self.map_names():
            head = split(name.pop(0))
            self['id'  ] += [head.pop(0)]
            self['name'] += [head.pop(0)]
            for tag in ['OS','OX','GN','PE','SV']:
                try:
                    self[tag] += [next(filter(lambda x: x.startswith(tag+'='), name))[3:]]
                except:
                    self[tag] += ['']
                    
        self['accession'] = np.array([x.split('|')[1] for x in  self['id']], dtype=object)
        for tag in ['id','name','OS','OX','GN','PE','SV']:
            self[tag] = np.array(self[tag], dtype=object)
        
    def parse_profiles(self):
        self['profile'] = []
        for name in self.map_names():
            try:
                name = tuple(i[8:] for i in filter(lambda x: x.startswith('profile='), name))
                self['profile'] += [name[-1]]
            except:
                self['profile'] += [()]
        self['profile'] = np.array(self['profile'], dtype=object)

    def parse_query(self, prefix):
        self[prefix] = []
        for name in self.map_names():
            try:
                lp   = len(prefix)
                name = tuple(i[lp:] for i in filter(lambda x: x.startswith(prefix), name))
                self[prefix] += [name[-1]]
            except:
                self[prefix] += [()]
        self[prefix] = np.array(self[prefix], dtype=object)
