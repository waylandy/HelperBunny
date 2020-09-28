import os
import sys
from datetime import datetime
import tarfile
import requests
import sqlite3

import numpy as np
import networkx as nx

LIBHB    = '%s/.libhb/' % os.path.expanduser('~')

class NCBItaxdump:
    def __init__(self, path=LIBHB):
        self.path = path
        if not os.path.exists(self.path):
            os.makedirs(self.path)
            
        self.db = '%staxdump.sqlite' % path
        self.update()
        
    def download(self):
        url      = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
        filename = self.path+'taxdump.tar.gz'
        sys.stderr.write('Downloading taxdump.tar.gz from ncbi\n')
        request  = requests.get(url, stream=True)
        with open(filename, 'wb') as f:
            for chunk in request.iter_content(chunk_size=1024):
                f.write(chunk)
                f.flush()
        tar      = tarfile.open(filename, mode='r')
        contents = dict((f,self.path+f) for f in tar.getnames())
        sys.stderr.write('Extracting contents\n')
        tar.extractall(path=self.path)
        
    def make_db(self):
        """ uses the download function """
        
        parse    = lambda x: x.rstrip('\t|\n').split('\t|\t')
        itaxdump = lambda x: (parse(i) for i in open(self.path+x, 'r'))
        
        conn = sqlite3.connect(self.db)
        c    = conn.cursor()
        
        sys.stderr.write('Creating sqlite table: nodes\n')
        c.execute('''CREATE TABLE nodes (
            "tax_id"                         INTEGER,
            "parent tax_id"                  INTEGER,
            "rank"                           TEXT,
            "embl code"                      TEXT,
            "division id"                    TEXT,
            "inherited div flag"             INTEGER,
            "genetic code id"                TEXT,
            "inherited GC flag"              INTEGER,
            "mitochondrial genetic code id"  TEXT,
            "inherited MGC flag"             INTEGER,
            "GenBank hidden flag"            INTEGER,
            "hidden subtree root flag"       TEXT
            )''')
        for row in itaxdump('nodes.dmp'):
            c.execute('INSERT INTO nodes VALUES (?,?,?,?,?,?,?,?,?,?,?,?)', row[:12]) # row 13 contains comments
        c.execute('CREATE UNIQUE INDEX nodes_index on nodes (tax_id)')
        conn.commit()

        sys.stderr.write('Creating sqlite table: names\n')
        c.execute('''CREATE TABLE names (
            "tax_id"                         INTEGER,
            "name_txt"                       TEXT,
            "unique name"                    TEXT,
            "name class"                     TEXT
            )''')
        for row in itaxdump('names.dmp'):
            c.execute('INSERT INTO names VALUES (?,?,?,?)', row)
        conn.commit()
    
    def clean(self):
        sys.stderr.write('Cleaning up extracted files\n')
        for file in ['citations.dmp','delnodes.dmp','division.dmp','gc.prt','gencode.dmp',
                     'merged.dmp','names.dmp','nodes.dmp','readme.txt']:
            os.remove(self.path+file)
            
    def update(self, force=False):
        if os.path.isfile(self.db):
            last   = os.path.getmtime(self.db)
            now    = datetime.now().timestamp()
            days   = (now-last)/86400 # seconds in a day
            update = days > 20
        else:
            update = True

        if update or force:
            if os.path.isfile(self.db):
                os.rename(self.db, self.db+'.bak')
            self.download()
            self.make_db()
            self.clean()
        
    def node(self, tax_id):
        conn = sqlite3.connect(self.db)
        return taxdumpNode(tax_id, conn)

class Species:
    def __init__(self, tax, species):
        self.species = np.array(species, dtype=object)
        self.lineage = np.array([tax.lineage(sp) for sp in species], dtype=object)
        self.lindict = np.array([{} if i==None else dict(j for j in  i if j[0] not in ('clade','no rank')) for i in self.lineage], dtype=object)
    
    def __repr__(self):
        x = len(self.species), len(set(self.species)), sum(self.lineage==None)
        return '<Species: %s sequences in %s species (%s unknown)>' % x
    
    def __getitem__ (self, key):
        # if key is taxonomic rank, returns that rank from all lineages
        # if key is name, returns all lineages containing that name
        rank = set(j[0] for i in self.lineage if i!=None for j in i)
        if key in rank and key not in ('clade','no rank'):
            return np.array([i[key] if key in i else None for i in self.lindict], dtype=object)
        else:
            return np.array([False if i==None else key in (j[1] for j in i) for i in self.lineage], dtype=bool)
    
    def taxa(self, mode='species'):
        # count by species or by sequences
        if mode=='species':
            x = np.unique(self.species, return_index=True)[1]
            x = self.lineage[x]
        elif mode=='sequence':
            x = self.lineage
        else:
            raise Exception('Mode not recognized. Choose either "species" or "sequence"')
            
        c = {}
        for l in x:
            if l == None:
                continue
            for t in l:
                t = tuple(t)
                if t not in c:
                    c[t] = 0
                c[t] += 1
        return sorted((c[i], i) for i in c)[::-1]
    
class TaxonomyDatabase:
    def __init__(self, db=LIBHB):
        self.db = db
        
    def __repr__(self):
        try:
            self.tree
            d = os.path.getmtime(self.db+'taxdump.sqlite')
            d = datetime.fromtimestamp(d).strftime('%Y-%m-%d')
            return '<TaxonomyDatabase: NCBI taxdump (%s)' % d
        except:
            return '<TaxonomyDatabase: no database loaded>'
        
    def load(self):
        NCBItaxdump(path=self.db)
        
        def ifetch(cursor):
            while True:
                f = cursor.fetchone()
                if f == None:
                    break
                yield f
                
        conn = sqlite3.connect(self.db+'/taxdump.sqlite')
        c    = conn.cursor()

        self.tree = nx.DiGraph() # this object requires about 3gb ram 

        c.execute("""SELECT "tax_id", "parent tax_id", "rank" FROM nodes """ )
        for node, parent, rank in ifetch(c):
            self.tree.add_edge(node,parent)
            self.tree.nodes[node]['rank'] = rank
        c.execute("""SELECT "tax_id", "name_txt" FROM names WHERE "name class" = "scientific name" """ )
        for node, name in ifetch(c):
            self.tree.nodes[node]['name'] = name
        self.name2ox = dict((self.tree.nodes[i]['name'],i) for i in self.tree.nodes)
        return self
    
    def lineage(self, ind, format='list', v=False):
        try:
            if ind in self.name2ox:
                ind = self.name2ox[ind]
            else:
                ind = int(ind)

            if format=='list':
                out = []
                for i in nx.dfs_successors(self.tree, ind):
                    out += [[self.tree.nodes[i]['rank'], self.tree.nodes[i]['name']]]
                return out 

            if format=='dict':
                out = {'no rank':[]}
                for i in nx.dfs_successors(self.tree, ind):
                    rank = self.tree.nodes[i]['rank']
                    name = self.tree.nodes[i]['name']
                    if rank=='no rank':
                        out[rank] += [name] 
                    else:
                        out[rank]  = name
                return out

        except:
            if v:
                sys.stderr.write('FAILED QUERY: %s\n' % ind)
            return None

    def taxonomy(self, species):
        return Species(self, species)



