import os
import sys
from datetime import datetime
import tarfile
import requests
import sqlite3

import networkx as nx

#tax  = taxdump().tree()

LIBHB    = '%s/.libhb/' % os.path.expanduser('~')

class taxdump:
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
            c.execute('INSERT INTO nodes VALUES (?,?,?,?,?,?,?,?,?,?,?,?)', row)
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
            update = True if days>20 else False
        else:
            update = True

        if update or force:
            self.download()
            self.make_db()
            self.clean()
        
    def node(self, tax_id):
        conn = sqlite3.connect(self.db)
        return taxdumpNode(tax_id, conn)
    
    def tree(self):
        return taxdumpTree(self.db)

class taxdumpNode(dict):
    def __init__(self, node, conn):
        self.conn      = conn
        self.c         = conn.cursor()
        self['tax_id'] = node
        self.get_node()
        self.get_name()
    
    def lineage(self):
        yield self
        query = self.parent()
        while True:
            if query['tax_id']==1:
                return
            else:
                query = query.parent()
                yield query
    
    def parent(self):
        return TaxDumpNode(self['parent'], conn=self.conn)
    
    def children(self):
        self.c.execute("""SELECT tax_id
                          FROM nodes 
                          WHERE "parent tax_id"=?
                       """, [self['tax_id']])
        for result in self.c.fetchall():
            yield TaxDumpNode(result[0], conn=self.conn)
    
    def get_node(self):
        self.c.execute("""SELECT "parent tax_id", "rank" 
                          FROM nodes 
                          WHERE tax_id=?
                       """, [self['tax_id']])
        self['parent'], self['rank'] = self.c.fetchone()
    
    def get_name(self):
        self.c.execute("""SELECT "name_txt"
                          FROM names 
                          WHERE tax_id=? AND "name class"="scientific name"
                       """, [self['tax_id']])
        self['name'] = self.c.fetchone()[0]

class taxdumpTree:
    def __init__(self, db):
        def ifetch(cursor):
            while True:
                f = cursor.fetchone()
                if f == None:
                    break
                yield f
                
        conn = sqlite3.connect(db)
        c    = conn.cursor()

        self.tree = nx.DiGraph() # this object requires about 3gb ram 
        c.execute("""SELECT "tax_id", "parent tax_id", "rank" FROM nodes """ )
        for node, parent, rank in ifetch(c):
            self.tree.add_edge(node,parent)
            self.tree.nodes[node]['rank'] = rank
        c.execute("""SELECT "tax_id", "name_txt" FROM names WHERE "name class" = "scientific name" """ )
        for node, name in ifetch(c):
            self.tree.nodes[node]['name'] = name

    def lineage(self, ind, format='list'):
        try:
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
            sys.stderr.write('FAILED QUERY: %s\n' % ind)
            return None



