import os
import sys
from datetime import datetime
import tarfile
import requests
import sqlite3

import numpy as np
import networkx as nx

LIBHB    = '%s/.libhb/' % os.path.expanduser('~')

"""
taxdump = TaxonomyDatabase()
q = np.array(['Boops boops',942147,'environmental samples','1964143',None,''])
t = taxdump.get_taxonomy(q)
"""

class TaxonomyDatabase:
    def __init__(self, path=LIBHB, update=28):
        self.path = path
        self.db_file = f'{path}/taxdump.sqlite'
        self.update = update
        self.db_date = 'None'

    def __repr__(self):
        return f'< TaxonomyDatabase: db_date={self.db_date} >'

    @staticmethod
    def download(dest_file):
        url  = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
        request  = requests.get(url, stream=True)
        with open(dest_file, 'wb') as w:
            for chunk in request.iter_content(chunk_size=1024):
                w.write(chunk)
                w.flush()

    @staticmethod
    def build(taxdump_file, taxdump_db, update=0):
        if os.path.isfile(taxdump_db):
            last_written = os.path.getmtime(taxdump_db)
            now = datetime.now().timestamp()
            delta = int((now - last_written) / 86400) # seconds / day
            if delta >= update:
                os.rename(taxdump_db, taxdump_db+'_old')
            else:
                return
            
        sys.stderr.write('Downloading taxdump\n')
        TaxonomyDatabase.download(taxdump_file)

        conn = sqlite3.connect(taxdump_db)
        cur = conn.cursor()

        tar = tarfile.open(taxdump_file, "r:gz")
        for member in tar.getmembers():
            handle = tar.extractfile(member)
            
            if member.name == 'nodes.dmp': # nodes.dmp (load first 3 cols)
                cur.execute('''CREATE TABLE nodes (
                    "tax_id"         INTEGER,
                    "parent tax_id"  INTEGER,
                    "rank"           TEXT
                    )''')
                for *row, _ in map(lambda x: x.decode().split('\t|\t', 3), handle):
                    cur.execute('INSERT INTO nodes VALUES (?,?,?)', row)
                conn.commit()
                sys.stderr.write('| creating table : nodes\n')
            
            if member.name == 'names.dmp': # names.dmp (load first 4 cols, keep first 2)
                cur.execute('''CREATE TABLE names (
                    "tax_id"         INTEGER,
                    "name_txt"       TEXT
                    )''')
                for row in map(lambda x: x.decode().rstrip('\t|\n').split('\t|\t'), handle):
                    if row[3] == 'scientific name':
                        cur.execute('INSERT INTO names VALUES (?,?)', row[:2])
                conn.commit()
                sys.stderr.write('| creating table : names\n')
            
            if member.name == 'merged.dmp': # merged.dmp (load all cols)
                cur.execute('''CREATE TABLE merged (
                    "old_tax_id"     INTEGER,
                    "new_tax_id"     INTEGER
                    )''')
                for row in map(lambda x: x.decode().rstrip('\t|\n').split('\t|\t'), handle):
                    cur.execute('INSERT INTO merged VALUES (?,?)', row[:2])
                conn.commit()
                sys.stderr.write('| creating table : merged\n')

        conn.close()

    def connect(self):
        taxdump_file = f'{self.path}/taxdump.tar.gz'
        TaxonomyDatabase.build(taxdump_file, self.db_file, update=self.update)
        self.db_date = datetime.fromtimestamp(os.path.getmtime(taxdump_file)).strftime('%m-%d-%Y')
        return sqlite3.connect(self.db_file)
    
    def get_taxonomy(self, query):
        conn = self.connect()
        cur = conn.cursor()

        # convert old to new ids
        cur.execute("""SELECT "old_tax_id", "new_tax_id" FROM merged""")
        merged = dict(cur.fetchall())

        # convert names to ids (if unique)
        cur.execute("""SELECT "name_txt", "tax_id" FROM names GROUP BY "name_txt" HAVING COUNT("tax_id") == 1""")
        name2id = dict(cur.fetchall())

        # convert ids to names
        cur.execute("""SELECT "tax_id", "name_txt" FROM names""")
        id2name = dict(cur.fetchall())
        
        # get query node ids
        size = len(query)
        query_ids = np.zeros(size, dtype=int)
        for n, i in enumerate(query):
            if type(i) == str:
                if i.isnumeric():
                    i = int(i)
                elif i in name2id:
                    query_ids[n] = name2id[i]
            if type(i) == int:
                if i in merged:
                    i = merged[i]
                if i in id2name:
                    query_ids[n] = i
                    
        del merged
        del name2id

        # load nodes
        cur.execute("""SELECT "tax_id", "parent tax_id", "rank" FROM nodes """)
        node_arr = dict(zip(["tax_id", "parent tax_id", "rank"], zip(*cur.fetchall())))
        node_arr["tax_id"       ] = np.array(node_arr["tax_id"       ], dtype=int)
        node_arr["parent tax_id"] = np.array(node_arr["parent tax_id"], dtype=int)
        node_arr["rank"         ] = np.array(node_arr["rank"         ], dtype=object)
        
        assert 0 not in node_arr['tax_id']
        assert 0 not in node_arr['parent tax_id']

        # mask ancestors
        n_nodes = 0
        mask = np.isin(node_arr['tax_id'], query_ids, assume_unique=False)
        while n_nodes != mask.sum():
            # print(f'|\no {n_nodes}')
            n_nodes = mask.sum()
            ancestors = node_arr["parent tax_id"][mask]
            mask += np.isin(node_arr['tax_id'], ancestors, assume_unique=False)

        node_arr["tax_id"       ] = node_arr["tax_id"       ][mask]
        node_arr["parent tax_id"] = node_arr["parent tax_id"][mask]
        id2rank = dict(zip(node_arr["tax_id"], node_arr["rank"][mask]))
        
        # build tree
        cladogram = nx.DiGraph()
        cladogram.add_edges_from(np.array([node_arr["tax_id"], node_arr["parent tax_id"]]).T)
        for i in cladogram.nodes():
            cladogram.nodes[i]['name'] = id2name[i]
            cladogram.nodes[i]['rank'] = id2rank[i]
            
        del node_arr
        del id2name
        del id2rank
            
        return Taxonomy(query_ids, cladogram)

class Taxonomy:
    def __init__(self, query_ids, cladogram):
        self.query_ids = query_ids
        self.cladogram = cladogram

        _unique_ids = {}
        for n, i in enumerate(self.query_ids):
            if i not in _unique_ids:
                _unique_ids[i] = []
            _unique_ids[i] += [n]
        
        self._unique_ids = np.array(list(_unique_ids.keys()), dtype=int)
        self._unique_indices = np.array(list(_unique_ids.values()), dtype=object)
        
        self._unique_ranks = []
        self._unique_lineages = []
        
        self._ranks = set()
        for i in self._unique_ids:
            ids = self.get_successors(i)[::-1]
            ranks = {}
            lineage = []

            for i in ids:
                node = self.cladogram.nodes[i]
                rank, name = node['rank'], node['name']
                if rank not in ('clade','no rank'):
                    self._ranks.add(rank)
                    ranks[rank] = name
                lineage += [name]
                
            self._unique_ranks += [ranks]
            self._unique_lineages += [tuple(lineage)]
            
        self._unique_ranks = np.array(self._unique_ranks, dtype=object)
        self._unique_lineages = np.array(self._unique_lineages, dtype=object)
    
    def __repr__(self):
        n_taxa = len(set(filter(lambda x: x!=0, self.query_ids)))
        n_nodes = self.cladogram.number_of_nodes()
        n_seqs = self.query_ids.size
        n_unmapped = (self.query_ids == 0).sum()
        return f'<Taxonomy: {n_taxa} taxa; {n_nodes} nodes; {n_seqs} queries; {n_unmapped} unmapped>'
    
    def __getitem__ (self, key):
        if key in self._ranks:
            arr = np.full(self.query_ids.size, '', dtype=object)
            for ndx, rank in zip(self._unique_indices, self._unique_ranks):
                if key in rank:
                    arr[ndx] = rank[key]
            return arr
        else:
            arr = np.full(self.query_ids.size, False)
            for ndx, lineage in zip(self._unique_indices, self._unique_lineages):
                arr[ndx] = key in lineage
            return arr
    
    def get_successors(self, query_id):
        if query_id != 0:
            return [query_id] + [i[0] for i in nx.dfs_successors(self.cladogram, query_id).values()]
        else:
            return []
    
    def count_rank(self, rank):
        if rank in self._ranks:
            ranks, counts = np.unique(self[rank], return_counts=True)
            sort = np.argsort(counts)[::-1]
            return np.array([counts[sort], ranks[sort]], dtype=object).T



