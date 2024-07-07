import os
import sys
import tarfile
import sqlite3
import requests
from datetime import datetime
from itertools import groupby

import numpy as np
import networkx as nx

LIBHB = '%s/.libhb/' % os.path.expanduser('~')

"""
taxdump = TaxonomyDatabase()
q = np.array(['Boops boops',942147,'environmental samples','1964143',None,''])
t = taxdump.get_taxonomy(q)
"""

class TaxonomyDatabase:
    def __init__(self, path=LIBHB, update=1):
        self.taxdump_gz = f'{path}/taxdump.tar.gz'
        self.taxdump_db = f'{path}/taxdump.db'

        if not os.path.isdir(path):
            os.makedirs(path)
        
        if os.path.isfile(self.taxdump_db):
            last_written = os.path.getmtime(self.taxdump_db)
            elapsed = (datetime.now().timestamp() - last_written) // 86400 # seconds / day
            if (elapsed != None) and (elapsed >= update):
                os.rename(self.taxdump_db, f'{self.taxdump_db}_old')

        if not os.path.exists(self.taxdump_db):
            self.download()
            self.build()
        
    def download(self):
        url  = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
        request  = requests.get(url, stream=True)
        sys.stderr.write('downloading taxdump\n')
        with open(self.taxdump_gz, 'wb') as w:
            for chunk in request.iter_content(chunk_size=1024):
                w.write(chunk)
                w.flush()
    
    def build(self):
        conn = sqlite3.connect(self.taxdump_db)
        cur = conn.cursor()

        tar = tarfile.open(self.taxdump_gz, "r:gz")
        for member in tar.getmembers():
            handle = tar.extractfile(member)
            
            if member.name == 'nodes.dmp': # nodes.dmp
                cur.execute('''CREATE TABLE nodes (
                    "tax_id"         INTEGER,
                    "parent tax_id"  INTEGER,
                    "rank"           TEXT
                    )''')
                for *row, _ in map(lambda x: x.decode().split('\t|\t', 3), handle):
                    cur.execute('INSERT INTO nodes VALUES (?,?,?)', row)
                conn.commit()
                sys.stderr.write('| creating table : nodes\n')
            
            if member.name == 'merged.dmp': # merged.dmp
                cur.execute('''CREATE TABLE merged (
                    "old_tax_id"     INTEGER,
                    "new_tax_id"     INTEGER
                    )''')
                for row in map(lambda x: x.decode().rstrip('\t|\n').split('\t|\t'), handle):
                    cur.execute('INSERT INTO merged VALUES (?,?)', row[:2])
                conn.commit()
                sys.stderr.write('| creating table : merged\n')
            
            if member.name == 'names.dmp': # names.dmp
                cur.execute('''CREATE TABLE names (
                    "tax_id"         INTEGER,
                    "name_txt"       TEXT,
                    "name_class"     TEXT
                    )''')
                for row in map(lambda x: x.decode().rstrip('\t|\n').split('\t|\t'), handle):
                    cur.execute('INSERT INTO names VALUES (?,?,?)', (*row[:2], row[3]))
                sys.stderr.write('| creating table : names\n')
                conn.commit()
        
        conn.close()

    def get_taxonomy(self, query):
        conn = sqlite3.connect(self.taxdump_db)
        cur = conn.cursor()
        
        # old to new ids
        cur.execute("""SELECT "old_tax_id", "new_tax_id" FROM merged""")
        merged = dict(cur.fetchall())
        
        # names to ids (if unique)
        name2id = {}
        cur.execute("""SELECT "name_txt", "tax_id" FROM names WHERE "name_class" == 'scientific name' GROUP BY "name_txt" HAVING COUNT("tax_id") == 1""")
        name2id['scientific'] = dict(cur.fetchall())
        cur.execute("""SELECT "name_txt", "tax_id" FROM names WHERE "name_class" == 'synonym' GROUP BY "name_txt" HAVING COUNT("tax_id") == 1""")
        name2id['synonym'] = dict(cur.fetchall())
        # cur.execute("""SELECT "name_txt", "tax_id" FROM names WHERE "name_class" == 'authority' GROUP BY "name_txt" HAVING COUNT("tax_id") == 1""")
        # name2id['authority'] = dict(cur.fetchall())
        
        # ids to names
        cur.execute("""SELECT "tax_id", "name_txt" FROM names WHERE "name_class" == 'scientific name'""")
        id2name = dict(cur.fetchall())
        
        taxids = np.zeros(len(query), dtype=int)
        for n, i in enumerate(query.tolist() if type(query) == np.ndarray else query):
            if type(i) == str:
                if i.strip().isdigit():
                    i = int(i)
                elif i in name2id['scientific']:
                    i = name2id['scientific'][i]
                elif i in name2id['synonym']:
                    i = name2id['synonym'][i]
                # elif i in name2id['authority']:
                #     i = name2id['authority'][i]
            if type(i) == int:
                if i in merged:
                    i = merged[i]
                if i in id2name:
                    taxids[n] = i
                
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
        mask = np.isin(node_arr['tax_id'], taxids, assume_unique=False)
        while n_nodes != mask.sum():
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

        return Taxonomy(taxids, cladogram)

class Taxonomy:
    def __init__(self, taxids, cladogram):
        self.taxids = taxids
        self.cladogram = cladogram

        groups = groupby(sorted(i[::-1] for i in enumerate(self.taxids)), lambda x: x[0])
        _unique_ids = {k: [i[1] for i in g] for k, g in groups}
        
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
        n_taxa = len(set(filter(lambda x: x!=0, self.taxids)))
        n_nodes = self.cladogram.number_of_nodes()
        n_seqs = self.taxids.size
        n_unmapped = (self.taxids == 0).sum()
        return f'<Taxonomy: {n_taxa} taxa; {n_nodes} nodes; {n_seqs} queries; {n_unmapped} unmapped>'
    
    def __getitem__ (self, key):
        if key in self._ranks:
            arr = np.full(self.taxids.size, '', dtype=object)
            for ndx, rank in zip(self._unique_indices, self._unique_ranks):
                if key in rank:
                    arr[ndx] = rank[key]
            return arr
        else:
            arr = np.full(self.taxids.size, False)
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

