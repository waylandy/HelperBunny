import networkx as nx


class Taxonomy:
    def __init__(self, nodes, names):
        # nodes = 'nodes.dmp'
        # names = 'names.dmp' 
        self.taxonomy = nx.DiGraph()
        with open(nodes) as r:
            for l in r:
                node, parent, rank = (i.strip() for i in l.split('|', 3)[:3])
                node, parent       = int(node), int(parent)
                self.taxonomy.add_edge(node,parent)
                self.taxonomy.nodes[node]['rank'] = rank
        with open(names) as r:
            for l in r:
                node, name, unique_name, name_class = (i.strip() for i in l.split('|')[:4])
                if name_class != 'scientific name':
                    continue
                self.taxonomy.nodes[int(node)]['name'] = name
        return
    
    def lineage(self, ind):
        out = {'no rank':[]}
        for i in nx.dfs_successors(self.taxonomy, ind):
            rank = self.taxonomy.nodes[i]['rank']
            name = self.taxonomy.nodes[i]['name']
            if rank=='no rank':
                out[rank] += [name] 
            else:
                out[rank]  = name
        return out 


