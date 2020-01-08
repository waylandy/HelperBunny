import pandas as pd

class CDD:
    def __init__(self, cdd_file):
        self.cdd_file = cdd_file
        self._build()
        return

    def __getitem__(self, key):
        return self.df.xs(key, level='Query')

    def _build(self):
        it, i = iter(open(self.cdd_file)), -1
        while not next(it).startswith('Query'): i+=1
        self.df = pd.read_csv(self.cdd_file, sep='\t', skiprows=i)
        self.df['mid']   =  (self.df['From'].values + self.df['To'].values)/2
        index, query     = list(zip(*[i.split('#')[1].split(' - ') for i in self.df['Query']]))
        self.df['Query'] = query
        self.df['Index'] = list(map(int, index))
        self.df.sort_values(by=['Index', 'mid'], ascending=[1, 1], inplace=True)
        self.df.set_index(['Query', 'Accession'], inplace=True)
        self.df = self.df[['From', 'To', 'Short name', 'Superfamily', 'Hit type', 'PSSM-ID', 'E-Value', 'Incomplete']]
