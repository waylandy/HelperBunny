

CODON = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
         'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
         'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
         'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
         'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
         'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
         'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
         'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
         'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
         'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
         'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
         'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
         'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
         'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
         'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
         'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

RES3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
           'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
           'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
           'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

RES1to3 = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
           'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
           'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
           'A': 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

PHI = ('C','N','CA','C') # n-1, n, n, n
      
PSI = ('N','CA','C','N') # n, n, n, n+1

CHI = {
    'CHI1':{
         'ARG ': ('N' , 'CA', 'CB' , 'CG' ),
         'ASN ': ('N' , 'CA', 'CB' , 'CG' ),
         'ASP ': ('N' , 'CA', 'CB' , 'CG' ),
         'CYS ': ('N' , 'CA', 'CB' , 'SG' ),
         'GLN ': ('N' , 'CA', 'CB' , 'CG' ),
         'GLU ': ('N' , 'CA', 'CB' , 'CG' ),
         'HIS ': ('N' , 'CA', 'CB' , 'CG' ),
         'ILE ': ('N' , 'CA', 'CB' , 'CG1'),
         'LEU ': ('N' , 'CA', 'CB' , 'CG' ),
         'LYS ': ('N' , 'CA', 'CB' , 'CG' ),
         'MET ': ('N' , 'CA', 'CB' , 'CG' ),
         'PHE ': ('N' , 'CA', 'CB' , 'CG' ),
         'PRO ': ('N' , 'CA', 'CB' , 'CG' ),
         'SER ': ('N' , 'CA', 'CB' , 'OG' ),
         'THR ': ('N' , 'CA', 'CB' , 'OG1'),
         'TRP ': ('N' , 'CA', 'CB' , 'CG' ),
         'TYR ': ('N' , 'CA', 'CB' , 'CG' ),
         'VAL ': ('N' , 'CA', 'CB' , 'CG1')
    },
    'CHI2':{
         'ARG ': ('CA', 'CB', 'CG' , 'CD' ),
         'ASN ': ('CA', 'CB', 'CG' , 'OD1'),
         'ASP ': ('CA', 'CB', 'CG' , 'OD1'),
         'GLN ': ('CA', 'CB', 'CG' , 'CD' ),
         'GLU ': ('CA', 'CB', 'CG' , 'CD' ),
         'HIS ': ('CA', 'CB', 'CG' , 'ND1'),
         'ILE ': ('CA', 'CB', 'CG1', 'CD' ),
         'LEU ': ('CA', 'CB', 'CG' , 'CD1'),
         'LYS ': ('CA', 'CB', 'CG' , 'CD' ),
         'MET ': ('CA', 'CB', 'CG' , 'SD' ),
         'PHE ': ('CA', 'CB', 'CG' , 'CD1'),
         'PRO ': ('CA', 'CB', 'CG' , 'CD' ),
         'TRP ': ('CA', 'CB', 'CG' , 'CD1'),
         'TYR ': ('CA', 'CB', 'CG' , 'CD1')
    },
    'CHI3':{
         'ARG ': ('CB', 'CG', 'CD' , 'NE' ),
         'GLN ': ('CB', 'CG', 'CD' , 'OE1'),
         'GLU ': ('CB', 'CG', 'CD' , 'OE1'),
         'LYS ': ('CB', 'CG', 'CD' , 'CE' ),
         'MET ': ('CB', 'CG', 'SD' , 'CE' )
    },
    'CHI4':{
         'ARG ': ('CG', 'CD', 'NE' , 'CZ' ), 
         'LYS ': ('CG', 'CD', 'CE' , 'NZ' )
    },
    'CHI5':
        {'ARG ': ('CD', 'NE', 'CZ' , 'Nh1')
    },
}
#define k_B 0.001982923700 // Boltzmann’s constant in kcal/mol K
#define k_B 0.0083144621 // Boltzmann’s constant kJ/mol-K
