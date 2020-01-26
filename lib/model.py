from modeller import *
from modeller.automodel import *

# align.ali 
"""
>P1;name1
structureX:file.pdb:88:A:92:A::::
PSRIA*
>P1;name2
sequence:name::::::::
MA-AG*
"""

env = environ()
for model in ['name2']:
    hm = automodel(env, alnfile='align.ali',
                   knowns='name1', sequence=model,
                   assess_methods=(assess.DOPE,assess.GA341))
    hm.starting_model = 1
    hm.ending_model   = 3
    hm.make()

