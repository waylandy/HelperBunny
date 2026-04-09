import os
import shutil
import random
import datetime
from subprocess import call, PIPE, STDOUT, DEVNULL

from Bio.Seq import Seq

"""
# https://manpages.ubuntu.com/manpages/jammy/man1/exonerate.1.html

# git clone https://github.com/nathanweeks/exonerate
# cd exonerate && git checkout v2.4.0 && ./configure
# make && make check && sudo make install

import exonerate

sequence1 = "GAVTTFVALYEYESRSETDLSFKKGERLQIVNNTEGDWWLAHSLSTGQTGYPSNYVAPSDS"
sequence2 = (
    "gcgggcccgctggccGGTggagtgaccacctttgtggccctctatgactatgagtctaggacggagacagaccgtc"
    "ctcaagaaaggcgagcggctccagattgtcaacaacactaggaagtggatgtcagccagacctggttcacattcat"
    "tgaatgtgaaccggtctggctgacatccaccttcctcgtgttgttgacaatctggagccgctcgcctttctatgag"
    "gacaggtctgtctccgtcctagactcatagtcatagagggccacaaaggtggtcactccACCggccagcgggcccg"
    )

output = exonerate.get_pairwise_alignment(sequence1, sequence2)
display(output)

"""

def get_pairwise_alignment(
    sequence1, 
    sequence2,
    cores = 1,
    model = "protein2dna",
    exonerate = "exonerate",
    ):

    ### staging
    get_random_string = lambda x: ''.join(random.choice("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") for _ in range(x))
    get_datetime_string = lambda: datetime.datetime.now().strftime("%y%m%d-%H%M%S")
    temp_name = f'{get_datetime_string()}-{get_random_string(8)}'

    os.makedirs(temp_name, exist_ok=True)
    query_fasta = f"{temp_name}/query.fa"
    target_fasta = f"{temp_name}/target.fa"
    output_map = f"{temp_name}/output.map"

    with open(query_fasta, "w") as w:
        w.write(f">query\n{sequence1}\n")
    with open(target_fasta, "w") as w:
        w.write(f">target\n{sequence2}\n")

    ### subprocess
    params = {
        "--query"              : query_fasta,
        "--target"             : target_fasta,
        "--showalignment"      : "False",
        "--showsugar"          : "False",
        "--showcigar"          : "False",
        "--showvulgar"         : "False",
        "--showquerygff"       : "False",
        "--showtargetgff"      : "False",
        "--ryo"                : "\nSTART\n%S %qt %tt\n{%Pqs,%Pqb,%Pqe,%Pts,%Ptb,%Pte,%Pl\n}END\n",
        "--model"              : model,
        "--alignmentwidth"     : "9999",
        "--forwardcoordinates" : "True",
        "--cores"              : str(cores),
        }
    cmd = [exonerate] + [j for i in params.items() for j in i]
    with open(output_map, 'wb') as wb:
        call(cmd, stdout=wb, stderr=DEVNULL)
    
    ### parse output : map
    state, output = 0, []
    with open(output_map, "r") as r:
        for line in map(lambda x: x.rstrip(), r):
            if line == "START":
                state = 1
            elif state == 1:
                state = 2
                sugar, transitions = [int(i) if i.isnumeric() else i for i in line.split(" ")], []
            elif line == "END":
                state = 1

                qi, qs, qe, qS, ti, ts, te, tS, S, qt, tt = sugar
                start1, end1 = 1 + min(qs, qe), max(qs, qe)
                start2, end2 = 1 + min(ts, te), max(ts, te)
                
                aligned1, aligned2 = [], []
                for Pqs, Pqb, Pqe, Pts, Ptb, Pte, Pl in transitions:
                    if Pqs != Pts:
                        aligned1 += [f" {Pqs} " if qt!="DNA" else f"{Pqs:-<3}"]
                        aligned2 += [f" {Pts} " if tt!="DNA" else f"{Pts:-<3}"]
                
                segment1 = "".join(aligned1).translate(str.maketrans("","","- "))
                assert sequence1[start1-1: end1] == Seq(segment1).reverse_complement() if qS == "-" else segment1
                segment2 = "".join(aligned2).translate(str.maketrans("","","- "))
                assert sequence2[start2-1: end2] == Seq(segment2).reverse_complement() if tS == "-" else segment2
                
                output += [{
                    "sequence1" : sequence1,
                    "sequence2" : sequence2,
                    "alignment1": ' '.join(aligned1),
                    "alignment2": ' '.join(aligned2),
                    "start1"    : start1,
                    "end1"      : end1,
                    "strand1"   : qS,
                    "type1"     : qt.lower(),
                    "start2"    : start2,
                    "end2"      : end2,
                    "strand2"   : tS,
                    "type2"     : tt.lower(),
                    "score"     : S,
                }]
            
            elif state == 2:
                transitions += [[int(i) if i.isnumeric() else i for i in line.split(",")]]
    shutil.rmtree(temp_name)
    
    ### format output : row oriented
    return output

