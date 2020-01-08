import subprocess

"""
This definately need to be abstracted later...

"""



def tmalign(pdb_a, pdb_b, bin='TMalign'):
    """
    Captures the output of TM-align; tested on TM-align (Version 20170708)
    Link: https://zhanglab.ccmb.med.umich.edu/TM-align/
    Use the fortran version!
    """

    result={}
    aln=False
    for l in subprocess.check_output([bin, pdb_a, pdb_b]).split('\n'):
        if   'Length of Chain_1' in l:
            result['Length of Chain_1'] = int(l.split()[3])
        elif 'Length of Chain_2' in l:
            result['Length of Chain_2'] = int(l.split()[3])
        elif 'Aligned length='   in l:
            for a in l.split(', '):
                a = a.split('=')
                if   'Alig'   in a[0]:
                    result['Aligned length'] = a[1]
                elif 'RMSD'   in a[0]:
                    result['RMSD'] = a[1]
                elif 'Seq_ID' in a[0]:
                    result['Seq_ID'] = a[2]
        elif 'normalized by length of Chain_1)' in l:
            result['TM-score1'] = float(l.split()[1])
        elif 'normalized by length of Chain_2)' in l:
            result['TM-score2'] = float(l.split()[1])
        elif '(":" denotes' in l:
            aln=1
        elif   aln==1:
            result['alignment']=[l]
            aln+=1
        elif aln==2:
            result['alignment'].append(l)
            aln+=1
        elif aln==3:
            result['alignment'].append(l)
            aln+=1
    return result
