#!/bin/bash

up='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/'
today=`date +%y%m%d`

if [ ! -d up$today ]; then
  mkdir up$today; 
fi
cd up$today

wget ${up}Reference_Proteomes_*.tar.gz
ref=`echo Reference_Proteomes_*.tar.gz`
tar -xvf $ref --wildcards --no-anchored README '*[^_DNA].fasta.gz'
#rm $ref

#zcat Eukaryota/*.fasta.gz > up_euk_$today.fa
#zcat Archaea/*.fasta.gz   > up_arc_$today.fa
#zcat Bacteria/*.fasta.gz  > up_bac_$today.fa
#zcat Viruses/*.fasta.gz   > up_vir_$today.fa

zcat Eukaryota/*/*.fasta.gz > up_euk_$today.fa
zcat Archaea/*/*.fasta.gz   > up_arc_$today.fa
zcat Bacteria/*/*.fasta.gz  > up_bac_$today.fa
zcat Viruses/*/*.fasta.gz   > up_vir_$today.fa

# rm {Archaea,Bacteria,Eukaryota,Viruses} -r

# makeblastdb -in uniprot_sprot.fasta -out uniprot_sprot -type prot -parse_seqids
