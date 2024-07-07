#!/bin/bash

today=`date +%y%m%d`

if [ ! -d nr$today ]; then
  mkdir nr$today;
fi
cd nr$today

wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

zcat nr.gz | cut -d $'\x01' -f 1 > nr$today.fa
