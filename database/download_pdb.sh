#!/bin/bash

# download structures:
rsync -rlpt -v -z --delete --port=33444 \
rsync.rcsb.org::ftp_data/structures/divided/pdb/ ./pdb

# biounits:
# ftp://ftp.wwpdb.org/pub/pdb/data/biounit/PDB/divided

#############################################################

today=`date +%y%m%d`

rm metadata/pdb_seqres.txt
wget 'ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt'
mv pdb_seqres.txt metadata/pdb_seqres.txt

rm metadata/entries.idx
wget 'ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx'
mv entries.idx metadata/entries.idx

rm metadata/pdbtosp.txt
wget 'https://www.uniprot.org/docs/pdbtosp.txt'
mv pdbtosp.txt metadata/pdbtosp.txt
