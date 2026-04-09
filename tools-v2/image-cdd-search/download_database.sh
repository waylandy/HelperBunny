#!/bin/bash

mkdir database
cd database

wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/README

mkdir ./db
wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz && tar -xzf Cdd_LE.tar.gz -C ./db && rm -f Cdd_LE.tar.gz
wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_NCBI_LE && tar -xzf Cdd_NCBI_LE -C ./db && rm -f Cdd_NCBI_LE
wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cog_LE.tar.gz && tar -xzf Cog_LE.tar.gz -C ./db && rm -f Cog_LE.tar.gz
wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Kog_LE.tar.gz && tar -xzf Kog_LE.tar.gz -C ./db && rm -f Kog_LE.tar.gz
wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Prk_LE.tar.gz && tar -xzf Prk_LE.tar.gz -C ./db && rm -f Prk_LE.tar.gz
wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Smart_LE.tar.gz && tar -xzf Smart_LE.tar.gz -C ./db && rm -f Smart_LE.tar.gz
wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Tigr_LE.tar.gz && tar -xzf Tigr_LE.tar.gz -C ./db && rm -f Tigr_LE.tar.gz

mkdir ./data
wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.info -O ./data/cdd.info
wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz -O ./data/cddid.tbl.gz && gzip -d ./data/cddid.tbl.gz
wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdtrack.txt -O ./data/cdtrack.txt
wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links -O ./data/family_superfamily_links
wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot.dat.gz -O ./data/cddannot.dat.gz && gzip -d ./data/cddannot.dat.gz
wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot_generic.dat.gz -O ./data/cddannot_generic.dat.gz && gzip -d ./data/cddannot_generic.dat.gz
wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/bitscore_specific.txt -O ./data/bitscore_specific.txt

cd ..

