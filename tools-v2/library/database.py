import os
import sys
import gzip
import tarfile
import datetime
import requests
from itertools import groupby

"""
import datetime
import database

get_date_string = lambda: datetime.datetime.now().strftime("%y%m%d")

database.download_nr(f"db/nr.{get_date_string()}")
database.download_uniref100(f"db/uniref100.{get_date_string()}")
database.download_reference_proteomes(f"db/reference_proteomes.{get_date_string()}")
database.download_mgy_proteins(f"db/mgy_proteins.{get_date_string()}")
database.download_pdb_seqres(f"db/pdb_seqres.{get_date_string()}")

"""

######### utils

def get_date_string():
    return datetime.datetime.now().strftime("%y%m%d")

def read_url(url):
    response = requests.get(url, stream=False)
    assert 200 == response.status_code
    return response.content

def download_url(url, filename, buffer_size=2048):
    response = requests.get(url, stream=True)
    assert 200 == response.status_code
    total_size = int(response.headers.get("content-length", 0))
    downloaded_size = 0
    with open(filename, "wb") as wb:
        for buffer in response.iter_content(buffer_size):
            wb.write(buffer)
            downloaded_size += len(buffer)
            sys.stdout.write(f'{100 * (downloaded_size / total_size):.2f}%\r')
    assert downloaded_size == total_size
    sys.stdout.write("\n")

def read_fasta(file):
    is_header = lambda x: x.startswith('>')
    compress  = lambda x: ''.join(_.strip() for _ in x)
    reader    = iter(groupby(open(file), is_header))
    reader    = iter(groupby(open(file), is_header)) if next(reader)[0] else reader
    for key, group in reader:
        if key:
            for header in group:
                header = header[1:].strip()
        else:
            sequence = compress(group)
            if sequence != '':
                yield header, sequence

######### database : nr

def download_nr(dirname):

    os.makedirs(dirname, exist_ok=True)
    gz_file = f"{dirname}/nr.gz"
    download_url(f"https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz", gz_file, buffer_size=2048)

    fasta_file = f"{dirname}/nr.fasta"
    n_sequences = 0
    with open(fasta_file, 'wb') as wb:
        for line in gzip.open(gz_file, 'rb'):
            wb.write(line)
            if line.startswith(b'>'):
                n_sequences += 1
            if n_sequences % 1000 == 0:
                sys.stdout.write(f'n={n_sequences}\r')
    sys.stdout.write(f'n={n_sequences}\n')
    return fasta_file

# def parse_header_nr(header):
#     accession, header = header.split(' ', 1)
#     if ' [' not in header:
#         return {'accession': accession, 'description': header, 'species': ''}
#     else:
#         description, header = header.split(' [', 1)
#         species = header[::-1].split(']', 1)[1][::-1]
#         return {'accession': accession, 'description': description, 'species': species}

######### database : uniref100

def download_uniref100(dirname):

    os.makedirs(dirname, exist_ok=True)
    gz_file = f"{dirname}/uniref100.fasta.gz"
    download_url(f"https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz", gz_file, buffer_size=2048)

    fasta_file = f"{dirname}/uniref100.fasta"
    n_sequences = 0
    with open(fasta_file, 'wb') as wb:
        for line in gzip.open(gz_file, 'rb'):
            wb.write(line)
            if line.startswith(b'>'):
                n_sequences += 1
            if n_sequences % 1000 == 0:
                sys.stdout.write(f'n={n_sequences}\r')
    sys.stdout.write(f'n={n_sequences}\n')
    return fasta_file

# def parse_header_reference_proteomes(header):
#     edges = [0] + sorted(header.index(i) for i in (' OS=', ' OX=', ' GN=', ' PE=', ' SV=') if i in header) + [len(header)]
#     header, *annotations = [header[i[0]:i[1]].strip() for i in zip(edges[:-1], edges[1:])]
#     header, proteinname = header.split(' ', 1)
#     db, uniqueidentifier, entryname = header.split('|')
#     output = {'db': db, 'UniqueID': uniqueidentifier, 'EntryName': entryname, 'ProteinName': proteinname}
#     output.update(dict(i.split('=', 1) for i in annotations))
#     return output

######### database : reference_proteomes

def download_reference_proteomes(dirname):
    
    os.makedirs(dirname, exist_ok=True)
    url_base = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/"
    version_date = read_url(url_base).decode().split('Reference_Proteomes_', 1)[1].split('.tar.gz', 1)[0]
    tar_gz_file = f"{dirname}/Reference_Proteomes_{version_date}.tar.gz"
    download_url(f"{url_base}/Reference_Proteomes_{version_date}.tar.gz", tar_gz_file, buffer_size=2048)

    tar = tarfile.open(tar_gz_file, mode='r:gz')
    handle = filter(lambda x: ('_DNA' not in x.name) and x.name.endswith('.fasta.gz'), tar)
    fasta_file = f"{dirname}/Reference_Proteomes_{version_date}.fasta"
    n_species, n_sequences = 0, 0
    with open(fasta_file, 'wb') as wb:
        for tarinfo in handle:
            n_species += 1
            for line in gzip.open(tar.extractfile(tarinfo), 'rb'):
                if line.startswith(b'>'):
                    n_sequences += 1
                wb.write(line)
            sys.stdout.write(f'sp={n_species}; n={n_sequences}\r')
    sys.stdout.write(f'sp={n_species}; n={n_sequences}\n')
    return fasta_file

######### database : mgnify

def download_mgy_proteins(dirname):

    os.makedirs(dirname, exist_ok=True)
    site_dir = read_url("https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/current_release/")
    parts = sorted(map(int,set(filter(lambda x: x.isnumeric(), map(lambda x: x.split('.fa.gz')[0], site_dir.decode().split('mgy_proteins_'))))))
    urls = [f"https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/current_release/mgy_proteins_{i}.fa.gz" for i in parts]
    gz_files = [f"{dirname}/mgy_proteins_{i}.fa.gz" for i in parts]
    for n, (url, gz_file) in enumerate(zip(urls, gz_files), 1):
        sys.stdout.write(f"        {n:>2} / {len(urls)}\r")
        download_url(url, gz_file, buffer_size=2048)
    
    fasta_file = f"{dirname}/mgy_proteins.fa"
    n_sequences = 0
    with open(fasta_file, 'wb') as wb:
        for n_parts, gz_file in enumerate(gz_files, 1):
            for line in gzip.open(gz_file, 'rb'):
                wb.write(line)
                if line.startswith(b'>'):
                    n_sequences += 1
                if n_sequences % 1000 == 0:
                    sys.stdout.write(f'part={n_parts}; n={n_sequences}\r')
    sys.stdout.write(f'part={n_parts}; n={n_sequences}\n')
    return fasta_file

######### database : pdb_seqres

def download_pdb_seqres(dirname):

    os.makedirs(dirname, exist_ok=True)
    gz_file = f"{dirname}/pdb_seqres.fasta.gz"
    download_url(f"https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz", gz_file, buffer_size=2048)

    fasta_file = f"{dirname}/pdb_seqres.fasta"
    n_sequences = 0
    with open(fasta_file, 'wb') as wb:
        for line in gzip.open(gz_file, 'rb'):
            wb.write(line)
            if line.startswith(b'>'):
                n_sequences += 1
            if n_sequences % 1000 == 0:
                sys.stdout.write(f'n={n_sequences}\r')
    sys.stdout.write(f'n={n_sequences}\n')
    return fasta_file


