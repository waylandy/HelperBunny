import os
import sys
import polars as pl

from src.hmmer import (
    get_hmm_hits_parquets,
    get_mapped_hits_with_context,
    )


INPUT_FASTA   = os.environ.get("INPUT_FASTA")
PROFILE_HMM   = os.environ.get("PROFILE_HMM")
OUTPUT_DIR    = os.environ.get("OUTPUT_DIR")
BUFFER_SIZE   = int(os.environ.get("BUFFER_SIZE", 300000))
E_VALUE       = float(os.environ.get("E_VALUE", 1e-2))
DATABASE_SIZE = int(os.environ.get("DATABASE_SIZE", 0))
N_THREADS     = int(os.environ.get("N_THREADS", 3))

assert INPUT_FASTA is not None
assert PROFILE_HMM is not None
DATABASE_SIZE = None if DATABASE_SIZE == 0 else DATABASE_SIZE

HITS_PARQUET   = f"{OUTPUT_DIR}/hits.parquet"
MAPPED_PARQUET = f"{OUTPUT_DIR}/mapped.parquet"
ALIGNMENTS_DIR = f"{OUTPUT_DIR}/mapped"
os.makedirs(ALIGNMENTS_DIR, exist_ok=True)

sys.stdout.write(f"\n")
sys.stdout.write(f"INPUT_FASTA   = {INPUT_FASTA}\n")
sys.stdout.write(f"PROFILE_HMM   = {PROFILE_HMM}\n")
sys.stdout.write(f"OUTPUT_DIR    = {OUTPUT_DIR}\n")
sys.stdout.write(f"BUFFER_SIZE   = {BUFFER_SIZE}\n")
sys.stdout.write(f"E_VALUE       = {E_VALUE}\n")
sys.stdout.write(f"DATABASE_SIZE = {DATABASE_SIZE}\n")
sys.stdout.write(f"N_THREADS     = {N_THREADS}\n")
sys.stdout.write(f"\n")

###

get_hmm_hits_parquets(
    INPUT_FASTA,
    PROFILE_HMM,
    HITS_PARQUET,
    fasta_buffer_size = BUFFER_SIZE,
    parquet_partition_size = 5000,
    e_value = E_VALUE,
    db_size = DATABASE_SIZE,
    cpus = N_THREADS,
    )

get_mapped_hits_with_context(
    HITS_PARQUET, 
    extend_by_gaps=False, 
    extend_terms=(0, 0),
    ).write_parquet(MAPPED_PARQUET)

pl.read_parquet(MAPPED_PARQUET)
