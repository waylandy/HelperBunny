import os
import sys

from utils import (
    iter_fasta,
    validate_a2m,
    get_couplings_npz,
    )

INPUT_A2M  = os.environ.get("INPUT_A2M", None)
OUTPUT_NPZ = os.environ.get("OUTPUT_NPZ", None)
DEVICE     = int(os.environ.get("DEVICE", 0))
THREADS    = int(os.environ.get("THREADS", 1))
NUMITER    = int(os.environ.get("NUMITER", 50))
EPSILON    = float(os.environ.get("EPSILON", 0.01))
LASTK      = int(os.environ.get("LASTK", 5))
IDTHRES    = float(os.environ.get("IDTHRES", 0.8))
LFACTOR    = float(os.environ.get("LFACTOR", 0.2))

assert INPUT_A2M is not None
OUTPUT_NPZ = f"{INPUT_A2M}.npz" if OUTPUT_NPZ is None else OUTPUT_NPZ

sys.stdout.write(f"\n")
sys.stdout.write(f"INPUT_A2M  = {INPUT_A2M}\n")
sys.stdout.write(f"OUTPUT_NPZ = {OUTPUT_NPZ}\n")
sys.stdout.write(f"DEVICE     = {DEVICE}\n")
sys.stdout.write(f"THREADS    = {THREADS}\n")
sys.stdout.write(f"NUMITER    = {NUMITER}\n")
sys.stdout.write(f"EPSILON    = {EPSILON}\n")
sys.stdout.write(f"LASTK      = {LASTK}\n")
sys.stdout.write(f"IDTHRES    = {IDTHRES}\n")
sys.stdout.write(f"LFACTOR    = {LFACTOR}\n")
sys.stdout.write(f"\n")

###

n_positions, n_sequences = validate_a2m(INPUT_A2M)
a3m_strings = [sequence for header, sequence in iter_fasta(INPUT_A2M)]

output_npz = get_couplings_npz(
    a3m_strings,
    OUTPUT_NPZ,
    device = DEVICE,
    threads = THREADS,
    numiter = NUMITER,
    epsilon = EPSILON,
    lastk = LASTK,
    idthres = IDTHRES,
    lfactor = LFACTOR,
    )
