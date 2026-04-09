import os
import sys

import numpy as np
import polars as pl

from utils import (
    validate_a2m,
    get_polars_from_fasta,
    a2m_to_aln,
    get_coupling_scores_from_a3m_string,
    )

INPUT_A2M   = os.environ.get("INPUT_A2M", None)
POTTS_NPZ   = os.environ.get("POTTS_NPZ", None)
OUTPUT_FILE = os.environ.get("OUTPUT_FILE", "output_scores")

assert INPUT_A2M is not None
assert POTTS_NPZ is not None

sys.stdout.write(f"\n")
sys.stdout.write(f"INPUT_A2M   = {INPUT_A2M}\n")
sys.stdout.write(f"POTTS_NPZ   = {POTTS_NPZ}\n")
sys.stdout.write(f"OUTPUT_FILE = {OUTPUT_FILE}\n")
sys.stdout.write(f"\n")

###

potts = np.load(POTTS_NPZ)
couplings, states = potts['couplings'], potts['states']

n_positions, n_sequences = validate_a2m(INPUT_A2M)
table = get_polars_from_fasta(INPUT_A2M).to_pandas().rename({"sequence": "a2m"}, axis=1)

scores = []
for sequence in table["a2m"].values:
    scores_ = get_coupling_scores_from_a3m_string(a2m_to_aln(sequence), couplings, states)
    scores += [scores_.sum()]

table["couplings"] = np.array(scores, dtype=np.float32)
table = pl.DataFrame(table)
table.write_parquet(OUTPUT_FILE)
