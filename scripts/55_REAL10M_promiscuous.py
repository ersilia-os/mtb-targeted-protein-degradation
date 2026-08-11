#!/usr/bin/env python3
"""
Pocket-level promiscuity screen over the surrogate model's (script 29/37, LazyQSAR RF, `bin_01`)
predicted probabilities for the Enamine REAL 9.56M library, one prediction per pocket across all
276 pockets - no protein-level collapse (a compound's 276 raw pocket scores are used directly,
not aggregated to one value per protein first).

Per pocket, raw probabilities are rank-transformed to a percentile in [0, 1]
(scipy.stats.rankdata(method="min"), so a tied block at the floor maps to exactly 0.0 - see
discussion in this script's originating conversation on why method="min" was chosen over
"average" tie-handling). Per compound, percentiles-of-percentiles are then computed across the
276 pockets: p1, p10, p25, p50, p75, p90, p99.

Threshold: p75 > 0.90, i.e. the compound ranks in the top 10%
of predicted probability in AT LEAST 25% of the 276 pockets. This is a deliberately pocket-level
(not protein-level) definition of "systematic" - a compound can pass by scoring consistently well
across many pockets of just 1-2 heavily-represented proteins (e.g. P9WFS9 alone contributes 40 of
the 276 pockets), so this is a looser, more permissive notion of promiscuity than "hits many
distinct targets" (contrast with scripts 52-54, which require >=2 of a small set of DISTINCT
GENES). Whether pocket-level or protein-level promiscuity is more meaningful is left to the user's
interpretation of the output.

Data path note: scripts 32/33 expect per-pocket inference arrays under
processed/unidock_docking/inference_probs/, but on this machine the downloaded data instead lives
under output/unidock_docking/inference_probs/ (mirrors the output/ vs processed/ drift already
present in src/docking_utils.py's SMILES_PATHS). This script reads directly from the output/
location.

Alongside the p75 filter, each filtered compound is also annotated with n_distinct_proteins_top10:
the number of distinct proteins (out of 21) for which the compound ranks in the top 10% in AT
LEAST ONE of that protein's own pockets. This does not affect the filter itself - it's a reported
diagnostic so pocket-level promiscuity (many pockets, possibly redundant structures of few
proteins) can be distinguished at a glance from genuine multi-target breadth.

Outputs (output/55_identify_promiscuous_enamine_REAL/):
  - promiscuous_hits.csv: id, smiles, p1/p10/p25/p50/p75/p90/p99, n_distinct_proteins_top10
  - promiscuous_indices.pkl: a set of row indices (positions into success_mols.pkl) of the
    filtered compounds, for reuse by downstream scripts without recomputing the full matrix

Usage:
    python 55_REAL10M_promiscuous.py
"""
import os
import pickle
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from docking_utils import lookup_smiles

PATH_TO_INFERENCE_PROBS = os.path.join(ROOT, "output", "unidock_docking", "inference_probs")
OUTPUT_DIR = os.path.join(ROOT, "output", "55_identify_promiscuous_enamine_REAL")
os.makedirs(OUTPUT_DIR, exist_ok=True)

PERCENTILE_LEVELS = [1, 10, 25, 50, 75, 90, 99]
FILTER_LEVEL = 75
FILTER_THRESHOLD = 0.90  # top-10% percentile in >=25% of pockets

PROTEIN_RE = re.compile(r"^[^_]+_(P9W[A-Z0-9]+)_")


def protein_of(fname):
    return PROTEIN_RE.match(fname).group(1)


def main():
    npz_files = sorted(f for f in os.listdir(PATH_TO_INFERENCE_PROBS) if f.endswith(".npz"))
    print(f"Found {len(npz_files)} pocket files")

    with open(os.path.join(PATH_TO_INFERENCE_PROBS, "success_mols.pkl"), "rb") as f_in:
        success_mols = np.array(pickle.load(f_in))

    n = len(success_mols)
    P = np.empty((n, len(npz_files)), dtype=np.float32)

    proteins = sorted(set(protein_of(f) for f in npz_files))
    protein_top10 = {p: np.zeros(n, dtype=bool) for p in proteins}
    print(f"{len(proteins)} distinct proteins across the {len(npz_files)} pockets")

    for j, f in enumerate(npz_files):
        arr = np.load(os.path.join(PATH_TO_INFERENCE_PROBS, f))["arr_0"]
        assert len(arr) == n, f"{f}: expected {n} compounds, got {len(arr)}"
        ranks = rankdata(arr, method="min")
        pct = (ranks - 1) / (n - 1)
        P[:, j] = pct
        protein_top10[protein_of(f)] |= (pct > FILTER_THRESHOLD)
        if (j + 1) % 50 == 0:
            print(f"  processed {j + 1}/{len(npz_files)} pockets...")

    print(f"\nPercentile matrix shape: {P.shape}  (n_compounds x n_pockets)")
    print(f"per-pocket min/max (should be 0.0 / 1.0 for all): "
          f"min={P.min(axis=0).max():.4f}, max={P.max(axis=0).min():.4f}")

    row_pctls = np.percentile(P, PERCENTILE_LEVELS, axis=1)  # shape (7, n_compounds)
    del P

    n_distinct_proteins = np.zeros(n, dtype=np.int8)
    for arr in protein_top10.values():
        n_distinct_proteins += arr.astype(np.int8)
    del protein_top10

    p75 = row_pctls[PERCENTILE_LEVELS.index(FILTER_LEVEL)]
    mask = p75 > FILTER_THRESHOLD
    n_hits = int(mask.sum())
    print(f"\nFilter: p{FILTER_LEVEL} > {FILTER_THRESHOLD} "
          f"(top-10% percentile in >={100 * (1 - FILTER_LEVEL / 100):.0f}% of the "
          f"{len(npz_files)} pockets)")
    print(f"Compounds passing filter: {n_hits} / {n} ({100 * n_hits / n:.4f}%)")

    idx_hits = np.where(mask)[0]
    hit_ids = success_mols[idx_hits]

    result = pd.DataFrame({"id": hit_ids})
    for k, level in enumerate(PERCENTILE_LEVELS):
        result[f"p{level}"] = row_pctls[k, idx_hits]
    result["n_distinct_proteins_top10"] = n_distinct_proteins[idx_hits]
    result = result.sort_values("p90", ascending=False).reset_index(drop=True)

    print(f"\nLooking up SMILES for {len(hit_ids)} compounds...")
    smiles_map = lookup_smiles(hit_ids.tolist(), "REAL_ROUND1")
    result["smiles"] = result["id"].map(smiles_map)
    n_missing_smiles = result["smiles"].isna().sum()
    if n_missing_smiles:
        print(f"  Warning: {n_missing_smiles} compounds missing a SMILES match")

    cols = ["id", "smiles"] + [f"p{level}" for level in PERCENTILE_LEVELS] + ["n_distinct_proteins_top10"]
    result = result[cols]

    out_path = os.path.join(OUTPUT_DIR, "promiscuous_hits.csv")
    result.to_csv(out_path, index=False)
    print(f"\nSaved {len(result)} compounds to {out_path}")

    print(f"\nDistinct-protein breadth (out of {len(proteins)}) among the {n_hits} filtered compounds:")
    print(result["n_distinct_proteins_top10"].describe())

    indices_path = os.path.join(OUTPUT_DIR, "promiscuous_indices.pkl")
    with open(indices_path, "wb") as f_out:
        pickle.dump(set(idx_hits.tolist()), f_out)
    print(f"Saved {len(idx_hits)} promiscuous compound indices (positions in success_mols.pkl) "
          f"to {indices_path}")


if __name__ == "__main__":
    main()
