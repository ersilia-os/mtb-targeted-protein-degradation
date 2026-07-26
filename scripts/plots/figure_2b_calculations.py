"""
Randomly samples the Enamine REAL 10B library for figure_2b: 1000 compounds per chunk (994
chunks), merged and re-sampled down to 100k, then run through the same physicochemical property
calculation as figure_2a_calculations.py.

Resumable: one CSV per chunk (processed/figure_2b_REAL10B_samples/{chunk}.csv), skips chunks
already sampled. Each chunk's ~90MB SMILES mapping is deleted after use to keep tmp/ bounded
across 994 chunks (herbert's root disk has little headroom).

Usage:
    python figure_2b_calculations.py
"""
import argparse
import glob
import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

import pandas as pd

from default import RANDOM_SEED
from docking_utils import compute_properties
from screening_10b_utils import download_file, list_smiles_chunks

ROOT = os.path.join(root, "..", "..")
TMP_DIR = os.path.join(ROOT, "tmp")
SAMPLES_DIR = os.path.join(ROOT, "processed", "figure_2b_REAL10B_samples")
PLOTS_DIR = os.path.join(ROOT, "plots", "figure_2")
os.makedirs(TMP_DIR, exist_ok=True)
os.makedirs(SAMPLES_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

N_PER_CHUNK = 1000
N_FINAL = 100_000


def banner(title):
    line = "=" * (len(title) + 10)
    print(line)
    print(f"==== {title} ====")
    print(line)


def sample_chunk(chunk, chunk_index):
    out_path = os.path.join(SAMPLES_DIR, f"{chunk}.csv")
    if os.path.isfile(out_path):
        return out_path

    ids_path = os.path.join(TMP_DIR, f"{chunk}_SMILES_IDs.tsv.zip")
    download_file(ids_path)
    try:
        df = pd.read_csv(ids_path, sep="\t")
        sample = df.sample(n=min(N_PER_CHUNK, len(df)), random_state=RANDOM_SEED + chunk_index)
        sample[["id", "smiles"]].to_csv(out_path, index=False)
    finally:
        os.remove(ids_path)
    return out_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-chunks", type=int, default=None,
                         help="Limit to the first N chunks (for a quick dry run).")
    args = parser.parse_args()

    banner("Discovering REAL 10B chunks")
    chunks = list_smiles_chunks()
    if args.max_chunks is not None:
        chunks = chunks[:args.max_chunks]
    print(f"{len(chunks)} chunk(s) to process")

    banner("Sampling 1000 compounds per chunk")
    n_processed, n_skipped = 0, 0
    for i, chunk in enumerate(chunks):
        out_path = os.path.join(SAMPLES_DIR, f"{chunk}.csv")
        if os.path.isfile(out_path):
            n_skipped += 1
            continue
        print(f"[{i + 1}/{len(chunks)}] {chunk}")
        sample_chunk(chunk, i)
        n_processed += 1
    print(f"Processed {n_processed}, skipped (already done) {n_skipped}")

    banner("Merging per-chunk samples")
    sample_files = sorted(glob.glob(os.path.join(SAMPLES_DIR, "*.csv")))
    merged = pd.concat([pd.read_csv(f) for f in sample_files], ignore_index=True)
    print(f"Merged pool: {len(merged):,} compounds from {len(sample_files)} chunk file(s)")

    final = merged.sample(n=min(N_FINAL, len(merged)), random_state=RANDOM_SEED)
    print(f"Final random sample: {len(final):,} compounds")

    banner("Computing properties")
    smiles_dict = dict(zip(final["id"], final["smiles"]))
    props = compute_properties(smiles_dict)
    output_path = os.path.join(PLOTS_DIR, "figure_2b_REAL10B.csv")
    props.to_csv(output_path)
    print(f"Saved to {output_path}")


if __name__ == "__main__":
    main()
