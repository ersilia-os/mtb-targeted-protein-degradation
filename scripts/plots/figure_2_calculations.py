"""
Computes all data feeding figure_2_plot.py: physicochemical properties for HL/REAL 10M/REAL 10B
(panel B) and per-pocket docking-score percentiles for all three libraries (panel C).

HL and REAL 10M (compute_hl_and_real10m): random 100k samples, physicochemical properties via
docking_utils.compute_properties(). REAL 10B (compute_real10b): the Enamine REAL 10B library is
too large to hold locally, so it's sampled 1000 compounds per chunk (994 chunks), merged and
re-sampled down to 100k, then run through the same property calculation. Resumable: one CSV per
chunk (processed/figure_2b_REAL10B_samples/{chunk}.csv), skips chunks already sampled - each
chunk's ~90MB SMILES mapping is deleted after use to keep tmp/ bounded across 994 chunks (herbert's
root disk has little headroom).

Docking percentiles (compute_docking_percentiles): p1/p0.1 (top-1,000/top-100 out of each
pocket's screened compound set) per pocket, for all three libraries. Stores precomputed
percentiles rather than raw score arrays - 276 pockets x ~100k compounds per library would make
for an unnecessarily large file. REAL 10M percentiles use each pocket's prioritized top-100k
"active" set only (docking_utils.load_real_positive_scores), not the shared ~12,958-compound
background sample also present in that library's docking output - matching HL's and REAL 10B's
fully-screened sets.

Usage:
    python figure_2_calculations.py [--max-chunks N]
"""
import argparse
import glob
import json
import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

import numpy as np
import pandas as pd

from default import RANDOM_SEED
from docking_utils import LIBRARIES, compute_properties, load_real_positive_scores, load_scores, sample_prescreened_smiles
from screening_10b_utils import download_file, list_smiles_chunks

ROOT = os.path.join(root, "..", "..")
TMP_DIR = os.path.join(ROOT, "tmp")
SAMPLES_DIR = os.path.join(ROOT, "processed", "figure_2b_REAL10B_samples")
plots_dir = os.path.join(ROOT, "plots", "figure_2")
os.makedirs(TMP_DIR, exist_ok=True)
os.makedirs(SAMPLES_DIR, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)

HL_SAMPLE_N = 100_000
REAL10M_SAMPLE_N = 100_000
N_PER_CHUNK = 1000
N_FINAL = 100_000


def banner(title):
    line = "=" * (len(title) + 10)
    print(line)
    print(f"==== {title} ====")
    print(line)


def compute_hl_and_real10m():
    banner("HL (Enamine Hit Locator 100K) - random 100k sample")
    hl_smiles = sample_prescreened_smiles("DL", exclude_ids=set(), n=HL_SAMPLE_N, seed=RANDOM_SEED)
    print(f"Sampled {len(hl_smiles):,} compounds")
    hl_props = compute_properties(hl_smiles)
    output_path = os.path.join(plots_dir, "figure_2a_HL.csv")
    hl_props.to_csv(output_path)
    print(f"Saved to {output_path}")

    banner("REAL 10M (Enamine REAL 9.56M) - random 100k sample")
    real10m_smiles = sample_prescreened_smiles("REAL_ROUND1", exclude_ids=set(), n=REAL10M_SAMPLE_N, seed=RANDOM_SEED)
    print(f"Sampled {len(real10m_smiles):,} compounds")
    real10m_props = compute_properties(real10m_smiles)
    output_path = os.path.join(plots_dir, "figure_2a_REAL10M.csv")
    real10m_props.to_csv(output_path)
    print(f"Saved to {output_path}")


def sample_real10b_chunk(chunk, chunk_index):
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


def compute_real10b(max_chunks=None):
    banner("Discovering REAL 10B chunks")
    chunks = list_smiles_chunks()
    if max_chunks is not None:
        chunks = chunks[:max_chunks]
    print(f"{len(chunks)} chunk(s) to process")

    banner("Sampling 1000 compounds per chunk")
    n_processed, n_skipped = 0, 0
    for i, chunk in enumerate(chunks):
        out_path = os.path.join(SAMPLES_DIR, f"{chunk}.csv")
        if os.path.isfile(out_path):
            n_skipped += 1
            continue
        print(f"[{i + 1}/{len(chunks)}] {chunk}")
        sample_real10b_chunk(chunk, i)
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
    output_path = os.path.join(plots_dir, "figure_2b_REAL10B.csv")
    props.to_csv(output_path)
    print(f"Saved to {output_path}")


def pocket_percentiles(pocket, scores, library, uniprot_to_gene):
    uniprot_ac = pocket.split("_model_")[0].split("_")[-1]
    return {
        "pocket": pocket,
        "uniprot_ac": uniprot_ac,
        "gene": uniprot_to_gene.get(uniprot_ac, "unknown"),
        "library": library,
        "n": len(scores),
        "p0_1": np.percentile(scores, 0.1),
        "p1": np.percentile(scores, 1),
    }


def compute_docking_percentiles():
    banner("Loading gene mapping from figure 1's color mapping")
    with open(os.path.join(ROOT, "plots", "figure_1", "color_mapping.json")) as f:
        uniprot_to_gene = json.load(f)["uniprot_to_gene"]

    rows = []

    banner("HL (output/unidock_docking) - percentiles per pocket")
    for pocket in sorted(os.listdir(LIBRARIES["DL"])):
        report_path = os.path.join(LIBRARIES["DL"], pocket, "report.csv")
        if not os.path.isfile(report_path):
            print(f"  Warning: report not found for pocket '{pocket}', skipping.")
            continue
        scores = load_scores(report_path).values
        rows.append(pocket_percentiles(pocket, scores, "HL", uniprot_to_gene))

    banner("REAL 10M (output/unidock_REAL_docking, prioritized 100k only) - percentiles per pocket")
    for pocket in sorted(os.listdir(LIBRARIES["DL"])):
        scores = load_real_positive_scores(pocket).values
        if len(scores) == 0:
            print(f"  Warning: no REAL 10M positive scores for pocket '{pocket}', skipping.")
            continue
        rows.append(pocket_percentiles(pocket, scores, "REAL 10M", uniprot_to_gene))

    banner("REAL 10B (output/unidock_REAL_docking_2) - percentiles per pocket")
    for pocket in sorted(os.listdir(LIBRARIES["REAL"])):
        report_path = os.path.join(LIBRARIES["REAL"], pocket, "report.csv")
        if not os.path.isfile(report_path):
            print(f"  Warning: report not found for pocket '{pocket}', skipping.")
            continue
        scores = load_scores(report_path).values
        rows.append(pocket_percentiles(pocket, scores, "REAL 10B", uniprot_to_gene))

    df = pd.DataFrame(rows)
    output_path = os.path.join(plots_dir, "figure_2c_docking_percentiles.csv")
    df.to_csv(output_path, index=False)
    print(f"Saved {len(df)} row(s) to {output_path}")
    for library, group in df.groupby("library"):
        print(f"  {library}: {len(group)} pocket(s)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-chunks", type=int, default=None,
                         help="Limit REAL 10B sampling to the first N chunks (for a quick dry run).")
    args = parser.parse_args()

    compute_hl_and_real10m()
    compute_real10b(max_chunks=args.max_chunks)
    compute_docking_percentiles()


if __name__ == "__main__":
    main()
