#!/usr/bin/env python3
"""
Build docking rank matrices for selectivity analysis.

For the provided tRNA synthetases, constructs:
  Matrix 1 — compounds × reference pockets (target genes only)
  Matrix 2 — compounds × all pockets (all other genes)

Scores are rank-transformed per column: rank 0 = best binder (lowest score).

Usage:
    python 48_docking_hits.py --trna pheS --lib REAL
    python 48_docking_hits.py --trna pheS,aspS,lysS,alaS --lib REAL
"""

import argparse
import os
import sys

import pandas as pd
from tqdm import tqdm

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

LIBRARIES = {
    "DL":   os.path.join(ROOT, "output", "unidock_docking",        "docking_results"),
    "REAL": os.path.join(ROOT, "output", "unidock_REAL_docking_2", "docking_results"),
}

REFERENCE_POCKET_CSV = os.path.join(ROOT, "output", "reference_pocket.csv")


def load_gene_map():
    path = os.path.join(ROOT, "data", "mtb_trna_synthetases_bosch_2021_fig5.csv")
    df = pd.read_csv(path)
    return {row["gene_name_in_bosch_2021"]: row["uniprot_ac"] for _, row in df.iterrows()}


def load_reference_pockets():
    if not os.path.isfile(REFERENCE_POCKET_CSV):
        print(f"Error: {REFERENCE_POCKET_CSV} not found.")
        print("Run script 47 and manually populate output/reference_pocket.csv")
        print("  (columns: gene_name, pocket_name)")
        sys.exit(1)
    df = pd.read_csv(REFERENCE_POCKET_CSV)
    return dict(zip(df["gene_name"], df["pocket_name"]))


def load_scores(pocket_path):
    """Return a pd.Series of scores indexed by compound id."""
    df = pd.read_csv(pocket_path)
    return df.set_index("compound")["score"]


def build_matrix(pocket_map, results_dir, label=""):
    """
    Build a raw-score DataFrame from a dict {column_label: pocket_name}.
    Rows = compounds, columns = column_labels. Missing cells → NaN.
    """
    series = {}
    for col_label, pocket_name in tqdm(pocket_map.items(), desc=label, unit="pocket"):
        report = os.path.join(results_dir, pocket_name, "report.csv")
        if not os.path.isfile(report):
            print(f"  Warning: report not found for pocket '{pocket_name}', skipping.")
            continue
        series[col_label] = load_scores(report)
    if not series:
        return pd.DataFrame()
    return pd.concat(series, axis=1)


def rank_transform(df):
    """
    Per-column ascending rank: rank 0 = lowest score (best binder).
    Returns Int32 DataFrame (nullable integer to preserve NaN).
    """
    ranked = df.rank(method="min", ascending=True, na_option="keep") - 1
    return ranked.astype("Int32")


def print_summary(name, score_df, rank_df):
    n_compounds = len(score_df)
    n_pockets = len(score_df.columns)
    n_nan = int(score_df.isna().sum().sum())
    print(f"  {name}")
    print(f"    Compounds : {n_compounds:,}")
    print(f"    Pockets   : {n_pockets}")
    print(f"    NaN cells : {n_nan:,}")


def main():
    parser = argparse.ArgumentParser(
        description="Build docking rank matrices for selectivity analysis."
    )
    parser.add_argument("--trna", required=True,
                        help="Comma-separated gene name(s) (e.g. pheS or pheS,aspS)")
    parser.add_argument("--lib", choices=["DL", "REAL"], required=True,
                        help="Compound library: DL or REAL")
    args = parser.parse_args()

    genes = [g.strip() for g in args.trna.split(",")]
    results_dir = LIBRARIES[args.lib]

    gene_map = load_gene_map()
    ref_pockets = load_reference_pockets()

    # Validate genes
    invalid = [g for g in genes if g not in gene_map]
    if invalid:
        print(f"Error: unknown gene(s): {', '.join(invalid)}. Available genes:")
        print("  " + ", ".join(sorted(gene_map)))
        sys.exit(1)

    missing_ref = [g for g in genes if g not in ref_pockets]
    if missing_ref:
        print(f"Error: no reference pocket defined for: {', '.join(missing_ref)}")
        print(f"Add entries to {REFERENCE_POCKET_CSV}")
        sys.exit(1)

    target_uniprot_acs = {gene_map[g] for g in genes}

    # --- Matrix 1: target reference pockets ---
    target_pocket_map = {g: ref_pockets[g] for g in genes}

    # --- Matrix 2: all pockets from non-target proteins ---
    all_pockets = sorted(
        p for p in os.listdir(results_dir)
        if os.path.isdir(os.path.join(results_dir, p))
    )
    non_target_pocket_map = {
        p: p for p in all_pockets
        if not any(ac in p for ac in target_uniprot_acs)
    }

    print(f"Library : {args.lib}")
    print(f"Targets : {', '.join(genes)}")
    print(f"Building Matrix 1 ({len(target_pocket_map)} reference pocket(s))...")
    scores1 = build_matrix(target_pocket_map, results_dir, label="Matrix 1")

    print(f"Building Matrix 2 ({len(non_target_pocket_map)} non-target pocket(s))...")
    scores2 = build_matrix(non_target_pocket_map, results_dir, label="Matrix 2")

    ranks1 = rank_transform(scores1)
    ranks2 = rank_transform(scores2)

    print("\nMatrix summary (raw scores):")
    print_summary("Matrix 1 — target reference pockets", scores1, ranks1)
    print_summary("Matrix 2 — non-target pockets", scores2, ranks2)


if __name__ == "__main__":
    main()
