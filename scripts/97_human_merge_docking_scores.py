#!/usr/bin/env python3
"""
Merges script 96's docking results (389 pockets x 1,095 compounds = 425,955 endpoints) with
script 93's pocket annotation (InterPro domain labels + classified AlphaFill ligand evidence) and
script 70's compound info (SMILES, source), into one long table -- one row per (pocket, compound).

Usage:
    python 97_human_merge_docking_scores.py
"""
import glob
import os

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

DOCKING_RESULTS_DIR = os.path.join(ROOT, "output", "96_human_docking", "docking_results")
MERGED_POCKET_DATA_CSV = os.path.join(ROOT, "output", "93_human_merge_pocket_annotations", "merged_pocket_data.csv")
FILTERED_HITS_CSV = os.path.join(ROOT, "output", "70_filtering", "filtered_hits.csv")

OUTPUT_DIR = os.path.join(ROOT, "output", "97_human_merge_docking_scores")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_docking_scores():
    """Long table: Uniprot AC, Pocket number, compound_id, score -- one row per (pocket, compound),
    read directly from script 96's per-pocket report.csv files."""
    report_paths = sorted(glob.glob(os.path.join(DOCKING_RESULTS_DIR, "*", "*", "report.csv")))
    rows = []
    for path in report_paths:
        pocket_number = int(os.path.basename(os.path.dirname(path)))
        uniprot_ac = os.path.basename(os.path.dirname(os.path.dirname(path)))
        df = pd.read_csv(path).rename(columns={"compound": "compound_id"})
        df.insert(0, "Pocket number", pocket_number)
        df.insert(0, "Uniprot AC", uniprot_ac)
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


def main():
    print("Loading docking scores from all pockets...")
    scores = load_docking_scores()
    n_pockets = scores.groupby(["Uniprot AC", "Pocket number"]).ngroups
    print(f"  {len(scores):,} (pocket, compound) rows from {n_pockets} pockets")

    pocket_annotation = pd.read_csv(MERGED_POCKET_DATA_CSV)
    merged = scores.merge(pocket_annotation, on=["Uniprot AC", "Pocket number"],
                           how="left", validate="many_to_one")
    assert len(merged) == len(scores), "pocket-annotation merge changed row count"

    compounds = pd.read_csv(FILTERED_HITS_CSV)[["compound_id", "smiles", "source"]]
    merged = merged.merge(compounds, on="compound_id", how="left", validate="many_to_one")
    assert len(merged) == len(scores), "compound-info merge changed row count"

    out_path = os.path.join(OUTPUT_DIR, "docking_scores.csv")
    merged.to_csv(out_path, index=False)

    n_expected = 389 * 1095
    print(f"\nSaved {len(merged):,} rows x {len(merged.columns)} columns -> {out_path}")
    print(f"Expected {n_expected:,} rows (389 pockets x 1,095 compounds): "
          f"{'MATCH' if len(merged) == n_expected else 'MISMATCH'}")
    print(f"\nScore summary:\n{merged['score'].describe()}")
    print(f"\nMissing scores: {merged['score'].isna().sum()}")


if __name__ == "__main__":
    main()
