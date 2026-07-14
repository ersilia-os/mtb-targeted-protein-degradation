#!/usr/bin/env python3
"""
Joins scripts 62 (physchem), 66 (docking scores/ranks) and 67 (Ersilia model outputs) into one
wide table for the 2,923 aggregated hits, then applies a 7-rule sequential filter: QED > 0.5,
300 < MW < 500, NSPS in [10, 40], no PAINS (eos2xeq has_pains flag), all 3 eos42ez cytotoxicity
endpoints < 0.3, top-80% (lowest, i.e. most permeable) eos5jv3 mycomembrane_permeation, and
(>=2 CAT pockets with rank <= 10,000 AND score <= -10) OR (>=1 NON-CAT pocket with rank <= 10,000
AND score <= -8). Rank is vs. the ~99,105-compound REAL round-2 background (script 66) and
normalizes for pocket-specific score-scale differences that a raw-score-only cutoff would miss.

Usage:
    python 70_filtering.py
"""
import os

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

AGGREGATED_HITS_CSV = os.path.join(ROOT, "output", "62_aggregate_hits", "aggregated_hits.csv")
MERGED_DOCKING_SCORES_CSV = os.path.join(ROOT, "output", "66_merge_docking_scores", "merged_docking_scores.csv")
ERSILIA_DIR = os.path.join(ROOT, "output", "67_ersilia_characterization")
SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")

OUTPUT_DIR = os.path.join(ROOT, "output", "70_filtering")
os.makedirs(OUTPUT_DIR, exist_ok=True)

ERSILIA_MODELS = [
    "eos12x7.csv",
    "eos5jv3.csv",
    "eos2xeq.csv",
    "eos42ez.csv",
]

NSPS_RANGE = (10, 40)
MW_RANGE = (300, 500)
CYTOTOXICITY_COLS = ["cytotoxicity_hepg2", "cytotoxicity_hskmc", "cytotoxicity_imr90"]
EOS5JV3_TOP_PCT = 0.80
SOURCE_RENAME = {
    "cat_promiscuous": "CAT-promiscuous",
    "cat_selective": "CAT-selective",
    "noncat_promiscuous": "NONCAT-promiscuous",
    "noncat_top100_10m": "NONCAT-top100-10M",
    "noncat_top100_10b": "NONCAT-top100-10B",
}
CAT_SCORE_THRESHOLD = -10
CAT_RANK_THRESHOLD = 10000
CAT_MIN_COUNT = 2
NONCAT_SCORE_THRESHOLD = -8
NONCAT_RANK_THRESHOLD = 10000
NONCAT_MIN_COUNT = 1


def load_pockets_by_site_type():
    pockets = pd.read_csv(SELECTED_POCKETS_CSV)
    return {site_type: group["pocket_name"].tolist() for site_type, group in pockets.groupby("site_type")}


def build_merged_table():
    df = pd.read_csv(AGGREGATED_HITS_CSV)
    n = len(df)

    df["source"] = df["source"].apply(lambda s: ";".join(SOURCE_RENAME[tok] for tok in s.split(";")))

    docking = pd.read_csv(MERGED_DOCKING_SCORES_CSV).drop(columns=["smiles", "source"])
    df = df.merge(docking, on="compound_id", how="left")
    assert len(df) == n, "merged_docking_scores.csv merge changed row count"

    for filename in ERSILIA_MODELS:
        model = pd.read_csv(os.path.join(ERSILIA_DIR, filename)).drop(columns=["key"])
        df = df.merge(model, left_on="smiles", right_on="input", how="left").drop(columns=["input"])
        assert len(df) == n, f"{filename} merge changed row count"

    return df


def apply_filters(df):
    n = len(df)
    print(f"Start: {n:,}")

    mask = df["QED"] > 0.5
    print(f"QED > 0.5: {mask.sum():,} ({100 * mask.sum() / n:.1f}%)")

    mask &= df["MW"].between(*MW_RANGE, inclusive="neither")
    print(f"300 < MW < 500: {mask.sum():,} ({100 * mask.sum() / n:.1f}%)")

    mask &= df["nsps_score"].between(*NSPS_RANGE, inclusive="both")
    print(f"NSPS in [10, 40]: {mask.sum():,} ({100 * mask.sum() / n:.1f}%)")

    mask &= df["has_pains"] == 0
    print(f"No PAINS: {mask.sum():,} ({100 * mask.sum() / n:.1f}%)")

    mask &= (df[CYTOTOXICITY_COLS] < 0.3).all(axis=1)
    print(f"Cytotoxicity (all 3 < 0.3): {mask.sum():,} ({100 * mask.sum() / n:.1f}%)")

    eos5jv3_cutoff = df["mycomembrane_permeation"].quantile(EOS5JV3_TOP_PCT)
    mask &= df["mycomembrane_permeation"] <= eos5jv3_cutoff
    print(f"eos5jv3 top-80% (<= {eos5jv3_cutoff:.3f}): {mask.sum():,} ({100 * mask.sum() / n:.1f}%)")

    pockets_by_site = load_pockets_by_site_type()
    cat_pockets = pockets_by_site["CAT"]
    noncat_pockets = pockets_by_site["NON-CAT"]

    cat_qualify = pd.concat(
        [(df[f"{p}_rank"] <= CAT_RANK_THRESHOLD) & (df[p] <= CAT_SCORE_THRESHOLD) for p in cat_pockets], axis=1)
    noncat_qualify = pd.concat(
        [(df[f"{p}_rank"] <= NONCAT_RANK_THRESHOLD) & (df[p] <= NONCAT_SCORE_THRESHOLD) for p in noncat_pockets], axis=1)

    cat_hits = cat_qualify.sum(axis=1)
    noncat_hits = noncat_qualify.sum(axis=1)
    mask &= (cat_hits >= CAT_MIN_COUNT) | (noncat_hits >= NONCAT_MIN_COUNT)
    print(f">=2 CAT pockets (rank<=10000, score<=-10) OR >=1 NON-CAT pocket (rank<=10000, score<=-8): "
          f"{mask.sum():,} ({100 * mask.sum() / n:.1f}%)")

    return df[mask]


def main():
    df = build_merged_table()
    merged_path = os.path.join(OUTPUT_DIR, "merged_all_results.csv")
    df.to_csv(merged_path, index=False)
    print(f"Saved {len(df):,} rows x {len(df.columns)} columns to {merged_path}\n")

    filtered = apply_filters(df)
    filtered_path = os.path.join(OUTPUT_DIR, "filtered_hits.csv")
    filtered.to_csv(filtered_path, index=False)
    print(f"\nSaved {len(filtered):,} rows x {len(filtered.columns)} columns to {filtered_path}")


if __name__ == "__main__":
    main()
