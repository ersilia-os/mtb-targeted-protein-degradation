#!/usr/bin/env python3
"""
Merge script 65's 12 per-pocket results.csv files (mean docking score across 5 replicates) into
one wide table: compound_id, smiles,
<pocket_1>, <pocket_1>_rank, <pocket_2>, <pocket_2>_rank, ... - score and rank vs. the Enamine REAL
round-2 library (~99,105 compounds) for that pocket.

Usage:
    python 66_merge_docking_scores.py
"""
import os

import numpy as np
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

AGGREGATED_HITS_CSV = os.path.join(ROOT, "output", "62_aggregate_hits", "aggregated_hits.csv")
SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")
DOCKING_RESULTS_DIR = os.path.join(ROOT, "output", "65_aggregated_docking", "docking_results")

REAL_ROUND2_DIR = os.path.join(ROOT, "output", "unidock_REAL_docking_2", "docking_results")
MULTIMER_DIR = os.path.join(ROOT, "output", "50_unidock_docking_multimers")
DIMER_POCKET = "7K98_pocket_6"

OUTPUT_DIR = os.path.join(ROOT, "output", "66_merge_docking_scores")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_report_scores(path):
    """{compound_id: score} from a report.csv, or {} if missing."""
    if not os.path.isfile(path):
        return {}
    df = pd.read_csv(path)
    df = df[df["score"].notna() & (df["score"] != "")]
    return dict(zip(df["compound"], df["score"].astype(float)))


def load_round2_background(pockets):
    """{pocket_name: sorted np.array of scores} - the Enamine REAL round-2 library (~99,105
    compounds), used exclusively as the ranking background (unlike load_old_scores' broader
    round-1/REAL10B-selective fallback). Dimer pocket's round-2-equivalent scores live in the
    separate multimer docking output (script 50), not REAL_ROUND2_DIR."""
    background = {}
    for pocket_name in pockets:
        base = MULTIMER_DIR if pocket_name == DIMER_POCKET else REAL_ROUND2_DIR
        scores = load_report_scores(os.path.join(base, pocket_name, "report.csv"))
        background[pocket_name] = np.sort(np.array(list(scores.values())))
    return background


def main():
    result = pd.read_csv(AGGREGATED_HITS_CSV, usecols=["compound_id", "smiles", "source"])
    n_compounds = len(result)
    print(f"Base compound list: {n_compounds:,} compounds")

    props = pd.read_csv(AGGREGATED_HITS_CSV, usecols=["QED", "MW"])
    n_qed = (props["QED"] > 0.5).sum()
    n_mw = ((props["MW"] > 300) & (props["MW"] < 500)).sum()
    print(f"QED > 0.5: {n_qed:,} / {n_compounds:,} ({100 * n_qed / n_compounds:.1f}%)")
    print(f"MW in (300, 500): {n_mw:,} / {n_compounds:,} ({100 * n_mw / n_compounds:.1f}%)")

    pockets = pd.read_csv(SELECTED_POCKETS_CSV)["pocket_name"].tolist()
    print(f"Pockets: {len(pockets)}")

    background = load_round2_background(pockets)

    std_cols = []
    for pocket_name in pockets:
        results_path = os.path.join(DOCKING_RESULTS_DIR, pocket_name, "results.csv")
        report = pd.read_csv(results_path, usecols=["compound", "score", "score_std"])
        std_col = f"{pocket_name}__std"
        report = report.rename(columns={"compound": "compound_id", "score": pocket_name, "score_std": std_col})
        report[pocket_name] = report[pocket_name].astype(float)
        std_cols.append(std_col)

        result = result.merge(report, on="compound_id", how="left")
        assert len(result) == n_compounds, f"{pocket_name}: merge changed row count (duplicate compound_ids?)"

        n_missing = result[pocket_name].isna().sum()
        if n_missing:
            print(f"  Warning: {pocket_name} missing {n_missing} score(s)")

        scores = result[pocket_name]
        ranks = np.searchsorted(background[pocket_name], scores.fillna(np.inf).to_numpy(), side="left")
        ranks = pd.array(ranks, dtype="Int64")
        ranks[scores.isna().to_numpy()] = pd.NA
        result[f"{pocket_name}_rank"] = ranks

    total_endpoints = n_compounds * len(pockets)
    all_scores = result[pockets].to_numpy()
    all_stds = result[std_cols].to_numpy()

    n_le_10 = (all_scores <= -10).sum()
    n_le_11 = (all_scores <= -11).sum()
    n_le_12 = (all_scores <= -12).sum()
    n_std_le_01 = (all_stds <= 0.1).sum()
    n_std_le_05 = (all_stds <= 0.5).sum()
    n_std_le_10 = (all_stds <= 1.0).sum()

    print(f"\nTotal endpoints (compounds x pockets): {total_endpoints:,}")
    print(f"endpoints with score <=-10: {n_le_10:,} ({100 * n_le_10 / total_endpoints:.1f}%)")
    print(f"endpoints with score <=-11: {n_le_11:,} ({100 * n_le_11 / total_endpoints:.1f}%)")
    print(f"endpoints with score <=-12: {n_le_12:,} ({100 * n_le_12 / total_endpoints:.1f}%)")
    print(f"endpoints with std <=0.1: {n_std_le_01:,} ({100 * n_std_le_01 / total_endpoints:.1f}%)")
    print(f"endpoints with std <=0.5: {n_std_le_05:,} ({100 * n_std_le_05 / total_endpoints:.1f}%)")
    print(f"endpoints with std <=1.0: {n_std_le_10:,} ({100 * n_std_le_10 / total_endpoints:.1f}%)")

    result = result.drop(columns=std_cols)

    out_path = os.path.join(OUTPUT_DIR, "merged_docking_scores.csv")
    result.to_csv(out_path, index=False)
    print(f"\nSaved {len(result):,} rows x {len(result.columns)} columns to {out_path}")


if __name__ == "__main__":
    main()
