#!/usr/bin/env python3
"""
Merge script 65's 12 per-pocket report.csv files into one wide table: compound_id, smiles,
<pocket_1>, <pocket_1>_rank, <pocket_2>, <pocket_2>_rank, ... - score and rank vs. the Enamine REAL
round-2 library (~99,105 compounds) for that pocket. Also reports reproducibility against whatever
old/precalculated scores already existed for the same (compound, pocket) pairs, from prior docking
rounds (round-1, round-2, the dimer's multimer round, and the NON-CAT REAL10B-selective round).

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
REAL_ROUND1_DIR = os.path.join(ROOT, "output", "unidock_REAL_docking", "docking_results")
REAL10B_SELECTIVE_DIR = os.path.join(ROOT, "output", "60_unidock_docking_noncat_selective_10B", "docking_results")
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


def load_old_scores(pockets):
    """{(compound_id, pocket_name): score} from whatever historical docking round already
    scored that pair - round-2/multimer/round-1/REAL10B-selective, first found wins."""
    old_scores = {}
    for pocket_name in pockets:
        if pocket_name == DIMER_POCKET:
            for cid, s in load_report_scores(os.path.join(MULTIMER_DIR, pocket_name, "report.csv")).items():
                old_scores[(cid, pocket_name)] = s
            continue
        for base in (REAL_ROUND2_DIR, REAL_ROUND1_DIR, REAL10B_SELECTIVE_DIR):
            for cid, s in load_report_scores(os.path.join(base, pocket_name, "report.csv")).items():
                old_scores.setdefault((cid, pocket_name), s)
    return old_scores


def main():
    result = pd.read_csv(AGGREGATED_HITS_CSV, usecols=["compound_id", "smiles", "source"])
    n_compounds = len(result)
    print(f"Base compound list: {n_compounds:,} compounds")

    pockets = pd.read_csv(SELECTED_POCKETS_CSV)["pocket_name"].tolist()
    print(f"Pockets: {len(pockets)}")

    background = load_round2_background(pockets)

    for pocket_name in pockets:
        report_path = os.path.join(DOCKING_RESULTS_DIR, pocket_name, "report.csv")
        report = pd.read_csv(report_path)
        report = report.rename(columns={"compound": "compound_id", "score": pocket_name})
        report[pocket_name] = report[pocket_name].astype(float)

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

    out_path = os.path.join(OUTPUT_DIR, "merged_docking_scores.csv")
    result.to_csv(out_path, index=False)
    print(f"\nSaved {len(result):,} rows x {len(result.columns)} columns to {out_path}")

    print("\nChecking reproducibility against precalculated scores from prior docking rounds...")
    old_scores = load_old_scores(pockets)
    diffs = []
    n_excluded = 0
    for _, row in result.iterrows():
        for pocket_name in pockets:
            old_score = old_scores.get((row["compound_id"], pocket_name))
            if old_score is None:
                continue
            new_score = row[pocket_name]
            if old_score > 0 or new_score > 0:  # failed/clashing pose artifact, not a real energy
                n_excluded += 1
                continue
            diffs.append(abs(new_score - old_score))

    diffs = pd.Series(diffs)
    print(f"  {len(diffs) + n_excluded:,} (compound, pocket) pairs had a precalculated score to compare against")
    print(f"  {n_excluded:,} excluded (failed/clashing-pose artifact - nonsensical positive energy on either side)")
    print("\nreproducibility to previous runs")
    print(f"<0.1: {100 * (diffs <= 0.1).sum() / len(diffs):.1f}%")
    print(f"<0.5: {100 * (diffs <= 0.5).sum() / len(diffs):.1f}%")
    print(f"<1.0: {100 * (diffs <= 1.0).sum() / len(diffs):.1f}%")


if __name__ == "__main__":
    main()
