#!/usr/bin/env python3
"""
Merge docking scores from the 1st (10M/9.56M-derived) and 2nd (10B-derived) Enamine
REAL docking campaigns for every NON-CAT pocket (output/selected_pockets.csv) into one
ranked compound table per pocket, then select the top 10,000 compounds per pocket while
capping chemical redundancy (max MAX_SYN occurrences of any given synthon among the
selected compounds). Finally, combine every pocket's top-10k selection into one file.

Both campaigns store one report.csv (compound,score) per pocket. For a given pocket,
compounds present in both round-1 and round-2 keep their best (lowest) score and are
labeled "Both" instead of "10M" or "10B". One NON-CAT pocket (7K98_pocket_6, docked
against a real multimeric structure) has no round-1 data, since round-1 docking predates
the multimer-detection effort — for that pocket, every row is source "10B" and the pool
is round-2's ~99,105 compounds instead of the ~213k combined pool.

Usage:
    python 54_merge_scores_select_hits.py
"""

import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

import pandas as pd

from docking_utils import LIBRARIES, REAL_ROUND1_RESULTS_DIR, load_scores, lookup_smiles

SELECTED_POCKETS_CSV = os.path.join(root, "..", "output", "selected_pockets.csv")
MULTIMER_DOCKING_DIR = os.path.join(root, "..", "output", "50_unidock_docking_multimers")

MAX_SYN = 3
TARGET_N = 10_000

output_dir = os.path.join(root, "..", "output", "54_merge_scores_select_hits")
os.makedirs(output_dir, exist_ok=True)


def load_noncat_pockets():
    df = pd.read_csv(SELECTED_POCKETS_CSV)
    df = df[df["site_type"] == "NON-CAT"]
    return list(zip(df["gene_name"], df["pocket_name"]))


def resolve_round2_dir(pocket_name):
    """Same "_model_" heuristic as scripts/53_NONCAT_hit_overlap.py: monomeric pockets
    live under the shared REAL round-2 results dir, the multimeric one (7K98_pocket_6)
    under its own script-50 output tree."""
    return LIBRARIES["REAL"] if "_model_" in pocket_name else MULTIMER_DOCKING_DIR


def get_synthons(compound_id):
    return compound_id.split("____")[1:]


def select_diverse_top_n(df, max_syn=MAX_SYN, target_n=TARGET_N):
    """df must already be sorted ascending by docking_score. Returns the full df (all
    rows), with observed_synthons/select computed throughout; select stops being set
    to 1 once target_n compounds have been selected, but observed_synthons keeps
    accumulating over every remaining row."""
    all_seen_counts = {}
    selected_counts = {}
    observed_col, select_col = [], []
    n_selected = 0
    n_considered_to_reach_target = None

    for i, compound_id in enumerate(df["compound_id"]):
        synthons = get_synthons(compound_id)
        observed = [all_seen_counts.get(s, 0) for s in synthons]
        can_select = n_selected < target_n and all(selected_counts.get(s, 0) < max_syn for s in synthons)

        observed_col.append(";".join(map(str, observed)))
        select_col.append(1 if can_select else 0)

        for s in synthons:
            all_seen_counts[s] = all_seen_counts.get(s, 0) + 1

        if can_select:
            for s in synthons:
                selected_counts[s] = selected_counts.get(s, 0) + 1
            n_selected += 1
            if n_selected == target_n:
                n_considered_to_reach_target = i + 1

    print(f"  Considered {n_considered_to_reach_target:,} molecules to select {target_n:,} compounds (MAX_SYN={max_syn}).")

    out_df = df.copy()
    out_df["observed_synthons"] = observed_col
    out_df["select"] = select_col
    return out_df


def process_pocket(pocket_name):
    round1_report = os.path.join(REAL_ROUND1_RESULTS_DIR, pocket_name, "report.csv")
    round2_report = os.path.join(resolve_round2_dir(pocket_name), pocket_name, "report.csv")

    s1 = load_scores(round1_report) if os.path.isfile(round1_report) else None
    if s1 is not None:
        print(f"  Round-1 (10M): {len(s1):,} compounds")
    else:
        print("  Round-1 (10M): not available for this pocket (multimer), skipping merge.")

    s2 = load_scores(round2_report)
    print(f"  Round-2 (10B): {len(s2):,} compounds")

    if s1 is None:
        df = pd.DataFrame({"compound_id": s2.index, "source": "10B", "docking_score": s2.values})
    else:
        only_1 = s1.index.difference(s2.index)
        only_2 = s2.index.difference(s1.index)
        both = s1.index.intersection(s2.index)
        print(f"  10M-only: {len(only_1):,}  10B-only: {len(only_2):,}  Both: {len(both):,}")

        rows = []
        rows.append(pd.DataFrame({"compound_id": only_1, "source": "10M", "docking_score": s1[only_1].values}))
        rows.append(pd.DataFrame({"compound_id": only_2, "source": "10B", "docking_score": s2[only_2].values}))
        best_both = pd.concat([s1[both], s2[both]], axis=1).min(axis=1)
        rows.append(pd.DataFrame({"compound_id": both, "source": "Both", "docking_score": best_both.values}))
        df = pd.concat(rows, ignore_index=True)

    real_ids = df.loc[df["source"].isin(["10B", "Both"]), "compound_id"].tolist()
    real_smiles = lookup_smiles(real_ids, "REAL")

    round1_only_ids = df.loc[df["source"] == "10M", "compound_id"].tolist()
    round1_smiles = lookup_smiles(round1_only_ids, "REAL_ROUND1") if round1_only_ids else {}

    smiles_map = {**round1_smiles, **real_smiles}
    df["smiles"] = df["compound_id"].map(smiles_map)

    df = df.sort_values("docking_score", ascending=True).reset_index(drop=True)
    df = df[["compound_id", "smiles", "source", "docking_score"]]

    print(f"  Total rows: {len(df):,}  {dict(df['source'].value_counts())}")
    print(f"  Best (min) score: {df['docking_score'].min()}  Worst (max) score: {df['docking_score'].max()}")

    df = select_diverse_top_n(df)

    out_path = os.path.join(output_dir, f"{pocket_name}_merged_docking_scores.csv")
    df.to_csv(out_path, index=False)
    print(f"  Saved: {out_path}")

    return df


def main():
    pockets = load_noncat_pockets()
    print(f"NON-CAT pockets: {len(pockets)}")

    combined_rows = []
    for gene_name, pocket_name in pockets:
        print(f"\n=== {gene_name} / {pocket_name} ===")
        df = process_pocket(pocket_name)
        selected = df[df["select"] == 1].copy()
        selected["pocket_name"] = pocket_name
        combined_rows.append(selected)

    combined = pd.concat(combined_rows, ignore_index=True)
    combined_path = os.path.join(output_dir, "NONCAT_top10k_selected_combined.csv")
    combined.to_csv(combined_path, index=False)
    print(f"\nCombined top-10k selections across {len(pockets)} pockets: {len(combined):,} rows")
    print(f"Saved: {combined_path}")

    n_unique = combined["compound_id"].nunique()
    pockets_per_compound = combined["compound_id"].value_counts()
    n_duplicated = (pockets_per_compound > 1).sum()
    print(f"\nUnique compounds: {n_unique:,} / {len(combined):,} rows ({n_duplicated:,} compounds hit more than one pocket)")
    print("Compounds by number of pockets hit:")
    print(pockets_per_compound.value_counts().sort_index().rename_axis("n_pockets").rename("n_compounds"))


if __name__ == "__main__":
    main()
