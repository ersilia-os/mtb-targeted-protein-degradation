#!/usr/bin/env python3
"""
Multi-target hit-overlap analysis restricted to the 4 catalytic (CAT) pockets - pheS, aspS, lysS,
alaS (pheT has no CAT entry) - against the second Enamine REAL screening
(output/unidock_REAL_docking_2, ~99,105 compounds).

Uses the shared helpers in src/docking_utils.py. Reference pockets are loaded directly from
output/reference_pocket_catalytic.csv (output/reference_pocket.csv no longer exists, having been
split into output/reference_pocket_catalytic.csv / _noncatalytic.csv).

For each of top-100, top-1,000, top-10,000:
  - an UpSet plot of hit overlap across the 4 genes (docking_utils.save_upset)
  - a printed observed-vs-expected-by-chance tally of compounds hitting 2, 3, or 4 targets
  - a CSV of every compound hitting >=2 of the 4 targets within that top-N, with per-gene
    score/rank, number of targets, SMILES, and physchem properties (MW, cLogP, TPSA, HBD, HBA,
    RotBonds, AromaticRings, QED, is_pains)

Then, for the >=2-target compounds from the top-1,000 CSV vs. a random 10,000-compound background
sample from the same REAL library (docking_utils.sample_prescreened_smiles, seed=RANDOM_SEED):
  - a physchem profiling comparison (docking_utils.plot_profiling)
  - a t-SNE projection of ECFP4 fingerprints (docking_utils.plot_tsne), background in gray,
    selected in purple

Also, per-gene raw docking score boxplots (docking_utils.plot_score_boxplots) of this pocket's
pre-screened (round-2 REAL, ~99,105 compounds) vs. selected (>=2-target, top-1,000) scores,
alongside two external reference distributions for the same 4 pockets: the Enamine Hit Locator
100k library (output/unidock_docking, the first screening round) and the Enamine REAL 9.56M
screening's "negative set" (output/unidock_REAL_docking, ~12,958 compounds per pocket - everything
beyond each pocket's top-100,000 prioritized compounds from that round).

Usage:
    python 52_CAT_hit_overlap.py
"""
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from default import RANDOM_SEED
from docking_utils import (
    LIBRARIES, build_matrix, compute_properties, load_real_negative_scores, lookup_smiles,
    plot_profiling, plot_score_boxplots, plot_tsne, sample_prescreened_smiles, save_upset,
)

CAT_POCKETS_CSV = os.path.join(ROOT, "output", "reference_pocket_catalytic.csv")
OUTPUT_DIR = os.path.join(ROOT, "output", "52_CAT_promiscuous")

LIB = "REAL"
TOP_NS = [100, 1_000, 10_000]
PROFILING_TOP_N = 1_000
PROFILING_MIN_TARGETS = 2
BG_SAMPLE_SIZE = 20_000
BOXPLOT_GENE_ORDER = ["pheS", "aspS", "lysS", "alaS"]

PROP_COLUMNS = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotBonds", "AromaticRings", "QED", "is_pains"]


def load_cat_pockets():
    """gene_name -> pocket_name for the 4 catalytic pockets."""
    df = pd.read_csv(CAT_POCKETS_CSV)
    return dict(zip(df["gene_name"], df["pocket_name"]))


def build_multi_target_table(genes, gene_top, rank_of, scores1, top_n):
    """
    Builds a DataFrame of every compound hitting >=2 of the 4 genes' top-top_n sets, with per-gene
    score/rank (NaN if not in that gene's top-top_n), n_targets, SMILES and physchem properties.
    Also returns the raw observed-vs-expected-by-chance tally (k -> compound count) for printing.
    """
    top_sets = {g: set(gene_top[g][:top_n]) for g in genes}
    all_in_top = set.union(*top_sets.values())
    n_targets_of = {cid: sum(1 for g in genes if cid in top_sets[g]) for cid in all_in_top}
    tally = Counter(n_targets_of.values())

    multi_ids = [cid for cid in all_in_top if n_targets_of[cid] >= 2]
    smiles_map = lookup_smiles(multi_ids, LIB)

    rows = []
    for cid in multi_ids:
        row = {"compound": cid, "smiles": smiles_map.get(cid, ""), "n_targets": n_targets_of[cid]}
        for g in genes:
            in_top = cid in top_sets[g]
            row[f"score_{g}"] = scores1.loc[cid, g] if in_top else np.nan
            row[f"rank_{g}"] = rank_of[g].get(cid, np.nan) if in_top else np.nan
        rows.append(row)
    multi_df = pd.DataFrame(rows)

    props = compute_properties(smiles_map)
    multi_df = multi_df.merge(props, left_on="compound", right_index=True, how="left")
    return multi_df, tally


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    pocket_map = load_cat_pockets()
    genes = sorted(pocket_map)
    trna_tag = "_".join(genes)
    print(f"Catalytic pockets: {pocket_map}")

    scores1 = build_matrix(pocket_map, LIBRARIES[LIB], label="Loading CAT pocket scores")
    M = len(scores1)

    gene_top = {g: scores1[g].dropna().sort_values(ascending=True).index.tolist() for g in genes}
    rank_of = {g: {cid: i + 1 for i, cid in enumerate(gene_top[g])} for g in genes}

    multi_target_dfs = {}
    for top_n in TOP_NS:
        print(f"\n--- Top-{top_n} ---")
        save_upset(gene_top, top_n, OUTPUT_DIR, LIB, f"{trna_tag}_CAT")

        multi_df, tally = build_multi_target_table(genes, gene_top, rank_of, scores1, top_n)
        for k in range(len(genes), 1, -1):
            expected = (top_n ** k) / (M ** (k - 1))
            print(f"  {k}/{len(genes)} targets: {tally.get(k, 0):,} compound(s)  (expected {expected:.1f})")

        out_csv = os.path.join(OUTPUT_DIR, f"{trna_tag}_{LIB}_CAT_multi_target_top{top_n}.csv")
        multi_df = multi_df.sort_values("n_targets", ascending=False)
        multi_df.to_csv(out_csv, index=False)
        print(f"  Saved: {out_csv}")
        multi_target_dfs[top_n] = multi_df

    sel_df = multi_target_dfs[PROFILING_TOP_N]
    sel_df = sel_df[sel_df["n_targets"] >= PROFILING_MIN_TARGETS]
    sel_ids = set(sel_df["compound"])

    print("\n--- Score boxplots ---")
    dl_scores = build_matrix(pocket_map, LIBRARIES["DL"], label="Loading DL (100k) scores")
    real_negative_scores = {g: load_real_negative_scores(pocket_map[g]) for g in genes}
    scores_path = os.path.join(OUTPUT_DIR, f"{trna_tag}_{LIB}_CAT_scores.png")
    plot_score_boxplots(scores1, sel_ids, BOXPLOT_GENE_ORDER, scores_path, dl_scores=dl_scores, real_negative_scores=real_negative_scores)
    print(f"Saved: {scores_path}")

    print("\n--- Physchem profiling ---")
    sel_props = sel_df.set_index("compound")[PROP_COLUMNS]
    sel_smiles = dict(zip(sel_df["compound"], sel_df["smiles"]))

    bg_smiles = sample_prescreened_smiles(LIB, sel_ids, BG_SAMPLE_SIZE, RANDOM_SEED)
    bg_props = compute_properties(bg_smiles)

    profiling_path = os.path.join(OUTPUT_DIR, f"{trna_tag}_{LIB}_CAT_profiling.png")
    plot_profiling(sel_props, bg_props, profiling_path)
    print(f"Saved: {profiling_path}")

    print("\n--- t-SNE (ECFP4) ---")
    tsne_path = os.path.join(OUTPUT_DIR, f"{trna_tag}_{LIB}_CAT_tsne.png")
    plot_tsne(sel_smiles, bg_smiles, tsne_path, seed=RANDOM_SEED)
    print(f"Saved: {tsne_path}")


if __name__ == "__main__":
    main()
