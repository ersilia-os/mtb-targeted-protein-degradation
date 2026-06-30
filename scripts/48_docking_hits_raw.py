#!/usr/bin/env python3
"""
Raw docking hit overlap analysis for tRNA synthetase targets.

For the provided genes, builds Matrix 1 (compounds × reference pockets),
then reports:
  - Shared-hit tables + UpSet plots at top 100 and top 1,000
  - Multi-target binder counts (compounds hitting k/N genes at top 100 and 1,000)
  - CSV of compounds with ≥2 targets at top 1,000 with per-gene score and rank

Usage:
    python 48_docking_hits_raw.py --trna pheS,aspS,lysS,alaS --lib REAL
"""

import argparse
import os
import sys
import warnings
from collections import Counter
from itertools import combinations

import matplotlib
matplotlib.use("Agg")
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*converter.*already registered.*")
warnings.filterwarnings("ignore", category=FutureWarning, module="upsetplot")

import matplotlib.pyplot as plt
import pandas as pd
from upsetplot import from_contents, plot as upset_plot

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from default import RANDOM_SEED
from docking_utils import (
    LIBRARIES,
    build_matrix,
    compute_properties,
    load_gene_map,
    load_reference_pockets,
    lookup_smiles,
    plot_profiling,
    plot_score_boxplots,
    sample_background_smiles,
)

BG_SAMPLE_SIZE = 25_000
UPSET_THRESHOLDS = [100, 1_000]


def save_upset(gene_top, top_n, output_dir, lib, trna_tag):
    """UpSet plot of top-N hit overlap across target genes."""
    top_sets = {g: set(ids[:top_n]) for g, ids in gene_top.items()}
    data = from_contents(top_sets)
    upset_plot(data)
    plt.suptitle(f"{lib} — top {top_n:,}")
    fname = f"{trna_tag}_{lib}_upset_top{top_n}.png"
    path = os.path.join(output_dir, fname)
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved: {path}")


def main():
    parser = argparse.ArgumentParser(
        description="Raw docking hit overlap analysis across tRNA synthetase targets."
    )
    parser.add_argument("--trna", required=True,
                        help="Comma-separated gene name(s) (e.g. pheS,aspS,lysS,alaS)")
    parser.add_argument("--lib", choices=["DL", "REAL"], required=True,
                        help="Compound library: DL or REAL")
    args = parser.parse_args()

    genes = [g.strip() for g in args.trna.split(",")]
    results_dir = LIBRARIES[args.lib]

    gene_map = load_gene_map()
    ref_pockets = load_reference_pockets()

    invalid = [g for g in genes if g not in gene_map]
    if invalid:
        print(f"Error: unknown gene(s): {', '.join(invalid)}. Available genes:")
        print("  " + ", ".join(sorted(gene_map)))
        sys.exit(1)

    missing_ref = [g for g in genes if g not in ref_pockets]
    if missing_ref:
        print(f"Error: no reference pocket defined for: {', '.join(missing_ref)}")
        sys.exit(1)

    output_dir = os.path.join(ROOT, "output", "48_docking_hits")
    os.makedirs(output_dir, exist_ok=True)
    trna_tag = "_".join(genes)

    print(f"Library : {args.lib}")
    print(f"Targets : {', '.join(genes)}")

    target_pocket_map = {g: ref_pockets[g] for g in genes}
    print(f"Building Matrix 1 ({len(target_pocket_map)} reference pocket(s))...")
    scores1 = build_matrix(target_pocket_map, results_dir, label="Matrix 1")

    gene_top = {g: scores1[g].sort_values().index.tolist() for g in genes if g in scores1.columns}

    if len(gene_top) < 2:
        print("Need at least two genes with docking results for overlap analysis.")
        sys.exit(0)

    # --- Raw hit overlap ---
    print("\n--- Raw hit overlap ---")
    for top_n in UPSET_THRESHOLDS:
        top_sets = {g: set(gene_top[g][:top_n]) for g in gene_top}
        print(f"\nShared hits — top {top_n:,}:")
        for r in range(2, len(gene_top) + 1):
            for combo in combinations(sorted(gene_top), r):
                shared = len(set.intersection(*(top_sets[g] for g in combo)))
                print(f"  {' & '.join(combo)}: {shared:,}")
        save_upset(gene_top, top_n, output_dir, args.lib, trna_tag)
    print()

    # --- Multi-target binder counts ---
    M = len(scores1)
    print("--- Multi-target binders ---")
    for top_n in [100, 1_000]:
        top_sets = {g: set(gene_top[g][:top_n]) for g in gene_top}
        all_in_top = set.union(*top_sets.values())
        tally = Counter(
            sum(1 for g in gene_top if cid in top_sets[g])
            for cid in all_in_top
        )
        print(f"Top {top_n:,}:")
        for k in range(len(gene_top), 1, -1):
            expected = (top_n ** k) / (M ** (k - 1))
            print(f"  {k}/{len(gene_top)} targets: {tally.get(k, 0):,} compound(s)  (expected {expected:.1f})")
    print()

    # --- Multi-target CSV report (≥2 targets at top 1,000) ---
    top1k_sets = {g: set(gene_top[g][:1_000]) for g in gene_top}
    multi_ids = [
        cid for cid in set.union(*top1k_sets.values())
        if sum(1 for g in gene_top if cid in top1k_sets[g]) >= 2
    ]
    if multi_ids:
        gene_ranks = {g: scores1[g].rank(ascending=True, method="first") - 1 for g in gene_top}
        rows = []
        for cid in multi_ids:
            row = {"compound": cid,
                   "n_targets": sum(1 for g in gene_top if cid in top1k_sets[g])}
            for g in genes:
                if g in gene_top:
                    row[f"score_{g}"] = scores1.loc[cid, g]
                    row[f"rank_{g}"]  = int(gene_ranks[g][cid])
                else:
                    row[f"score_{g}"] = float("nan")
                    row[f"rank_{g}"]  = pd.NA
            rows.append(row)
        multi_df = pd.DataFrame(rows).set_index("compound")
        rank_cols = [f"rank_{g}" for g in genes]
        multi_df["_rank_sum"] = multi_df[rank_cols].sum(axis=1)
        multi_df = multi_df.sort_values(["n_targets", "_rank_sum"], ascending=[False, True])
        multi_df = multi_df.drop(columns=["_rank_sum"])
        sel_smiles = lookup_smiles(multi_df.index.tolist(), args.lib)
        multi_df.insert(0, "smiles", multi_df.index.map(sel_smiles))
        print(f"  Computing properties for {len(sel_smiles):,} selected...")
        sel_props = compute_properties(sel_smiles)
        multi_df = multi_df.join(sel_props, how="left")
        multi_path = os.path.join(output_dir, f"{trna_tag}_{args.lib}_multi_target_top1000.csv")
        multi_df.to_csv(multi_path)
        print(f"Multi-target report: {len(multi_df):,} compounds → {multi_path}")

        scores_path = os.path.join(output_dir, f"{trna_tag}_{args.lib}_scores.png")
        plot_score_boxplots(scores1, set(multi_ids), genes, scores_path,
                            sel_label="Multi-target ≥2 targets, top 1,000")
        print(f"  Saved: {scores_path}")

        print("\nPhysicochemical profiling...")
        bg_smiles  = sample_background_smiles(args.lib, set(multi_ids), BG_SAMPLE_SIZE, RANDOM_SEED)
        print(f"  Background: {len(bg_smiles):,} compounds...")
        bg_props  = compute_properties(bg_smiles)
        profiling_path = os.path.join(output_dir, f"{trna_tag}_{args.lib}_profiling.png")
        plot_profiling(sel_props, bg_props, profiling_path)
        print(f"  Saved profiling figure to {profiling_path}")


if __name__ == "__main__":
    main()
