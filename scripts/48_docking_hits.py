#!/usr/bin/env python3
"""
For each gene, identify the reference pocket (highest P2Rank probability),
then count docking hits from the specified library at score cutoffs -15 to -9.

Usage:
    python 48_docking_hits.py --trna pheS --lib REAL
    python 48_docking_hits.py --trna pheS,aspS --lib DL
"""

import argparse
import os
import sys
from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from upsetplot import from_contents, plot as upset_plot

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

LIBRARIES = {
    "DL":   os.path.join(ROOT, "output", "unidock_docking",        "docking_results"),
    "REAL": os.path.join(ROOT, "output", "unidock_REAL_docking_2", "docking_results"),
}

CUTOFFS = list(range(-15, -8))  # -15, -14, ..., -9


def load_gene_map():
    path = os.path.join(ROOT, "data", "mtb_trna_synthetases_bosch_2021_fig5.csv")
    df = pd.read_csv(path)
    return {row["gene_name_in_bosch_2021"]: row["uniprot_ac"] for _, row in df.iterrows()}


def load_pocket_data():
    path = os.path.join(ROOT, "output", "pocket_detection_data.csv")
    df = pd.read_csv(path)
    probs = {}
    for _, row in df.iterrows():
        key = f"{row['File name'].replace('.pdb', '')}_pocket_{row['Pocket number']}"
        probs[key] = row["Pocket probability"]
    return probs


def save_upset(genes, gene_compounds, top_n, output_dir):
    top_sets = {g: set(gene_compounds[g][:top_n]) for g in genes}
    data = from_contents(top_sets)
    upset_plot(data)
    fname = "_".join(sorted(genes)) + f"_top{top_n}.png"
    path = os.path.join(output_dir, fname)
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {path}")


def reference_pocket(uniprot_ac, pocket_probs, results_dir):
    candidates = [p for p in os.listdir(results_dir) if uniprot_ac in p]
    if not candidates:
        return None
    return max(candidates, key=lambda p: pocket_probs.get(p, 0))


def main():
    parser = argparse.ArgumentParser(description="Count docking hits for the reference pocket of each gene.")
    parser.add_argument("--trna", required=True,
                        help="Comma-separated gene name(s) (e.g. pheS or pheS,aspS)")
    parser.add_argument("--lib", choices=["DL", "REAL"], required=True,
                        help="Compound library: DL (Enamine DL-HLL-100) or REAL (Enamine REAL 10B)")
    parser.add_argument("--silent-cutoffs", action="store_true",
                        help="Suppress the per-cutoff hit table")
    parser.add_argument("--top", default=None,
                        help="Comma-separated integers: report shared hits across all gene combinations at each top-N threshold (e.g. 100,500,1000)")
    args = parser.parse_args()

    genes = [g.strip() for g in args.trna.split(",")]
    cutoffs = CUTOFFS
    top_values = [int(x.strip()) for x in args.top.split(",")] if args.top else []

    gene_map = load_gene_map()
    pocket_probs = load_pocket_data()
    results_dir = LIBRARIES[args.lib]

    invalid = [g for g in genes if g not in gene_map]
    if invalid:
        print(f"Error: unknown gene(s): {', '.join(invalid)}. Available genes:")
        print("  " + ", ".join(sorted(gene_map)))
        sys.exit(1)

    gene_compounds = {}  # gene -> sorted compound array (best first)

    for gene in genes:
        uniprot_ac = gene_map[gene]
        ref = reference_pocket(uniprot_ac, pocket_probs, results_dir)
        if ref is None:
            print(f"[{gene}] No docking results found in {args.lib}.")
            print()
            continue

        prob = pocket_probs.get(ref, None)
        report_path = os.path.join(results_dir, ref, "report.csv")
        df_report = pd.read_csv(report_path).sort_values("score")
        scores = df_report["score"].values
        n_total = len(scores)

        print(f"Gene     : {gene} ({uniprot_ac})")
        print(f"Pocket   : {ref}  (prob: {round(prob, 3) if prob is not None else 'N/A'})")
        print(f"Library  : {args.lib}  ({n_total} compounds)")
        if not args.silent_cutoffs:
            print(f"{'Cutoff':<10} {'Hits':>8} {'%':>8}")
            print("-" * 30)
            for cutoff in cutoffs:
                hits = int(np.sum(scores < cutoff))
                pct = 100 * hits / n_total if n_total > 0 else 0.0
                print(f"{cutoff:<10} {hits:>8} {pct:>7.1f}%")
        print()

        gene_compounds[gene] = df_report["compound"].values  # already sorted by score

    if top_values and len(gene_compounds) < 2:
        print("Nothing to share: provide at least two genes to compute shared hits.")
    elif top_values:
        available = [g for g in genes if g in gene_compounds]
        output_dir = os.path.join(ROOT, "output", "hit_overlap", args.lib)
        os.makedirs(output_dir, exist_ok=True)
        for top_n in top_values:
            print(f"Shared hits — top {top_n}:")
            top_sets = {g: set(gene_compounds[g][:top_n]) for g in available}
            for r in range(2, len(available) + 1):
                for combo in combinations(available, r):
                    shared = len(set.intersection(*(top_sets[g] for g in combo)))
                    label = " & ".join(combo)
                    print(f"  {label}: {shared}")
            print()
            save_upset(available, gene_compounds, top_n, output_dir)


if __name__ == "__main__":
    main()
