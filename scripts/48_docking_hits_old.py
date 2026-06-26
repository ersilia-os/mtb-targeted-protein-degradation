#!/usr/bin/env python3
"""
For each gene, identify the reference pocket (highest P2Rank probability),
then count docking hits from the specified library at score cutoffs -15 to -9.

Usage:
    python 48_docking_hits.py --trna pheS --lib REAL
    python 48_docking_hits.py --trna pheS,aspS --lib DL
    python 48_docking_hits.py --trna pheS,aspS,lysS,alaS --lib REAL --top 1000 --profile
"""

import argparse
import os
import sys
import warnings
from collections import Counter
from itertools import combinations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, mannwhitneyu
from tqdm import tqdm
from rdkit.Chem import MolFromSmiles, rdMolDescriptors
from rdkit.Chem.QED import qed as calc_qed
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
warnings.filterwarnings("ignore", category=FutureWarning, module="upsetplot")
from upsetplot import from_contents, plot as upset_plot

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

LIBRARIES = {
    "DL":   os.path.join(ROOT, "output", "unidock_docking",        "docking_results"),
    "REAL": os.path.join(ROOT, "output", "unidock_REAL_docking_2", "docking_results"),
}

COMPOUNDS_PATH = os.path.join(
    ROOT, "output", "unidock_REAL_docking", "inference_10B", "clustered_compounds.csv"
)

CUTOFFS = list(range(-15, -8))  # -15, -14, ..., -9
RANDOM_SEED = 42
BG_SAMPLE_SIZE = 5_000

PROP_COLUMNS = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotBonds", "AromaticRings", "QED"]


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


def save_upset(genes, gene_compounds, top_n, output_dir, lib):
    top_sets = {g: set(gene_compounds[g][:top_n]) for g in genes}
    data = from_contents(top_sets)
    upset_plot(data)
    plt.suptitle(f"{lib} — top {top_n}")
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


def compute_properties(smiles_dict):
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    catalog = FilterCatalog(params)

    records = []
    for cid, smi in tqdm(smiles_dict.items(), desc="Computing properties"):
        mol = MolFromSmiles(smi)
        if mol is None:
            continue
        records.append({
            "id": cid,
            "MW": rdMolDescriptors.CalcExactMolWt(mol),
            "cLogP": rdMolDescriptors.CalcCrippenDescriptors(mol)[0],
            "TPSA": rdMolDescriptors.CalcTPSA(mol),
            "HBD": rdMolDescriptors.CalcNumHBD(mol),
            "HBA": rdMolDescriptors.CalcNumHBA(mol),
            "RotBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "QED": calc_qed(mol),
            "is_pains": catalog.HasMatch(mol),
        })
    return pd.DataFrame(records).set_index("id")


def run_profiling(genes, gene_compounds, top_n, props_df, lib, output_dir):
    import stylia

    n_genes = len(genes)
    group_labels = {0: "Background"} | {i: f"{i} target{'s' if i > 1 else ''}" for i in range(1, n_genes + 1)}

    top_sets = {g: set(gene_compounds[g][:top_n]) for g in genes}
    target_count = Counter()
    for compounds in top_sets.values():
        for c in compounds:
            target_count[c] += 1

    props = props_df.copy()
    props["group"] = props.index.map(lambda x: target_count.get(x, 0))

    group_order = list(range(n_genes, -1, -1))  # [4,3,2,1,0]
    pal = stylia.CategoricalPalette("npg")
    colors = pal.get(len(group_order))
    color_map = {g: colors[i] for i, g in enumerate(group_order)}

    # --- Distribution plots (2×4 grid, one panel per property) ---
    fig, axes = plt.subplots(2, 4, figsize=(12, 6))
    axes = axes.flatten()
    for idx, prop in enumerate(PROP_COLUMNS):
        ax = axes[idx]
        for g in group_order:
            data = props[props["group"] == g][prop].dropna().values
            if len(data) < 2:
                continue
            kde = gaussian_kde(data)
            x = np.linspace(data.min(), data.max(), 300)
            ax.plot(x, kde(x), color=color_map[g], label=group_labels[g], lw=1.5)
        ax.set_xlabel(prop)
        ax.set_ylabel("Density")
        ax.legend(fontsize=7)
    fig.tight_layout()
    dist_path = os.path.join(output_dir, f"{lib}_top{top_n}_distributions.png")
    fig.savefig(dist_path, dpi=150, bbox_inches="tight")
    plt.close("all")

    # --- PAINS bar chart ---
    fig, ax = plt.subplots(figsize=(6, 4))
    bar_labels, bar_pcts, bar_colors = [], [], []
    for g in group_order:
        subset = props[props["group"] == g]
        if len(subset) == 0:
            continue
        bar_labels.append(group_labels[g])
        bar_pcts.append(100 * subset["is_pains"].sum() / len(subset))
        bar_colors.append(color_map[g])
    ax.bar(range(len(bar_labels)), bar_pcts, color=bar_colors)
    ax.set_xticks(range(len(bar_labels)))
    ax.set_xticklabels(bar_labels, rotation=45, ha="right")
    ax.set_ylabel("PAINS (%)")
    fig.tight_layout()
    pains_path = os.path.join(output_dir, f"{lib}_top{top_n}_pains.png")
    fig.savefig(pains_path, dpi=150, bbox_inches="tight")
    plt.close("all")

    # --- Stats table (Mann-Whitney U, one-sided: hits > background) ---
    bg_data = props[props["group"] == 0]
    stats_records = []
    for prop in PROP_COLUMNS:
        bg_vals = bg_data[prop].dropna().values
        for g in range(1, n_genes + 1):
            hit_vals = props[props["group"] == g][prop].dropna().values
            if len(hit_vals) < 2:
                continue
            stat, pval = mannwhitneyu(hit_vals, bg_vals, alternative="greater")
            stats_records.append({
                "property": prop,
                "group": group_labels[g],
                "n": len(hit_vals),
                "median_hits": float(np.median(hit_vals)),
                "median_bg": float(np.median(bg_vals)),
                "U": stat,
                "p_value": pval,
            })
    stats_df = pd.DataFrame(stats_records)
    stats_path = os.path.join(output_dir, f"{lib}_top{top_n}_stats.csv")
    stats_df.to_csv(stats_path, index=False)

    print(f"\nPhysicochemical profiling — top {top_n}:")
    print(f"  {'Property':<15} {'Group':<12} {'N':>6} {'Med(hits)':>10} {'Med(bg)':>10} {'p-value':>12}")
    print("  " + "-" * 68)
    for _, row in stats_df.iterrows():
        print(f"  {row['property']:<15} {row['group']:<12} {int(row['n']):>6} "
              f"{row['median_hits']:>10.3f} {row['median_bg']:>10.3f} {row['p_value']:>12.4f}")

    print(f"\n  Saved: {dist_path}")
    print(f"  Saved: {pains_path}")
    print(f"  Saved: {stats_path}")


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
    parser.add_argument("--profile", action="store_true",
                        help="Run physicochemical profiling of hit groups (requires --top and --lib REAL)")
    args = parser.parse_args()

    if args.profile and not args.top:
        print("Error: --profile requires --top.")
        sys.exit(1)
    if args.profile and args.lib != "REAL":
        print("Error: --profile is only supported for --lib REAL.")
        sys.exit(1)

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
            target_count = Counter()
            for compounds in top_sets.values():
                for c in compounds:
                    target_count[c] += 1
            tally = Counter(target_count.values())
            print(f"  Target coverage summary (top {top_n}):")
            for n in range(len(available), 0, -1):
                print(f"    {n} target(s): {tally.get(n, 0)} compound(s)")
            print()
            save_upset(available, gene_compounds, top_n, output_dir, args.lib)

        if args.profile:
            print("\nLoading compound library for physicochemical profiling...")
            compounds_df = pd.read_csv(COMPOUNDS_PATH)
            all_ids = set(compounds_df["id"])

            # Union of all hits across every requested top-N
            hit_ids = set()
            for top_n in top_values:
                for g in available:
                    hit_ids.update(gene_compounds[g][:top_n])
            hit_ids &= all_ids

            # Random sample of BG_SAMPLE_SIZE from the non-hit remainder
            non_hit_ids = list(all_ids - hit_ids)
            rng = np.random.default_rng(RANDOM_SEED)
            sample_size = min(BG_SAMPLE_SIZE, len(non_hit_ids))
            bg_ids = set(rng.choice(non_hit_ids, size=sample_size, replace=False))

            profiling_ids = hit_ids | bg_ids
            subset = compounds_df[compounds_df["id"].isin(profiling_ids)]
            smiles_dict = dict(zip(subset["id"], subset["smiles"]))

            print(f"  {len(hit_ids)} hit compounds + {len(bg_ids)} background samples")
            props_df = compute_properties(smiles_dict)

            profile_dir = os.path.join(ROOT, "output", "physicochemical_profiling")
            os.makedirs(profile_dir, exist_ok=True)
            for top_n in top_values:
                run_profiling(available, gene_compounds, top_n, props_df, args.lib, profile_dir)


if __name__ == "__main__":
    main()
