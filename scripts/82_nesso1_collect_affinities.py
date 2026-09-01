#!/usr/bin/env python3
"""
Reads script 80's affinity results CSV (one row per gene_name x compound) and joins in smiles
(from compounds.csv), writing one enriched, long-format CSV. Revised for the per-gene design (see
script 78's docstring): no more pocket-level broadcast/dedup bookkeeping (explode-to-pockets,
shared_structure_with) since there's no pocket concept left at all -- Nesso-1's results are
already exactly at the granularity they were computed at (21 genes), nothing to fan back out to.

Stays raw/lossless (one row per gene x compound result) -- does not aggregate further, that's a
downstream concern. Also plots, per gene with data so far, an affinity_pred_value distribution
(x-axis in IC50 uM, not raw log10) and an affinity_pred_value vs affinity_probability_binary
scatter annotated with Pearson/Spearman correlation, as a quick model-behavior sanity check (same
convention as Boltz-2's script 75).

Usage:
    python 82_nesso1_collect_affinities.py [--out-subdir nesso_results]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from docking_utils import _safe_import_stylia  # noqa: E402

DOCKING_DIR = os.path.join(ROOT, "output", "80_nesso1_docking")

PROTEIN_SEQUENCES_CSV = os.path.join(ROOT, "output", "78_nesso1_prepare_inputs", "protein_sequences.csv")
COMPOUNDS_CSV = os.path.join(ROOT, "output", "71_boltz2_prepare_inputs", "compounds.csv")

OUTPUT_DIR = os.path.join(ROOT, "output", "82_nesso1_collect_affinities")
PLOTS_DIR = os.path.join(OUTPUT_DIR, "plots")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

N_COMPOUNDS = 1095

AFFINITY_FIELDS = ["affinity_pred_value", "affinity_pred_value1", "affinity_pred_value2",
                   "affinity_logits_binary", "affinity_probability_binary"]
ENTROPY_FIELDS = ["entropy_pp", "entropy_pl", "entropy_ll",
                   "entropy_crop_pp", "entropy_crop_pl", "entropy_crop_ll"]
RESULT_COLUMNS = (
    ["gene_name", "uniprot_ac", "compound_id", "smiles"]
    + AFFINITY_FIELDS + ENTROPY_FIELDS
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-subdir", default="nesso_results",
                         help="Output subdirectory script 80 wrote to, e.g. for a scoped/smoke-test "
                              "run (default: nesso_results)")
    return parser.parse_args()


def local_csv_path(out_subdir):
    path = os.path.join(DOCKING_DIR, f"{out_subdir}_affinity_results.csv")
    if not os.path.isfile(path):
        sys.exit(f"No affinity results CSV found at {path}; has script 80 written it yet?")
    return path


def merge_results(raw_df, proteins, compounds):
    df = raw_df.merge(proteins[["gene_name", "uniprot_ac"]], on="gene_name", how="left")
    df = df.merge(compounds[["compound_id", "smiles"]], on="compound_id", how="left")

    unmatched_genes = sorted(df.loc[df["uniprot_ac"].isna(), "gene_name"].unique())
    unmatched_compounds = sorted(df.loc[df["smiles"].isna(), "compound_id"].unique())

    df = (df[RESULT_COLUMNS]
          .sort_values(["gene_name", "compound_id"])
          .reset_index(drop=True))
    return df, unmatched_genes, unmatched_compounds


def report_unmatched(unmatched_genes, unmatched_compounds):
    """Never dropped, just flagged -- rows with no uniprot_ac or smiles match are kept with those
    fields NaN."""
    if unmatched_genes:
        print(f"WARNING: {len(unmatched_genes)} gene_name(s) with no uniprot_ac match: {unmatched_genes}")
    if unmatched_compounds:
        print(f"WARNING: {len(unmatched_compounds)} compound_id(s) with no smiles match "
              f"(showing up to 5): {unmatched_compounds[:5]}")


def report_progress(df, proteins):
    counts = df.groupby("gene_name").size().reindex(proteins["gene_name"], fill_value=0)
    progress = proteins.set_index("gene_name").assign(n_rows=counts)
    progress["affinities"] = progress["n_rows"].astype(str) + f"/{N_COMPOUNDS}"
    progress["pct"] = (100 * progress["n_rows"] / N_COMPOUNDS).round(1)
    progress = progress[["uniprot_ac", "n_rows", "affinities", "pct"]].sort_values("pct")

    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(progress.to_string())

    zero = progress[progress["n_rows"] == 0]
    if len(zero):
        print(f"\nGenes with 0 results so far: {zero.index.tolist()}")

    n_genes_with_data = int((progress["n_rows"] > 0).sum())
    n_unique_compounds = df["compound_id"].nunique()
    print(f"\nTotal: {len(df):,} rows, {n_unique_compounds}/{N_COMPOUNDS} unique compounds seen, "
          f"{n_genes_with_data}/{len(proteins)} genes with any data")


def _um_ticks(vmin, vmax):
    """Integer log10(IC50 uM) tick positions spanning [vmin, vmax], labeled in IC50 uM (e.g.
    0.1 uM, 1 uM, 10 uM) instead of raw log10 values -- easier to read at a glance."""
    lo, hi = int(np.floor(vmin)), int(np.ceil(vmax))
    ticks = list(range(lo, hi + 1))
    labels = [f"{10 ** t:g} uM" for t in ticks]
    return ticks, labels


def plot_gene_affinities(df):
    """One PNG per gene with data: affinity_pred_value distribution (x-axis in IC50 uM) +
    affinity_pred_value vs affinity_probability_binary scatter annotated with Pearson/Spearman
    correlation. Format: slide | Style: ersilia -- change with stylia.set_format() /
    stylia.set_style()."""
    stylia = _safe_import_stylia()
    stylia.set_format("slide")
    stylia.set_style("ersilia")
    nc = stylia.NamedColors()

    n_written = 0
    correlations = []
    for gene_name, group in df.groupby("gene_name"):
        n_before = len(group)
        group = group.dropna(subset=["affinity_pred_value", "affinity_probability_binary"])
        n_dropped = n_before - len(group)
        if n_dropped:
            print(f"  WARNING: {gene_name}: dropping {n_dropped} row(s) with missing affinity values from plot")
        if group.empty:
            continue

        fig, axs = stylia.create_figure(1, 2)

        ax = axs.next()
        ax.hist(group["affinity_pred_value"], color=nc.blue)
        ticks, tick_labels = _um_ticks(group["affinity_pred_value"].min(), group["affinity_pred_value"].max())
        ax.set_xticks(ticks)
        ax.set_xticklabels(tick_labels)
        stylia.label(ax, xlabel="Predicted IC50", ylabel="Count",
                     title=f"{gene_name}: affinity distribution (n={len(group):,})", abc="A")

        r, _ = pearsonr(group["affinity_probability_binary"], group["affinity_pred_value"])
        rho, _ = spearmanr(group["affinity_probability_binary"], group["affinity_pred_value"])
        correlations.append({
            "gene_name": gene_name,
            "n": len(group),
            "pearson_r": r,
            "spearman_rho": rho,
        })

        ax = axs.next()
        ax.scatter(group["affinity_probability_binary"], group["affinity_pred_value"], color=nc.purple)
        ax.text(0.03, 0.97, f"Pearson r = {r:.2f}\nSpearman ρ = {rho:.2f}",
                transform=ax.transAxes, ha="left", va="top")
        stylia.label(ax, xlabel="P(binder)", ylabel="affinity_pred_value (log10 IC50, uM)",
                     title=f"{gene_name}: affinity vs P(binder)", abc="B")

        out_path = os.path.join(PLOTS_DIR, f"{gene_name}.png")
        stylia.save_figure(out_path)
        print(f"  {gene_name}: {len(group):,} points -> {out_path}")
        n_written += 1

    print(f"\nWrote {n_written} affinity plot(s) to {PLOTS_DIR}")

    corr_df = pd.DataFrame(correlations).sort_values("gene_name").reset_index(drop=True)
    corr_path = os.path.join(OUTPUT_DIR, "affinity_probability_correlations.csv")
    corr_df.to_csv(corr_path, index=False)
    print(f"\nSaved affinity_pred_value vs P(binder) correlations (Pearson/Spearman) per gene to {corr_path}")
    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(corr_df.to_string(index=False))


def main():
    args = parse_args()
    local_csv = local_csv_path(args.out_subdir)

    raw_df = pd.read_csv(local_csv)
    proteins = pd.read_csv(PROTEIN_SEQUENCES_CSV)
    compounds = pd.read_csv(COMPOUNDS_CSV)

    df, unmatched_genes, unmatched_compounds = merge_results(raw_df, proteins, compounds)
    report_unmatched(unmatched_genes, unmatched_compounds)
    print()
    report_progress(df, proteins)

    out_path = os.path.join(OUTPUT_DIR, "affinity_results.csv")
    df.to_csv(out_path, index=False)
    print(f"\nSaved {len(df):,} rows x {len(df.columns)} columns to {out_path}")

    print()
    plot_gene_affinities(df)


if __name__ == "__main__":
    main()
