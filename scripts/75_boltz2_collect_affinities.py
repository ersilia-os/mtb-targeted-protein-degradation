#!/usr/bin/env python3
"""
Reads script 73's affinity results CSV (now written directly to local disk -- no longer requires
an rsync pull from nebula), joins in gene_name/site_type (from pocket_sequences.csv) and smiles
(from compounds.csv), and writes one enriched, long-format CSV for downstream consumers (e.g. the
molecule-auditing explorer). Safe to re-run anytime during the multi-day run -- always reflects
script 73's latest snapshot on disk, no append/destructive logic. Stays raw/lossless (one row per
pocket x compound result): does not merge pheS/pheT into pheST and does not aggregate to a
best-score-per-gene-x-site-type -- that's a downstream concern (see
filtered_hits_explorer/prepare_audit_input.py, which does it for docking scores). Also plots, per
pocket with data so far, an affinity_pred_value distribution (x-axis in IC50 uM, not raw log10) and
an affinity_pred_value vs affinity_probability_binary scatter annotated with Pearson/Spearman
correlation, as a quick model-behavior sanity check.

Usage:
    python 75_boltz2_collect_affinities.py [--out-subdir boltz_results]
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

DOCKING_DIR = os.path.join(ROOT, "output", "73_boltz2_docking")

SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")
POCKET_SEQUENCES_CSV = os.path.join(ROOT, "output", "71_boltz2_prepare_inputs", "pocket_sequences.csv")
COMPOUNDS_CSV = os.path.join(ROOT, "output", "71_boltz2_prepare_inputs", "compounds.csv")

OUTPUT_DIR = os.path.join(ROOT, "output", "75_boltz2_collect_affinities")
PLOTS_DIR = os.path.join(OUTPUT_DIR, "plots")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

N_COMPOUNDS = 1095

RESULT_COLUMNS = [
    "gene_name", "site_type", "pocket_name", "compound_id", "smiles",
    "affinity_pred_value", "affinity_probability_binary",
    "affinity_pred_value1", "affinity_probability_binary1",
    "affinity_pred_value2", "affinity_probability_binary2",
    "confidence_score", "ptm", "iptm", "ligand_iptm", "complex_plddt",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-subdir", default="boltz_results",
                         help="Output subdirectory script 73 wrote to, e.g. for a scoped/smoke-test "
                              "run (default: boltz_results)")
    return parser.parse_args()


def local_csv_path(out_subdir):
    path = os.path.join(DOCKING_DIR, f"{out_subdir}_affinity_results.csv")
    if not os.path.isfile(path):
        sys.exit(f"No affinity results CSV found at {path}; has script 73 written it yet?")
    return path


def load_raw_results(path):
    return pd.read_csv(path)


def load_pocket_sequences():
    return pd.read_csv(POCKET_SEQUENCES_CSV)


def load_compounds():
    return pd.read_csv(COMPOUNDS_CSV)


def sanity_check_pocket_lists(selected_pockets, pocket_sequences):
    """Script 71 now covers all 12 pockets, including the 7K98_pocket_6 dimer (2-chain pheS+pheT
    complex), so no diff between the two pocket lists is expected at all."""
    diff = set(selected_pockets["pocket_name"]) - set(pocket_sequences["pocket_name"])
    if diff:
        print(f"WARNING: unexpected pockets in selected_pockets.csv but not pocket_sequences.csv: {diff}")
    else:
        print("Pocket lists check out: all pockets in selected_pockets.csv are present in pocket_sequences.csv.")


def merge_results(raw_df, pocket_sequences, compounds):
    df = raw_df.merge(pocket_sequences[["pocket_name", "gene_name", "site_type"]], on="pocket_name", how="left")
    df = df.merge(compounds[["compound_id", "smiles"]], on="compound_id", how="left")

    unmatched_pockets = sorted(df.loc[df["gene_name"].isna(), "pocket_name"].unique())
    unmatched_compounds = sorted(df.loc[df["smiles"].isna(), "compound_id"].unique())

    df = (df[RESULT_COLUMNS]
          .sort_values(["gene_name", "site_type", "pocket_name", "compound_id"])
          .reset_index(drop=True))
    return df, unmatched_pockets, unmatched_compounds


def report_unmatched(unmatched_pockets, unmatched_compounds):
    """Never dropped, just flagged -- rows with no gene/site_type or smiles match are kept with
    those fields NaN."""
    if unmatched_pockets:
        print(f"WARNING: {len(unmatched_pockets)} pocket_name(s) with no gene_name/site_type match: "
              f"{unmatched_pockets}")
    if unmatched_compounds:
        print(f"WARNING: {len(unmatched_compounds)} compound_id(s) with no smiles match "
              f"(showing up to 5): {unmatched_compounds[:5]}")


def report_progress(df, pocket_sequences):
    counts = df.groupby("pocket_name").size().reindex(pocket_sequences["pocket_name"], fill_value=0)
    progress = pocket_sequences.set_index("pocket_name").assign(n_rows=counts)
    progress["affinities"] = progress["n_rows"].astype(str) + f"/{N_COMPOUNDS}"
    progress["pct"] = (100 * progress["n_rows"] / N_COMPOUNDS).round(1)
    progress = progress[["gene_name", "site_type", "n_rows", "affinities", "pct"]].sort_values("pct")

    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(progress.to_string())

    zero = progress[progress["n_rows"] == 0]
    if len(zero):
        print(f"\nPockets with 0 results so far: {zero.index.tolist()}")

    n_pockets_with_data = int((progress["n_rows"] > 0).sum())
    n_unique_compounds = df["compound_id"].nunique()
    print(f"\nTotal: {len(df):,} rows, {n_unique_compounds}/{N_COMPOUNDS} unique compounds seen, "
          f"{n_pockets_with_data}/{len(pocket_sequences)} pockets with any data")


def _um_ticks(vmin, vmax):
    """Integer log10(IC50 uM) tick positions spanning [vmin, vmax], labeled in IC50 uM (e.g.
    0.1 uM, 1 uM, 10 uM) instead of raw log10 values -- easier to read at a glance."""
    lo, hi = int(np.floor(vmin)), int(np.ceil(vmax))
    ticks = list(range(lo, hi + 1))
    labels = [f"{10 ** t:g} uM" for t in ticks]
    return ticks, labels


def plot_pocket_affinities(df):
    """One PNG per pocket with data: affinity_pred_value distribution (x-axis in IC50 uM) +
    affinity_pred_value vs affinity_probability_binary scatter annotated with Pearson/Spearman
    correlation. Format: slide | Style: ersilia -- change with stylia.set_format() /
    stylia.set_style()."""
    stylia = _safe_import_stylia()
    stylia.set_format("slide")
    stylia.set_style("ersilia")
    nc = stylia.NamedColors()

    n_written = 0
    correlations = []
    for pocket_name, group in df.groupby("pocket_name"):
        n_before = len(group)
        group = group.dropna(subset=["affinity_pred_value", "affinity_probability_binary"])
        n_dropped = n_before - len(group)
        if n_dropped:
            print(f"  WARNING: {pocket_name}: dropping {n_dropped} row(s) with missing affinity values from plot")
        if group.empty:
            continue

        fig, axs = stylia.create_figure(1, 2)

        ax = axs.next()
        ax.hist(group["affinity_pred_value"], color=nc.blue)
        ticks, tick_labels = _um_ticks(group["affinity_pred_value"].min(), group["affinity_pred_value"].max())
        ax.set_xticks(ticks)
        ax.set_xticklabels(tick_labels)
        stylia.label(ax, xlabel="Predicted IC50", ylabel="Count",
                     title=f"{pocket_name}: affinity distribution (n={len(group):,})", abc="A")

        r, _ = pearsonr(group["affinity_probability_binary"], group["affinity_pred_value"])
        rho, _ = spearmanr(group["affinity_probability_binary"], group["affinity_pred_value"])
        correlations.append({
            "pocket_name": pocket_name,
            "gene_name": group["gene_name"].iloc[0],
            "site_type": group["site_type"].iloc[0],
            "n": len(group),
            "pearson_r": r,
            "spearman_rho": rho,
        })

        ax = axs.next()
        ax.scatter(group["affinity_probability_binary"], group["affinity_pred_value"], color=nc.purple)
        ax.text(0.03, 0.97, f"Pearson r = {r:.2f}\nSpearman ρ = {rho:.2f}",
                transform=ax.transAxes, ha="left", va="top")
        stylia.label(ax, xlabel="P(binder)", ylabel="affinity_pred_value (log10 IC50, uM)",
                     title=f"{pocket_name}: affinity vs P(binder)", abc="B")

        out_path = os.path.join(PLOTS_DIR, f"{pocket_name}.png")
        stylia.save_figure(out_path)
        print(f"  {pocket_name}: {len(group):,} points -> {out_path}")
        n_written += 1

    print(f"\nWrote {n_written} affinity plot(s) to {PLOTS_DIR}")

    corr_df = pd.DataFrame(correlations).sort_values("pocket_name").reset_index(drop=True)
    corr_path = os.path.join(OUTPUT_DIR, "affinity_probability_correlations.csv")
    corr_df.to_csv(corr_path, index=False)
    print(f"\nSaved affinity_pred_value vs P(binder) correlations (Pearson/Spearman) per pocket to {corr_path}")
    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(corr_df.to_string(index=False))


def main():
    args = parse_args()
    local_csv = local_csv_path(args.out_subdir)

    raw_df = load_raw_results(local_csv)
    pocket_sequences = load_pocket_sequences()
    compounds = load_compounds()
    selected_pockets = pd.read_csv(SELECTED_POCKETS_CSV)

    sanity_check_pocket_lists(selected_pockets, pocket_sequences)

    df, unmatched_pockets, unmatched_compounds = merge_results(raw_df, pocket_sequences, compounds)
    report_unmatched(unmatched_pockets, unmatched_compounds)
    print()
    report_progress(df, pocket_sequences)

    out_path = os.path.join(OUTPUT_DIR, "affinity_results.csv")
    df.to_csv(out_path, index=False)
    print(f"\nSaved {len(df):,} rows x {len(df.columns)} columns to {out_path}")

    print()
    plot_pocket_affinities(df)


if __name__ == "__main__":
    main()
