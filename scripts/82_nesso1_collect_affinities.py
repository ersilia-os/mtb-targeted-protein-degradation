#!/usr/bin/env python3
"""
Reads script 80's affinity results CSV (one row per structure_id x compound), broadcasts each
unique structure's result back out to every selected_pockets.csv pocket_name that maps to it
(script 78 deduplicated pockets sharing an identical structure/sequence -- e.g. lysS's CAT and
NON-CAT pockets both sit on the same AlphaFold3 model -- since Nesso-1 has no pocket-conditioning
to tell them apart), and writes one enriched, long-format CSV with the same
gene_name/site_type/pocket_name/compound_id/smiles shape as Boltz-2's
output/75_boltz2_collect_affinities/affinity_results.csv. A `shared_structure_with` column names
any sibling pocket_name(s) that received the identical underlying Nesso-1 result, so the
deduplication is visible in the data, never silent.

Stays raw/lossless (one row per pocket x compound result): does not merge pheS/pheT into pheST
and does not aggregate to a best-score-per-gene-x-site-type -- that's a downstream concern.
Also plots, per structure with data so far, an affinity_pred_value distribution (x-axis in IC50
uM, not raw log10) and an affinity_pred_value vs affinity_probability_binary scatter annotated
with Pearson/Spearman correlation, as a quick model-behavior sanity check (same convention as
script 75, one plot per structure rather than per pocket_name so lysS's two duplicate pockets
don't produce two identical plots).

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

SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")
STRUCTURE_SEQUENCES_CSV = os.path.join(ROOT, "output", "78_nesso1_prepare_inputs", "structure_sequences.csv")
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
    ["gene_name", "site_type", "pocket_name", "structure_id", "compound_id", "smiles"]
    + AFFINITY_FIELDS + ENTROPY_FIELDS
    + ["used_trimmed_dimer_fallback", "shared_structure_with"]
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


def explode_to_pockets(structures):
    """One row per (pocket_name, structure_id), where structures['pocket_names'] is the
    space-separated list of every original selected_pockets.csv pocket_name mapped to this
    unique structure (script 78). `shared_structure_with` names the other pocket_name(s) sharing
    that same structure, empty string if none."""
    rows = []
    for _, srow in structures.iterrows():
        names = srow["pocket_names"].split()
        for pocket_name in names:
            siblings = [n for n in names if n != pocket_name]
            rows.append({
                "pocket_name": pocket_name,
                "structure_id": srow["structure_id"],
                "shared_structure_with": " ".join(siblings),
            })
    return pd.DataFrame(rows)


def build_structure_labels(structures, selected_pockets):
    """{structure_id: display label} for plot titles, e.g. 'lysS (CAT/NON-CAT)
    [alphafold3_..._pocket_1, alphafold3_..._pocket_2]'."""
    sp = selected_pockets.set_index("pocket_name")
    labels = {}
    for _, srow in structures.iterrows():
        names = srow["pocket_names"].split()
        site_types = sorted({sp.loc[n, "site_type"] for n in names if n in sp.index})
        labels[srow["structure_id"]] = f"{srow['gene_name']} ({'/'.join(site_types)}) [{', '.join(names)}]"
    return labels


def merge_results(raw_df, structures, selected_pockets, compounds):
    """Broadcasts each real result in raw_df out to every pocket_name sharing its structure_id.
    Starts FROM raw_df (not from the full pocket list), so -- same raw/lossless convention as
    Boltz-2's script 75 -- a pocket with zero results so far produces zero rows, not a
    placeholder NaN row."""
    pocket_map = explode_to_pockets(structures)
    df = raw_df.merge(pocket_map, on="structure_id", how="left")
    df = df.merge(selected_pockets[["pocket_name", "gene_name", "site_type"]], on="pocket_name", how="left")
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


def report_progress(df, structures):
    pocket_names = [n for names in structures["pocket_names"].str.split() for n in names]
    counts = df.groupby("pocket_name").size().reindex(pocket_names, fill_value=0)
    progress = pd.DataFrame({"pocket_name": pocket_names, "n_rows": counts.values})
    progress["affinities"] = progress["n_rows"].astype(str) + f"/{N_COMPOUNDS}"
    progress["pct"] = (100 * progress["n_rows"] / N_COMPOUNDS).round(1)
    progress = progress.sort_values("pct")

    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(progress.to_string(index=False))

    zero = progress[progress["n_rows"] == 0]
    if len(zero):
        print(f"\nPockets with 0 results so far: {zero['pocket_name'].tolist()}")

    n_pockets_with_data = int((progress["n_rows"] > 0).sum())
    n_unique_compounds = df["compound_id"].nunique()
    print(f"\nTotal: {len(df):,} rows, {n_unique_compounds}/{N_COMPOUNDS} unique compounds seen, "
          f"{n_pockets_with_data}/{len(pocket_names)} pockets with any data")


def _um_ticks(vmin, vmax):
    """Integer log10(IC50 uM) tick positions spanning [vmin, vmax], labeled in IC50 uM (e.g.
    0.1 uM, 1 uM, 10 uM) instead of raw log10 values -- easier to read at a glance."""
    lo, hi = int(np.floor(vmin)), int(np.ceil(vmax))
    ticks = list(range(lo, hi + 1))
    labels = [f"{10 ** t:g} uM" for t in ticks]
    return ticks, labels


def plot_structure_affinities(df, structure_labels):
    """One PNG per structure_id with data (not per pocket_name -- lysS's two duplicate pockets
    would otherwise produce two identical plots): affinity_pred_value distribution (x-axis in
    IC50 uM) + affinity_pred_value vs affinity_probability_binary scatter annotated with
    Pearson/Spearman correlation. Format: slide | Style: ersilia -- change with
    stylia.set_format() / stylia.set_style()."""
    stylia = _safe_import_stylia()
    stylia.set_format("slide")
    stylia.set_style("ersilia")
    nc = stylia.NamedColors()

    n_written = 0
    correlations = []
    for structure_id, group in df.groupby("structure_id"):
        n_before = len(group)
        group = group.dropna(subset=["affinity_pred_value", "affinity_probability_binary"])
        n_dropped = n_before - len(group)
        if n_dropped:
            print(f"  WARNING: {structure_id}: dropping {n_dropped} row(s) with missing affinity values from plot")
        if group.empty:
            continue

        fig, axs = stylia.create_figure(1, 2)

        ax = axs.next()
        ax.hist(group["affinity_pred_value"], color=nc.blue)
        ticks, tick_labels = _um_ticks(group["affinity_pred_value"].min(), group["affinity_pred_value"].max())
        ax.set_xticks(ticks)
        ax.set_xticklabels(tick_labels)
        stylia.label(ax, xlabel="Predicted IC50", ylabel="Count",
                     title=f"{structure_labels[structure_id]}: affinity distribution (n={len(group):,})", abc="A")

        r, _ = pearsonr(group["affinity_probability_binary"], group["affinity_pred_value"])
        rho, _ = spearmanr(group["affinity_probability_binary"], group["affinity_pred_value"])
        correlations.append({
            "structure_id": structure_id,
            "pocket_names": group["pocket_name"].drop_duplicates().tolist(),
            "n": len(group),
            "pearson_r": r,
            "spearman_rho": rho,
        })

        ax = axs.next()
        ax.scatter(group["affinity_probability_binary"], group["affinity_pred_value"], color=nc.purple)
        ax.text(0.03, 0.97, f"Pearson r = {r:.2f}\nSpearman ρ = {rho:.2f}",
                transform=ax.transAxes, ha="left", va="top")
        stylia.label(ax, xlabel="P(binder)", ylabel="affinity_pred_value (log10 IC50, uM)",
                     title=f"{structure_labels[structure_id]}: affinity vs P(binder)", abc="B")

        out_path = os.path.join(PLOTS_DIR, f"{structure_id}.png")
        stylia.save_figure(out_path)
        print(f"  {structure_id}: {len(group):,} points -> {out_path}")
        n_written += 1

    print(f"\nWrote {n_written} affinity plot(s) to {PLOTS_DIR}")

    corr_df = pd.DataFrame(correlations).sort_values("structure_id").reset_index(drop=True)
    corr_path = os.path.join(OUTPUT_DIR, "affinity_probability_correlations.csv")
    corr_df.to_csv(corr_path, index=False)
    print(f"\nSaved affinity_pred_value vs P(binder) correlations (Pearson/Spearman) per structure to {corr_path}")
    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(corr_df.to_string(index=False))


def main():
    args = parse_args()
    local_csv = local_csv_path(args.out_subdir)

    raw_df = pd.read_csv(local_csv)
    structures = pd.read_csv(STRUCTURE_SEQUENCES_CSV)
    compounds = pd.read_csv(COMPOUNDS_CSV)
    selected_pockets = pd.read_csv(SELECTED_POCKETS_CSV)

    df, unmatched_pockets, unmatched_compounds = merge_results(raw_df, structures, selected_pockets, compounds)
    report_unmatched(unmatched_pockets, unmatched_compounds)
    print()
    report_progress(df, structures)

    out_path = os.path.join(OUTPUT_DIR, "affinity_results.csv")
    df.to_csv(out_path, index=False)
    print(f"\nSaved {len(df):,} rows x {len(df.columns)} columns to {out_path}")

    n_dedup_rows = df["shared_structure_with"].astype(bool).sum()
    if n_dedup_rows:
        print(f"({n_dedup_rows:,} row(s) share their underlying Nesso-1 result with a sibling "
              f"pocket_name on the same structure -- see 'shared_structure_with'.)")

    structure_labels = build_structure_labels(structures, selected_pockets)
    print()
    plot_structure_affinities(df, structure_labels)


if __name__ == "__main__":
    main()
