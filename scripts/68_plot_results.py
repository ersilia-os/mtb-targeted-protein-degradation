#!/usr/bin/env python3
"""
Plots/reports script 67's Ersilia characterization outputs for the aggregated hit set, plus docking
score boxplots for the 12 curated pockets: script 65/66's aggregated cross-docking set against Hit
Locator / REAL round 1 / REAL round 2 references, mirroring script 61's design. A missing input
(e.g. eos42ez.csv, not yet run via script 67) prints a warning and is skipped, not fatal.

Usage:
    python 68_plot_results.py
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from docking_utils import LIBRARIES, load_real_positive_scores, load_scores, plot_score_boxplots_multi

ERSILIA_DIR = os.path.join(ROOT, "output", "67_ersilia_characterization")
SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")
AGGREGATED_DOCKING_DIR = os.path.join(ROOT, "output", "65_aggregated_docking", "docking_results")
MERGED_DOCKING_SCORES_CSV = os.path.join(ROOT, "output", "66_merge_docking_scores", "merged_docking_scores.csv")
MULTIMER_DOCKING_DIR = os.path.join(ROOT, "output", "50_unidock_docking_multimers")
DIMER_POCKET = "7K98_pocket_6"
OUTPUT_DIR = os.path.join(ROOT, "output", "68_plot_results")
os.makedirs(OUTPUT_DIR, exist_ok=True)

TOP_N = 10
DOCKING_SCORE_THRESHOLDS = [-12, -11, -10]
RANK_THRESHOLDS = [100, 500, 1000]
FLAG_COLS = ["has_pains", "has_brenk", "is_sim_known_ab", "nitrofuran_motif",
             "fluoroquinolone_motif", "carbepenem_motif", "betalactam_motif"]
HIST_MODELS = [
    ("eos12x7.csv", ["sps_score", "nsps_score"], "eos12x7_distributions.png", "eos12x7 (molecular complexity)"),
    ("eos5jv3.csv", ["mycomembrane_permeation"], "eos5jv3_distribution.png", "eos5jv3 (permeability)"),
    ("eos42ez.csv", ["cytotoxicity_hepg2", "cytotoxicity_hskmc", "cytotoxicity_imr90"], "eos42ez_distributions.png", "eos42ez (cytotoxicity)"),
]


def _safe_import_stylia():
    import matplotlib as mpl
    os.makedirs(mpl.get_cachedir(), exist_ok=True)
    import stylia
    return stylia


def plot_histograms(df, columns, out_path):
    stylia = _safe_import_stylia()
    stylia.set_format("slide")
    stylia.set_style("ersilia")
    palette_src = stylia.NamedColors()
    palette = [palette_src.blue, palette_src.purple, palette_src.mint, palette_src.orange]

    fig, axs = stylia.create_figure(1, len(columns), width=0.5 * len(columns), height=0.5)
    for col, color in zip(columns, palette):
        ax = axs.next()
        ax.hist(df[col].dropna(), bins=30, color=color)
        stylia.label(ax, xlabel=col, ylabel="Count")
    stylia.save_figure(out_path)


def load_report(directory, pocket):
    path = os.path.join(directory, pocket, "report.csv")
    return load_scores(path) if os.path.isfile(path) else pd.Series(dtype=float)


def load_real2_report(pocket):
    """REAL round 2's report for `pocket` - the dimer pocket's equivalent (docked against the same
    ~99k round-2 set, see scripts/50_unidock_docking_multimers.py) lives under MULTIMER_DOCKING_DIR
    instead of LIBRARIES["REAL"], same distinction script 54 already makes."""
    directory = MULTIMER_DOCKING_DIR if pocket == DIMER_POCKET else LIBRARIES["REAL"]
    return load_report(directory, pocket)


def build_present(loader, pockets):
    """{pocket: scores} for pockets where loader(pocket) returns non-empty scores."""
    return {p: s for p in pockets if len(s := loader(p))}


def load_pocket_ranks():
    """{pocket: pd.Series(rank, indexed by compound_id)} from script 66's merged_docking_scores.csv
    - rank vs. the Enamine REAL round-2 (~99,105) background library, not vs. our own 2,923-compound
    aggregated set, so values can go well beyond 2,923 (up to ~100k)."""
    df = pd.read_csv(MERGED_DOCKING_SCORES_CSV, index_col="compound_id")
    rank_cols = [c for c in df.columns if c.endswith("_rank")]
    return {c.removesuffix("_rank"): df[c].dropna() for c in rank_cols}


def load_protein_pockets():
    """{protein: [pocket_name, ...]} for all 12 curated pockets, pheS+pheT merged into "pheST" -
    same merge convention as scripts/54_NONCAT_promiscuous.py's load_noncat_targets(), extended
    here to include CAT pockets too (not NON-CAT only)."""
    df = pd.read_csv(SELECTED_POCKETS_CSV)
    df["gene_name"] = df["gene_name"].replace({"pheS": "pheST", "pheT": "pheST"})
    return df.groupby("gene_name")["pocket_name"].apply(list).to_dict()


def summarize_by_n_targets(hit_sets):
    """'1 target: X, 2 targets: Y, ...' - how many compounds hit exactly k of the proteins in
    hit_sets, for k from 1 to len(hit_sets)."""
    from collections import Counter
    all_ids = set().union(*hit_sets.values())
    counts = Counter(sum(cid in s for s in hit_sets.values()) for cid in all_ids)
    n = len(hit_sets)
    return ", ".join(f"{k} target{'s' if k != 1 else ''}: {counts.get(k, 0)}" for k in range(1, n + 1))


def annotate_degree_groups(ax_bars, upset):
    """Draws a horizontal bracket line (same height for every group, based on the tallest bar
    overall) + total-count label above each contiguous group of bars sharing the same degree
    (number of proteins in that intersection) - sort_by="degree" in save_protein_upset guarantees
    same-degree bars are contiguous, so a single pass suffices."""
    sizes = upset.intersections.to_numpy()
    degrees = [sum(idx) for idx in upset.intersections.index]
    bars = ax_bars.patches

    y_line = max(b.get_height() for b in bars) * 1.1
    label_y = y_line * 1.05

    i = 0
    while i < len(bars):
        j = i
        while j < len(bars) and degrees[j] == degrees[i]:
            j += 1
        x_start = bars[i].get_x()
        x_end = bars[j - 1].get_x() + bars[j - 1].get_width()
        degree = degrees[i]
        total = int(sizes[i:j].sum())
        ax_bars.plot([x_start, x_end], [y_line, y_line], color="black", lw=1, clip_on=False)
        ax_bars.text((x_start + x_end) / 2, label_y, f"{degree} target{'s' if degree != 1 else ''}; n={total}",
                     ha="center", va="bottom", fontsize=8, clip_on=False)
        i = j
    ax_bars.set_ylim(top=max(ax_bars.get_ylim()[1], label_y * 1.15))


def save_protein_upset(hit_sets, title, out_path):
    """Plain matplotlib + upsetplot, no stylia - matches docking_utils.save_upset and script 54's
    save_target_upset, since upsetplot.plot() doesn't compose with stylia's figure API."""
    import warnings
    import matplotlib.pyplot as plt
    from upsetplot import UpSet, from_contents
    warnings.filterwarnings("ignore", category=FutureWarning, module="upsetplot")
    warnings.filterwarnings("ignore", category=UserWarning, message=".*tight_layout.*")

    upset = UpSet(from_contents(hit_sets), sort_by="degree")
    axes = upset.plot()
    annotate_degree_groups(axes["intersections"], upset)

    plt.suptitle(f"{title}\n{summarize_by_n_targets(hit_sets)}")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close("all")


def plot_ersilia_distributions():
    for i, (filename, columns, out_name, label) in enumerate(HIST_MODELS):
        if i:
            print()
        print(f"--- {label} ---")
        path = os.path.join(ERSILIA_DIR, filename)
        if not os.path.isfile(path):
            print(f"Warning: {path} not found - skipping.")
            continue
        out_path = os.path.join(OUTPUT_DIR, out_name)
        plot_histograms(pd.read_csv(path), columns, out_path)
        print(f"Saved: {out_path}")


def report_eos2xeq_flags():
    print("--- eos2xeq (motif/structural-alert flags) ---")
    path = os.path.join(ERSILIA_DIR, "eos2xeq.csv")
    if not os.path.isfile(path):
        print(f"Warning: {path} not found - skipping.")
        return
    df = pd.read_csv(path)
    n = len(df)
    print(f"eos2xeq motif/structural-alert flags (n={n:,}):")
    for col in FLAG_COLS:
        count = int(df[col].sum())
        print(f"  {col:<22}: {count:>5,} ({100 * count / n:5.1f}%)")


def plot_docking_score_boxplots(pockets, agg_scores):
    print("--- Docking score boxplots (12 curated pockets) ---")
    score_sources = [
        ("Hit Locator", "mint", build_present(lambda p: load_report(LIBRARIES["DL"], p), pockets), True),
        ("REAL 1 - preselected", "orange", build_present(load_real_positive_scores, pockets), False),
        ("REAL 2 - preselected", "yellow", build_present(load_real2_report, pockets), True),
        ("All selected", "purple", agg_scores, True),
        ("Top-10 selected", "pink", {p: s.sort_values().head(TOP_N) for p, s in agg_scores.items()}, False),
    ]
    print("Data availability:")
    for label, _, scores, _ in score_sources:
        print(f"  {label}: {len(scores)}/{len(pockets)} pockets")

    out_path = os.path.join(OUTPUT_DIR, "aggregated_docking_score_boxplots.png")
    plot_score_boxplots_multi(pockets, score_sources, out_path, xtick_rotation=30)
    print(f"Saved: {out_path}")


def plot_protein_upsets(agg_scores):
    protein_pockets = load_protein_pockets()

    print("--- Protein-level docking score UpSet plots ---")
    for threshold in DOCKING_SCORE_THRESHOLDS:
        hit_sets = {
            protein: set().union(*(
                set(agg_scores[p][agg_scores[p] <= threshold].index)
                for p in pockets_for_protein if p in agg_scores
            ))
            for protein, pockets_for_protein in protein_pockets.items()
        }
        print(f"  score <= {threshold}: " + ", ".join(f"{p}={len(s)}" for p, s in hit_sets.items()))
        out_path = os.path.join(OUTPUT_DIR, f"upset_score_{abs(threshold)}.png")
        save_protein_upset(hit_sets, f"Docking score ≤ {threshold}", out_path)
        print(f"  Saved: {out_path}")

    print("--- Protein-level docking rank UpSet plots ---")
    pocket_ranks = load_pocket_ranks()
    for rank in RANK_THRESHOLDS:
        hit_sets = {
            protein: set().union(*(
                set(pocket_ranks[p][pocket_ranks[p] <= rank].index)
                for p in pockets_for_protein if p in pocket_ranks
            ))
            for protein, pockets_for_protein in protein_pockets.items()
        }
        print(f"  rank <= {rank}: " + ", ".join(f"{p}={len(s)}" for p, s in hit_sets.items()))
        out_path = os.path.join(OUTPUT_DIR, f"upset_rank_{rank}.png")
        save_protein_upset(hit_sets, f"Rank ≤ {rank} vs. REAL round-2 (~99,105)", out_path)
        print(f"  Saved: {out_path}")


def main():
    plot_ersilia_distributions()
    print()
    report_eos2xeq_flags()
    print()
    pockets = pd.read_csv(SELECTED_POCKETS_CSV)["pocket_name"].tolist()
    agg_scores = build_present(lambda p: load_report(AGGREGATED_DOCKING_DIR, p), pockets)
    plot_docking_score_boxplots(pockets, agg_scores)
    print()
    plot_protein_upsets(agg_scores)


if __name__ == "__main__":
    main()
