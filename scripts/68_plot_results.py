#!/usr/bin/env python3
"""
Plots/reports script 67's Ersilia characterization outputs, plus docking score boxplots and
protein-level UpSet plots for the 12 curated pockets. A missing input is skipped with a warning,
not fatal.

Usage:
    python 68_plot_results.py [--annotate] [--only-upset]
"""
import argparse
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
DOCKING_SCORE_THRESHOLDS = [-12, -11, -10, -9]
RANK_THRESHOLDS = [100, 500, 1000]
FLAG_COLS = ["has_pains", "has_brenk", "is_sim_known_ab", "nitrofuran_motif",
             "fluoroquinolone_motif", "carbepenem_motif", "betalactam_motif"]
HIST_MODELS = [
    ("eos12x7.csv", ["sps_score", "nsps_score"], "eos12x7_distributions.png", "eos12x7 (molecular complexity)"),
    ("eos5jv3.csv", ["mycomembrane_permeation"], "eos5jv3_distribution.png", "eos5jv3 (permeability)"),
    ("eos42ez.csv", ["cytotoxicity_hepg2", "cytotoxicity_hskmc", "cytotoxicity_imr90"], "eos42ez_distributions.png", "eos42ez (cytotoxicity)"),
]
PERMEABILITY_MODELS = [
    ("eos3ujl.csv", "permeability_probability", "eos3ujl"),
    ("eos8d8a.csv", "mycpermcheck", "eos8d8a"),
    ("eos1lb5.csv", "perm_proba_lepori_mtb", "eos1lb5_lepori_mtb"),
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


def load_report(directory, pocket, filename="report.csv"):
    path = os.path.join(directory, pocket, filename)
    return load_scores(path) if os.path.isfile(path) else pd.Series(dtype=float)


def load_real2_report(pocket):
    """REAL 10B's report for `pocket` - the dimer pocket lives under MULTIMER_DOCKING_DIR instead."""
    directory = MULTIMER_DOCKING_DIR if pocket == DIMER_POCKET else LIBRARIES["REAL"]
    return load_report(directory, pocket)


def build_present(loader, pockets):
    """{pocket: scores} for pockets where loader(pocket) returns non-empty scores."""
    return {p: s for p in pockets if len(s := loader(p))}


def load_pocket_ranks():
    """{pocket: pd.Series(rank, indexed by compound_id)} - rank vs. the ~100k-compound REAL
    10B background library, not vs. our own aggregated set."""
    df = pd.read_csv(MERGED_DOCKING_SCORES_CSV, index_col="compound_id")
    rank_cols = [c for c in df.columns if c.endswith("_rank")]
    return {c.removesuffix("_rank"): df[c].dropna() for c in rank_cols}


def load_protein_pockets(site_type=None):
    """{protein: [pocket_name, ...]}, pheS+pheT merged into "pheST". Optionally filter to one
    site_type ("CAT" or "NON-CAT")."""
    df = pd.read_csv(SELECTED_POCKETS_CSV)
    if site_type is not None:
        df = df[df["site_type"] == site_type]
    df["gene_name"] = df["gene_name"].replace({"pheS": "pheST", "pheT": "pheST"})
    return df.groupby("gene_name")["pocket_name"].apply(list).to_dict()


def summarize_by_n_targets(hit_sets):
    """'1 target: X, 2 targets: Y, ...' plus a single- vs multi-target proportion line."""
    from collections import Counter
    all_ids = set().union(*hit_sets.values())
    counts = Counter(sum(cid in s for s in hit_sets.values()) for cid in all_ids)
    n = len(hit_sets)
    breakdown = ", ".join(f"{k} target{'s' if k != 1 else ''}: {counts.get(k, 0)}" for k in range(1, n + 1))
    total = sum(counts.values())
    if not total:
        return breakdown
    single_pct = 100 * counts.get(1, 0) / total
    return f"{breakdown}\nSingle-target: {single_pct:.1f}%, Multi-target: {100 - single_pct:.1f}%"


def tighten_totals_gap(axes, fig, buffer=0.025):
    """upsetplot reserves a wide blank column between the totals-bar panel and the matrix/
    category-label panel (sized to fit the label text with a generous margin). Shrinks it to
    just fit the rendered category-label text (+ `buffer`, figure-fraction) by shifting the
    matrix + intersections panels left."""
    totals_ax = axes.get("totals")
    matrix_ax = axes.get("matrix")
    if totals_ax is None or matrix_ax is None:
        return
    labels = matrix_ax.get_yticklabels()
    if not labels:
        return
    renderer = fig.canvas.get_renderer()
    fig_width_px = fig.get_window_extent(renderer=renderer).width
    text_width_frac = max(t.get_window_extent(renderer=renderer).width for t in labels) / fig_width_px

    new_x0 = totals_ax.get_position().x1 + text_width_frac + buffer
    shift = matrix_ax.get_position().x0 - new_x0
    if shift <= 0:
        return
    for name in ("matrix", "intersections"):
        ax = axes.get(name)
        if ax is None:
            continue
        pos = ax.get_position()
        ax.set_position([pos.x0 - shift, pos.y0, pos.width + shift, pos.height])


def style_upset_by_degree(upset):
    """Colors both the intersection bars and matrix dots by degree (# proteins in the
    intersection) - must be called before upset.plot(). Returns {degree: color} for the legend."""
    stylia = _safe_import_stylia()
    ac = stylia.ArticleColors()
    palette = [ac.crimson, ac.tangerine, ac.amber, ac.lime, ac.turquoise, ac.cobalt, ac.periwinkle, ac.orchid]

    degrees = sorted(set(sum(idx) for idx in upset.intersections.index))
    degree_colors = {d: palette[(d - 1) % len(palette)] for d in degrees}
    for d, color in degree_colors.items():
        upset.style_subsets(min_degree=d, max_degree=d, facecolor=color)
    return degree_colors


def add_degree_legend(ax_bars, degree_colors):
    """Legend mapping each degree color to '{d} target(s)', placed to the right of the bars."""
    import matplotlib.patches as mpatches
    handles = [
        mpatches.Patch(color=color, label=f"{d} target{'s' if d != 1 else ''}")
        for d, color in degree_colors.items()
    ]
    ax_bars.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1), borderaxespad=0, fontsize=8)


def annotate_degree_groups(ax_bars, upset):
    """Draws a bracket + total-count label above each contiguous group of same-degree bars."""
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
        ax_bars.text((x_start + x_end) / 2, label_y, f"{degree} target{'s' if degree != 1 else ''}\nn={total}",
                     ha="center", va="bottom", fontsize=8, clip_on=False)
        i = j
    ax_bars.set_ylim(top=max(ax_bars.get_ylim()[1], label_y * 1.15))


def save_protein_upset(hit_sets, title, out_path, annotate=False):
    """upsetplot manages its own figure/axes, so this can't go through stylia's create_figure()/
    label() - but stylia's slide/article rcParams (fonts, sizes, black text) are applied first so
    the plot still matches the rest of the report."""
    import warnings
    import matplotlib.pyplot as plt
    from upsetplot import UpSet, from_contents
    warnings.filterwarnings("ignore", category=FutureWarning, module="upsetplot")
    warnings.filterwarnings("ignore", category=UserWarning, message=".*tight_layout.*")

    stylia = _safe_import_stylia()
    stylia.set_format("slide")
    stylia.set_style("article")

    upset = UpSet(from_contents(hit_sets), sort_by="degree")
    degree_colors = style_upset_by_degree(upset)
    axes = upset.plot()
    fig = plt.gcf()
    fig.canvas.draw()  # force a layout pass so text extents below are accurate, not placeholders
    tighten_totals_gap(axes, fig)
    add_degree_legend(axes["intersections"], degree_colors)
    if annotate:
        annotate_degree_groups(axes["intersections"], upset)
        title = f"{title}\n{summarize_by_n_targets(hit_sets)}"

    plt.suptitle(title)
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


def load_permeability_scores():
    """Merge eos3ujl, eos8d8a and eos1lb5's perm_proba_lepori_mtb column on `key` into one
    DataFrame (key, eos3ujl, eos8d8a, eos1lb5_lepori_mtb), or None if any input is missing."""
    merged = None
    for filename, column, label in PERMEABILITY_MODELS:
        path = os.path.join(ERSILIA_DIR, filename)
        if not os.path.isfile(path):
            print(f"Warning: {path} not found - skipping permeability comparison.")
            return None
        df = pd.read_csv(path, usecols=["key", column]).rename(columns={column: label})
        merged = df if merged is None else merged.merge(df, on="key", how="inner")
    return merged


def plot_permeability_comparison():
    print("--- Permeability model comparison (eos3ujl, eos8d8a, eos1lb5_lepori_mtb) ---")
    merged = load_permeability_scores()
    if merged is None:
        return

    stylia = _safe_import_stylia()
    stylia.set_format("slide")
    stylia.set_style("ersilia")
    nc = stylia.NamedColors()

    labels = [label for _, _, label in PERMEABILITY_MODELS]
    fig, axs = stylia.create_figure(2, 2)

    ax = axs.next()
    for label, color in zip(labels, [nc.blue, nc.purple, nc.mint]):
        ax.hist(merged[label], bins=30, color=color, alpha=0.5, label=label)
    ax.legend()
    stylia.label(ax, xlabel="Permeability probability", ylabel="Count", title=f"Distributions (n={len(merged):,})")

    for x_label, y_label in [(labels[0], labels[1]), (labels[0], labels[2]), (labels[1], labels[2])]:
        ax = axs.next()
        r = merged[x_label].corr(merged[y_label], method="spearman")
        ax.scatter(merged[x_label], merged[y_label], color=nc.blue, s=3, alpha=0.3)
        stylia.label(ax, xlabel=x_label, ylabel=y_label, title=f"n={len(merged):,}, r={r:.2f}")

    out_path = os.path.join(OUTPUT_DIR, "permeability_model_comparison.png")
    stylia.save_figure(out_path)
    print(f"Saved: {out_path}")


def report_nsps():
    print("--- eos12x7 (NSPS drug-likeness range) ---")
    path = os.path.join(ERSILIA_DIR, "eos12x7.csv")
    if not os.path.isfile(path):
        print(f"Warning: {path} not found - skipping.")
        return
    df = pd.read_csv(path)
    n = len(df)
    count = int(((df["nsps_score"] >= 10) & (df["nsps_score"] <= 40)).sum())
    print(f"NSPS in [10, 40]: {count:,} / {n:,} ({100 * count / n:.1f}%)")


def report_cytotoxicity():
    print("--- eos42ez (cytotoxicity, all 3 endpoints < 0.3) ---")
    path = os.path.join(ERSILIA_DIR, "eos42ez.csv")
    if not os.path.isfile(path):
        print(f"Warning: {path} not found - skipping.")
        return
    df = pd.read_csv(path)
    n = len(df)
    cols = ["cytotoxicity_hepg2", "cytotoxicity_hskmc", "cytotoxicity_imr90"]
    count = int((df[cols] < 0.3).all(axis=1).sum())
    print(f"All 3 cytotoxicity scores < 0.3: {count:,} / {n:,} ({100 * count / n:.1f}%)")


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
        ("REAL 10M - preselected", "orange", build_present(load_real_positive_scores, pockets), False),
        ("REAL 10B - preselected", "yellow", build_present(load_real2_report, pockets), True),
        ("All selected", "purple", agg_scores, True),
        ("Top-10 selected", "pink", {p: s.sort_values().head(TOP_N) for p, s in agg_scores.items()}, False),
    ]
    print("Data availability:")
    for label, _, scores, _ in score_sources:
        print(f"  {label}: {len(scores)}/{len(pockets)} pockets")

    out_path = os.path.join(OUTPUT_DIR, "aggregated_docking_score_boxplots.png")
    plot_score_boxplots_multi(pockets, score_sources, out_path, xtick_rotation=30)
    print(f"Saved: {out_path}")


def _score_upsets(agg_scores, protein_pockets, suffix, label, annotate=False):
    for threshold in DOCKING_SCORE_THRESHOLDS:
        hit_sets = {
            protein: set().union(*(
                set(agg_scores[p][agg_scores[p] <= threshold].index)
                for p in pockets_for_protein if p in agg_scores
            ))
            for protein, pockets_for_protein in protein_pockets.items()
        }
        print(f"  score <= {threshold}: " + ", ".join(f"{p}={len(s)}" for p, s in hit_sets.items()))
        out_path = os.path.join(OUTPUT_DIR, f"upset_score_{abs(threshold)}{suffix}.png")
        title_label = f"\n{label.strip()}" if label else ""
        save_protein_upset(hit_sets, f"Docking score ≤ {threshold}{title_label}", out_path, annotate=annotate)
        print(f"  Saved: {out_path}")


def _rank_upsets(pocket_ranks, protein_pockets, suffix, label, annotate=False):
    for rank in RANK_THRESHOLDS:
        hit_sets = {
            protein: set().union(*(
                set(pocket_ranks[p][pocket_ranks[p] <= rank].index)
                for p in pockets_for_protein if p in pocket_ranks
            ))
            for protein, pockets_for_protein in protein_pockets.items()
        }
        print(f"  rank <= {rank}: " + ", ".join(f"{p}={len(s)}" for p, s in hit_sets.items()))
        out_path = os.path.join(OUTPUT_DIR, f"upset_rank_{rank}{suffix}.png")
        title_label = f"\n{label.strip()}" if label else ""
        save_protein_upset(hit_sets, f"Rank ≤ {rank} in REAL 10B (~100k){title_label}", out_path, annotate=annotate)
        print(f"  Saved: {out_path}")


def plot_protein_upsets(agg_scores, annotate=False):
    """Runs the score/rank UpSet plots for all pockets, then again for CAT-only and NON-CAT-only."""
    pocket_ranks = load_pocket_ranks()
    variants = [("", None), ("_CAT", "CAT"), ("_NONCAT", "NON-CAT")]

    for suffix, site_type in variants:
        label = f" ({site_type} pockets only)" if site_type else ""
        protein_pockets = load_protein_pockets(site_type)

        print(f"--- Protein-level docking score UpSet plots{label} ---")
        _score_upsets(agg_scores, protein_pockets, suffix, label, annotate=annotate)

        print(f"--- Protein-level docking rank UpSet plots{label} ---")
        _rank_upsets(pocket_ranks, protein_pockets, suffix, label, annotate=annotate)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--annotate", action="store_true",
                         help="Draw the per-degree bracket + count label on UpSet plots (off by default).")
    parser.add_argument("--only-upset", action="store_true",
                         help="Skip everything except the protein-level UpSet plots.")
    return parser.parse_args()


def main():
    args = parse_args()
    if not args.only_upset:
        plot_ersilia_distributions()
        print()
        plot_permeability_comparison()
        print()
        report_nsps()
        print()
        report_cytotoxicity()
        print()
        report_eos2xeq_flags()
        print()
    pockets = pd.read_csv(SELECTED_POCKETS_CSV)["pocket_name"].tolist()
    agg_scores = build_present(lambda p: load_report(AGGREGATED_DOCKING_DIR, p, filename="results.csv"), pockets)
    if not args.only_upset:
        plot_docking_score_boxplots(pockets, agg_scores)
        print()
    plot_protein_upsets(agg_scores, annotate=args.annotate)


if __name__ == "__main__":
    main()
