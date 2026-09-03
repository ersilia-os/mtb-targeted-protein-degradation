"""
Figure 1.5 (main figure, sitting between Figure 1 and Figure 2): top row is panel a (65%
width) - the PocketVec t-SNE (3 reference canonical pockets highlighted by gene color,
bold-lettered callout badges pinned to fixed corners) - and panel b (35% width, currently
blank - reserved placeholder, per request). The Enamine-library-docking t-SNEs (HLL, REAL
10B) that used to sit alongside PocketVec here have moved out to their own supplementary
figure, FigSupp/tSNE_Enamine.py. Panel c is the 19-row pocket-score strip figure across all
276 detected pockets (ported from FigSupp/figure_supp_pocket.py - P2Rank info, domain
annotation, docking-score summaries, aaRS-class strips), spanning the full width below.
Panels a/b and c are promoted from supplementary to main-figure status; the 2 original
source scripts were retired once this one was confirmed working.

Built on figure_1_plot.py's panel-per-file architecture: each panel is its own standalone
figure, saved as its own PDF (Fig_1_5a.pdf, Fig_1_5b.pdf, Fig_1_5c.pdf) at an EXACT physical
size read from output/plots/figure_1_5/panel_layout.csv (columns: panel, x, delta_x, y,
delta_y, padding, in cm), then pasted onto one positioned master PDF (Fig_1_5_full.pdf, via
pypdf) and flattened to a PNG (via pymupdf). panel_layout.csv is a first-draft starting
point, meant to be tuned iteratively by rendering and looking, not solved analytically up
front.

Panel c's pocket-score data prep (residue-count PDB parsing across ~178 structures +
per-pocket docking report reads across 3 libraries x 276 pockets) is slow, so it's cached
to figure_1_5_pocket_scores.csv and only recomputed with --rerun - same convention as
figure_1_plot.py's own --rerun for its slow PyMOL rendering step.

Usage:
    python figure_1_5_plot.py [--rerun] [--subpanels a,b,c]
"""
import argparse
import json
import math
import os
import sys

os.environ["QT_QPA_PLATFORM"] = "offscreen"

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymupdf
import stylia
from pypdf import PdfReader, PdfWriter, Transformation
from stylia.config import get_fg_color
from stylia.figure.figure import stylize

from default import RANDOM_SEED
from pocket_group_plots import calculate_density, density_to_sizes, load_canonical_pocket_labels
from pocketvec_tsne_embedding import compute_pocketvec_embedding

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

output_dir = os.path.join(root, "..", "..", "output")
plots_dir = os.path.join(root, "..", "..", "output", "plots", "figure_1_5")
os.makedirs(plots_dir, exist_ok=True)

PANEL_LETTERS = ["a", "b", "c"]
panel_layout_path = os.path.join(plots_dir, "panel_layout.csv")

# ===========================================================================
# Panel a data + helper: t-SNE comparison
# (ported from FigSupp/reference_pockets_tsne_comparison.py)
# ===========================================================================

TARGET_CLUSTERS = ["gltS_cluster1", "tyrS_cluster1", "ileS_cluster1"]
# Mathtext ($\bf{...}$) bolds only the letter, not "Pocket" or the gene name, per request.
CLUSTER_LABELS = {
    "ileS_cluster1": r"Pocket $\bf{A}$ (ileS)",
    "tyrS_cluster1": r"Pocket $\bf{B}$ (tyrS)",
    "gltS_cluster1": r"Pocket $\bf{C}$ (gltS)",
}

with open(os.path.join(output_dir, "plots", "figure_1", "color_mapping.json")) as f:
    GENE_TO_COLOR = json.load(f)["gene_to_color"]


def cluster_color(cluster):
    gene = cluster.split("_cluster")[0]
    return GENE_TO_COLOR[gene]


POINT_MIN_SIZE, POINT_MAX_SIZE = 5, 150  # wider density-size range than pocket_group_plots's default 10-70

# Fixed corner per cluster, in axes-fraction coordinates - same as a legend pinned to a
# location ("upper left", etc.) rather than placed relative to the data: (x, y, ha, va).
LABEL_CORNERS = {
    "ileS_cluster1": (0.03, 0.97, "left", "top"),      # top-left
    "tyrS_cluster1": (0.97, 0.97, "right", "top"),     # top-right
    "gltS_cluster1": (0.03, 0.03, "left", "bottom"),   # bottom-left
}

# HLL docking / REAL 10B docking t-SNEs moved out to FigSupp/tSNE_Enamine.py (per request) -
# panel a is now just the single PocketVec embedding, so it carries the reference-pocket
# callout badges itself instead of only the (now-removed) HLL panel.
EMBEDDING = ("PocketVec", lambda: compute_pocketvec_embedding(output_dir, RANDOM_SEED))


def plot_tsne_panel(ax, coords, canonical, title, annotate):
    nc = stylia.NamedColors()
    canonical = np.array(canonical)
    is_target = np.isin(canonical, TARGET_CLUSTERS)

    # "Other" points: light gray, sized by local 2D density (same convention as
    # pocketvec_tsne.py's standalone density plot, widened min/max range per request); the 3
    # reference clusters keep their own fixed color on top so they stay identifiable
    # regardless of local density.
    sizes = density_to_sizes(calculate_density(coords), min_size=POINT_MIN_SIZE, max_size=POINT_MAX_SIZE)
    light_gray = nc.get("silver", lighten=0.5)

    # Explicit z-order throughout (background < connector line < each cluster's own dots <
    # caption badge) so the stacking is guaranteed regardless of call order, per request: the
    # line must sit above the background but below the dots of the cluster it points to.
    ax.scatter(coords[~is_target, 0], coords[~is_target, 1], color=light_gray, s=sizes[~is_target], zorder=1)
    for cluster in TARGET_CLUSTERS:
        mask = canonical == cluster
        if not mask.any():
            continue
        ax.scatter(coords[mask, 0], coords[mask, 1], color=cluster_color(cluster), s=sizes[mask], zorder=3)

    # Fix the view limits now (data range + a bit of headroom), *before* placing labels, so
    # the caption boxes never cross the panel border.
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()
    xlim = (x_min - 0.1 * (x_max - x_min), x_max + 0.1 * (x_max - x_min))
    ylim = (y_min - 0.15 * (y_max - y_min), y_max + 0.15 * (y_max - y_min))
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if annotate:
        for cluster in TARGET_CLUSTERS:
            mask = canonical == cluster
            if not mask.any():
                continue
            # One ax.annotate() call (not split into 2 artists) so matplotlib auto-shrinks the
            # connector line to stop at the badge's own edge instead of running into the text
            # underneath it. zorder=2 puts the whole thing (line + badge) above the background
            # (zorder=1) but below this cluster's own dots (zorder=3) - the badge itself never
            # actually overlaps those dots anyway, since it's pinned to a far corner.
            cx, cy = coords[mask, 0].mean(), coords[mask, 1].mean()
            lx, ly, ha, va = LABEL_CORNERS[cluster]
            ax.annotate(CLUSTER_LABELS[cluster], xy=(cx, cy), xycoords="data",
                        xytext=(lx, ly), textcoords="axes fraction",
                        ha=ha, va=va, fontsize=stylia.FONTSIZE, color="black", zorder=2,
                        bbox=dict(facecolor=cluster_color(cluster), edgecolor="black",
                                  alpha=0.6, boxstyle="square,pad=0.3"),
                        # linewidth explicit: stylia's article style sets rcParam
                        # patch.linewidth=0, which silently zeroes a FancyArrowPatch's stroke
                        # (the connector line) unless overridden here.
                        arrowprops=dict(arrowstyle="-", color="black", linewidth=stylia.LINEWIDTH))
    stylia.label(ax, xlabel="t-SNE 1", ylabel="t-SNE 2", title=title)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_box_aspect(1)


# ===========================================================================
# Panel c data + helpers: pocket-score strips across all 276 detected pockets
# (ported from FigSupp/figure_supp_pocket.py)
# ===========================================================================

# A pocket counts as "catalytic" at catalytic_confidence >= 3 (either strong direct-PDB or
# strong AlphaFill ligand evidence for the Catalytic Domain (ATP/ligase) label - see
# scripts/77_pocket_annotation/10_assemble_final_table.py's confidence scale, 0-4), vs. the
# rest (weak/no evidence or no catalytic label at all, confidence 0-2) - a threshold the
# user specified directly, not chosen here.
CATALYTIC_CONFIDENCE_MIN = 3

# The other 3 domain rows (tRNA binding, Editing, Anticodon binding) have no ligand-evidence
# confidence score like catalytic_confidence - only a binary "this InterPro domain label is
# present for this pocket" (curated_labels). Present there means label present AND
# catalytic_confidence < CATALYTIC_CONFIDENCE_MIN - mutually exclusive with the Catalytic row,
# since a pocket can carry both a catalytic and a non-catalytic label at once.
CURATED_LABEL_COLUMNS = {
    "is_trna_binding": "tRNA Binding Domain",
    "is_editing": "Editing Domain",
    "is_anticodon_binding": "Anticodon Binding Domain",
}

# aaRS class per gene, ported from figure_1_plot.py's own AARS_CLASS_LABELS - gatA/gatB
# aren't canonical aaRS classes, so they fall to "Other" via the .get() default.
AARS_CLASS_LABELS = {
    "alaS": "Class II", "argS": "Class I", "aspS": "Class II", "cysS1": "Class I",
    "gltS": "Class I", "glyS": "Class II", "hisS": "Class II", "ileS": "Class I",
    "leuS": "Class I", "lysS": "Class II", "metS": "Class I", "pheS": "Class II",
    "pheT": "Class II", "proS": "Class II", "serS": "Class II", "thrS": "Class II",
    "trpS": "Class I", "tyrS": "Class I", "valS": "Class I",
}

# Gene name + color per pocket's own protein, straight from figure_1_calculations.py's own
# output (output/plots/figure_1/color_mapping.json) - "following the coloring scheme from
# figure 1" per request, rather than recomputing a separate palette here.
with open(os.path.join(output_dir, "plots", "figure_1", "color_mapping.json")) as f:
    figure_1_mappings = json.load(f)
uniprot_to_gene = figure_1_mappings["uniprot_to_gene"]
gene_to_color = figure_1_mappings["gene_to_color"]

# Number of pocket residues within RESIDUE_DISTANCE_RADII (6A, 12A) of each pocket's own
# centroid - ported from scripts/91_human_detect_pockets.py's own parse_residue_geometry/
# residues_within_radius (pLDDT/B-factor dropped, not needed here). Structure files are
# cached per (Uniprot AC, File name) since many pockets share the same structure.
RESIDUE_DISTANCE_RADII = [6.0, 12.0]
STRUCTURES_DIR = os.path.join(output_dir, "aligned_relaxed_structures")


def parse_residue_geometry(pdb_file):
    """Residue number -> (x, y, z), from each residue's CA atom."""
    res_geom = {}
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                resnum = int(line[22:26])
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                res_geom[resnum] = (x, y, z)
    return res_geom


_res_geom_cache = {}


def get_residue_geometry(uniprot_ac, file_name):
    key = (uniprot_ac, file_name)
    if key not in _res_geom_cache:
        _res_geom_cache[key] = parse_residue_geometry(os.path.join(STRUCTURES_DIR, uniprot_ac, file_name))
    return _res_geom_cache[key]


# Per-pocket docking score summary (median + 1st percentile) against 3 compound libraries:
# "hll" = the Enamine Hit Locator Library 100k
# (data/Enamine/Enamine_Hit_Locator_Library_100K_Set_*.smiles, ~100k compounds/pocket,
# docked in output/unidock_docking/docking_results/); "real10m" = the 9.56M Enamine REAL
# library's own round-1 docking (~113k compounds/pocket,
# output/unidock_REAL_docking/docking_results/); "real10b" = the ~100k compounds selected
# from the 10B Enamine REAL library for final validation docking
# (output/unidock_REAL_docking_2/docking_results/). Each pocket's own report.csv (columns:
# compound, score) lives at <docking_dir>/<file_stem>_pocket_<pocket_number>/report.csv.
DOCKING_LIBRARIES = {
    "hll": os.path.join(output_dir, "unidock_docking", "docking_results"),
    "real10m": os.path.join(output_dir, "unidock_REAL_docking", "docking_results"),
    "real10b": os.path.join(output_dir, "unidock_REAL_docking_2", "docking_results"),
}

scores_path = os.path.join(plots_dir, "figure_1_5_pocket_scores.csv")


def _compute_pocket_scores():
    pockets = pd.read_csv(os.path.join(output_dir, "pocket_detection_data.csv"))
    pockets = pockets.sort_values("Pocket probability", ascending=False).reset_index(drop=True)

    interpro = pd.read_csv(os.path.join(output_dir, "77_pocket_annotation", "pocket_detection_interpro_updated.csv"),
                            keep_default_na=False)
    pockets = pockets.merge(
        interpro[["Uniprot AC", "File name", "Pocket number", "catalytic_confidence", "curated_labels",
                  "has_direct_ligand_evidence", "alphafill_ligands"]],
        on=["Uniprot AC", "File name", "Pocket number"], how="left",
    )

    scores = pockets[["Uniprot AC", "File name", "Pocket number", "Pocket probability", "Pocket score", "catalytic_confidence"]].copy()
    scores.insert(0, "pocket_rank", range(1, len(scores) + 1))
    scores = scores.rename(columns={"Pocket probability": "pocket_probability", "Pocket score": "pocket_score"})
    scores["is_catalytic"] = scores["catalytic_confidence"] >= CATALYTIC_CONFIDENCE_MIN
    has_label = {
        col: pockets["curated_labels"].apply(lambda labels, label=label: label in labels.split("|"))
        for col, label in CURATED_LABEL_COLUMNS.items()
    }
    for col, mask in has_label.items():
        scores[col] = mask & ~scores["is_catalytic"]

    # 3-way ligand evidence, per request: "direct" (has_direct_ligand_evidence == "yes", a real
    # ligand seen in an experimental PDB structure of this protein) takes priority over
    # "alphafill" (alphafill_ligands non-empty - a ligand only transplanted in from a homolog,
    # weaker evidence), else "none".
    has_direct_ligand = pockets["has_direct_ligand_evidence"].str.lower() == "yes"
    has_alphafill_ligand = pockets["alphafill_ligands"] != ""
    scores["ligand_evidence"] = pd.Series("none", index=pockets.index)
    scores.loc[has_alphafill_ligand, "ligand_evidence"] = "alphafill"
    scores.loc[has_direct_ligand, "ligand_evidence"] = "direct"

    scores["gene"] = scores["Uniprot AC"].map(uniprot_to_gene)
    scores["aars_class"] = scores["gene"].map(lambda g: AARS_CLASS_LABELS.get(g, "Other"))

    for radius in RESIDUE_DISTANCE_RADII:
        counts = []
        for _, row in pockets.iterrows():
            res_geom = get_residue_geometry(row["Uniprot AC"], row["File name"])
            cx, cy, cz = (float(v) for v in row["Pocket centroid coordinate (x y z)"].split())
            counts.append(sum(1 for x, y, z in res_geom.values() if math.dist((x, y, z), (cx, cy, cz)) <= radius))
        scores[f"n_residues_{int(radius)}A"] = counts
        print(f"N residues within {radius:.0f}A: min={min(counts)}, max={max(counts)}")

    for lib_key, docking_dir in DOCKING_LIBRARIES.items():
        medians, perc1s = [], []
        for _, row in pockets.iterrows():
            pocket_dir = row["File name"].replace(".pdb", "") + "_pocket_" + str(row["Pocket number"])
            report_scores = pd.read_csv(os.path.join(docking_dir, pocket_dir, "report.csv"), usecols=["score"])["score"]
            medians.append(report_scores.median())
            perc1s.append(report_scores.quantile(0.01))
        scores[f"median_{lib_key}"] = medians
        scores[f"perc1_{lib_key}"] = perc1s
        print(f"{lib_key} median: min={min(medians):.2f}, max={max(medians):.2f}; "
              f"perc1: min={min(perc1s):.2f}, max={max(perc1s):.2f}")

    return scores


def prepare_panel_c_data(rerun=False):
    """Loads (or recomputes + caches) pocket_scores, then derives every colormap Normalize
    that depends on its data range - both as module globals, since every plot_*_row function
    below references them by bare name (unchanged from figure_supp_pocket.py)."""
    global pocket_scores
    global P2RANK_NORM, P2RANK_SCORE_NORM, RESIDUE_COUNT_NORMS, CATALYTIC_CONFIDENCE_NORM, DOCKING_NORMS

    if not rerun and os.path.exists(scores_path):
        pocket_scores = pd.read_csv(scores_path)
        print(f"Loaded cached pocket scores from {scores_path}")
    else:
        pocket_scores = _compute_pocket_scores()
        pocket_scores.to_csv(scores_path, index=False)
        print(f"Saved {len(pocket_scores):,} row(s) to {scores_path}")

    P2RANK_NORM = mcolors.Normalize(vmin=0, vmax=1.0, clip=True)
    P2RANK_SCORE_NORM = mcolors.Normalize(vmin=0, vmax=pocket_scores["pocket_score"].max(), clip=True)
    RESIDUE_COUNT_NORMS = {
        radius: mcolors.Normalize(vmin=0, vmax=pocket_scores[f"n_residues_{int(radius)}A"].max(), clip=True)
        for radius in RESIDUE_DISTANCE_RADII
    }
    CATALYTIC_CONFIDENCE_NORM = mcolors.Normalize(vmin=0, vmax=4, clip=True)
    DOCKING_NORMS = {
        column: mcolors.Normalize(vmin=pocket_scores[column].quantile(0.10), vmax=pocket_scores[column].quantile(0.90), clip=True)
        for lib_key in DOCKING_LIBRARIES
        for column in (f"median_{lib_key}", f"perc1_{lib_key}")
    }


# White-anchored version of a stylia FadingColormap: its own pale end is a light tint
# (e.g. "plum"'s pale lilac, "crimson"'s blush), not true white - this instead fades from
# true white at 0 up to that colormap's own full-saturation endpoint (cmap(1.0)) at 1.
def white_to_color_cmap(base_cmap):
    return mcolors.LinearSegmentedColormap.from_list("white_to_color", ["white", base_cmap(1.0)])


# Reverse direction of white_to_color_cmap - full color at 0, true white at 1. Docking
# scores are lower-is-better, so pairing this with Normalize(vmin=best/lowest,
# vmax=worst/highest) puts full color on the best score and white on the worst, instead of
# the other way around.
def color_to_white_cmap(base_cmap):
    return mcolors.LinearSegmentedColormap.from_list("color_to_white", [base_cmap(1.0), "white"])


# Gradient for the P2Rank probability row - true white at 0 up to full "plum" (ersilia's
# own FadingColormap preset) at 1.0, no cutoff - a plain continuous gradient over the
# pocket's full probability range, per request.
P2RANK_CMAP = white_to_color_cmap(stylia.FadingColormap("plum", transformation=None).cmap)

# Gradient for the P2Rank raw-score row (pocket_score, unbounded - unlike pocket_probability
# it isn't capped at 1.0). Same "plum" preset as the probability row above, per request.
# vmax is the data's own max, not a fixed cap - so the pocket with the single highest raw
# score gets the full deep-plum color, not a flattened-out plateau.
P2RANK_SCORE_CMAP = stylia.FadingColormap("plum", transformation=None).cmap

# Gradient for the 2 residue-count rows (n_residues_6A, n_residues_12A) - same "plum"
# colormap as the other 2 P2Rank rows above, per request ("colored with the same color").
RESIDUE_COUNT_CMAP = stylia.FadingColormap("plum", transformation=None).cmap

# Gradient for the continuous catalytic_confidence row (0-4, see CATALYTIC_CONFIDENCE_MIN
# above) - true white at 0 up to full "crimson" (same red family as the binarized
# "Catalytic" domain strip) at 4, per request.
CATALYTIC_CONFIDENCE_CMAP = white_to_color_cmap(stylia.FadingColormap("crimson", transformation=None).cmap)

# The 4 domain bands: column in pocket_scores -> band title -> its own "present" color,
# against a shared white "absent" background. Catalytic uses a ligand-evidence confidence
# threshold (CATALYTIC_CONFIDENCE_MIN); the other 3 use plain InterPro label presence,
# mutually exclusive with Catalytic - see CURATED_LABEL_COLUMNS above.
DOMAIN_STRIP_COLUMNS = [
    ("is_catalytic", "Catalytic", "crimson"),
    ("is_trna_binding", "tRNA binding", "cobalt"),
    ("is_editing", "Editing", "amber"),
    ("is_anticodon_binding", "Anticodon binding", "lime"),
]

# 3-way ligand evidence row color, per request: black = direct (real ligand seen in an
# experimental PDB structure), gray = alphafill (only transplanted in from a homolog),
# white = none.
LIGAND_EVIDENCE_COLORS = {"direct": "black", "alphafill": "gray", "none": "white"}

# Each domain-row "present" mark's width, in pocket_rank units (1 unit = 1 pocket's own
# rank spacing) - 2x a pocket's natural 1-unit spacing, for better legibility.
DOMAIN_ROW_BAR_WIDTH = 2.0

# The 21 individual gene rows collapsed into 3, per request - Class I / Class II (from
# AARS_CLASS_LABELS) / Other (gatA, gatB).
AARS_CLASSES = ["Class I", "Class II", "Other"]

# Docking-score section (median + 1st percentile, 3 libraries, see DOCKING_LIBRARIES
# above) - red ("crimson") for HL, turquoise for REAL 10M, yellow ("amber") for REAL 10B,
# per request. Each row's vmin/vmax is its own column's 10th/90th percentile (see
# DOCKING_NORMS in prepare_panel_c_data), not its true min/max: 1-2 outlier pockets per
# column sit far past the rest, and a true-min/max scale would compress the other 274
# pockets into a barely-distinguishable sliver of the gradient - a call the user made
# directly (tried p1/p99, then p5/p95, settled on p10/p90), not chosen here.
_nc = stylia.NamedColors()
DOCKING_CMAPS = {
    "hll": color_to_white_cmap(stylia.FadingColormap("crimson", transformation=None).cmap),
    "real10m": color_to_white_cmap(stylia.FadingColormap("turquoise", transformation=None).cmap),
    "real10b": mcolors.LinearSegmentedColormap.from_list("amber_to_white", [_nc.amber, "white"]),
}

# The 6 docking-score rows, in order: column in pocket_scores -> row label -> library key
# (column indexes DOCKING_NORMS, library key indexes DOCKING_CMAPS). REAL 10M sits in the
# middle, per request.
DOCKING_ROWS = [
    ("median_hll", "med. HL", "hll"),
    ("perc1_hll", "perc1 HL", "hll"),
    ("median_real10m", "med. REAL 10M", "real10m"),
    ("perc1_real10m", "perc1 REAL 10M", "real10m"),
    ("median_real10b", "med. REAL 10B", "real10b"),
    ("perc1_real10b", "perc1 REAL 10B", "real10b"),
]

# Vertical gap between the 4 row blocks (P2Rank info / domain strips / docking scores /
# aaRS-class strips) - each block's own rows sit flush against each other (hspace=0, no
# gap). This is hspace on the OUTER 4-row gridspec (not a single flat gridspec across all
# rows), so the same fraction reads as a much bigger absolute gap.
ROWS_HSPACE = 0.2

# Fixed absolute left margin, wide enough to fit the longest row label ("Anticodon
# binding") at its own fixed point size.
ROW_LABEL_LEFT_MARGIN_IN = 0.75

# P2Rank block: pocket_probability + pocket_score + 2 residue-count rows +
# plot_ligand_evidence_row.
N_P2RANK_ROWS = 2 + len(RESIDUE_DISTANCE_RADII) + 1
# Domain block: DOMAIN_STRIP_COLUMNS (4) + plot_catalytic_confidence_row (1).
N_DOMAIN_ROWS = len(DOMAIN_STRIP_COLUMNS) + 1


def plot_pocket_scores_row(ax):
    """1st of 6 stacked rows (alongside plot_pocket_score_row and the 4
    plot_domain_strip_row bands) - a gradient strip, not a distribution: one thin cell per
    pocket (same pocket_rank x-order as plot_domain_strip_row), colored by
    P2RANK_CMAP/P2RANK_NORM directly from its own probability value (true white at 0 up to
    deep plum at 1.0, continuous - no cutoff)."""
    colors = P2RANK_CMAP(P2RANK_NORM(pocket_scores["pocket_probability"]))
    ax.bar(pocket_scores["pocket_rank"], 1, width=DOMAIN_ROW_BAR_WIDTH, color=colors, edgecolor="none", zorder=2)
    ax.set_xlim(pocket_scores["pocket_rank"].min(), pocket_scores["pocket_rank"].max())
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    stylia.label(ax, xlabel="", ylabel="P2Rank prob.")
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.label.set_ha("right")
    ax.yaxis.label.set_va("center")


def plot_pocket_score_row(ax):
    """2nd of 6 stacked rows - same gradient-strip treatment as plot_pocket_scores_row, but
    for the raw P2Rank pocket_score instead of pocket_probability (P2RANK_SCORE_CMAP/
    P2RANK_SCORE_NORM: pale plum at 0 -> deep plum at the data's own max, no fixed cap)."""
    colors = P2RANK_SCORE_CMAP(P2RANK_SCORE_NORM(pocket_scores["pocket_score"]))
    ax.bar(pocket_scores["pocket_rank"], 1, width=DOMAIN_ROW_BAR_WIDTH, color=colors, edgecolor="none", zorder=2)
    ax.set_xlim(pocket_scores["pocket_rank"].min(), pocket_scores["pocket_rank"].max())
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    stylia.label(ax, xlabel="", ylabel="P2Rank score")
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.label.set_ha("right")
    ax.yaxis.label.set_va("center")


def plot_residue_count_row(ax, radius):
    """3rd/4th of 6 stacked rows - same gradient-strip treatment as plot_pocket_scores_row/
    plot_pocket_score_row, but for the number of residues within `radius` (6 or 12,
    RESIDUE_DISTANCE_RADII) Angstrom of the pocket's own centroid (RESIDUE_COUNT_CMAP/
    RESIDUE_COUNT_NORMS[radius])."""
    column = f"n_residues_{int(radius)}A"
    colors = RESIDUE_COUNT_CMAP(RESIDUE_COUNT_NORMS[radius](pocket_scores[column]))
    ax.bar(pocket_scores["pocket_rank"], 1, width=DOMAIN_ROW_BAR_WIDTH, color=colors, edgecolor="none", zorder=2)
    ax.set_xlim(pocket_scores["pocket_rank"].min(), pocket_scores["pocket_rank"].max())
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    stylia.label(ax, xlabel="", ylabel=f"Res. at {radius:.0f}A")
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.label.set_ha("right")
    ax.yaxis.label.set_va("center")


def plot_domain_strip_row(ax, column, title, color_name):
    """One of the 4 domain-strip rows - a thin colored cell per pocket (all 276, same
    pocket_rank x-order - sorted by P2Rank probability descending, rank 1 leftmost), in
    this band's own color where `column` is True, white otherwise."""
    nc = stylia.NamedColors()
    present_color = getattr(nc, color_name)
    colors = [present_color if v else nc.white for v in pocket_scores[column]]
    ax.bar(pocket_scores["pocket_rank"], 1, width=DOMAIN_ROW_BAR_WIDTH, color=colors, edgecolor="none", zorder=2)
    ax.set_xlim(pocket_scores["pocket_rank"].min(), pocket_scores["pocket_rank"].max())
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    stylia.label(ax, xlabel="", ylabel=title)
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.label.set_ha("right")
    ax.yaxis.label.set_va("center")


def plot_catalytic_confidence_row(ax):
    """Continuous companion to the binarized "Catalytic" domain strip - same gradient-strip
    treatment as plot_pocket_scores_row/plot_pocket_score_row, but for catalytic_confidence
    (0-4) via CATALYTIC_CONFIDENCE_CMAP/CATALYTIC_CONFIDENCE_NORM (true white at 0 up to
    full crimson at 4)."""
    colors = CATALYTIC_CONFIDENCE_CMAP(CATALYTIC_CONFIDENCE_NORM(pocket_scores["catalytic_confidence"]))
    ax.bar(pocket_scores["pocket_rank"], 1, width=DOMAIN_ROW_BAR_WIDTH, color=colors, edgecolor="none", zorder=2)
    ax.set_xlim(pocket_scores["pocket_rank"].min(), pocket_scores["pocket_rank"].max())
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    stylia.label(ax, xlabel="", ylabel="Catalytic conf.")
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.label.set_ha("right")
    ax.yaxis.label.set_va("center")


def plot_ligand_evidence_row(ax):
    """Last row of the P2Rank info block - 3-way ligand_evidence (LIGAND_EVIDENCE_COLORS:
    black "direct" / gray "alphafill" / white "none")."""
    colors = [LIGAND_EVIDENCE_COLORS[v] for v in pocket_scores["ligand_evidence"]]
    ax.bar(pocket_scores["pocket_rank"], 1, width=DOMAIN_ROW_BAR_WIDTH, color=colors, edgecolor="none", zorder=2)
    ax.set_xlim(pocket_scores["pocket_rank"].min(), pocket_scores["pocket_rank"].max())
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    stylia.label(ax, xlabel="", ylabel="Ligand evidence")
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.label.set_ha("right")
    ax.yaxis.label.set_va("center")


def plot_class_strip_row(ax, aars_class):
    """One of the 3 collapsed aaRS-class strips (AARS_CLASSES) - a thin colored cell per
    pocket belonging to this class, still in that pocket's own PROTEIN color (gene_to_color,
    same as figure_1_plot.py), not one flat per-class color - white for every pocket outside
    this class."""
    nc = stylia.NamedColors()
    colors = [
        gene_to_color[gene] if cls == aars_class else nc.white
        for gene, cls in zip(pocket_scores["gene"], pocket_scores["aars_class"])
    ]
    ax.bar(pocket_scores["pocket_rank"], 1, width=DOMAIN_ROW_BAR_WIDTH, color=colors, edgecolor="none", zorder=2)
    ax.set_xlim(pocket_scores["pocket_rank"].min(), pocket_scores["pocket_rank"].max())
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    stylia.label(ax, xlabel="", ylabel=aars_class)
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.label.set_ha("right")
    ax.yaxis.label.set_va("center")


def plot_docking_score_row(ax, column, label, lib_key):
    """One of the 6 docking-score rows (DOCKING_ROWS) - same gradient-strip treatment as
    the other continuous rows, but full color = BEST (lowest) score, true white = WORST
    (highest) score, via DOCKING_CMAPS[lib_key] (shared per library) and DOCKING_NORMS
    [column] (this row's own, independent of every other row)."""
    colors = DOCKING_CMAPS[lib_key](DOCKING_NORMS[column](pocket_scores[column]))
    ax.bar(pocket_scores["pocket_rank"], 1, width=DOMAIN_ROW_BAR_WIDTH, color=colors, edgecolor="none", zorder=2)
    ax.set_xlim(pocket_scores["pocket_rank"].min(), pocket_scores["pocket_rank"].max())
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    stylia.label(ax, xlabel="", ylabel=label)
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.label.set_ha("right")
    ax.yaxis.label.set_va("center")


# ===========================================================================
# Panel-saving scaffolding (ported from figure_1_plot.py)
# ===========================================================================

PANEL_LABEL_MARGIN = 0.02


def add_panel_label(fig, letter):
    """Bold panel letter at the top-left of the FIGURE (page), fixed regardless of padding."""
    fig.text(PANEL_LABEL_MARGIN, 1 - PANEL_LABEL_MARGIN, letter, fontweight="bold",
              fontsize=stylia.FONTSIZE_BIG, color=get_fg_color(), ha="left", va="top",
              transform=fig.transFigure)


MAX_WIDTH_IN = stylia.SIZE  # Nature two-column guideline (~7.09in/18cm)


def load_panel_sizes(panels):
    if not os.path.exists(panel_layout_path):
        raise FileNotFoundError(
            f"{panel_layout_path} not found. Expected columns: panel, x, delta_x, y, delta_y (cm), "
            "one row per panel: " + ", ".join(PANEL_LETTERS))
    df = pd.read_csv(panel_layout_path).set_index("panel")
    missing = [p for p in panels if p not in df.index]
    if missing:
        raise ValueError(f"{panel_layout_path} is missing row(s) for panel(s): {missing}")
    return {p: (df.loc[p, "delta_x"] / 2.54, df.loc[p, "delta_y"] / 2.54) for p in panels}


def load_panel_padding(panels):
    df = pd.read_csv(panel_layout_path).set_index("panel")
    if "padding" not in df.columns:
        return {p: 0.0 for p in panels}
    missing = [p for p in panels if p not in df.index]
    if missing:
        raise ValueError(f"{panel_layout_path} is missing row(s) for panel(s): {missing}")
    return {p: df.loc[p, "padding"] for p in panels}


def apply_padding(fig, padding):
    scale = 1 - padding
    for ax in fig.axes:
        pos = ax.get_position()
        x0 = 0.5 + (pos.x0 - 0.5) * scale
        x1 = 0.5 + (pos.x1 - 0.5) * scale
        y0 = 0.5 + (pos.y0 - 0.5) * scale
        y1 = 0.5 + (pos.y1 - 0.5) * scale
        ax.set_position([x0, y0, x1 - x0, y1 - y0])


def save_panel(fig, letter, use_tight_layout=True, tight_pad=1.08, tight_w_pad=None,
               subplots_adjust=None, padding=0.0):
    """Saves at exactly panel_layout.csv's delta_x/delta_y (no bbox_inches="tight")."""
    fig.set_tight_layout(False)
    if use_tight_layout:
        plt.tight_layout(pad=tight_pad, w_pad=tight_w_pad)
    else:
        fig.subplots_adjust(**(subplots_adjust or dict(left=0.01, right=0.99, top=0.99, bottom=0.01)))
    apply_padding(fig, padding)
    add_panel_label(fig, letter)
    output_path = os.path.join(plots_dir, f"Fig_1_5{letter}.pdf")
    plt.savefig(output_path, dpi=600, transparent=False)
    plt.close(fig)
    print(f"Saved Fig_1_5{letter}.pdf")


CM_TO_PT = 72 / 2.54


def merge_panels():
    """Pastes each Fig_1_5{letter}.pdf onto one positioned master PDF (Fig_1_5_full.pdf) per
    panel_layout.csv's x/y - pure translation (no rescaling, true vector via pypdf), since
    each panel already saved at its exact target size."""
    df = pd.read_csv(panel_layout_path).set_index("panel")
    missing = [p for p in PANEL_LETTERS if p not in df.index]
    if missing:
        raise ValueError(f"{panel_layout_path} is missing row(s) for panel(s): {missing}")

    total_width_cm = max(df.loc[p, "x"] + df.loc[p, "delta_x"] for p in PANEL_LETTERS)
    total_height_cm = max(df.loc[p, "y"] + df.loc[p, "delta_y"] for p in PANEL_LETTERS)

    writer = PdfWriter()
    master_page = writer.add_blank_page(width=total_width_cm * CM_TO_PT, height=total_height_cm * CM_TO_PT)

    for p in PANEL_LETTERS:
        panel_path = os.path.join(plots_dir, f"Fig_1_5{p}.pdf")
        if not os.path.exists(panel_path):
            print(f"  Skipping {p}: {panel_path} not found - Fig_1_5_full.pdf will be missing it.")
            continue
        panel_page = PdfReader(panel_path).pages[0]
        x_top_cm, y_top_cm, delta_y_cm = df.loc[p, "x"], df.loc[p, "y"], df.loc[p, "delta_y"]
        x_pt = x_top_cm * CM_TO_PT
        y_bottom_pt = (total_height_cm - y_top_cm - delta_y_cm) * CM_TO_PT
        master_page.merge_transformed_page(panel_page, Transformation().translate(tx=x_pt, ty=y_bottom_pt))
        print(f"  Placed Fig_1_5{p}.pdf at x={x_top_cm}cm, y={y_top_cm}cm (top-left)")

    output_path = os.path.join(plots_dir, "Fig_1_5_full.pdf")
    with open(output_path, "wb") as f:
        writer.write(f)
    print(f"Saved merged master figure ({total_width_cm:.2f} x {total_height_cm:.2f} cm) to {output_path}")

    png_path = os.path.join(plots_dir, "Fig_1_5_full.png")
    pdf_doc = pymupdf.open(output_path)
    zoom = 600 / 72  # pymupdf's Pixmap render defaults to 72dpi
    pix = pdf_doc[0].get_pixmap(matrix=pymupdf.Matrix(zoom, zoom))
    pix.save(png_path)
    pdf_doc.close()
    print(f"Saved {png_path}")


# ===========================================================================
# Panel builders
# ===========================================================================

def build_panel_a(size, padding):
    title, embed_fn = EMBEDDING
    keys, coords = embed_fn()
    canonical = load_canonical_pocket_labels(output_dir, keys)

    fig, ax = plt.subplots(figsize=size)
    fig.patch.set_facecolor("white")
    stylize(ax)
    plot_tsne_panel(ax, coords, canonical, title, annotate=True)
    save_panel(fig, "a", use_tight_layout=False,
               subplots_adjust=dict(left=0.15, right=0.85, top=0.85, bottom=0.15), padding=padding)


def build_panel_b(size, padding):
    """Reserved placeholder, top row next to panel a (35% width) - left blank for now,
    per request. Still gets a panel letter/border-free blank page so it participates in
    merge_panels() and can be filled in later without touching the layout."""
    fig, ax = plt.subplots(figsize=size)
    fig.patch.set_facecolor("white")
    ax.set_axis_off()
    save_panel(fig, "b", use_tight_layout=False,
               subplots_adjust=dict(left=0.01, right=0.99, top=0.99, bottom=0.01), padding=padding)


def build_panel_c(size, padding, rerun=False):
    prepare_panel_c_data(rerun=rerun)

    fig = plt.figure(figsize=size)
    fig.patch.set_facecolor("white")
    # 4 blocks, in order: P2Rank info, domain strips, docking scores, aaRS-class strips -
    # each with its own inner subgridspec at hspace=0 (no gap between rows within a block),
    # while ROWS_HSPACE is kept only between the blocks themselves.
    block_sizes = [N_P2RANK_ROWS, N_DOMAIN_ROWS, len(DOCKING_ROWS), len(AARS_CLASSES)]
    outer_gs = fig.add_gridspec(4, 1, height_ratios=block_sizes, hspace=ROWS_HSPACE)

    p2rank_gs = outer_gs[0, 0].subgridspec(block_sizes[0], 1, hspace=0)
    plot_pocket_scores_row(stylize(fig.add_subplot(p2rank_gs[0, 0])))
    plot_pocket_score_row(stylize(fig.add_subplot(p2rank_gs[1, 0])))
    for i, radius in enumerate(RESIDUE_DISTANCE_RADII):
        plot_residue_count_row(stylize(fig.add_subplot(p2rank_gs[2 + i, 0])), radius)
    plot_ligand_evidence_row(stylize(fig.add_subplot(p2rank_gs[2 + len(RESIDUE_DISTANCE_RADII), 0])))

    domain_gs = outer_gs[1, 0].subgridspec(block_sizes[1], 1, hspace=0)
    row_idx = 0
    for column, title, color_name in DOMAIN_STRIP_COLUMNS:
        ax = stylize(fig.add_subplot(domain_gs[row_idx, 0]))
        plot_domain_strip_row(ax, column=column, title=title, color_name=color_name)
        row_idx += 1
        if column == "is_catalytic":
            plot_catalytic_confidence_row(stylize(fig.add_subplot(domain_gs[row_idx, 0])))
            row_idx += 1

    docking_gs = outer_gs[2, 0].subgridspec(block_sizes[2], 1, hspace=0)
    for i, (column, label, lib_key) in enumerate(DOCKING_ROWS):
        ax = stylize(fig.add_subplot(docking_gs[i, 0]))
        plot_docking_score_row(ax, column=column, label=label, lib_key=lib_key)

    class_gs = outer_gs[3, 0].subgridspec(block_sizes[3], 1, hspace=0)
    for i, aars_class in enumerate(AARS_CLASSES):
        ax = stylize(fig.add_subplot(class_gs[i, 0]))
        plot_class_strip_row(ax, aars_class)

    save_panel(fig, "c", use_tight_layout=False,
               subplots_adjust=dict(left=ROW_LABEL_LEFT_MARGIN_IN / size[0], right=0.995, top=0.93, bottom=0.01),
               padding=padding)


def main(rerun=False, subpanels=None):
    if subpanels is None:
        subpanels = PANEL_LETTERS
    sizes = load_panel_sizes(subpanels)
    paddings = load_panel_padding(subpanels)

    for letter in subpanels:
        if sizes[letter][0] > MAX_WIDTH_IN:
            print(f"WARNING: panel '{letter}' delta_x is {sizes[letter][0] * 2.54:.2f}cm, exceeding "
                  f"the {MAX_WIDTH_IN * 2.54:.2f}cm Nature two-column guideline by "
                  f"{(sizes[letter][0] - MAX_WIDTH_IN) * 2.54:.2f}cm.")

    if "a" in subpanels:
        build_panel_a(sizes["a"], paddings["a"])

    if "b" in subpanels:
        build_panel_b(sizes["b"], paddings["b"])

    if "c" in subpanels:
        build_panel_c(sizes["c"], paddings["c"], rerun=rerun)

    merge_panels()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rerun", action="store_true", default=False,
                         help="Force panel c's pocket-score data recomputation instead of reusing the cached CSV")
    parser.add_argument("--subpanels", type=str, default=None,
                         help="Comma-separated subset of panels to (re)generate (e.g. 'a' or 'c'), from "
                              f"{{{','.join(PANEL_LETTERS)}}}. Default: all. Fig_1_5_full.pdf is always "
                              "re-merged from whatever Fig_1_5{letter}.pdf files currently exist on disk.")
    args = parser.parse_args()
    if args.subpanels is None:
        subpanels = PANEL_LETTERS
    else:
        subpanels = [p.strip() for p in args.subpanels.split(",")]
        invalid = [p for p in subpanels if p not in PANEL_LETTERS]
        if invalid:
            parser.error(f"Unknown panel(s) {invalid}, must be from {PANEL_LETTERS}")
    main(rerun=args.rerun, subpanels=subpanels)
