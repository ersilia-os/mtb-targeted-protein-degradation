"""
Figure 1.5 (main figure, sitting between Figure 1 and Figure 2): top row is panel a
(square) - the PocketVec t-SNE (3 reference canonical pockets highlighted by gene color,
bold-lettered callout badges pinned to fixed corners) - and panel b (remaining row width,
currently blank - reserved placeholder, per request). The Enamine-library-docking t-SNEs
(HLL, REAL 9.92B) that used to sit alongside PocketVec here have moved out to their own
supplementary figure, FigSupp/tSNE_Enamine.py. Panel c is the 17-row pocket-score strip
figure across all 276 detected pockets (ported from FigSupp/figure_supp_pocket.py - P2Rank
info + 6 physicochemical descriptors ported in from FigSupp/pocket_physchem_molmap.py,
domain annotation, aaRS-class strips; the docking-score summary rows that used to sit
between domain and aaRS-class, and the 2 residue-count-within-radius rows, have both been
removed per request), spanning the full width below. Panels a/b and c are promoted from
supplementary to main-figure status; the 2 original source scripts were retired once this
one was confirmed working.

Built on figure_1_plot.py's panel-per-file architecture: each panel is its own standalone
figure, saved as its own PDF (Fig_1_5a.pdf, Fig_1_5b.pdf, Fig_1_5c.pdf) at an EXACT physical
size read from output/plots/figure_1_5/panel_layout.csv (columns: panel, x, delta_x, y,
delta_y, padding, in cm), then pasted onto one positioned master PDF (Fig_1_5_full.pdf, via
pypdf) and flattened to a PNG (via pymupdf). panel_layout.csv is a first-draft starting
point, meant to be tuned iteratively by rendering and looking, not solved analytically up
front.

Panel c's pocket-score data prep (residue-identity PDB parsing across ~178 structures, for
the physicochemical descriptors) is slow, so it's cached to figure_1_5_pocket_scores.csv
and only recomputed with --rerun - same convention as figure_1_plot.py's own --rerun for
its slow PyMOL rendering step.

Usage:
    python figure_1_5_plot.py [--rerun] [--subpanels a,b,c]
"""
import argparse
import json
import os
import pickle
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
from pocket_group_plots import calculate_density, density_to_sizes, load_canonical_pocket_labels, make_key
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

# Order is A, B, C (ileS, tyrS, gltS) - also drives panel b's top-to-bottom block order.
TARGET_CLUSTERS = ["ileS_cluster1", "gltS_cluster1", "tyrS_cluster1"]
# Mathtext ($\bf{...}$) bolds only the letter, not "Pocket" or the gene name, per request.
# A=ileS, B=gltS, C=tyrS (per explicit correction - not alphabetical by gene name).
CLUSTER_LABELS = {
    "ileS_cluster1": "Pocket $\\bf{A}$\n(ileS)",
    "gltS_cluster1": "Pocket $\\bf{B}$\n(gltS)",
    "tyrS_cluster1": "Pocket $\\bf{C}$\n(tyrS)",
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
    "ileS_cluster1": (0.03, 0.925, "left", "center"),    # A: top-left
    "gltS_cluster1": (0.97, 0.925, "right", "center"),   # B: top-right
    "tyrS_cluster1": (0.03, 0.075, "left", "center"),    # C: bottom-left
}

# HLL docking / REAL 10B docking t-SNEs moved out to FigSupp/tSNE_Enamine.py (per request) -
# panel a is now just the single PocketVec embedding, so it carries the reference-pocket
# callout badges itself instead of only the (now-removed) HLL panel.
EMBEDDING = ("PocketVec descriptors", lambda: compute_pocketvec_embedding(output_dir, RANDOM_SEED))


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
            # CLUSTER_LABELS is the 2-line ("Pocket A" / "(ileS)") version used by panel b's
            # gallery labels; panel a's own corner badges stay single-line, per request.
            badge_text = CLUSTER_LABELS[cluster].replace("\n", " ")
            ax.annotate(badge_text, xy=(cx, cy), xycoords="data",
                        xytext=(lx, ly), textcoords="axes fraction",
                        ha=ha, va=va,
                        fontsize=stylia.FONTSIZE, color="black", zorder=2,
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

# Structure files (for residue-name resolution, see PHYSCHEM_ROWS below) are cached per
# (Uniprot AC, File name) since many pockets share the same structure.
STRUCTURES_DIR = os.path.join(output_dir, "aligned_relaxed_structures")

# 6 physicochemical descriptors ported from FigSupp/pocket_physchem_molmap.py (pocket size +
# 4 composition fractions + average hydropathy + average B-factor/pLDDT; that script's other
# 2 descriptors, P2Rank score/probability, are already their own rows above). Same
# hand-rolled PDB-line parsing as parse_residue_geometry below (not Biopython, to avoid a
# new dependency) - residue *name* is read from the same ATOM/CA line as the coordinates.
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "CYS"}
AROMATIC = {"PHE", "TYR", "TRP"}
POSITIVE = {"LYS", "ARG"}
NEGATIVE = {"ASP", "GLU"}
# Kyte & Doolittle, J Mol Biol 157:105-132 (1982)
KYTE_DOOLITTLE = {
    "ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5, "GLN": -3.5, "GLU": -3.5,
    "GLY": -0.4, "HIS": -3.2, "ILE": 4.5, "LEU": 3.8, "LYS": -3.9, "MET": 1.9, "PHE": 2.8,
    "PRO": -1.6, "SER": -0.8, "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2,
}
PHYSCHEM_ROWS = [
    ("pocket_size", "Pocket size"),
    ("frac_hydrophobic", "Hydroph."),
    ("frac_aromatic", "Aromatic"),
    ("frac_positive", "Positive"),
    ("frac_negative", "Negative"),
    ("avg_hydropathy", "Hydropathy"),
]
# Panel b's physchem gallery shows all 7 non-P2Rank descriptors from pocket_physchem_molmap.py
# (PHYSCHEM_ROWS's 6 plus avg_bfactor, which stays excluded from panel c per earlier request).
PANEL_B_PHYSCHEM_ROWS = PHYSCHEM_ROWS + [("avg_bfactor", "Avg. B-factor/pLDDT")]


def parse_residue_geometry(pdb_file):
    """Residue number -> (x, y, z, three-letter residue name), from each residue's CA atom."""
    res_geom = {}
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                resnum = int(line[22:26])
                resname = line[17:20].strip()
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                res_geom[resnum] = (x, y, z, resname)
    return res_geom


_res_geom_cache = {}


def get_residue_geometry(uniprot_ac, file_name):
    key = (uniprot_ac, file_name)
    if key not in _res_geom_cache:
        _res_geom_cache[key] = parse_residue_geometry(os.path.join(STRUCTURES_DIR, uniprot_ac, file_name))
    return _res_geom_cache[key]


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

    # 6 physicochemical descriptors (see PHYSCHEM_ROWS above) - resolved from the pocket's own
    # residue list (Pocket residues (chain_resn), e.g. "A_196"), keyed by residue number only
    # (matching parse_residue_geometry/get_residue_geometry's own chain-agnostic convention,
    # which already assumes single-chain structures elsewhere in this function). avg_bfactor
    # is a 7th descriptor, computed the same way as pocket_physchem_molmap.py's own
    # compute_descriptors (mean of the "B-factors" column, already aligned 1:1 with "Pocket
    # residues (chain_resn)" tokens in pocket_detection_data.csv - no PDB re-parsing needed).
    # It's kept out of PHYSCHEM_ROWS/panel c (removed there per earlier request) but used by
    # panel b's physchem gallery, via PANEL_B_PHYSCHEM_ROWS below.
    physchem = {name: [] for name, _ in PHYSCHEM_ROWS}
    bfactors_avg = []
    for _, row in pockets.iterrows():
        res_geom = get_residue_geometry(row["Uniprot AC"], row["File name"])
        resnames = []
        for token in row["Pocket residues (chain_resn)"].split():
            resnum = int(token.split("_", 1)[1])
            entry = res_geom.get(resnum)
            if entry is not None and entry[3] in KYTE_DOOLITTLE:
                resnames.append(entry[3])
        n = len(resnames)
        frac = lambda s: sum(r in s for r in resnames) / n
        physchem["pocket_size"].append(n)
        physchem["frac_hydrophobic"].append(frac(HYDROPHOBIC))
        physchem["frac_aromatic"].append(frac(AROMATIC))
        physchem["frac_positive"].append(frac(POSITIVE))
        physchem["frac_negative"].append(frac(NEGATIVE))
        physchem["avg_hydropathy"].append(np.mean([KYTE_DOOLITTLE[r] for r in resnames]))
        bfactors_avg.append(np.mean([float(b) for b in row["B-factors"].split()]))
    for column, _ in PHYSCHEM_ROWS:
        scores[column] = physchem[column]
    scores["avg_bfactor"] = bfactors_avg

    return scores


def prepare_panel_c_data(rerun=False):
    """Loads (or recomputes + caches) pocket_scores, then derives every colormap Normalize
    that depends on its data range - both as module globals, since every plot_*_row function
    below references them by bare name (unchanged from figure_supp_pocket.py)."""
    global pocket_scores
    global P2RANK_NORM, P2RANK_SCORE_NORM, CATALYTIC_CONFIDENCE_NORM, PHYSCHEM_NORMS

    if not rerun and os.path.exists(scores_path):
        pocket_scores = pd.read_csv(scores_path)
        print(f"Loaded cached pocket scores from {scores_path}")
    else:
        pocket_scores = _compute_pocket_scores()
        pocket_scores.to_csv(scores_path, index=False)
        print(f"Saved {len(pocket_scores):,} row(s) to {scores_path}")

    # Every Normalize below is fully data-adaptive (vmin/vmax = this row's own observed
    # min/max, not a fixed theoretical bound like 0-1 for a fraction or 0-4 for
    # catalytic_confidence) - per request, so each row uses its full color range across the
    # pockets actually observed, however narrow that range happens to be.
    P2RANK_NORM = mcolors.Normalize(vmin=pocket_scores["pocket_probability"].min(),
                                     vmax=pocket_scores["pocket_probability"].max(), clip=True)
    P2RANK_SCORE_NORM = mcolors.Normalize(vmin=pocket_scores["pocket_score"].min(),
                                           vmax=pocket_scores["pocket_score"].max(), clip=True)
    CATALYTIC_CONFIDENCE_NORM = mcolors.Normalize(vmin=pocket_scores["catalytic_confidence"].min(),
                                                   vmax=pocket_scores["catalytic_confidence"].max(), clip=True)
    # Built from PANEL_B_PHYSCHEM_ROWS (7) rather than PHYSCHEM_ROWS (6) so avg_bfactor also
    # gets a Normalize entry, for panel b's physchem gallery - panel c only ever looks up the
    # 6 PHYSCHEM_ROWS keys, so this superset dict doesn't change its behavior.
    PHYSCHEM_NORMS = {
        column: mcolors.Normalize(vmin=pocket_scores[column].min(), vmax=pocket_scores[column].max(), clip=True)
        for column, _ in PANEL_B_PHYSCHEM_ROWS
    }


# P2Rank probability/score rows - stylia's "ersilia"-preset "plum" FadingColormap (pale
# lilac -> deep plum), per request.
P2RANK_CMAP = stylia.FadingColormap("plum", transformation=None).cmap
P2RANK_SCORE_CMAP = stylia.FadingColormap("plum", transformation=None).cmap

# Physicochemical rows - orange, anchored to gatA's own canonical gene color (gene_to_color,
# same source as every other gene color in this figure) rather than a stylia preset (none of
# FadingColormap's presets - crimson/cobalt/turquoise/orchid/lime/plum - are orange). Built
# the same way stylia's own presets fade: a pale tint of the hue -> the hue itself (not true
# white), per request.
_gatA_orange = np.array(mcolors.to_rgb(gene_to_color["gatA"]))
_gatA_orange_pale = 0.92 * np.array([1.0, 1.0, 1.0]) + 0.08 * _gatA_orange
SECTION1_CMAP = mcolors.LinearSegmentedColormap.from_list("gatA_orange", [_gatA_orange_pale, _gatA_orange])

# Panel b's PocketVec heatmap - same pale-tint -> hue construction, anchored to pheS's own
# canonical gene color, per request (replacing the earlier reversed stylia "plum" preset).
_pheS_blue = np.array(mcolors.to_rgb(gene_to_color["pheS"]))
_pheS_blue_pale = 0.92 * np.array([1.0, 1.0, 1.0]) + 0.08 * _pheS_blue
POCKETVEC_CMAP = mcolors.LinearSegmentedColormap.from_list("pheS_blue", [_pheS_blue_pale, _pheS_blue])

# Gradient for the continuous catalytic_confidence row - stylia's own native "crimson"
# FadingColormap (same red family as the binarized "Catalytic" domain strip), per request,
# over the pocket's own observed confidence range (CATALYTIC_CONFIDENCE_NORM above).
CATALYTIC_CONFIDENCE_CMAP = stylia.FadingColormap("crimson", transformation=None).cmap

# The 4 domain bands: column in pocket_scores -> band title -> its own "present" color,
# against a shared white "absent" background. Catalytic uses a ligand-evidence confidence
# threshold (CATALYTIC_CONFIDENCE_MIN); the other 3 use plain InterPro label presence,
# mutually exclusive with Catalytic - see CURATED_LABEL_COLUMNS above.
DOMAIN_STRIP_COLUMNS = [
    ("is_catalytic", "Catalytic", "crimson"),
    ("is_trna_binding", "tRNA binding", "cobalt"),
    ("is_editing", "Editing", "amber"),
    ("is_anticodon_binding", "Anticodon", "lime"),
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

# Vertical gap between the 3 row blocks (P2Rank info / domain strips / aaRS-class strips) -
# each block's own rows sit flush against each other (hspace=0, no gap). This is hspace on
# the OUTER 3-row gridspec (not a single flat gridspec across all rows), so the same
# fraction reads as a much bigger absolute gap.
ROWS_HSPACE = 0.2

# Fixed absolute left margin, comfortably wider than the longest row label ("Catalytic
# conf.", ~0.44in at FONTSIZE_SMALL) so the text has visible breathing room rather than
# sitting flush against the page edge. Panel a's own left margin (build_panel_a) is widened
# to match this same absolute page position, rather than shrinking this one to the limit.
ROW_LABEL_LEFT_MARGIN_IN = 1.1

# P2Rank block: pocket_probability + pocket_score + 6 physchem rows + plot_ligand_evidence_row.
N_P2RANK_ROWS = 2 + len(PHYSCHEM_ROWS) + 1
# Domain block: DOMAIN_STRIP_COLUMNS (4) + plot_catalytic_confidence_row (1).
N_DOMAIN_ROWS = len(DOMAIN_STRIP_COLUMNS) + 1


def plot_pocket_scores_row(ax):
    """1st of 6 stacked rows (alongside plot_pocket_score_row and the 4
    plot_domain_strip_row bands) - a gradient strip, not a distribution: one thin cell per
    pocket (same pocket_rank x-order as plot_domain_strip_row), colored by
    P2RANK_CMAP/P2RANK_NORM directly from its own probability value (pale lilac at this
    row's own observed min up to deep plum at its own max, continuous - no cutoff)."""
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
    ax.yaxis.label.set_fontsize(stylia.FONTSIZE_SMALL)
    # Default axes.labelpad (~4pt) would otherwise push the label further left than its own
    # text width accounts for, clipping the longest labels against the page edge given how
    # tight ROW_LABEL_LEFT_MARGIN_IN is tuned.
    ax.yaxis.labelpad = 6


def plot_pocket_score_row(ax):
    """2nd of 6 stacked rows - same gradient-strip treatment as plot_pocket_scores_row, but
    for the raw P2Rank pocket_score instead of pocket_probability (P2RANK_SCORE_CMAP/
    P2RANK_SCORE_NORM: pale lilac to deep plum, over the data's own min-max range, no
    fixed cap)."""
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
    ax.yaxis.label.set_fontsize(stylia.FONTSIZE_SMALL)
    # Default axes.labelpad (~4pt) would otherwise push the label further left than its own
    # text width accounts for, clipping the longest labels against the page edge given how
    # tight ROW_LABEL_LEFT_MARGIN_IN is tuned.
    ax.yaxis.labelpad = 6


def plot_physchem_row(ax, column, label):
    """One of the 6 PHYSCHEM_ROWS - same gradient-strip treatment as plot_pocket_scores_row/
    plot_pocket_score_row, but for a physicochemical descriptor (SECTION1_CMAP/
    PHYSCHEM_NORMS[column]; see PHYSCHEM_NORMS in prepare_panel_c_data for each column's own
    vmin/vmax)."""
    colors = SECTION1_CMAP(PHYSCHEM_NORMS[column](pocket_scores[column]))
    ax.bar(pocket_scores["pocket_rank"], 1, width=DOMAIN_ROW_BAR_WIDTH, color=colors, edgecolor="none", zorder=2)
    ax.set_xlim(pocket_scores["pocket_rank"].min(), pocket_scores["pocket_rank"].max())
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    stylia.label(ax, xlabel="", ylabel=label)
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.label.set_ha("right")
    ax.yaxis.label.set_va("center")
    ax.yaxis.label.set_fontsize(stylia.FONTSIZE_SMALL)
    # Default axes.labelpad (~4pt) would otherwise push the label further left than its own
    # text width accounts for, clipping the longest labels against the page edge given how
    # tight ROW_LABEL_LEFT_MARGIN_IN is tuned.
    ax.yaxis.labelpad = 6


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
    ax.yaxis.label.set_fontsize(stylia.FONTSIZE_SMALL)
    # Default axes.labelpad (~4pt) would otherwise push the label further left than its own
    # text width accounts for, clipping the longest labels against the page edge given how
    # tight ROW_LABEL_LEFT_MARGIN_IN is tuned.
    ax.yaxis.labelpad = 6


def plot_catalytic_confidence_row(ax):
    """Continuous companion to the binarized "Catalytic" domain strip - same gradient-strip
    treatment as plot_pocket_scores_row/plot_pocket_score_row, but for catalytic_confidence
    (0-4) via CATALYTIC_CONFIDENCE_CMAP/CATALYTIC_CONFIDENCE_NORM (pale crimson at this
    row's own observed min up to full crimson at its own max)."""
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
    ax.yaxis.label.set_fontsize(stylia.FONTSIZE_SMALL)
    # Default axes.labelpad (~4pt) would otherwise push the label further left than its own
    # text width accounts for, clipping the longest labels against the page edge given how
    # tight ROW_LABEL_LEFT_MARGIN_IN is tuned.
    ax.yaxis.labelpad = 6


def plot_ligand_evidence_row(ax):
    """Last row of the P2Rank info block - 3-way ligand_evidence (LIGAND_EVIDENCE_COLORS:
    black "direct" / gray "alphafill" / white "none")."""
    colors = [LIGAND_EVIDENCE_COLORS[v] for v in pocket_scores["ligand_evidence"]]
    ax.bar(pocket_scores["pocket_rank"], 1, width=DOMAIN_ROW_BAR_WIDTH, color=colors, edgecolor="none", zorder=2)
    ax.set_xlim(pocket_scores["pocket_rank"].min(), pocket_scores["pocket_rank"].max())
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    stylia.label(ax, xlabel="", ylabel="Lig. evidence")
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.label.set_ha("right")
    ax.yaxis.label.set_va("center")
    ax.yaxis.label.set_fontsize(stylia.FONTSIZE_SMALL)
    # Default axes.labelpad (~4pt) would otherwise push the label further left than its own
    # text width accounts for, clipping the longest labels against the page edge given how
    # tight ROW_LABEL_LEFT_MARGIN_IN is tuned.
    ax.yaxis.labelpad = 6


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
    ax.yaxis.label.set_fontsize(stylia.FONTSIZE_SMALL)
    # Default axes.labelpad (~4pt) would otherwise push the label further left than its own
    # text width accounts for, clipping the longest labels against the page edge given how
    # tight ROW_LABEL_LEFT_MARGIN_IN is tuned.
    ax.yaxis.labelpad = 6


# ===========================================================================
# Panel-saving scaffolding (ported from figure_1_plot.py)
# ===========================================================================

PANEL_LABEL_MARGIN_IN = 0.08  # absolute inches, not a fraction of this panel's own size


def add_panel_label(fig, letter):
    """Bold panel letter at the top-left of the FIGURE (page), fixed regardless of padding.
    Uses a fixed ABSOLUTE margin (PANEL_LABEL_MARGIN_IN), not a fraction of this panel's
    own width/height - a fractional margin gives every panel a different physical inset
    (e.g. panel c, much wider and shorter than a/b, would get a much larger horizontal
    inset and smaller vertical one), which reads as misaligned once panels of different
    sizes are merged onto one page. A fixed absolute inset keeps every panel's letter the
    same physical distance from its own top-left corner, so panels sharing a left edge
    (a/c, both x=0) or a top edge (a/b, both y=0) get letters that line up cleanly."""
    width_in, height_in = fig.get_size_inches()
    x = PANEL_LABEL_MARGIN_IN / width_in
    y = 1 - PANEL_LABEL_MARGIN_IN / height_in
    fig.text(x, y, letter, fontweight="bold",
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
BOTTOM_MARGIN_CM = 0.3  # blank space below the lowest panel, see merge_panels()


def merge_panels():
    """Pastes each Fig_1_5{letter}.pdf onto one positioned master PDF (Fig_1_5_full.pdf) per
    panel_layout.csv's x/y - pure translation (no rescaling, true vector via pypdf), since
    each panel already saved at its exact target size."""
    df = pd.read_csv(panel_layout_path).set_index("panel")
    missing = [p for p in PANEL_LETTERS if p not in df.index]
    if missing:
        raise ValueError(f"{panel_layout_path} is missing row(s) for panel(s): {missing}")

    total_width_cm = max(df.loc[p, "x"] + df.loc[p, "delta_x"] for p in PANEL_LETTERS)
    # BOTTOM_MARGIN_CM padding below the lowest panel - every panel's own y is measured from
    # the top, so simply growing the canvas taller (without touching any panel's own y)
    # pushes it away from the now-taller page's bottom edge, leaving blank space there.
    total_height_cm = max(df.loc[p, "y"] + df.loc[p, "delta_y"] for p in PANEL_LETTERS) + BOTTOM_MARGIN_CM

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
    # right=0.97 (not the symmetric 0.85) shifts the square further right within its own
    # page, closing most of the horizontal gap to panel b - box_aspect(1) still keeps it
    # square, just re-centered within this now heavily off-center margin box (the box is
    # height-constrained by top/bottom, so the box's own left edge only moves a third as
    # fast as this "left" value - box_left = 3*(left+right) - 2.1, in cm). left=0.28 puts
    # the box's left edge at the same absolute page position as panel c's own plot box (see
    # ROW_LABEL_LEFT_MARGIN_IN) - panel a's own left margin is widened to match c, rather
    # than shrinking c's comfortable label margin down to the limit.
    save_panel(fig, "a", use_tight_layout=False,
               subplots_adjust=dict(left=0.28, right=0.97, top=0.85, bottom=0.15), padding=padding)


def build_panel_b(size, padding, rerun=False):
    """PocketVec 128-dim rank-fingerprint heatmap gallery, one thin row per individual
    pocket structure belonging to the 3 reference canonical pockets highlighted in panel a
    (TARGET_CLUSTERS: ileS_cluster1, gltS_cluster1, tyrS_cluster1 - ~30 rows total, grouped
    into 3 blocks in that order: A, B, C). Same probe column order and fixed vmin/vmax=1/131
    (the rank range including dummy/penalty ranks) as FigSupp/pocketvec_molmap.py's own
    gallery_by_canonical_pocket - reuses that script's cached probe_order.pkl directly for
    exact between-figure consistency, rather than recomputing the correlation-clustering
    ordering here. POCKETVEC_CMAP (pheS-anchored, reversed) colors it - low rank/strong hit
    renders as the salient deep color.

    Each block has 4 side-by-side pieces (BLOCK_WIDTH_RATIOS): a text label (CLUSTER_LABELS,
    right-aligned so it sits right up against the color swatch), a solid color swatch
    (cluster_color, its own axes rather than a spine on the heatmap - a distinct visual
    element, not touching the heatmap), the PocketVec heatmap, and - right next to it - all 7
    non-P2Rank physicochemical descriptors from pocket_physchem_molmap.py (PANEL_B_PHYSCHEM_ROWS
    - PHYSCHEM_ROWS's 6 plus avg_bfactor, which panel c excludes per earlier request), one
    column each, for that same pocket - reusing panel c's own PHYSCHEM_NORMS/SECTION1_CMAP so a
    given color means the same thing in both panels. subplots_adjust's own right= (0.95, not
    0.99) reserves a small blank margin at the panel's right edge for other content later -
    chosen so the content's right edge lands at the same absolute page position as panel c's
    own (right=0.9667 there, over c's full 18cm vs. b's 12cm width: 0.9667*18 == 6 + 0.95*12,
    both landing 0.6cm short of the page's own right edge)."""
    prepare_panel_c_data(rerun=rerun)
    physchem = pocket_scores.copy()
    physchem["key"] = [make_key(fn, pn) for fn, pn in zip(physchem["File name"], physchem["Pocket number"])]
    key_to_physchem = physchem.set_index("key")

    fps = pickle.load(open(os.path.join(output_dir, "pocketvec_RUN", "fps_rank.pkl"), "rb"))
    keys = sorted(fps)
    key_to_canonical = dict(zip(keys, load_canonical_pocket_labels(output_dir, keys)))

    order = pickle.load(open(
        os.path.join(output_dir, "plots", "FigSupp", "pocketvec_molmap", "probe_order.pkl"), "rb"))

    cmap = POCKETVEC_CMAP.reversed()
    vmin, vmax = 1, 131

    blocks = [(cluster, sorted(k for k in keys if key_to_canonical[k] == cluster)) for cluster in TARGET_CLUSTERS]
    block_sizes = [len(members) for _, members in blocks]

    fig = plt.figure(figsize=size)
    fig.patch.set_facecolor("white")
    outer_gs = fig.add_gridspec(len(blocks), 1, height_ratios=block_sizes, hspace=0.15)

    # LABEL_W narrowed now that CLUSTER_LABELS are 2 short lines ("Pocket A" / "(ileS)")
    # rather than one long line - the label text is right-aligned within this column, so a
    # too-wide column just leaves dead blank space at the panel's left edge. The reserved
    # right-hand margin lives in subplots_adjust (below), not in these ratios, so it stays in
    # sync with panel c's own margin regardless of block width tweaks here.
    LABEL_W, GAP_W, SWATCH_W, HEATMAP_W, PHYSCHEM_W = 0.09, 0.015, 0.02, 0.65, 0.18
    BLOCK_WIDTH_RATIOS = [LABEL_W, GAP_W, SWATCH_W, GAP_W, HEATMAP_W, GAP_W * 2, PHYSCHEM_W]

    for i, (cluster, members) in enumerate(blocks):
        color = cluster_color(cluster)
        block_gs = outer_gs[i, 0].subgridspec(1, 7, width_ratios=BLOCK_WIDTH_RATIOS, wspace=0)

        label_ax = fig.add_subplot(block_gs[0, 0])
        label_ax.axis("off")
        label_ax.text(1.0, 0.5, CLUSTER_LABELS[cluster], transform=label_ax.transAxes,
                      ha="right", va="center", fontsize=stylia.FONTSIZE_SMALL, color="black")

        swatch_ax = fig.add_subplot(block_gs[0, 2])
        swatch_ax.set_facecolor(color)
        swatch_ax.set_xticks([])
        swatch_ax.set_yticks([])
        for spine in swatch_ax.spines.values():
            spine.set_visible(False)

        heatmap_gs = block_gs[0, 4].subgridspec(len(members), 1, hspace=0)
        physchem_gs = block_gs[0, 6].subgridspec(len(members), 1, hspace=0)
        for j, key in enumerate(members):
            ax = stylize(fig.add_subplot(heatmap_gs[j, 0]))
            row = np.asarray(fps[key])[order].reshape(1, -1)
            ax.imshow(row, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
            ax.set_xticks([])
            ax.set_yticks([])
            title = "PocketVec descriptors" if (i == 0 and j == 0) else ""
            stylia.label(ax, xlabel="", ylabel="", title=title)

            # Each of the 7 columns is normalized (and colored) by its own PHYSCHEM_NORMS
            # entry, so the RGBA image is built cell-by-cell rather than via a single
            # cmap/vmin/vmax imshow call (unlike the PocketVec heatmap, whose 128 columns
            # all share one rank scale).
            pc_row = key_to_physchem.loc[key]
            rgba = np.array([[SECTION1_CMAP(PHYSCHEM_NORMS[column](pc_row[column]))
                               for column, _ in PANEL_B_PHYSCHEM_ROWS]])
            pax = stylize(fig.add_subplot(physchem_gs[j, 0]))
            pax.imshow(rgba, aspect="auto")
            pax.set_xticks([])
            pax.set_yticks([])
            pc_title = "Physchem. properties" if (i == 0 and j == 0) else ""
            stylia.label(pax, xlabel="", ylabel="", title=pc_title)

    # top=0.85/bottom=0.15 matches panel a's own subplots_adjust exactly (both panels are
    # the same 6cm height), so panel b's gallery content is vertically aligned with panel
    # a's actual plot box, not just sharing the same overall page height.
    save_panel(fig, "b", use_tight_layout=False,
               subplots_adjust=dict(left=0.28, right=0.95, top=0.85, bottom=0.15), padding=padding)


def build_panel_c(size, padding, rerun=False):
    prepare_panel_c_data(rerun=rerun)

    fig = plt.figure(figsize=size)
    fig.patch.set_facecolor("white")
    # 3 blocks, in order: P2Rank info, domain strips, aaRS-class strips - each with its own
    # inner subgridspec at hspace=0 (no gap between rows within a block), while ROWS_HSPACE
    # is kept only between the blocks themselves. Docking-score rows removed per request.
    block_sizes = [N_P2RANK_ROWS, N_DOMAIN_ROWS, len(AARS_CLASSES)]
    outer_gs = fig.add_gridspec(3, 1, height_ratios=block_sizes, hspace=ROWS_HSPACE)

    p2rank_gs = outer_gs[0, 0].subgridspec(block_sizes[0], 1, hspace=0)
    plot_pocket_scores_row(stylize(fig.add_subplot(p2rank_gs[0, 0])))
    plot_pocket_score_row(stylize(fig.add_subplot(p2rank_gs[1, 0])))
    row_idx = 2
    for column, label in PHYSCHEM_ROWS:
        plot_physchem_row(stylize(fig.add_subplot(p2rank_gs[row_idx, 0])), column, label)
        row_idx += 1
    plot_ligand_evidence_row(stylize(fig.add_subplot(p2rank_gs[row_idx, 0])))

    domain_gs = outer_gs[1, 0].subgridspec(block_sizes[1], 1, hspace=0)
    row_idx = 0
    for column, title, color_name in DOMAIN_STRIP_COLUMNS:
        ax = stylize(fig.add_subplot(domain_gs[row_idx, 0]))
        plot_domain_strip_row(ax, column=column, title=title, color_name=color_name)
        row_idx += 1
        if column == "is_catalytic":
            plot_catalytic_confidence_row(stylize(fig.add_subplot(domain_gs[row_idx, 0])))
            row_idx += 1

    class_gs = outer_gs[2, 0].subgridspec(block_sizes[2], 1, hspace=0)
    for i, aars_class in enumerate(AARS_CLASSES):
        ax = stylize(fig.add_subplot(class_gs[i, 0]))
        plot_class_strip_row(ax, aars_class)

    # right=0.90 (not 0.995) reserves blank space at the panel's right edge for other content
    # later, matching panel b's own reserved margin.
    save_panel(fig, "c", use_tight_layout=False,
               subplots_adjust=dict(left=0.28, right=0.9667, top=0.88, bottom=0.12),
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
        build_panel_b(sizes["b"], paddings["b"], rerun=rerun)

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
