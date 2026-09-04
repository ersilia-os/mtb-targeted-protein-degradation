"""
Figure 4, panel a: PyMOL screenshot of PDB 7K98, the crystal structure of M. tuberculosis PheRS
(pheS/pheT, the project's Tier-1 target) in its "preaminoacylation" complex with cognate tRNA(Phe)
and the phenylalanyl-adenylate mimic F-AMS (ligand code W5Y). This structure is already used
downstream in the pipeline as 7K98_pocket_6 ("pocket at dimerization interface",
output/selected_pockets.csv, scripts 48-50/65/71-75) - panel 4a introduces the real experimental
structure behind that pocket.

Source: data/structures/pdbe_database/P9WFU3/P9WFU3_updated/7k98_updated.cif. Citation (from the
deposited file's own JRNL header, not invented): Michalska K, Jedrzejczak R, Wower J, Chang C,
Baragana B, Gilbert IH, Forte B, Joachimiak A. "Mycobacterium tuberculosis Phe-tRNA synthetase:
structural insights into tRNA recognition and aminoacylation." Nucleic Acids Res. 2021;49(9):5351.
DOI: 10.1093/nar/gkab272. Deposited 2020-09-28, resolution 2.19 A.

Shows the full deposited complex: PheRS alpha subunit (pheS, chains A/D), PheRS beta subunit
(pheT, chains B/E, one copy - chain B - also shown as an opaque surface), tRNA(Phe) (chains C/F),
and the native co-crystallized ligand W5Y (F-AMS, fuchsia sticks - not a docked pose, so
intentionally NOT the project's usual orange COLOR_LIGAND_DOCKED). Solvent (waters) and metallic
atoms (the bound Mg2+ ions), along with Cl-/glycerol, are not shown. Every custom color (pheS
cobalt, pheT amber, tRNA turquoise, ligand fuchsia) is sourced from stylia.ArticleColors() rather
than a native PyMOL color name or a hand-picked RGB, so the palette stays tied to the same source
as every other stylia-driven plot in the project - and pheS/pheT get distinct colors rather than
the project's own near-identical gene_to_color hexes for this pair (adjacent on the alphabetical
spectral palette), since distinguishing the two subunits is this panel's whole point.

Also shown: the 3 curated pheS/pheT pockets from output/selected_pockets.csv, as large spheres
(add_pocket_spheres()), one distinct article-palette color per pocket (POCKET_COLORS). That file
actually lists 4 pheS/pheT rows, but
alphafold3_P9WFU3_model_2_pocket_2 ("falls at interface") and 7K98_pocket_6 ("pocket at
dimerization interface") describe the same physical pocket, independently detected on an AF3
model and on this real structure - shown once, via 7K98_pocket_6's own native centroid (already
in this structure's own coordinate frame, unlike the other two pockets, which are P2Rank pockets
detected on separate monomeric models and need aligning onto "structure" first - see
_add_aligned_pocket_sphere()). User-confirmed dedup.

Panel size comes from output/plots/figure_4/panel_layout.csv (columns: panel, x, delta_x, y,
delta_y, padding - x/y are each panel's top-left placement, used by merge_panels() (called
automatically at the end of main(); formerly a separate figure_4_merge.py script, now folded into
this same file) to build one master PDF from the standalone per-panel PDFs; delta_x/delta_y (cm)
are each panel's exact saved page size; padding (0-1, see apply_padding()) is a per-dimension
linear shrink factor toward the panel's own center, adding deliberate blank margin around a
panel's content without affecting any font/tick/label size), same convention as figure_3_plot.py's
own per-panel-standalone-PDF setup (load_panel_padding/apply_padding ported near-verbatim from
there).

Panel b: 4x3 grid, one snapshot per output/selected_pockets.csv pocket (structure surface + best
docked ligand pose, restricted to script 62's ~2,923-compound aggregated/prioritized pool -
AGGREGATED_HITS_CSV, characterized with Ersilia ADMET models in script 67 - rather than each
pocket's full raw report.csv, so the showcased hit is one that's actually part of the pipeline's
own prioritized compound set, not just whatever scored lowest in isolation for that one pocket
(user-confirmed). See plot_pocket_grid_panel/render_pocket_snapshot for the full recipe and its
own POCKET_VIEWS per-pocket hand-tuned camera views). Still being tuned one pocket at a time.

Panel c: 2 protein-level (pheST/aspS/lysS/alaS) UpSet plots built directly from script 65's
aggregated docking scores (output/65_aggregated_docking/docking_results) - CAT pockets at docking
score <= CAT_SCORE_THRESHOLD (-11) and NON-CAT pockets at docking score <= NONCAT_SCORE_THRESHOLD
(-10). Computed here rather than reusing script 68's own upset_score_*.png renders (script 68 is a
numbered pipeline script and can't be imported anyway - its filename isn't a valid module
identifier) for two reasons: those renders are for script 68's own broader threshold sweep, and -
the actual reason this panel doesn't just composite them like panel a/b's PyMOL PNGs - each one is
independently auto-scaled by matplotlib/upsetplot to its own data, so a bar of the same height in
the two doesn't represent the same intersection size. _render_upset_figure() builds each as its
own standalone plt.figure() (upsetplot manages its own internal axes layout, so the two can't be
laid out jointly in one shared figure/subfigure - confirmed upset.plot(fig=...) requires a real
Figure, not a matplotlib SubFigure), then plot_upset_panel() reads back both figures'
"intersections" (vertical bars) and "totals" (horizontal per-protein bars) axis limits, takes the
max across the two, and forces both to that shared scale before rasterizing - only then are they
composited side by side the same render-to-PNG-then-imshow way as panel a/b.

Panels d and e (interchanged with the old panel c - user request, same content/code, only the
letter/position changed - see main(); later split from one combined panel into two standalone
ones, also per request - plot_sibling_robustness_panel/plot_affinity_panel): two structure-choice
robustness checks on the docking scores, both boxplots binned by docking score, whiskers at
p1/p99 (whis=(1, 99), not matplotlib's default 1.5-IQR rule - user request; min/max
(whis=(0, 100)) was tried too and reverted - p1/p99 preferred despite some very small bins (e.g.
n=4) collapsing their whiskers to zero length when the true min/max falls outside the 1st/99th
percentile threshold) with fliers/caps hidden (showfliers=False/showcaps=False - user request)
since they'd otherwise mark everything outside p1/p99 as an "outlier" dot, per-bin n printed to
console each run.

Panel d: does a compound's docking score in one structure of a pocket predict its score in that
same pocket's SIBLING structure(s) - independently-detected structures P2Rank found to occupy the
same physical site (output/77_pocket_annotation/pocket_clusters.csv's 6.14 A centroid dedup, see
scripts/77_pocket_annotation/09_cluster_pockets.py)? load_sibling_avg_scores(): for the 10 curated
pockets with at least one sibling (7K98's dimer pocket and pheT's own singleton-cluster pocket are
excluded - nothing to compare them against), single-run score (report_old.csv, pre-replicate),
boxplotted (SIBLING_AVG_BIN_EDGES + an open-ended -inf floor prepended at plot time - user
request: the leftmost "-13" bin catches everything more negative than -13, not just a fixed
window - only the loose/upper end at -8 is a real cutoff) against that SAME compound's average
score across the pocket's sibling structure(s) - sibling scores pulled from report_old.csv when
the sibling is itself curated, else from the round-2 library. ~43% of pairs get dropped because
the compound was never docked against a given sibling at all (not an invalid-pose issue -
confirmed via direct S3 size-check that local data isn't just unsynced, and traced to compounds
sourced from scripts 56/61's NON-CAT top100 selections, which only ever docked their picks against
the curated pocket, never broadly) - a real compound-library coverage gap, not a bug.

Panel e: Boltz-2-predicted IC50 (output/75_boltz2_collect_affinities, converted from its
log10(IC50 uM) affinity_pred_value), binned by docking score (output/70_filtering's
FILTERED_HITS_CSV), over the ~1,095 FILTERED hits x 12 pockets (load_docking_vs_affinity(),
13,140 pairs). A direct scatter of docking score vs. IC50 was tried first and rejected (too noisy
to read at this population size); a boxplot per BOXPLOT_BIN_EDGES 1-unit docking-score bin
instead, plus an open-ended -inf floor prepended at plot time (user request: the leftmost "-13"
bin catches everything more negative than -13, including the single pair below -14, not just a
fixed window - same bin structure as panel d's own SIBLING_AVG_BIN_EDGES) up to [-9,-8)
(stopping at the loosest bin the trend still rises through, since predicted IC50 flattens/reverses
past it) shows the same trend - does a stronger docking score correspond to a better (lower)
predicted IC50?

(An earlier version of this panel paired each bin's Boltz-2 box against a "Nesso" placeholder box
that just duplicated the Boltz-2 values verbatim, pending real Nesso-1 predictions. Real Nesso-1
results now exist, but were found to sit on a systematically different absolute IC50 scale from
Boltz-2's - not validated as cross-model comparable - so Nesso-1 was dropped from this panel
rather than plotted on the same absolute nM axis.)

Panel f: a PANEL_F_N_ROWS x PANEL_F_N_COLS (2x3, 6 total) grid of compound-structure slots, for a user-curated shortlist
(output/plots/figure_4/compounds/ids.csv, column cpd_id holding selection_table.csv source_key
values - not cpd_id values, despite the column's name - row order = slot order). Each slot is the
compound's 2D structure (RDKit, rdCoordGen layout, from
output/70_filtering/filtered_hits_explorer/selection_table.csv - the molecule-auditing skill's
own curated table) via OffsetImage/AnnotationBbox at a uniform zoom, inside a square (not
rounded) border, plus a small bordered 4-square grid of real per-protein data (pheST, aspS,
lysS, alaS) in each of the 4 corners (see CORNER_GRID_SPECS): top corners are catalytic/
non-catalytic docking score, bottom corners are catalytic/non-catalytic predicted affinity
(nM) - color intensity encodes strength (white = weak, saturated = strong), same gradient
convention as output/70_filtering/filtered_hits_explorer/viz_meta.json's own explorer. Any
slot beyond ids.csv's own row count is left empty, border still drawn.

Usage:
    python figure_4_plot.py [--rerun] [--subpanels a,b,c,d,e]
"""
import argparse
import io
import json
import os
import re
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

os.environ["QT_QPA_PLATFORM"] = "offscreen"

import tarfile
import tempfile

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymol
import stylia
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib.patches import FancyBboxPatch, Rectangle
from pymol import cmd
from pypdf import PdfReader, PdfWriter, Transformation
import pymupdf
from rdkit import Chem
from rdkit.Chem import rdCoordGen
from rdkit.Chem.Draw import rdMolDraw2D
from scipy.ndimage import gaussian_filter
from scipy.stats import pearsonr, rankdata, spearmanr
from stylia.config import get_fg_color
from stylia.figure.figure import stylize

from docking_utils import LIBRARIES, load_scores

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

ROOT = os.path.join(root, "..", "..")
plots_dir = os.path.join(ROOT, "output", "plots", "figure_4")
renders_dir = os.path.join(plots_dir, "renders")
os.makedirs(plots_dir, exist_ok=True)
os.makedirs(renders_dir, exist_ok=True)

PANEL_LETTERS = ["a", "b", "c", "d", "e", "f"]
panel_layout_path = os.path.join(plots_dir, "panel_layout.csv")

# Panels d/e (DOCKING_RESULTS_DIR_65 is also reused by panel c's _load_agg_scores)
DOCKING_RESULTS_DIR_65 = os.path.join(ROOT, "output", "65_aggregated_docking", "docking_results")
FILTERED_HITS_CSV = os.path.join(ROOT, "output", "70_filtering", "filtered_hits.csv")
AFFINITY_RESULTS_CSV = os.path.join(ROOT, "output", "75_boltz2_collect_affinities", "affinity_results.csv")
# Right subplot's boxplot bins (user-specified range) - 1-unit docking-score bins from [-15,-14)
# to [-9,-8) (the loosest bin the predicted-IC50 trend actually still rises through - it
# flattens/reverses past this point, per the percentile tables checked by hand before building
# this panel), pd.cut(..., right=False).
BOXPLOT_BIN_EDGES = list(range(-13, -7))

# Left subplot - sibling-structure comparison (see load_sibling_avg_scores). Pockets with no
# monomer-cluster siblings to compare against: 7K98 (dimer, own coordinate frame, never entered
# the 6.14 A centroid dedup below) and pheT's own curated pocket (singleton cluster,
# n_models_in_cluster=1) - both user-confirmed exclusions.
POCKET_CLUSTERS_CSV = os.path.join(ROOT, "output", "77_pocket_annotation", "pocket_clusters.csv")
# Round-2 REAL library (single-run) - the source for a sibling's score whenever that sibling
# isn't itself one of the 12 curated pockets (i.e. no report_old.csv of its own).
DOCKING_RESULTS_DIR_REAL2 = os.path.join(ROOT, "output", "unidock_REAL_docking_2", "docking_results")
SIBLING_EXCLUDED_POCKETS = {"7K98_pocket_6", "alphafold2_P9WFU1_model_0_pocket_1"}
# 1-unit docking-score bins over the user-specified [-14,-8) range.
SIBLING_AVG_BIN_EDGES = list(range(-13, -7))

# Panel b - the ~2,923 compounds prioritized by script 62's aggregation across the 5 upstream
# CAT/NON-CAT promiscuous/selective/top100 selections (scripts 52-61), later characterized with
# Ersilia ADMET models in script 67. Panel b's "best compound per pocket" is restricted to this
# set (see render_pocket_snapshot) rather than each pocket's full raw report.csv, so the
# showcased hit is one that's actually part of the pipeline's own prioritized/characterized
# pool - not just whatever scored lowest in isolation for that one pocket (user-confirmed).
AGGREGATED_HITS_CSV = os.path.join(ROOT, "output", "62_aggregate_hits", "aggregated_hits.csv")

# Panel c - built from script 65's aggregated docking scores directly (see module docstring),
# not from script 68's own renders (those auto-scale independently per plot and aren't comparable).
CAT_SCORE_THRESHOLD = -10
NONCAT_SCORE_THRESHOLD = -10
# upsetplot's own UpSet(..., element_size=...) kwarg (points/row-unit, default 32) - controls how
# large the chart geometry (bars, matrix rows) renders, independent of text size (see
# _render_upset_figure). Tuned down from the default so the panel can be small without the text
# shrinking along with it - retune alongside panel_layout.csv's "c" row if that row changes a lot.
UPSET_ELEMENT_SIZE = 12
# Blank gap (inches) between the CAT and NON-CAT images in plot_upset_panel - flush (0) read as
# too cramped (user-flagged); taken out of the shared display height (see disp_h_in) rather than
# added on top of it, so the pair still fits within panel_layout.csv's own "d" row width exactly.
UPSET_GAP_IN = 0.06

# Panel f - user-curated compound-structure row (see module docstring). ids.csv doesn't ship
# with the repo (it's a hand-picked shortlist, dropped in by the user); selection_table.csv is
# the molecule-auditing skill's own output, reused as-is (not regenerated here) purely for its
# cpd_id/source_key/smiles columns.
PANEL_F_CPD_IDS_CSV = os.path.join(plots_dir, "compounds", "ids.csv")
SELECTION_TABLE_CSV = os.path.join(ROOT, "output", "70_filtering", "filtered_hits_explorer",
                                    "selection_table.csv")
# Fixed slot count (user-confirmed) - independent of how many rows PANEL_F_CPD_IDS_CSV
# currently has; slots beyond the file's own row count are left empty rather than shrinking
# the grid, so the panel's layout doesn't reflow every time a compound is added.
PANEL_F_N = 6
# 2x3 grid (user-confirmed, cut down from 3x3/9 - fewer, larger compound cards) - 3 columns
# because panel_layout.csv's own "e" row is narrow (aligned under columns c/d), too narrow for
# PANEL_F_N slots side by side in a single row without overlapping.
PANEL_F_N_ROWS = 2
PANEL_F_N_COLS = 3
# Same wspace/hspace passed to add_gridspec() in plot_compound_cards_panel, and the same margins
# save_panel's own default subplots_adjust applies - factored out here because _mol_image_zoom()
# needs them to compute each slot's actual physical size analytically (see that function).
GRIDSPEC_SPACE = 0.15
GRIDSPEC_MARGIN = 0.01
# Mol image sizing (see _mol_image_zoom): fraction of the tighter slot dimension the image
# should fill - panel_layout.csv's panel-e row size is still being tuned by hand, so the zoom is
# derived from the actual slot size at render time (below) rather than a fixed value tuned for
# one specific panel size, which broke every time that row was resized.
MOL_IMAGE_FILL_FRAC = 0.9  # bigger (user request) - was 0.85, then 0.7 while the panel still had
# center grids stealing vertical space; now that those are gone (replaced by small corner grids,
# see CORNER_GRID_N below) the molecule can take up most of the card again.
# Border inset (axes-fraction units) - keeps the border from cutting into the molecule at the
# corners, and doubles as the corner grids' own inset from the card edge below; purely cosmetic,
# doesn't affect mol image sizing.
CARD_MARGIN = 0.04
# Thinner than stylia.LINEWIDTH (0.5) - the full-weight line read too heavy against these small
# slots (user-confirmed).
CARD_BORDER_LINEWIDTH = 0.3
# One CORNER_GRID_N-wide horizontal row of bordered square cells in each of the card's 4 corners
# (user request - was a vertical CORNER_GRID_N-tall stack, switched to horizontal + shrunk, see
# CORNER_CELL_FRAC) - replaces the earlier version's pair of wide center grids under the
# molecule. Each corner shows real per-protein data (user-confirmed, superseding this panel's
# earlier random-placeholder grids): one cell per CORNER_GRID_PROTEINS entry, colored via the
# same white(low)->color(high) linear gradient (_gradient_hex_rgb) and (white_at, red_at, hex)
# values as output/70_filtering/filtered_hits_explorer/viz_meta.json's own "value_gradients" - so
# these corners use the exact same color convention as that explorer for the same columns, at the
# same red/green/blue/orange 2x2 layout the user specified.
CORNER_GRID_N = 4
# One cell's edge length, as a fraction of the slot's own WIDTH - a plain x-fraction and
# y-fraction of the same value are only equal physical distances when the slot itself is square
# (true back when PANEL_F_N_ROWS == PANEL_F_N_COLS at 3x3; no longer true at 2x3, where slots are
# wider than tall). plot_compound_cards_panel converts this into a separate y-fraction, scaled by
# the slot's own width/height ratio (via _slot_size_in), so each corner cell renders as a true
# square regardless of the slot's aspect ratio. Shrunk twice on user request ("slightly
# smaller"): 0.09 -> 0.075 -> 0.065.
CORNER_CELL_FRAC = 0.065
# Cell order (left to right) within each corner's 4-cell row - same order as
# selection_table.csv's own columns.
CORNER_GRID_PROTEINS = ["pheST", "aspS", "lysS", "alaS"]
# corner -> (selection_table.csv column suffix, white_at, red_at, hex_color), one quadrant per
# metric: CAT/NON-CAT docking score (kcal/mol) and CAT/NON-CAT predicted affinity (nM).
CORNER_GRID_SPECS = {
    "top_left": ("cat", -10, -14, "#e63946"),          # red - CAT docking score
    "top_right": ("noncat", -8, -12, "#6bbf59"),        # green - NON-CAT docking score
    "bottom_left": ("cat, nM", 500, 50, "#0969da"),     # blue - CAT predicted affinity (nM)
    "bottom_right": ("noncat, nM", 1000, 100, "#f4845f"),  # orange - NON-CAT predicted affinity (nM)
}

SOURCE_CIF = os.path.join(ROOT, "data", "structures", "pdbe_database", "P9WFU3", "P9WFU3_updated",
                           "7k98_updated.cif")

# Every custom color below is sourced from stylia's own article palette (not a native PyMOL
# color name or a raw hand-picked RGB), so this render's palette stays tied to the same source
# as every stylia-driven plot in the project.
ac = stylia.ArticleColors()
# Same neutral structure-surface color as panel b's own COLOR_STRUCTURE_NEUTRAL below (moved up
# here so panel a can use it too) - project convention for de-emphasized structural context,
# reserving saturated color for whatever a panel actually wants to draw the eye to.
COLOR_STRUCTURE_NEUTRAL = [0.7804, 0.8275, 0.8667]
COLOR_PHES = ac.cobalt              # alpha subunit (chain A/D) - the catalytic chain this panel
                                     # is actually about, so it keeps a real color; pheT recedes
                                     # into COLOR_STRUCTURE_NEUTRAL below (user-confirmed: the
                                     # previous amber cartoon+surface combo "does not look good
                                     # enough to be published").
COLOR_PHET = COLOR_STRUCTURE_NEUTRAL  # beta subunit (chain B/E), cartoon AND surface both -
                                       # was ac.amber (yellow), replaced for a cleaner, more
                                       # muted look consistent with panel b's own neutral surfaces.
COLOR_TRNA = ac.turquoise   # replaces the native PyMOL "palecyan" used elsewhere in the project
                             # (scripts/77_pocket_annotation/11_build_pymol_sessions.py's
                             # COLOR_EXPERIMENTAL_TRNA) with the closest article-palette color.
# Was ac.fuchsia (closest article-palette match to the project's usual "native/experimental
# ligand" magenta, scripts/47b_reference_pocket_visualization.py's COLOR_LIGAND_PDB) - changed to
# periwinkle (user-confirmed) because this ligand (W5Y) sits directly inside the crimson CAT-
# pocket sphere (POCKET_COLORS), and fuchsia/crimson are both red-family hues that visually
# blended under the sphere's transparency, making the ligand's own sticks hard to distinguish
# from the sphere. Breaks the cross-project magenta-for-native-ligand convention for this one
# panel specifically, in exchange for actually being visible here.
COLOR_LIGAND = ac.periwinkle
COLOR_MG = "green"
# One distinct article-palette color per pocket (none reused from COLOR_PHES/PHET/TRNA/LIGAND
# above), keyed by pocket_name.
POCKET_COLORS = {
    "swissmodel_P9WFU3_model_0_pocket_1": ac.crimson,      # pheS, CAT
    "7K98_pocket_6": ac.lime,                              # pheS, NON-CAT (dimerization interface)
    "alphafold2_P9WFU1_model_0_pocket_1": ac.orchid,       # pheT, NON-CAT
}

# output/selected_pockets.csv's 3 curated pheS/pheT pockets (see module docstring dedup note),
# rendered as large spheres. 7K98_pocket_6 is already native to this structure's own coordinate
# frame (output/48_detect_pocket_multimers/); the other two are P2Rank pockets detected on
# separate monomeric models (output/pocket_detection_data.csv) and need aligning onto "structure"
# first - offset=0 for both P9WFU3 and P9WFU1 vs 7K98 in
# processed/pocket_annotation/chain_manifest.json, so pocket residue numbers correspond directly,
# no renumbering needed (unlike some other structure pairs in that manifest).
SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")
POCKET_DATA_CSV = os.path.join(ROOT, "output", "pocket_detection_data.csv")
MULTIMER_POCKET_DATA_CSV = os.path.join(ROOT, "output", "48_detect_pocket_multimers",
                                         "pocket_detection_data.csv")
ALIGNED_DIR = os.path.join(ROOT, "output", "aligned_relaxed_structures")
DETECTED_POCKETS_DIR = os.path.join(ROOT, "output", "detected_pockets")
# Extended beyond pheS/pheT (panel a's own pockets) to all 5 genes in output/selected_pockets.csv,
# for panel b's pocket grid below.
GENE_TO_UNIPROT = {"pheS": "P9WFU3", "pheT": "P9WFU1", "aspS": "P9WFW3", "lysS": "P9WFU9", "alaS": "P9WFW7"}
GENE_TO_TARGET_CHAIN = {"pheS": "A", "pheT": "B"}  # which "structure" chain to align onto
# alphafold3_P9WFU3_model_2_pocket_2 ("falls at interface") and 7K98_pocket_6 ("pocket at
# dimerization interface") describe the same physical pocket, independently detected on an AF3
# model and on this real structure - shown once, via 7K98_pocket_6's own native coordinates
# (user-confirmed dedup).
SKIP_POCKET_NAMES = {"alphafold3_P9WFU3_model_2_pocket_2"}

POCKET_SPHERE_RADIUS = 4.0
# Slightly transparent - fully opaque (0.0) read as too flat/solid; this keeps the sphere
# reading as a clean marker without completely hiding what's behind it.
POCKET_SPHERE_TRANSPARENCY = 0.15


def load_panel_sizes(panels):
    """panel -> (width_in, height_in), from output/plots/figure_4/panel_layout.csv - same
    convention as figure_3_plot.py's own load_panel_sizes."""
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
    """panel -> padding fraction (0-1) from panel_layout.csv's padding column - see
    apply_padding() for what this controls. Missing column defaults every panel to 0 (today's
    already-tuned zero-margin behavior, unchanged) - ported near-verbatim from
    figure_3_plot.py's own load_panel_padding."""
    df = pd.read_csv(panel_layout_path).set_index("panel")
    if "padding" not in df.columns:
        return {p: 0.0 for p in panels}
    missing = [p for p in panels if p not in df.index]
    if missing:
        raise ValueError(f"{panel_layout_path} is missing row(s) for panel(s): {missing}")
    return {p: df.loc[p, "padding"] for p in panels}


def apply_padding(fig, padding):
    """Scales every axes in `fig` toward the figure's own center by (1-padding) - see
    panel_layout.csv's padding column. Pure geometric Bbox transform, applied AFTER each panel's
    own tight-fit/clip-avoidance logic has already run (a no-op at padding=0, so it only ADDS
    margin on top of an already clip-safe baseline, never re-introduces clipping). Font/tick/
    label sizes are unaffected - matplotlib renders those at fixed absolute point sizes regardless
    of an axes' box size, so this only moves/shrinks the plotted geometry, not the text. Ported
    near-verbatim from figure_3_plot.py's own apply_padding."""
    scale = 1 - padding
    for ax in fig.axes:
        pos = ax.get_position()
        x0 = 0.5 + (pos.x0 - 0.5) * scale
        x1 = 0.5 + (pos.x1 - 0.5) * scale
        y0 = 0.5 + (pos.y0 - 0.5) * scale
        y1 = 0.5 + (pos.y1 - 0.5) * scale
        ax.set_position([x0, y0, x1 - x0, y1 - y0])


PANEL_LABEL_MARGIN = 0.02  # figure-fraction inset from the page corner - ported from figure_3_plot.py


def add_panel_label(fig, letter, x_margin=PANEL_LABEL_MARGIN):
    """Bold panel letter at the top-left of the FIGURE (page) itself, not the axes - stays at a
    fixed page corner regardless of panel_layout.csv's padding. Ported near-verbatim from
    figure_3_plot.py's own add_panel_label (already reused as-is by figure_1_plot.py and the
    former Figure 2, since retired; figure_4 was the only one of the four missing it). x_margin overrides the
    default inset - panels c/d/e use a smaller one (user request: "a bit to the left"), a/b keep
    the default."""
    fig.text(x_margin, 1 - PANEL_LABEL_MARGIN, letter, fontweight="bold",
             fontsize=stylia.FONTSIZE_BIG, color=get_fg_color(), ha="left", va="top",
             transform=fig.transFigure)


def save_panel(fig, letter, use_tight_layout=True, tight_pad=1.08, tight_w_pad=None,
               subplots_adjust=None, padding=0.0, letter_margin=PANEL_LABEL_MARGIN):
    # No bbox_inches="tight" - saves at exactly panel_layout.csv's delta_x/delta_y, same
    # convention as figure_3_plot.py's save_panel. use_tight_layout=False (panel b's grid, and
    # panels d/e below) skips plt.tight_layout() - it silently fails ("Axes not compatible with
    # tight_layout" warning; confirmed on panel e's log-scale axis too, not just panel
    # b's bbox_to_anchor corner badges) and produces IDENTICAL output regardless of tight_pad,
    # clipping axis-label text at the saved edge since nothing here uses bbox_inches="tight" to
    # grow the canvas to fit. subplots_adjust (dict of left/right/top/bottom/wspace fractions)
    # replaces that failed auto-layout with explicit, hand-tuned margins; falls back to the old
    # near-zero-margin default (fine for panel b's label-free grid) when not given.
    # stylia sets rcParams["figure.autolayout"]=True, which silently re-runs an automatic
    # tight_layout() (with default padding) at draw/savefig time on top of whatever layout is
    # set below - overriding subplots_adjust's hand-tuned near-zero margins with much larger
    # ones (confirmed: ~1% intended vs. ~6% actual on every side). Disabling it here restores
    # subplots_adjust/tight_layout(pad=...) as the actual, final layout.
    fig.set_tight_layout(False)
    if use_tight_layout:
        plt.tight_layout(pad=tight_pad, w_pad=tight_w_pad)
    else:
        fig.subplots_adjust(**(subplots_adjust or dict(left=0.01, right=0.99, top=0.99, bottom=0.01)))
    apply_padding(fig, padding)  # after tight_layout/subplots_adjust, so it isn't undone by them
    add_panel_label(fig, letter, x_margin=letter_margin)
    output_path = os.path.join(plots_dir, f"Fig_4{letter}.pdf")
    plt.savefig(output_path, dpi=600, transparent=False)
    plt.close(fig)
    print(f"Saved Fig_4{letter}.pdf")


CM_TO_PT = 72 / 2.54


def merge_panels():
    """Pastes each Fig_4{letter}.pdf onto one positioned master PDF (Fig_4_full.pdf) per
    panel_layout.csv's x/y - pure translation (no rescaling, true vector via pypdf), since
    each panel already saves at its exact target size. Folded in from the former
    figure_4_merge.py (now removed) so plot + merge happen in one script, matching
    figure_1_plot_v2.py/figure_2_plot_v2.py/figure_3_plot.py's own convention."""
    df = pd.read_csv(panel_layout_path).set_index("panel")
    missing = [p for p in PANEL_LETTERS if p not in df.index]
    if missing:
        raise ValueError(f"{panel_layout_path} is missing row(s) for panel(s): {missing}")

    total_width_cm = max(df.loc[p, "x"] + df.loc[p, "delta_x"] for p in PANEL_LETTERS)
    total_height_cm = max(df.loc[p, "y"] + df.loc[p, "delta_y"] for p in PANEL_LETTERS)

    writer = PdfWriter()
    master_page = writer.add_blank_page(width=total_width_cm * CM_TO_PT, height=total_height_cm * CM_TO_PT)

    for p in PANEL_LETTERS:
        panel_path = os.path.join(plots_dir, f"Fig_4{p}.pdf")
        if not os.path.exists(panel_path):
            print(f"  Skipping {p}: {panel_path} not found - Fig_4_full.pdf will be missing it.")
            continue
        panel_page = PdfReader(panel_path).pages[0]
        x_top_cm, y_top_cm, delta_y_cm = df.loc[p, "x"], df.loc[p, "y"], df.loc[p, "delta_y"]
        x_pt = x_top_cm * CM_TO_PT
        y_bottom_pt = (total_height_cm - y_top_cm - delta_y_cm) * CM_TO_PT
        master_page.merge_transformed_page(panel_page, Transformation().translate(tx=x_pt, ty=y_bottom_pt))
        print(f"  Placed Fig_4{p}.pdf at x={x_top_cm}cm, y={y_top_cm}cm (top-left)")

    output_path = os.path.join(plots_dir, "Fig_4_full.pdf")
    with open(output_path, "wb") as f:
        writer.write(f)
    print(f"Saved merged master figure ({total_width_cm:.2f} x {total_height_cm:.2f} cm) to {output_path}")

    # Flattened PNG alongside the vector PDF (user request), rendered at the same dpi=600
    # save_panel's own PDFs already use for their embedded raster content - pymupdf renders the
    # already-merged page directly (no external Poppler binary, unlike pdftoppm/pdf2image).
    png_path = os.path.join(plots_dir, "Fig_4_full.png")
    pdf_doc = pymupdf.open(output_path)
    zoom = 600 / 72  # pymupdf's Pixmap render defaults to 72dpi
    pix = pdf_doc[0].get_pixmap(matrix=pymupdf.Matrix(zoom, zoom))
    pix.save(png_path)
    pdf_doc.close()
    print(f"Saved {png_path}")


def autocrop_to_content(img, padding_frac=0.05, background_frac=0.98):
    """Crop a PyMOL render to its non-white bounding box (plus a small margin) - same helper as
    figure_1_plot.py's autocrop_to_content."""
    gray = img[..., :3].mean(axis=2)
    max_val = 1.0 if np.issubdtype(img.dtype, np.floating) else 255
    mask = gray < background_frac * max_val
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    if rows.size == 0 or cols.size == 0:
        return img

    h, w = img.shape[:2]
    y0, y1 = rows[0], rows[-1]
    x0, x1 = cols[0], cols[-1]
    pad_y = int((y1 - y0) * padding_frac)
    pad_x = int((x1 - x0) * padding_frac)
    y0 = max(y0 - pad_y, 0)
    y1 = min(y1 + pad_y, h - 1)
    x0 = max(x0 - pad_x, 0)
    x1 = min(x1 + pad_x, w - 1)
    return img[y0:y1 + 1, x0:x1 + 1]


def crop_to_aspect(img, target_aspect):
    """Center-crops img to exactly target_aspect (width/height) - a "cover" crop (like a CSS
    background-size: cover), trimming whichever dimension is oversized relative to the target
    box's own shape. Used so a render fills its panel edge-to-edge with no letterboxing, at the
    deliberate cost of trimming some content off the long axis, rather than adjustable="box"'s
    own "contain" behavior (shrinks to fit, leaving blank margin) - see plot_structure_panel."""
    h, w = img.shape[:2]
    cur_aspect = w / h
    if cur_aspect > target_aspect:
        new_w = int(h * target_aspect)
        x0 = (w - new_w) // 2
        return img[:, x0:x0 + new_w]
    else:
        new_h = int(w / target_aspect)
        y0 = (h - new_h) // 2
        return img[y0:y0 + new_h, :]


def _gradient_hex_rgb(v, white_at, red_at, hex_color):
    """Same white(low)->color(high) linear gradient as output/70_filtering/filtered_hits_explorer/
    make_visualizer.py's own gradientHex() JS function, ported to Python so panel f's corner grids
    (plot_compound_cards_panel) use the identical color convention as that explorer for the same
    columns. Returns an (r, g, b) tuple in [0, 1], usable directly as a matplotlib facecolor."""
    t = (v - white_at) / (red_at - white_at)
    t = max(0.0, min(1.0, t))
    c2 = tuple(int(hex_color[i:i + 2], 16) / 255 for i in (1, 3, 5))
    return tuple(1.0 + (c2[i] - 1.0) * t for i in range(3))


def _color_ligand(selection, carbon_color, color_name):
    """Color-by-element with a custom carbon color - same helper as figure_3_plot.py's
    _color_ligand / 47b's color_ligand."""
    cmd.util.cbag(selection)
    cmd.set_color(color_name, carbon_color)
    cmd.color(color_name, f"{selection} and elem C")


def _parse_pocket_residues(resi_str):
    """"D_148 D_154 ..." -> (chain, [148, 154, ...]) - same chain_resn format used throughout
    the project's pocket_detection_data.csv files. Assumes a single chain (true for every pocket
    used here)."""
    tokens = resi_str.split()
    chain = tokens[0].split("_")[0]
    resnums = [int(tok.split("_")[1]) for tok in tokens]
    return chain, resnums


def _add_native_pocket_sphere(obj_name, centroid_xyz, color_name):
    """7K98_pocket_6's own centroid, already in "structure"'s coordinate frame - no alignment."""
    cmd.pseudoatom(obj_name, pos=list(centroid_xyz), vdw=POCKET_SPHERE_RADIUS)
    # A freshly created pseudoatom defaults to some small marker representation (nonbonded
    # cross and/or nb_spheres dot) that stays on even after enabling spheres - hiding "nonbonded"
    # alone didn't clear it, so hide everything unconditionally first instead of guessing which
    # rep it actually is.
    cmd.hide("everything", obj_name)
    cmd.show("spheres", obj_name)
    cmd.color(color_name, obj_name)
    cmd.set("sphere_transparency", POCKET_SPHERE_TRANSPARENCY, obj_name)


def _add_aligned_pocket_sphere(obj_name, mobile_path, mobile_chain, resnums, centroid_xyz, target_chain, color_name):
    """Aligns a monomeric model's pocket (via its own resi CA atoms) onto "structure"'s
    corresponding chain, then keeps only a sphere at that pocket's centroid - now transformed
    into "structure"'s own coordinate frame."""
    resi_str = "+".join(str(r) for r in resnums)
    mobile_obj = f"_mobile_{obj_name}"
    cmd.load(mobile_path, mobile_obj)
    cmd.pseudoatom(f"_centroid_{obj_name}", pos=list(centroid_xyz), vdw=POCKET_SPHERE_RADIUS)
    cmd.create(mobile_obj, f"{mobile_obj} or _centroid_{obj_name}")
    cmd.delete(f"_centroid_{obj_name}")

    rms, n_atoms = cmd.align(
        f"{mobile_obj} and chain {mobile_chain} and resi {resi_str} and name CA and polymer.protein",
        f"structure and chain {target_chain} and resi {resi_str} and name CA and polymer.protein",
    )[:2]
    print(f"    Aligned {obj_name} onto structure chain {target_chain}: RMSD={rms:.2f} A over {n_atoms} atoms")

    cmd.create(obj_name, f"{mobile_obj} and resn PSD")
    cmd.delete(mobile_obj)
    # Same default-marker issue as _add_native_pocket_sphere - hide everything before showing spheres.
    cmd.hide("everything", obj_name)
    cmd.show("spheres", obj_name)
    cmd.color(color_name, obj_name)
    cmd.set("sphere_transparency", POCKET_SPHERE_TRANSPARENCY, obj_name)


def add_pocket_spheres():
    """Adds the 3 curated pheS/pheT pockets from output/selected_pockets.csv (see module
    docstring / SKIP_POCKET_NAMES) as large spheres, one distinct POCKET_COLORS color each."""
    for pocket_name, color in POCKET_COLORS.items():
        cmd.set_color(f"pocket_color_{pocket_name}", list(color))

    selected = pd.read_csv(SELECTED_POCKETS_CSV)
    selected = selected[selected["gene_name"].isin(GENE_TO_TARGET_CHAIN) & ~selected["pocket_name"].isin(SKIP_POCKET_NAMES)]

    multimer_data = pd.read_csv(MULTIMER_POCKET_DATA_CSV)
    monomer_data = pd.read_csv(POCKET_DATA_CSV)

    for _, row in selected.iterrows():
        gene, pocket_name = row["gene_name"], row["pocket_name"]
        obj_name = f"pocket_{pocket_name}"

        if "_model_" not in pocket_name:
            # Native multimer pocket (e.g. 7K98_pocket_6) - already in "structure"'s own frame.
            structure_name, pocket_number = pocket_name.rsplit("_pocket_", 1)
            prow = multimer_data[(multimer_data["File name"] == f"{structure_name}.pdb")
                                  & (multimer_data["Pocket number"] == int(pocket_number))].iloc[0]
            centroid = [float(v) for v in prow["Pocket centroid coordinate (x y z)"].split()]
            _add_native_pocket_sphere(obj_name, centroid, f"pocket_color_{pocket_name}")
            print(f"  {gene}: {pocket_name} - native centroid, no alignment needed")
        else:
            structure_name, pocket_number = pocket_name.rsplit("_pocket_", 1)
            uniprot_ac = GENE_TO_UNIPROT[gene]
            prow = monomer_data[(monomer_data["Uniprot AC"] == uniprot_ac)
                                 & (monomer_data["File name"] == f"{structure_name}.pdb")
                                 & (monomer_data["Pocket number"] == int(pocket_number))].iloc[0]
            centroid = [float(v) for v in prow["Pocket centroid coordinate (x y z)"].split()]
            mobile_chain, resnums = _parse_pocket_residues(prow["Pocket residues (chain_resn)"])
            mobile_path = os.path.join(ALIGNED_DIR, uniprot_ac, f"{structure_name}.pdb")
            _add_aligned_pocket_sphere(obj_name, mobile_path, mobile_chain, resnums, centroid,
                                        GENE_TO_TARGET_CHAIN[gene], f"pocket_color_{pocket_name}")


PNG_PATH = os.path.join(renders_dir, "figure_4a_7K98.png")
PSE_PATH = os.path.join(renders_dir, "figure_4a_7K98.pse")


def render_7k98(rerun=False):
    if os.path.exists(PNG_PATH) and not rerun:
        print(f"Reusing existing render: {PNG_PATH}")
        return PNG_PATH

    pymol.finish_launching(["pymol", "-cq"])
    cmd.reinitialize()

    cmd.load(SOURCE_CIF, "structure")
    cmd.hide("everything", "structure")

    # Only one copy of each half of the biological assembly is shown - chain A (pheS), chain B
    # (pheT), chain C (tRNA) - not the second D/E/F copy: all 3 pocket spheres below are only
    # meaningful relative to this A/B copy anyway (7K98_pocket_6's own residues are A+B, and the
    # two aligned pockets were placed onto chain A / chain B specifically), so D/E/F would just
    # be a redundant, unannotated duplicate competing for space in this already-tight framing
    # (user-confirmed trade-off against showing the full symmetric tetramer).
    cmd.set_color("phes_color", list(COLOR_PHES))
    cmd.set_color("phet_color", list(COLOR_PHET))
    cmd.show("cartoon", "structure and chain A and polymer.protein")
    cmd.color("phes_color", "structure and chain A and polymer.protein")
    cmd.show("cartoon", "structure and chain B and polymer.protein")
    cmd.color("phet_color", "structure and chain B and polymer.protein")

    # tRNA(Phe), chain C only (see above).
    cmd.set_color("trna_color", list(COLOR_TRNA))
    cmd.show("cartoon", "structure and chain C and polymer.nucleic")
    cmd.color("trna_color", "structure and chain C and polymer.nucleic")

    # Native co-crystallized ligand F-AMS (W5Y), chain A's copy only (see above).
    ligand_sel = "structure and resn W5Y and chain A"
    cmd.show("sticks", ligand_sel)
    _color_ligand(ligand_sel, COLOR_LIGAND, "ligC_w5y")

    # Chain B (pheT/beta, one of its two copies) also gets an opaque surface on top of its
    # cartoon, to show the subunit's molecular envelope. Computed on a standalone copy
    # (cmd.create), not directly on "structure and chain B" - PyMOL's surface calculation
    # otherwise still "sees" the rest of "structure" (chain A/C and beyond) as neighboring atoms,
    # which cuts a flat, artificially trimmed face where the selection boundary sits rather than
    # a true molecular envelope. The copy inherits chain B's already-applied phet_color (colored
    # above, before this point).
    cmd.create("chainB_surface", "structure and chain B and polymer.protein")
    # cmd.create carries over whatever reps were already on for these atoms (cartoon, from
    # earlier) - hidden here so it isn't drawn twice (once via "structure", once via this copy).
    cmd.hide("cartoon", "chainB_surface")
    cmd.show("surface", "chainB_surface")
    cmd.set("transparency", 0, "chainB_surface")

    # Solvent (waters) and metallic atoms (the bound Mg2+ ions) are never shown - along with
    # Cl-/glycerol, none of these crystallization/cryo additives are part of the rendered view.
    cmd.hide("everything", "structure and (solvent or metals)")

    # 3 curated pheS/pheT pockets (output/selected_pockets.csv), as large crimson spheres.
    add_pocket_spheres()

    # Standard project render settings, matching figure_1_plot.py / figure_3_plot.py /
    # 47b_reference_pocket_visualization.py's identical block.
    cmd.bg_color("white")
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_trace_fog", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_trace_gain", 0.02)
    cmd.set("spec_reflect", 0)
    cmd.set("specular", 0)
    cmd.set("antialias", 2)

    # Hand-tuned view (via PyMOL's own get_view), matching figure_1_plot.py's MY_VIEWS
    # per-structure convention.
    # cmd.set_view((
    #     0.394189388, -0.309228092, 0.865439832,
    #     -0.092989191, 0.923435509, 0.372304469,
    #     -0.914309084, -0.227235228, 0.335254282,
    #     -0.000425093, -0.000246342, -259.175048828,
    #     -6.897897720, 7.815885544, -67.098320007,
    #     159.797851562, 358.588134766, -20.000000000,
    # ))

    cmd.set_view ((
     0.688825428,   -0.414848149,    0.594485402,\
     0.114377744,    0.871984839,    0.475970566,\
    -0.715845942,   -0.259863168,    0.648097098,\
    -0.000063151,   -0.000413832, -271.099914551,\
    -4.010656357,    7.557247162,  -64.323608398,\
   171.754882812,  370.545410156,  -20.000000000 ))

    cmd.ray(2400, 2400)
    cmd.png(PNG_PATH, dpi=600)
    cmd.save(PSE_PATH)
    print(f"Saved render: {PNG_PATH}")
    print(f"Saved session: {PSE_PATH}")
    return PNG_PATH


# ===========================================================================
# Panel b (prototype): 4x3 grid, one snapshot per output/selected_pockets.csv pocket - structure
# surface + pocket centroid sphere + top-1 docked ligand pose (user-confirmed style, over a
# plain surface+sphere or a cartoon+residues alternative). Grid order matches
# selected_pockets.csv's own row order, filled row-major (row 1 = pheS's 3 pockets, row 2 =
# pheT + aspS's first 2, etc.) - not one row per gene, since row counts differ per gene (4/3/2/3)
# and wouldn't tile a plain 4x3 grid evenly. Built one pocket at a time (--preview-pocket) before
# assembling the full grid, mirroring how panel a's own recipe (view/colors/pockets) was tuned
# through several rounds of visual feedback rather than decided upfront.
#
# Each cell also gets a gene-colored circle badge (_corner_badge, top-left) + a CAT/NON-CAT
# label (bottom center, white alpha-blended box) - same conventions as figure_1_plot.py's
# render_structure_panel (colored-circle-as-legend-handle trick) and figure_3_plot.py's
# _corner_badge / docking-score text box, applied here to every cell (not just the ones that
# rendered - see plot_pocket_grid_panel). Colors come from figure_1's own gene_to_color
# (output/plots/figure_1/color_mapping.json), matching every other figure's gene coloring -
# not the target-group scheme this used to have before the user asked to match prior figures.
# ===========================================================================

GENE_COLOR_MAPPING_PATH = os.path.join(ROOT, "output", "plots", "figure_1", "color_mapping.json")
SITE_TYPE_LABELS = {"CAT": "Catalytic", "NON-CAT": "Non-catalytic"}


def _corner_badge(ax, dot_color, label, loc="upper left", anchor_point=(0, 1)):
    """Colored-circle-as-legend-handle badge - same trick as figure_1_plot.py's
    render_structure_panel / figure_3_plot.py's _corner_badge. anchor_point defaults to
    axes-fraction (ax.legend()'s own default bbox_transform), so it works whether or not this
    cell actually has an image."""
    # stylia only defines 3 fontsize/markersize tiers (SMALL/plain/BIG - no "MEDIUM") - FONTSIZE_BIG
    # (8pt) read fine but its matched MARKERSIZE_BIG (30, a 3x jump from plain MARKERSIZE=10) made
    # the circle balloon past the label. Settled on the plain middle tier for both - FONTSIZE
    # (6pt, "medium") is legible at this grid's small per-cell size without the BIG tier's circle
    # blowing up, and MARKERSIZE at its own plain value actually fits that text.
    handle = [Line2D([0], [0], marker="o", color="w", label=label,
                      markerfacecolor=dot_color, markeredgecolor="black", markeredgewidth=0.5,
                      markersize=stylia.MARKERSIZE * 0.6)]
    legend = ax.legend(handles=handle, loc=loc, bbox_to_anchor=anchor_point, frameon=True,
                        framealpha=0.85, edgecolor="black", fontsize=stylia.FONTSIZE,
                        handletextpad=0.35, borderpad=0.3, labelspacing=0)
    ax.add_artist(legend)


def _site_type_label(ax, site_type):
    """CAT/NON-CAT label, bottom center (with a small gap above the axes' own bottom border,
    not flush against it), white alpha-blended box - same text-box style as figure_3_plot.py's
    docking-score label (_draw_compound_pose_row). FONTSIZE (stylia's plain "medium" tier) to
    match _corner_badge's text size."""
    ax.text(0.5, 0.06, SITE_TYPE_LABELS[site_type], transform=ax.transAxes,
            ha="center", va="bottom", fontsize=stylia.FONTSIZE, color="black",
            bbox=dict(facecolor="white", edgecolor="black", alpha=0.85, boxstyle="square,pad=0.3"))

POCKET_GRID_N_ROWS = 4
POCKET_GRID_N_COLS = 3

# COLOR_STRUCTURE_NEUTRAL is defined above, alongside panel a's own colors (now shared by both -
# unlike panel a, this panel's surfaces aren't gene-colored, since here color is reserved for the
# target-group badge).
# Orange - the project's standard docked-pose color, verbatim from scripts/47b/51/figure_3_plot.py.
COLOR_LIGAND_DOCKED = [0xF5 / 255, 0xA6 / 255, 0x3A / 255]

# Same fixed Angstrom cube (centered on the ligand) as figure_3_plot.py's SHOWCASE_ZOOM_BOX, so
# every pocket snapshot in the grid is framed at the same physical scale regardless of protein size.
# 16.0 (figure_3's own value) put the camera essentially inside the pocket cavity here - close
# enough that the far cavity wall rendered as a solid black backface and the ligand washed out
# under it - backed off to give the surface room to actually curve around the ligand. Used as a
# fallback (auto orient + fixed-box zoom) for any pocket without its own file in VIEWS_DIR.
POCKET_SNAPSHOT_ZOOM_BOX = 28.0

# One user-dropped <pocket_name>.txt per pocket, each holding a raw PyMOL "set_view (...)" block
# (literally pasted from PyMOL's own get_view, backslash line continuations and the "### cut
# above here ###" trailer and all - _load_pocket_view() below just regexes out the 18 numbers,
# so the exact paste formatting doesn't matter). Same per-structure hand-tuned-view idea as
# figure_1_plot.py's MY_VIEWS, but file-based instead of a hardcoded dict, since each of the 12
# pockets sits on its own structure/coordinate frame and views get dropped in one at a time.
VIEWS_DIR = os.path.join(plots_dir, "views")
os.makedirs(VIEWS_DIR, exist_ok=True)


def _load_pocket_view(pocket_name):
    """Returns an 18-float cmd.set_view() tuple parsed out of VIEWS_DIR/<pocket_name>.txt, or
    None if that file doesn't exist yet OR is still an empty placeholder (every pocket got one
    pre-created via "create one empty txt file per pocket", ready for a view to be dropped in
    later) - both cases fall back to auto-orient + fixed-box zoom. Only raises if the file has
    SOME content that doesn't parse to exactly 18 numbers, a real formatting mistake."""
    view_path = os.path.join(VIEWS_DIR, f"{pocket_name}.txt")
    if not os.path.isfile(view_path):
        return None
    with open(view_path) as f:
        text = f.read()
    if not text.strip():
        return None
    numbers = [float(n) for n in re.findall(r"[-+]?\d+\.\d+", text)]
    if len(numbers) != 18:
        raise ValueError(f"{view_path} should contain exactly 18 numbers (a PyMOL set_view "
                          f"block), found {len(numbers)}")
    return tuple(numbers)


# Multimeric (experimental PDB, e.g. 7K98_pocket_6) pipeline paths - scripts 48-50, same
# convention as scripts/51_selected_pockets_visualization.py's own MULTIMER_STRUCTURES_DIR /
# MULTIMER_DOCKING_DIR.
MULTIMER_STRUCTURES_DIR = os.path.join(ROOT, "output", "48_detect_pocket_multimers", "stripped_structures")
MULTIMER_DOCKING_DIR = os.path.join(ROOT, "output", "50_unidock_docking_multimers")


def _resolve_pocket(gene_name, pocket_name):
    """(structure_path, docking_dir) for a pocket_name, branching on the monomeric
    ("_model_" in pocket_name - AF2/AF3/SwissModel/Chai1) vs multimeric (e.g. 7K98_pocket_6)
    convention - same resolution logic as scripts/51_selected_pockets_visualization.py's
    resolve_pocket()."""
    if "_model_" in pocket_name:
        structure_name, _pocket_number = pocket_name.rsplit("_pocket_", 1)
        uniprot_ac = GENE_TO_UNIPROT[gene_name]
        structure_path = os.path.join(ALIGNED_DIR, uniprot_ac, f"{structure_name}.pdb")
        docking_dir = os.path.join(LIBRARIES["REAL"], pocket_name)
    else:
        pdb_code, _pocket_number = pocket_name.rsplit("_pocket_", 1)
        structure_path = os.path.join(MULTIMER_STRUCTURES_DIR, f"{pdb_code}.pdb")
        docking_dir = os.path.join(MULTIMER_DOCKING_DIR, pocket_name)
    return structure_path, docking_dir


def _extract_ligand_sdf(docking_dir, compound_id):
    """Extracts docking/<compound_id>_out.sdf from docking_dir/docking.tar.gz to a temp file,
    returning its path (caller must remove it) - same tarfile idiom as figure_3_plot.py's
    render_showcase_pocket / scripts/51's load_ligand_objects."""
    tar_path = os.path.join(docking_dir, "docking.tar.gz")
    with tarfile.open(tar_path, "r|gz") as tf:
        for member in tf:
            if member.name == f"docking/{compound_id}_out.sdf":
                data = tf.extractfile(member).read()
                break
        else:
            raise FileNotFoundError(f"{compound_id}_out.sdf not found in {tar_path}")
    tmp_fd, tmp_path = tempfile.mkstemp(suffix=".sdf")
    os.close(tmp_fd)
    with open(tmp_path, "wb") as f:
        f.write(data)
    return tmp_path


def render_pocket_snapshot(gene_name, site_type, pocket_name, rerun=False):
    """Renders one output/selected_pockets.csv pocket: structure surface (neutral gray) + the
    pocket's best-scoring REAL round-2 docked pose (orange sticks, COLOR_LIGAND_DOCKED - matches
    every other docked-pose render project-wide). No pocket-centroid marker (dropped,
    user-confirmed - see comment below). Both monomeric and multimeric pockets are handled (see
    _resolve_pocket). Also saves a .pse session alongside the PNG, so an exact camera view
    worked out interactively in PyMOL can be copied back in via cmd.set_view() - same workflow
    as panel a's own hand-tuned MY_VIEWS-style view."""
    png_path = os.path.join(renders_dir, f"figure_4b_{pocket_name}.png")
    pse_path = os.path.join(renders_dir, f"figure_4b_{pocket_name}.pse")
    if os.path.exists(png_path) and not rerun:
        print(f"Reusing existing render: {png_path}")
        return png_path

    # Pocket centroid sphere dropped (user-confirmed) - pocket_N.pdb is the full cluster of pocket
    # residue atoms, not a single point, so a sphere there rendered as a huge, nearly-opaque blob
    # burying the ligand rather than a discreet marker. Just surface + ligand for now, matching
    # the user-supplied mockup exactly.
    structure_path, docking_dir = _resolve_pocket(gene_name, pocket_name)
    scores = load_scores(os.path.join(docking_dir, "report.csv"))
    # Restrict to script 62's aggregated/prioritized compound pool (see AGGREGATED_HITS_CSV)
    # before picking the best - every pocket's report.csv shares the same ~99k-compound scored
    # universe, of which ~1,524 are in this pool (confirmed non-empty for all 12 pockets).
    prioritized_ids = pd.read_csv(AGGREGATED_HITS_CSV, usecols=["compound_id"])["compound_id"]
    scores = scores.loc[scores.index.intersection(prioritized_ids)]
    best_compound = scores.idxmin()
    print(f"  {pocket_name}: best prioritized compound {best_compound}, score {scores.min():.3f}")
    sdf_path = _extract_ligand_sdf(docking_dir, best_compound)

    pymol.finish_launching(["pymol", "-cq"])
    cmd.reinitialize()

    cmd.load(structure_path, "structure")
    cmd.set_color("structure_color_neutral", COLOR_STRUCTURE_NEUTRAL)
    cmd.color("structure_color_neutral", "structure")
    cmd.hide("everything", "structure")
    cmd.show("surface", "structure")
    # 0.3, matching scripts/51_selected_pockets_visualization.py's own surface+ligand
    # convention - 0.1 left the surface nearly opaque, muting the ligand underneath it to a
    # washed-out gray instead of its own orange.
    cmd.set("transparency", 0.3, "structure")

    cmd.load(sdf_path, "ligand")
    os.remove(sdf_path)
    cmd.show("sticks", "ligand")
    cmd.hide("lines", "ligand")
    _color_ligand("ligand", COLOR_LIGAND_DOCKED, "ligC_docked")

    cmd.bg_color("white")
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_trace_fog", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_trace_gain", 0.02)
    cmd.set("spec_reflect", 0)
    cmd.set("specular", 0)
    cmd.set("antialias", 2)

    view = _load_pocket_view(pocket_name)
    if view is not None:
        cmd.set_view(view)
    else:
        cmd.orient("ligand")
        _zoom_to_fixed_box_b("ligand", box_size=POCKET_SNAPSHOT_ZOOM_BOX)

    cmd.ray(1200, 1200)
    cmd.png(png_path, dpi=600)
    cmd.save(pse_path)
    print(f"Saved render: {png_path}")
    print(f"Saved session: {pse_path}")
    return png_path


def _zoom_to_fixed_box_b(selection, box_size, box_name="_zoom_box_b"):
    """Same fixed-cube zoom helper as figure_3_plot.py's _zoom_to_fixed_box."""
    (xmin, ymin, zmin), (xmax, ymax, zmax) = cmd.get_extent(selection)
    cx, cy, cz = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0.5 * (zmin + zmax)
    h = box_size / 2.0
    cmd.delete(box_name)
    corners = [
        (cx - h, cy - h, cz - h), (cx - h, cy - h, cz + h),
        (cx - h, cy + h, cz - h), (cx - h, cy + h, cz + h),
        (cx + h, cy - h, cz - h), (cx + h, cy - h, cz + h),
        (cx + h, cy + h, cz - h), (cx + h, cy + h, cz + h),
    ]
    for i, (x, y, z) in enumerate(corners):
        cmd.pseudoatom(f"{box_name}_{i}", pos=[x, y, z])
    cmd.group(box_name, f"{box_name}_*")
    cmd.zoom(box_name, buffer=0.0, complete=1)
    cmd.delete(box_name)
    cmd.delete(f"{box_name}_*")


# Default "zoom in a bit" applied in the compositing step (not by touching each pocket's own
# PyMOL camera/view file) - a plain center-crop keeping this fraction of the width and height,
# cropped evenly off all 4 sides. Used whenever a pocket has no override in POCKET_ZOOMS.
POCKET_SNAPSHOT_ZOOM_FRAC = 0.8

# Per-pocket center-crop fraction override, hardcoded (not file-based, unlike POCKET_VIEWS) -
# fill in as each of the 12 pockets gets fine-tuned; a pocket not listed here just uses
# POCKET_SNAPSHOT_ZOOM_FRAC.
POCKET_ZOOMS = {}
POCKET_ZOOMS["swissmodel_P9WFU3_model_0_pocket_1"] = 0.6
POCKET_ZOOMS["alphafold2_P9WFU1_model_0_pocket_1"] = 0.6

# Each grid cell is 4:3 (user-confirmed - first tried 16:9, too wide to leave room for the
# corner/site-type labels) rather than the near-square cell the render used to be forced into
# (panel_layout.csv's "b" row delta_y set to 7cm to match, so POCKET_GRID_N_ROWS x
# POCKET_GRID_N_COLS cells at 7cm/POCKET_GRID_N_COLS wide come out exactly 4:3 tall) -
# crop_to_aspect below trims each render to that same ratio before display, same "cover" crop as
# plot_structure_panel's own panel a.
POCKET_GRID_CELL_ASPECT = 4 / 3


def _center_crop(img, frac):
    h, w = img.shape[:2]
    y0, y1 = int(h * (1 - frac) / 2), int(h * (1 + frac) / 2)
    x0, x1 = int(w * (1 - frac) / 2), int(w * (1 + frac) / 2)
    return img[y0:y1, x0:x1]


def plot_pocket_grid_panel(letter, size, rerun=False, padding=0.0):
    """Panel b: POCKET_GRID_N_ROWS x POCKET_GRID_N_COLS grid, one cell per output/
    selected_pockets.csv row, filled row-major in that file's own order (see module docstring).
    Every cell gets a gene-colored circle badge (_corner_badge) + CAT/NON-CAT label
    (_site_type_label), even ones that fail to render, so the grid stays informative about
    which pocket a blank cell belongs to. Both monomeric and multimeric pockets are handled
    (render_pocket_snapshot / _resolve_pocket) - the recipe is still being tuned one pocket at a
    time, so any pocket that fails to render for some other reason is left as a blank (but
    still labeled) cell with a printed warning, rather than crashing the whole panel."""
    selected = pd.read_csv(SELECTED_POCKETS_CSV)
    with open(GENE_COLOR_MAPPING_PATH) as f:
        gene_to_color = json.load(f)["gene_to_color"]

    fig = plt.figure(figsize=size)
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(POCKET_GRID_N_ROWS, POCKET_GRID_N_COLS, hspace=0.05, wspace=0.05)

    for i, row in selected.iterrows():
        ax = stylize(fig.add_subplot(gs[i // POCKET_GRID_N_COLS, i % POCKET_GRID_N_COLS]))
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_edgecolor("black")
            spine.set_linewidth(stylia.LINEWIDTH)
        stylia.label(ax, xlabel="", ylabel="")

        gene_name, site_type, pocket_name = row["gene_name"], row["site_type"], row["pocket_name"]
        _corner_badge(ax, gene_to_color[gene_name], gene_name)
        _site_type_label(ax, site_type)

        try:
            png_path = render_pocket_snapshot(gene_name, site_type, pocket_name, rerun=rerun)
        except Exception as e:
            print(f"  Warning: failed to render {pocket_name} ({gene_name}, {site_type}): {e}. Leaving blank.")
            continue
        img = _center_crop(mpimg.imread(png_path), POCKET_ZOOMS.get(pocket_name, POCKET_SNAPSHOT_ZOOM_FRAC))
        img = crop_to_aspect(img, POCKET_GRID_CELL_ASPECT)
        ax.imshow(img, interpolation="none")
        ax.set_aspect("equal", adjustable="datalim")

    save_panel(fig, letter, use_tight_layout=False, padding=padding)


def plot_structure_panel(letter, size, rerun=False, padding=0.0):
    png_path = render_7k98(rerun=rerun)

    fig, ax = plt.subplots(figsize=size)
    fig.patch.set_facecolor("white")
    stylize(ax)

    img = autocrop_to_content(mpimg.imread(png_path))
    # crop_to_aspect (not adjustable="box"/"datalim") - panel_layout.csv's own "a" row is now
    # much shorter/wider (delta_y 6->3cm) than render_7k98's hand-tuned camera view was framed
    # for. adjustable="box" (an earlier attempt) shrinks the axes box to the image's own (much
    # taller) aspect instead, which keeps every pixel of content but leaves it looking like a
    # small square plopped in the middle of a wide panel with blank margin either side
    # (user-flagged: "still looks squared"). Cropping the source image itself to the panel's own
    # aspect first - same "fill the whole cell" convention as figure_1_plot.py's own structure
    # panels - trims some of the (already fairly tight) top/bottom margin around the structure
    # instead, so the render fills panel a edge-to-edge.
    img = crop_to_aspect(img, size[0] / size[1])
    ax.imshow(img, interpolation="none")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis("off")
    stylia.label(ax, xlabel="", ylabel="")

    # Legend text for this panel's 3 pocket spheres, one label per POCKET_COLORS entry (lime,
    # crimson, orchid, in that left-to-right order per user request) - same box style as panel
    # b's own _site_type_label, but each box's own facecolor now matches its pocket sphere's
    # color (user request), instead of the plain white every other caption box in this figure
    # uses - so text stays legible, edgecolor/text stay black throughout.
    for x, label, color in zip((1 / 6, 0.5, 5 / 6),
                                ["Non catalytic (int.)", "Catalytic (pheS)", "Non catalytic (pheT)"],
                                (ac.lime, ac.crimson, ac.orchid)):
        ax.text(x, 0.05, label, transform=ax.transAxes,
                ha="center", va="bottom", fontsize=stylia.FONTSIZE, color="black",
                bbox=dict(facecolor=color, edgecolor="black", alpha=0.6, boxstyle="square,pad=0.3"))

    # use_tight_layout=False (same near-zero-margin path as panel b) so the render fills the
    # whole panel box instead of tight_layout's conservative padding around a blank/off axis.
    save_panel(fig, letter, use_tight_layout=False, padding=padding)


def load_replicate_robustness():
    """Pool script 65's 12 per-pocket results.csv files (compound, score, score_std - mean/std
    across 5 Uni-Dock replicate runs) into one long (compound_id, pocket_name, score, score_std)
    table over the ~2,923 aggregated hits x 12 pockets."""
    pockets = pd.read_csv(SELECTED_POCKETS_CSV)["pocket_name"].tolist()
    frames = []
    for pocket_name in pockets:
        df = pd.read_csv(os.path.join(DOCKING_RESULTS_DIR_65, pocket_name, "results.csv"),
                          usecols=["compound", "score", "score_std"])
        df = df.rename(columns={"compound": "compound_id"})
        df["pocket_name"] = pocket_name
        frames.append(df)
    return pd.concat(frames, ignore_index=True)[["compound_id", "pocket_name", "score", "score_std"]]


def load_docking_vs_affinity():
    """Mean docking score (output/70_filtering/filtered_hits.csv, one column per pocket) vs.
    Boltz-2-predicted IC50 in nM (output/75_boltz2_collect_affinities/affinity_results.csv's
    affinity_pred_value, which is log10(IC50 in uM) - see that script's own docstring), over the
    ~1,095 filtered hits x 12 pockets, joined on (compound_id, pocket_name)."""
    pockets = pd.read_csv(SELECTED_POCKETS_CSV)["pocket_name"].tolist()
    scores = pd.read_csv(FILTERED_HITS_CSV, usecols=["compound_id"] + pockets)
    scores = scores.melt(id_vars="compound_id", value_vars=pockets, var_name="pocket_name", value_name="score")

    affinity = pd.read_csv(AFFINITY_RESULTS_CSV, usecols=["compound_id", "pocket_name", "affinity_pred_value",
                                                           "affinity_probability_binary"])
    merged = scores.merge(affinity, on=["compound_id", "pocket_name"], how="inner")
    merged["ic50_nM"] = (10 ** merged["affinity_pred_value"]) * 1000
    return merged


def _parse_pocket_name(pocket_name):
    """"<pred_type>_<UniprotAC>_model_<N>_pocket_<K>" -> (uniprot_ac, file_name, pocket_number) -
    same parsing used to look up a pocket's row in POCKET_CLUSTERS_CSV."""
    m = re.match(r"^(.+)_([A-Z0-9]+)_model_(\d+)_pocket_(\d+)$", pocket_name)
    pred, ac, model, pnum = m.groups()
    return ac, f"{pred}_{ac}_model_{model}.pdb", int(pnum)


def _sibling_pockets(pocket_name, clusters):
    """Other pocket-instances in `pocket_name`'s own spatial cluster (POCKET_CLUSTERS_CSV, same
    6.14 A centroid dedup as output/77_pocket_annotation/09_cluster_pockets.py /
    figure_1_calculations.py) - i.e. structures independently detecting the same physical site."""
    ac, file_name, pnum = _parse_pocket_name(pocket_name)
    row = clusters[(clusters["Uniprot AC"] == ac) & (clusters["File name"] == file_name)
                   & (clusters["Pocket number"] == pnum)]
    cluster_id = row["spatial_cluster_id"].iloc[0]
    members = clusters[(clusters["Uniprot AC"] == ac) & (clusters["spatial_cluster_id"] == cluster_id)]
    names = [m["File name"].replace(".pdb", "") + f"_pocket_{m['Pocket number']}"
             for _, m in members.iterrows()]
    return [n for n in names if n != pocket_name]


def load_sibling_avg_scores():
    """For each of the 10 non-excluded curated pockets (SIBLING_EXCLUDED_POCKETS dropped - no
    sibling structures to compare against) x its ~2,923 aggregated hits, single-run score
    (report_old.csv, pre-replicate - "forget the replica" per this session), restricted only to
    score < SIBLING_AVG_BIN_EDGES[-1] (no lower floor - the leftmost bin at plot time is
    open-ended below, see plot_sibling_robustness_panel): the SAME compound's average docking
    score across that
    pocket's sibling structure(s) (_sibling_pockets - same physical site, different structural
    model). Sibling scores come from report_old.csv when the sibling is itself one of the 12
    curated pockets, else from the round-2 library (DOCKING_RESULTS_DIR_REAL2) - most siblings
    were never part of the curated set. Scores >= 0 ("invalid pose" - see module docstring) are
    excluded from the average, never from disk. A pair is dropped entirely if no sibling has a
    valid (< 0) score for that compound - confirmed this is overwhelmingly a real compound-library
    coverage gap (~43%, traced to compounds sourced from scripts 56/61's NON-CAT top100
    selections, which only ever docked their picks against the curated pocket, never broadly),
    not invalid-pose noise (only a handful of pairs are dropped for that reason)."""
    sel_all = pd.read_csv(SELECTED_POCKETS_CSV)
    sel = sel_all[~sel_all["pocket_name"].isin(SIBLING_EXCLUDED_POCKETS)]
    curated_pockets = set(sel_all["pocket_name"])
    clusters = pd.read_csv(POCKET_CLUSTERS_CSV)

    frames = []
    for pocket in sel["pocket_name"]:
        df = pd.read_csv(os.path.join(DOCKING_RESULTS_DIR_65, pocket, "report_old.csv"))
        df["pocket_name"] = pocket
        frames.append(df)
    all_df = pd.concat(frames, ignore_index=True)

    # No lower floor (user request) - the leftmost bin (SIBLING_AVG_BIN_EDGES[0]) is open-ended
    # below and catches everything more negative, not just a fixed window - only the loose
    # (upper) end is a real cutoff.
    hi = SIBLING_AVG_BIN_EDGES[-1]
    in_range = all_df[all_df["score"] < hi]

    sibling_score_cache = {}

    def sibling_scores(pocket):
        if pocket not in sibling_score_cache:
            if pocket in curated_pockets:
                path = os.path.join(DOCKING_RESULTS_DIR_65, pocket, "report_old.csv")
            else:
                path = os.path.join(DOCKING_RESULTS_DIR_REAL2, pocket, "report.csv")
            sibling_score_cache[pocket] = pd.read_csv(path).set_index("compound")["score"]
        return sibling_score_cache[pocket]

    rows = []
    n_no_compound_data, n_only_invalid_pose = 0, 0
    for _, r in in_range.iterrows():
        sibs = _sibling_pockets(r["pocket_name"], clusters)
        found_any, vals = False, []
        for s in sibs:
            scores = sibling_scores(s)
            if r["compound"] in scores.index:
                found_any = True
                v = scores.loc[r["compound"]]
                if v < 0:
                    vals.append(v)
        if not vals:
            n_only_invalid_pose += 1 if found_any else 0
            n_no_compound_data += 0 if found_any else 1
            continue
        rows.append(dict(compound=r["compound"], pocket_name=r["pocket_name"], score=r["score"],
                          sibling_avg=np.mean(vals)))

    result = pd.DataFrame(rows)
    print(f"Panel d: {len(all_df)} pairs (10 pockets), {len(in_range)} in (-inf,{hi}), "
          f"excluded {n_no_compound_data} (compound never docked against any sibling) + "
          f"{n_only_invalid_pose} (every sibling score was invalid-pose) -> {len(result)} analyzable")
    return result


def _density(x, y, bins=100, smoothing_sigma=2.0):
    """Fast (O(n + bins^2), not O(n^2)) per-point density: 2D histogram, Gaussian-smoothed (in
    bin units) to read as a continuous field rather than blocky per-bin cells, + bin-index
    lookup - avoids a gaussian_kde, which doesn't scale to this many points."""
    x, y = np.asarray(x), np.asarray(y)
    hist, xedges, yedges = np.histogram2d(x, y, bins=bins)
    hist = gaussian_filter(hist, sigma=smoothing_sigma)
    xidx = np.clip(np.digitize(x, xedges) - 1, 0, hist.shape[0] - 1)
    yidx = np.clip(np.digitize(y, yedges) - 1, 0, hist.shape[1] - 1)
    return hist[xidx, yidx]


SIZE_MIN, SIZE_MAX = 2, 50


def _density_order_and_size(x, y, log_x=False, log_y=False):
    """log_x/log_y=True bins density on log10(x)/log10(y) instead of the raw value - needed
    whenever the caller also sets ax.set_xscale("log")/ax.set_yscale("log"), otherwise the
    linear-space histogram bins can't resolve density among data that spans multiple orders of
    magnitude (e.g. IC50 in nM): the bulk of the points collapse into a couple of bins near zero,
    so the density stops tracking what the log-scale scatter actually looks dense. Returns
    (x, y, raw_frac, sizes) all reordered ascending by density, so the densest (biggest) points
    are drawn last, on top - otherwise sparse points plotted later in the array can sit over the
    dense core. raw_frac is the [0, 1] percentile-rank of density (uncapped - unlike
    _density_scatter's own color `frac`, which floors at low_cutoff so color doesn't fade to
    white; size has no such floor, so the sparsest points shrink to SIZE_MIN)."""
    x, y = np.asarray(x), np.asarray(y)
    density = _density(np.log10(x) if log_x else x, np.log10(y) if log_y else y)

    order = np.argsort(density)
    x, y, density = x[order], y[order], density[order]

    raw_frac = (rankdata(density, method="average") - 1) / (len(density) - 1)
    sizes = SIZE_MIN + raw_frac * (SIZE_MAX - SIZE_MIN)
    return x, y, raw_frac, sizes


def _density_scatter(ax, x, y, low_cutoff=0.35, log_x=False, log_y=False, zorder=2):
    x, y, raw_frac, sizes = _density_order_and_size(x, y, log_x=log_x, log_y=log_y)

    # frac (floored at low_cutoff) drives COLOR only, so color doesn't fade out to white -
    # separate from raw_frac's uncapped size scaling above.
    frac = low_cutoff + (1 - low_cutoff) * raw_frac

    cm = stylia.FadingColormap("turquoise")
    ax.scatter(x, y, s=sizes, color=cm.cmap(frac), linewidth=0, zorder=zorder)


def plot_sibling_robustness_panel(letter, size, padding=0.0):
    """Sibling-structure comparison - one boxplot per SIBLING_AVG_BIN_EDGES 1-unit
    docking-score bin, over load_sibling_avg_scores()'s ~10,006 analyzable pairs (10 curated
    pockets, single-run) - does a good score in one structure of a pocket predict a good score in
    that same pocket's sibling structure(s)? Standalone panel, split off from the former
    combined "robustness/affinity" panel (was its left subplot) per request - see
    plot_affinity_panel for the other half."""
    fig, ax = plt.subplots(figsize=size)
    stylize(ax)

    left_df = load_sibling_avg_scores()
    # Leftmost bin open-ended below (user request) - -inf prepended to SIBLING_AVG_BIN_EDGES so
    # the "-13" bin catches everything more negative than -13, not just a fixed window; every
    # other bin (and the upper/loose cutoff at the last edge) is unchanged. Tick labels are each
    # bin's own upper edge - same convention as plot_affinity_panel (e.g. the [-9,-8) bin reads
    # as "-8", not "-9,-8").
    left_cut_edges = [-np.inf] + SIBLING_AVG_BIN_EDGES
    left_bin_keys = [f"{lo},{hi}" for lo, hi in zip(left_cut_edges[:-1], left_cut_edges[1:])]
    left_tick_labels = [str(hi) for hi in left_cut_edges[1:]]
    left_df = left_df.assign(_bin=pd.cut(left_df["score"], bins=left_cut_edges,
                                          right=False, labels=left_bin_keys))
    left_groups = [left_df.loc[left_df["_bin"] == b, "sibling_avg"].to_numpy() for b in left_bin_keys]

    nc = stylia.NamedColors()
    # Same boxplot styling as plot_affinity_panel (translucent turquoise fill, opaque thin
    # edges/whiskers/caps/median, no fliers) - kept visually consistent across the two former
    # halves. No set_yscale("log") here (unlike plot_affinity_panel) - these are docking-score
    # kcal/mol values (roughly -12 to -6), not IC50, so a linear axis is the right one.
    line_kwargs = dict(linewidth=stylia.LINEWIDTH)
    left_positions = list(range(len(left_bin_keys)))
    ax.boxplot(left_groups, positions=left_positions, showfliers=False, showcaps=False,
               whis=(1, 99), patch_artist=True,
               boxprops=dict(facecolor=(*nc.turquoise, 0.7), edgecolor="black", **line_kwargs),
               whiskerprops=dict(color="black", **line_kwargs),
               medianprops=dict(color="black", **line_kwargs))
    ax.set_xticks(left_positions)
    ax.set_xticklabels(left_tick_labels)  # single upper-edge number, horizontal - matches
    # plot_affinity_panel's own convention, short enough not to need rotation.
    stylia.label(ax, xlabel="Docking score", ylabel="Avg. sibling score")
    ax.tick_params(axis="y", pad=1.5)
    # Integer-only y-ticks (user request) - the default locator picked half-integer steps
    # (-7.5, -10.0), which read oddly next to this axis's whole-kcal/mol values.
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%d"))
    # Extend (never shrink) the top of the view so -6 is always in range and gets its own tick
    # (user request) - the data itself doesn't reach that high, so autoscale alone stopped short.
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, max(ymax, -6))

    # Own left margin now that this axis no longer shares the former combined figure's wspace
    # gap with plot_affinity_panel - needs real room for the "Avg. sibling score" ylabel at this
    # panel's own narrower (half of the old combined) width.
    save_panel(fig, letter, use_tight_layout=False, padding=padding,
               subplots_adjust=dict(left=0.28, right=0.95, top=0.84, bottom=0.3),
               letter_margin=0.01)


def plot_affinity_panel(letter, size, padding=0.0):
    """Boltz-2-predicted IC50 (nM, log scale) binned by docking score - one boxplot per
    BOXPLOT_BIN_EDGES 1-unit bin, over load_docking_vs_affinity()'s 13,140 pairs - a direct
    scatter was tried first and rejected as too noisy to read at this population size (see
    module docstring). Standalone panel, split off from the former combined
    "robustness/affinity" panel (was its right subplot) per request - see
    plot_sibling_robustness_panel for the other half."""
    fig, ax = plt.subplots(figsize=size)
    stylize(ax)

    df = load_docking_vs_affinity()
    # Leftmost bin open-ended below (user request) - -inf prepended to BOXPLOT_BIN_EDGES so the
    # "-14" bin catches everything more negative than -14, not just a fixed window. Tick labels
    # are each bin's own upper edge (e.g. the [-9,-8) bin reads as "-8", not "-9,-8").
    cut_edges = [-np.inf] + BOXPLOT_BIN_EDGES
    bin_keys = [f"{lo},{hi}" for lo, hi in zip(cut_edges[:-1], cut_edges[1:])]
    tick_labels = [str(hi) for hi in cut_edges[1:]]
    df = df.assign(_bin=pd.cut(df["score"], bins=cut_edges, right=False, labels=bin_keys))

    groups, counts = [], []
    for bin_key in bin_keys:
        vals = df.loc[df["_bin"] == bin_key, "ic50_nM"].to_numpy()
        groups.append(vals)
        counts.append(len(vals))
    print("Panel e bin counts:", dict(zip(bin_keys, counts)))

    nc = stylia.NamedColors()
    positions = list(range(len(bin_keys)))
    # facecolor carries its own alpha (0.7, translucent fill); boxprops has no top-level "alpha"
    # key so it doesn't also apply to edgecolor (a Patch's own alpha, when set, overrides the
    # alpha channel of both face and edge uniformly) - edgecolor/whisker/cap/median lines stay
    # fully opaque. linewidth reduced to stylia.LINEWIDTH (matplotlib's boxplot rcParams default
    # is 1.0, thicker than this project's standard line weight). No explicit widths= override -
    # same default box width as plot_sibling_robustness_panel's own boxplot.
    line_kwargs = dict(linewidth=stylia.LINEWIDTH)
    box_kwargs = dict(showfliers=False, showcaps=False, whis=(1, 99), patch_artist=True,
                       whiskerprops=dict(color="black", **line_kwargs),
                       medianprops=dict(color="black", **line_kwargs))
    ax.boxplot(groups, positions=positions,
               boxprops=dict(facecolor=(*nc.turquoise, 0.7), edgecolor="black", **line_kwargs),
               **box_kwargs)
    ax.set_yscale("log")
    ax.set_xticks(positions)
    ax.set_xticklabels(tick_labels)  # single upper-edge number, horizontal (user request,
    # no ranges) - short enough not to need the earlier range-string's 45-degree rotation.

    stylia.label(ax, xlabel="Docking score", ylabel="IC50 (nM)")
    ax.tick_params(axis="y", pad=1.5)

    # Own left margin now that this axis no longer shares the former combined figure's wspace
    # gap with plot_sibling_robustness_panel - needs real room for the "IC50 (nM)" ylabel +
    # log-scale tick labels at this panel's own narrower (half of the old combined) width.
    # top=0.84 - gives this panel the same ~0.48cm absolute top inset that panel b currently has
    # at its own padding=0.12 (0.12 * b's 7cm height / 2), so this panel's content lines up with
    # panel b's despite the very different panel heights (3cm vs 7cm) - a uniform `padding`
    # fraction alone can't do this since apply_padding scales toward each panel's own center
    # using ITS OWN height. Set directly here (not via panel_layout.csv's padding column, which
    # stays 0 for this panel) so it doesn't also inflate the already-tuned left/bottom margins
    # above. right=0.95 adds a small right-hand margin (user request, carried over from the
    # former combined panel). letter_margin=0.01 (half PANEL_LABEL_MARGIN's default 0.02) nudges
    # the bold "e" a bit left, same as c/d/f (user request) - a/b keep the default.
    save_panel(fig, letter, use_tight_layout=False, padding=padding,
               subplots_adjust=dict(left=0.32, right=0.95, top=0.84, bottom=0.3),
               letter_margin=0.01)


def _load_protein_pockets(site_type):
    """{protein: [pocket_name, ...]} for `site_type` ("CAT"/"NON-CAT"), pheS+pheT merged into
    "pheST" - ported from script 68's load_protein_pockets (script 68 can't be imported, see
    module docstring), using .assign() instead of that script's direct column assignment on a
    filtered slice to avoid a SettingWithCopyWarning."""
    df = pd.read_csv(SELECTED_POCKETS_CSV)
    df = df[df["site_type"] == site_type]
    df = df.assign(gene_name=df["gene_name"].replace({"pheS": "pheST", "pheT": "pheST"}))
    return df.groupby("gene_name")["pocket_name"].apply(list).to_dict()


def _load_agg_scores(pockets):
    """{pocket: pd.Series(score, indexed by compound_id)} from script 65's aggregated docking
    results.csv, for pockets that have one - ported from script 68's build_present/load_report."""
    scores = {}
    for pocket in pockets:
        path = os.path.join(DOCKING_RESULTS_DIR_65, pocket, "results.csv")
        if os.path.isfile(path):
            scores[pocket] = load_scores(path)
    return scores


def _hit_sets_by_threshold(agg_scores, protein_pockets, threshold):
    """{protein: set(compound_id)} with docking score <= threshold, unioned across that
    protein's own pockets - ported from script 68's _score_upsets (same threshold logic, inlined
    there rather than its own named function)."""
    return {
        protein: set().union(*(
            set(agg_scores[p][agg_scores[p] <= threshold].index)
            for p in pockets_for_protein if p in agg_scores
        ))
        for protein, pockets_for_protein in protein_pockets.items()
    }


def _style_upset_by_degree(upset, degrees=None):
    """Colors both the intersection bars and matrix dots by degree (# proteins in the
    intersection) - ported from script 68's style_upset_by_degree, extended with an optional
    explicit `degrees` list. Must be called before upset.plot(). Returns {degree: color} for the
    legend.

    `degrees` lets a caller force degrees onto the returned color mapping (and thus the legend)
    beyond what's actually observed in this upset's own data - e.g. NON-CAT pockets never reach
    a 4-target intersection at its threshold (lysS has 0 hits), but its legend should still show
    "4 targets" in the same color CAT uses, since the two are shown side by side and share one
    color scheme. style_subsets() for a degree with no matching (non-empty) subset is a no-op,
    not an error. Defaults to the degrees actually present in `upset`'s own data."""
    ac = stylia.ArticleColors()
    palette = [ac.crimson, ac.tangerine, ac.amber, ac.lime, ac.turquoise, ac.cobalt, ac.periwinkle, ac.orchid]

    if degrees is None:
        degrees = sorted(set(sum(idx) for idx in upset.intersections.index))
    degree_colors = {d: palette[(d - 1) % len(palette)] for d in degrees}
    for d, color in degree_colors.items():
        upset.style_subsets(min_degree=d, max_degree=d, facecolor=color)
    return degree_colors


def _tighten_totals_gap(axes, fig, buffer=0.025):
    """upsetplot reserves a wide blank column between the totals-bar panel and the matrix/
    category-label panel (sized to fit the label text with a generous margin). Shrinks it to
    just fit the rendered category-label text (+ `buffer`, figure-fraction) by shifting the
    matrix + intersections panels left - ported from script 68's tighten_totals_gap."""
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


def _add_degree_legend(ax_bars, degree_colors):
    """Legend mapping each degree color to '{d} target(s)', placed to the right of the bars -
    ported from script 68's add_degree_legend, fontsize switched to stylia's own predefined
    tier (stylia.FONTSIZE) instead of a hardcoded point size."""
    import matplotlib.patches as mpatches
    handles = [
        mpatches.Patch(color=color, label=f"{d} target{'s' if d != 1 else ''}")
        for d, color in degree_colors.items()
    ]
    ax_bars.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1), borderaxespad=0,
                    fontsize=stylia.FONTSIZE)


def _render_upset_figure(hit_sets, element_size, title, show_legend=True, legend_degrees=None, show_yaxis=True):
    """Standalone UpSet figure (its own plt.figure() - upsetplot manages its own internal axes
    layout via fig.add_gridspec and can't be embedded in a shared gridspec/subfigure;
    upset.plot(fig=...) requires a real Figure, a matplotlib SubFigure doesn't work - it calls
    fig.get_figwidth(), which SubFigure doesn't have). Degree-colored, gap-tightened - same
    recipe as script 68's save_protein_upset, but returns the live Figure/axes instead of saving
    directly, so plot_upset_panel can force a shared intersections/totals scale across the CAT
    and NON-CAT figures before either is rasterized. No legend by default (dropped per user
    request).

    show_yaxis=False drops the "Number of compounds" axis's tick numbers only (ylabel text stays
    either way) - kept as an option, but plot_upset_panel now calls both CAT and NON-CAT with
    show_yaxis=True (user-confirmed): even though the two intersection axes share one forced
    scale (see plot_upset_panel), NON-CAT keeps its own tick numbers instead of relying on
    reading them off CAT's plot.

    element_size (points, upsetplot's own UpSet(..., element_size=...) kwarg - NOT a figsize)
    controls how large upsetplot renders each matrix row/bar-chart unit; it's the only thing
    that actually shrinks the chart geometry. Passing a smaller/larger figsize to plt.figure()
    here does nothing on its own - UpSet.plot(fig=...) always recomputes and overwrites the
    figure's own physical size from element_size and the category/subset count
    (upsetplot/plotting.py's make_grid: fig.set_figwidth/set_figheight), ignoring whatever
    figsize the figure was created with. Text (stylia's rcParams font sizes, in absolute points)
    is NOT governed by element_size, so shrinking element_size makes the chart geometry denser
    relative to fixed-size text - exactly what "make the panel smaller, keep the text size"
    means. (Rendering at a small figsize and then letting plot_upset_panel's later
    imshow-into-a-smaller-box shrink the whole raster, as this used to do, shrinks the text right
    along with the bars - that was the original bug; the follow-up mistake was then picking the
    smallest element_size that merely avoided label overlap, rather than checking legibility at
    actual print size - the panel must be viewed at ~true physical size, not zoomed in, before
    trusting it looks right.)"""
    import warnings
    import matplotlib.collections
    from upsetplot import UpSet, from_contents
    warnings.filterwarnings("ignore", category=FutureWarning, module="upsetplot")
    warnings.filterwarnings("ignore", category=UserWarning, message=".*tight_layout.*")

    upset = UpSet(from_contents(hit_sets), sort_by="degree", element_size=element_size)
    degree_colors = _style_upset_by_degree(upset, degrees=legend_degrees)
    fig = plt.figure()
    axes = upset.plot(fig=fig)
    # The connector lines between a subset's colored dots are drawn via ax.vlines(..., lw=2) deep
    # inside upsetplot's own plot_matrix (upsetplot/plotting.py) - not exposed as a UpSet(...)/
    # style_subsets() kwarg, so the only way to thin them is to find that LineCollection
    # afterwards and override its linewidth directly. Matched on type (LineCollection, not
    # PathCollection) so this doesn't also catch/shrink the dots themselves (ax.scatter's own
    # marker edges), which sit in the same axes.collections list.
    for collection in axes["matrix"].collections:
        if isinstance(collection, matplotlib.collections.LineCollection):
            collection.set_linewidth(1.0)
    fig.canvas.draw()  # force a layout pass so text extents below are accurate, not placeholders
    _tighten_totals_gap(axes, fig, buffer=0.05)
    # Re-styles upsetplot's own default gridlines to match stylia's own grid recipe (stylia's
    # stylize(), scripts/plots/figure_4_plot.py isn't itself using stylia.create_figure()/stylize()
    # here since upsetplot manages this axis, but the grid should still look like every other
    # stylia-driven plot in the project) - thin, light gray, semi-transparent, drawn behind data.
    axes["intersections"].set_axisbelow(True)
    axes["intersections"].grid(visible=True, axis="y", linewidth=stylia.LINEWIDTH * 0.5,
                                color="#DDDDDD", alpha=0.6)
    # upsetplot draws its own tick/category text at matplotlib's default size, not one of
    # stylia's tiers - visibly larger than the ylabel/title text just below and above, which ARE
    # explicitly set to stylia sizes. FONTSIZE_SMALL ("tick labels, annotations" per stylia's own
    # tier convention) matches how every other tick label in this project's figures is sized.
    axes["intersections"].tick_params(axis="y", labelsize=stylia.FONTSIZE_SMALL)
    axes["totals"].tick_params(axis="x", labelsize=stylia.FONTSIZE_SMALL)
    axes["matrix"].tick_params(axis="y", labelsize=stylia.FONTSIZE_SMALL)
    # ylabel shown on both figures (per user request); tick numbers gated by show_yaxis, which
    # plot_upset_panel now passes as True for both CAT and NON-CAT (user-confirmed). y=0.6 (up
    # from the default centered 0.5 - user request) is passed straight to the label Text's own
    # set_y, not read back via get_position()/set_label_coords first - that round-trip returned a
    # bogus pre-layout value here and blew up the saved figure's bbox_inches="tight" extent,
    # rendering almost the entire panel blank.
    axes["intersections"].set_ylabel("Number of\ncompounds", fontsize=stylia.FONTSIZE, y=0.6)
    if not show_yaxis:
        axes["intersections"].tick_params(axis="y", left=False, labelleft=False)
    # Totals bars' own count axis moved from bottom to top (user request) - tick_top() only
    # relocates the ticks/labels, not the axis spine itself, so the bottom spine (now sitting
    # where the ticks used to be, disconnected from them) is hidden and the top one shown instead.
    axes["totals"].xaxis.tick_top()
    axes["totals"].spines["bottom"].set_visible(False)
    axes["totals"].spines["top"].set_visible(True)
    # Centered over the matrix/intersections columns only, excluding the "totals" bars to their
    # left (fig.suptitle's own default x=0.5 centers over the WHOLE figure, which visibly skews
    # the title toward the totals bars - user-flagged). matrix and intersections share the same
    # x0/x1 (same category columns, one above the other), so either's get_position() gives the
    # right center; matrix is used since it's present regardless of show_yaxis.
    title_x = sum(axes["matrix"].get_position().intervalx) / 2.0
    # y lowered from suptitle's own default (~0.98, top of figure) toward the chart itself - user
    # request, fine with the title overlapping the plot area rather than sitting in its own
    # separate band above it.
    fig.suptitle(title, x=title_x, y=0.92, fontsize=stylia.FONTSIZE_BIG)
    if show_legend:
        _add_degree_legend(axes["intersections"], degree_colors)
    return fig, axes


def plot_upset_panel(letter, size, padding=0.0):
    """Panel c: CAT pockets at docking score <= CAT_SCORE_THRESHOLD vs. NON-CAT pockets at
    docking score <= NONCAT_SCORE_THRESHOLD, as 2 protein-level UpSet plots built directly from
    script 65's aggregated docking scores (see module docstring for why this isn't simply
    composited from script 68's own renders). Both the intersection-size (vertical bars) and
    per-protein-total (horizontal bars) axes are forced to a shared scale across the two figures
    before rasterizing, so bar heights are directly comparable between CAT and NON-CAT."""
    pockets = pd.read_csv(SELECTED_POCKETS_CSV)["pocket_name"].tolist()
    agg_scores = _load_agg_scores(pockets)

    cat_hits = _hit_sets_by_threshold(agg_scores, _load_protein_pockets("CAT"), CAT_SCORE_THRESHOLD)
    noncat_hits = _hit_sets_by_threshold(agg_scores, _load_protein_pockets("NON-CAT"), NONCAT_SCORE_THRESHOLD)
    print(f"Panel c - CAT (score <= {CAT_SCORE_THRESHOLD}): "
          + ", ".join(f"{p}={len(s)}" for p, s in cat_hits.items()))
    print(f"Panel c - NON-CAT (score <= {NONCAT_SCORE_THRESHOLD}): "
          + ", ".join(f"{p}={len(s)}" for p, s in noncat_hits.items()))

    # No legend on either figure (dropped per user request) - the degree-color mapping (still
    # applied to both, via _style_upset_by_degree) is otherwise identical between the two, so
    # colors stay consistent even without a legend spelling out what they mean.
    fig_cat, axes_cat = _render_upset_figure(
        cat_hits, UPSET_ELEMENT_SIZE, f"Catalytic", show_legend=False)
    fig_noncat, axes_noncat = _render_upset_figure(
        noncat_hits, UPSET_ELEMENT_SIZE, f"Non catalytic",
        show_legend=False, show_yaxis=True)

    # Intersection-size (vertical bars) axis is still forced to a shared scale across both
    # figures - each was auto-scaled independently by upsetplot to its own data, so without this
    # a bar of the same height in the two wouldn't represent the same intersection size (see
    # module docstring). The totals (horizontal bars) axis is intentionally left un-shared
    # (user-confirmed): forcing it to CAT's own much larger max (1793 vs. NON-CAT's 460) made
    # NON-CAT's own bars read as tiny against a mostly-empty scale - each now auto-scales to its
    # own data, at the cost of totals bar length no longer being directly comparable between the
    # two (unlike the intersection bars above, which stay comparable).
    shared_y = max(axes_cat["intersections"].get_ylim()[1], axes_noncat["intersections"].get_ylim()[1])
    for axes in (axes_cat, axes_noncat):
        axes["intersections"].set_ylim(0, shared_y)

    cat_png = os.path.join(renders_dir, "figure_4c_upset_CAT.png")
    noncat_png = os.path.join(renders_dir, "figure_4c_upset_NONCAT.png")
    for src_fig, png_path in ((fig_cat, cat_png), (fig_noncat, noncat_png)):
        src_fig.canvas.draw()
        # pad_inches above matplotlib's own 0.1in default - at this small a figsize (see
        # UPSET_ELEMENT_SIZE), bbox_inches="tight"'s own extent estimate for the rotated
        # "Number of compounds" ylabel runs a little short, clipping its last character/glyphs
        # without the extra margin.
        src_fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.15)
        plt.close(src_fig)

    # CAT (15 bars) and NON-CAT (7 bars) crop to different aspect ratios and both need to end up
    # exactly the same height. A GridSpec with width_ratios=aspects (an earlier attempt) only
    # allocates each column's width relative to the OTHER column, with no knowledge of the row's
    # actual available height - at panel c's own wide/short shape (see panel_layout.csv's "c"
    # row), the two images together don't need the full available width at that height, and the
    # unused width was leftover blank space split after EACH column (i.e. visibly BETWEEN the two
    # plots, not just at the panel's trailing edge). Computing each axes' rect directly (below)
    # sizes both images to the shared available height with zero gap between them, packing all
    # leftover width into a single trailing margin after NON-CAT instead.
    # padding_frac=0 (not autocrop_to_content's own 5% default) - the display height below
    # already fills the panel's full available height exactly (see disp_h_in), so any margin
    # baked into the cropped image itself shows up as literal blank space at the panel's top and
    # bottom edges rather than as harmless slack (user-flagged).
    imgs = [autocrop_to_content(mpimg.imread(p), padding_frac=0.0) for p in (noncat_png, cat_png)]
    aspects = [img.shape[1] / img.shape[0] for img in imgs]  # width/height, post-crop

    # Same margins previously passed to save_panel's subplots_adjust for this panel - baked in
    # here directly since the axes below are placed via add_axes (figure-fraction rects), which
    # fig.subplots_adjust() has no effect on. margin_right lowered from 0.995 to add a small
    # right-hand margin (user request), same as panels d/e.
    margin_left, margin_right, margin_top, margin_bottom = 0.005, 0.95, 0.99, 0.005
    fig_w_in, fig_h_in = size
    avail_w_in = fig_w_in * (margin_right - margin_left)
    avail_h_in = fig_h_in * (margin_top - margin_bottom)

    # Shared display height: the taller of the two constraints binds - either the full available
    # height fits both images plus UPSET_GAP_IN between them side by side (leftover width, panel
    # d's actual case), or the full available width is used up first (leftover height instead).
    disp_h_in = min(avail_h_in, (avail_w_in - UPSET_GAP_IN) / sum(aspects))
    widths_in = [disp_h_in * a for a in aspects]
    # Left-anchored (not centered) - any width leftover trails after the last (NON-CAT) image
    # instead of adding a leading margin before CAT, same "W" preference as the previous
    # GridSpec-based approach.
    x0_in = margin_left * fig_w_in
    y0_in = margin_bottom * fig_h_in + (avail_h_in - disp_h_in) / 2  # center any height leftover

    fig = plt.figure(figsize=size)
    fig.patch.set_facecolor("white")
    axs = []
    x_in = x0_in
    for img, w_in in zip(imgs, widths_in):
        ax = stylize(fig.add_axes([x_in / fig_w_in, y0_in / fig_h_in, w_in / fig_w_in, disp_h_in / fig_h_in]))
        ax.imshow(img, interpolation="none")
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis("off")
        axs.append(ax)
        x_in += w_in + UPSET_GAP_IN  # UPSET_GAP_IN blank strip before the next image
    stylia.label(axs[0], xlabel="", ylabel="")

    # letter_margin=0.01 nudges the bold "c" a bit left, along with d/e (user request).
    save_panel(fig, letter, use_tight_layout=False, padding=padding, letter_margin=0.01)


# Same 250x180 canvas proportions as the molecule-auditing skill's own svg_for(), scaled up 4x
# (1000x720) before rasterizing. The skill renders to an HTML <img> at native size, so 250x180 is
# plenty; save_panel here rasterizes the whole panel to a 600dpi PDF, and each card's mol cell is
# only ~1.5in wide (~900px at 600dpi) - at the skill's own 250px canvas that cell was upsampled
# and visibly blurry.
MOL_CANVAS_SCALE = 4
# RDKit's bondLineWidth is a fixed pixel count that does NOT auto-scale with canvas size, so this
# is set relative to MOL_CANVAS_SCALE (not a bare constant) to keep the skill's own bond-weight-
# to-molecule-size ratio as MOL_CANVAS_SCALE changes. The extra x1.75 on top of that 1:1 ratio is
# a deliberate departure from the skill's own weight (user request - thicker bonds).
MOL_BOND_LINE_WIDTH = MOL_CANVAS_SCALE * 1.75


def _mol_image(smiles, width=250 * MOL_CANVAS_SCALE, height=180 * MOL_CANVAS_SCALE):
    """Rasterizes one compound's 2D structure - the molecule-auditing skill's own svg_for()
    recipe (~/github/claude-ersilia-skills/skills/molecule-auditing/scripts/make_visualizer.py:
    rdCoordGen layout on a 250x180 canvas, bondLineWidth=1 at that native size), scaled up by
    MOL_CANVAS_SCALE (see above, also see MOL_BOND_LINE_WIDTH for a deliberate bond-weight
    departure from the skill's own 1:1 ratio) and through MolDraw2DCairo to a PNG buffer instead
    of MolDraw2DSVG, since this feeds a matplotlib OffsetImage rather than an HTML <img>.

    One deliberate deviation from the skill's own recipe: clearBackground=True (skill: False).
    An isolated single-molecule test (render at this panel's own tiny zoom, round-tripped through
    a save-to-PDF-then-rasterize cycle) showed no visible or measurable color difference between
    a transparent and an opaque canvas - an earlier finding to the contrary was a measurement
    artifact (accidentally sampling a nearby corner-grid legend swatch, which happens to use a
    similarly-hued blue/red, instead of the molecule's own N/O atoms). clearBackground=True is
    kept anyway since it's the simpler, more conservative choice for a raster that gets
    embedded+resized rather than rendered natively (no alpha-blend path to reason about at all),
    but it is not fixing a confirmed bug.

    Also autocrops to the drawn structure's own content bounding box (RDKit fits each molecule
    to the fixed 250x180 canvas preserving its own aspect ratio, so a long/shallow molecule like
    several of this panel's own compounds leaves wide dead margins above/below within that
    canvas - user-flagged, wanted the structure itself to fill the card the way the skill's own
    tightly-laid-out HTML cards do)."""
    mol = Chem.MolFromSmiles(smiles)
    rdCoordGen.AddCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    opts = drawer.drawOptions()
    opts.clearBackground = True
    opts.bondLineWidth = MOL_BOND_LINE_WIDTH
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    img = mpimg.imread(io.BytesIO(drawer.GetDrawingText()), format="png")
    return autocrop_to_content(img)


def _slot_size_in(fig_size, n_rows, n_cols):
    """Physical (width, height) in inches of one add_gridspec(n_rows, n_cols) slot within a
    fig_size (in inches) panel - computed analytically (not measured post-render) from
    add_gridspec's own wspace/hspace (GRIDSPEC_SPACE) and save_panel's default subplots_adjust
    margins (GRIDSPEC_MARGIN), using matplotlib's own equal-slot-with-fractional-gaps
    convention: total = slots + (slots-1) * space, in slot-width/-height units."""
    usable_w = fig_size[0] * (1 - 2 * GRIDSPEC_MARGIN)
    usable_h = fig_size[1] * (1 - 2 * GRIDSPEC_MARGIN)
    slot_w_in = usable_w / (n_cols + (n_cols - 1) * GRIDSPEC_SPACE)
    slot_h_in = usable_h / (n_rows + (n_rows - 1) * GRIDSPEC_SPACE)
    return slot_w_in, slot_h_in


def _mol_image_zoom(slot_w_in, slot_h_in, raster_w, raster_h):
    """OffsetImage zoom that makes a raster_w x raster_h raster fill MOL_IMAGE_FILL_FRAC of the
    tighter dimension of a slot_w_in x slot_h_in (inches) grid slot. OffsetImage's zoom->size
    mapping is dpi-independent (a zoom-of-1 image is nx/72 x ny/72 inches, verified empirically
    against matplotlib's own dpi_cor), so this needs no live draw/renderer to get right - deriving
    zoom this way (instead of a fixed constant) means it stays correct as panel_layout.csv's own
    "e" row size changes, which it does often while this panel's layout is still being tuned.
    raster_w/raster_h are the caller's own choice (not always _mol_image()'s fixed nominal canvas
    size) - see plot_compound_cards_panel, which passes the largest AUTOCROPPED raster among its
    own filled compounds, so every card's molecule stays honestly proportional to its real size
    rather than each cropped raster being independently fit to the slot (which would normalize
    away real size differences between compounds)."""
    return MOL_IMAGE_FILL_FRAC * min(slot_w_in * 72 / raster_w, slot_h_in * 72 / raster_h)


# Compound-number label style (user request: "bold Arial") - Arial is actually installed on this
# machine (~/.fonts/Arial.ttf, confirmed in matplotlib's own font manager), not silently
# substituted with a metric-compatible stand-in. fontsize moved from a hand-picked 7 to stylia's
# own BIG tier (8) for full stylia-tier consistency (user-confirmed) - _text_half_size_in below
# measures the real rendered glyph at whatever size this holds, so _label_xy_for_mol's
# anti-overlap search still works correctly at the new size with no other change needed.
LABEL_FONT = {"family": "Arial", "fontweight": "bold", "fontsize": stylia.FONTSIZE_BIG}
# _label_xy_for_mol searches outward from the molecule's own center, along its emptiest
# quadrant's own diagonal, for the CLOSEST point (smallest reach) whose own footprint is
# verified clear of both the structure's ink and the 4 corner grids - reach is a fraction of
# that quadrant's own half-extent (1.0 = the raster's own edge pixel). Min/max bound the search;
# step is its resolution. Replaces an earlier fixed reach (0.8, then 0.35 on user request) that
# only ever encoded "how far", never verified "is this specific point actually clear" - which a
# per-quadrant AVERAGE ink density can't guarantee for a quadrant that's mostly empty but still
# has, say, one carbonyl sitting inside it (user-flagged, with an annotated screenshot: the "6"
# label sat squarely on top of a C=O bond after LABEL_QUADRANT_REACH was pulled in to 0.35).
LABEL_REACH_MIN = 0.15
LABEL_REACH_MAX = 0.95
LABEL_REACH_STEP = 0.05
# Safety margin on the label's own MEASURED footprint (see _text_half_size_in) - also the knob
# for how much breathing room sits between a label and the structure, since _label_xy_for_mol
# keeps searching outward until this padded footprint (not the bare glyph) is ink-free. Bumped
# from 1.25 (labels sat right at the edge of the nearest bond - user-flagged, "too close").
LABEL_FOOTPRINT_MARGIN = 2.0


def _text_half_size_in(s):
    """Rendered (half_width, half_height) of `s` in LABEL_FONT, in inches - measured via a
    throwaway figure/renderer rather than estimated from a font-size heuristic. An earlier,
    hand-guessed constant (a flat 0.03 axes-fraction, same for both x and y) turned out to be
    under half the real glyph's height once actually measured this way (0.09in tall at
    LABEL_FONT's fontsize=7 vs the ~0.043in the guess implied at this panel's own slot size) -
    the "clear" check in _label_xy_for_mol was passing against a box less than half the size of
    what actually gets drawn, so it approved positions the real label still overlapped
    (user-flagged, with an annotated screenshot showing a label sitting on top of a C=O bond)."""
    fig = plt.figure()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    t = fig.text(0, 0, s, **LABEL_FONT)
    fig.canvas.draw()
    bbox = t.get_window_extent(renderer=renderer)
    plt.close(fig)
    return (LABEL_FOOTPRINT_MARGIN * bbox.width / fig.dpi / 2,
            LABEL_FOOTPRINT_MARGIN * bbox.height / fig.dpi / 2)


def _label_xy_for_mol(img, mol_zoom, slot_w_in, slot_h_in, corner_xy, corner_w, corner_h, label_text):
    """Axes-fraction (x, y) placement for a compound-number label next to (not on top of) a
    centered molecule - user request: numbered 1-6, "located depending on each molecule
    structure, non-overlapping with anything else".

    Splits the molecule's own (already autocropped) raster into 2x2 quadrants and picks the one
    with the lowest ink density (RDKit's diagonal-chain layout for most of this panel's own
    compounds leaves one or two quadrants nearly empty - e.g. the corner opposite the two ring
    clusters a chain runs between) as the label's general direction from center. Within that
    quadrant, walks outward from LABEL_REACH_MIN to LABEL_REACH_MAX (step LABEL_REACH_STEP) and
    takes the first (closest-to-structure) reach whose own MEASURED footprint (_text_half_size_in
    of label_text) doesn't overlap either the molecule's own ink or any of the 4 corner-grid
    rectangles - checked against the actual raster pixels, not just the quadrant's own average
    density, so a quadrant that's mostly empty but still has one bond/atom sitting inside it
    doesn't produce a label on top of that bond/atom. Falls back to LABEL_REACH_MAX (a
    guaranteed-empty corner, per the quadrant's own low average density) if no reach in the
    search range is fully clear."""
    h, w = img.shape[:2]
    # mpimg.imread returns floats in [0, 1] for a PNG (confirmed empirically - this array is
    # never uint8), not 0-255 - same max_val dodge as autocrop_to_content, which this raster
    # already went through once in _mol_image.
    max_val = 1.0 if np.issubdtype(img.dtype, np.floating) else 255
    gray = img[..., :3].mean(axis=2)
    ink = gray < 0.98 * max_val  # opaque white background (see _mol_image) is ~max_val
    mid_y, mid_x = h // 2, w // 2
    quadrants = {
        "tl": (slice(0, mid_y), slice(0, mid_x)),
        "tr": (slice(0, mid_y), slice(mid_x, w)),
        "bl": (slice(mid_y, h), slice(0, mid_x)),
        "br": (slice(mid_y, h), slice(mid_x, w)),
    }
    density = {k: ink[rows, cols].mean() for k, (rows, cols) in quadrants.items()}
    best = min(density, key=density.get)

    disp_w_frac = (mol_zoom * w / 72) / slot_w_in
    disp_h_frac = (mol_zoom * h / 72) / slot_h_in
    x_sign = 1 if "r" in best else -1
    y_sign = 1 if "t" in best else -1

    half_w_in, half_h_in = _text_half_size_in(label_text)
    half_x_frac = half_w_in / slot_w_in  # label footprint, in axes-fraction - x and y kept
    half_y_frac = half_h_in / slot_h_in  # separate since slot_w_in != slot_h_in in general
    half_px_x = max(1, round(half_x_frac / disp_w_frac * w))
    half_px_y = max(1, round(half_y_frac / disp_h_frac * h))

    def _xy_at(reach):
        return (0.5 + x_sign * reach * disp_w_frac / 2,
                0.5 + y_sign * reach * disp_h_frac / 2)

    def _clear_of_ink(reach):
        col = round(mid_x + x_sign * reach * mid_x)
        row = round(mid_y - y_sign * reach * mid_y)  # row 0 = top, so "t" moves row DOWN toward 0
        r0, r1 = max(0, row - half_px_y), min(h, row + half_px_y)
        c0, c1 = max(0, col - half_px_x), min(w, col + half_px_x)
        return not ink[r0:r1, c0:c1].any()

    def _clear_of_corners(x, y):
        for x0, y0 in corner_xy.values():
            if x0 - half_x_frac <= x <= x0 + corner_w + half_x_frac \
                    and y0 - half_y_frac <= y <= y0 + corner_h + half_y_frac:
                return False
        return True

    reach = LABEL_REACH_MIN
    while reach <= LABEL_REACH_MAX:
        x, y = _xy_at(reach)
        if _clear_of_ink(reach) and _clear_of_corners(x, y):
            return x, y
        reach += LABEL_REACH_STEP
    return _xy_at(LABEL_REACH_MAX)


def plot_compound_cards_panel(letter, size, padding=0.0):
    """Panel f: a PANEL_F_N_ROWS x PANEL_F_N_COLS grid of PANEL_F_N compound-structure slots, each
    the 2D structure inside a square (not rounded) border, plus a small bordered 4-square grid of
    real per-protein data in each of the 4 corners (see CORNER_GRID_SPECS, user-confirmed). Slots
    are filled in order from PANEL_F_CPD_IDS_CSV (see module docstring); any slot beyond that
    file's own row count is left empty, border still drawn.

    PANEL_F_CPD_IDS_CSV's column is named "cpd_id" but, as supplied, actually holds
    selection_table.csv's own source_key values (e.g. "s_22____28884674____59849" for CPD-897) -
    looked up against that table's source_key column rather than its cpd_id column
    (user-confirmed: keep the file's existing column/values, match on source_key)."""
    if not os.path.exists(PANEL_F_CPD_IDS_CSV):
        raise FileNotFoundError(
            f"{PANEL_F_CPD_IDS_CSV} not found. Expected a CSV with a single column 'cpd_id' "
            f"holding selection_table.csv source_key values, one row per compound, in the "
            "order the structures should be drawn.")
    source_keys = pd.read_csv(PANEL_F_CPD_IDS_CSV)["cpd_id"].tolist()
    if len(source_keys) > PANEL_F_N:
        raise ValueError(f"{PANEL_F_CPD_IDS_CSV} lists {len(source_keys)} compound(s), more than "
                          f"PANEL_F_N ({PANEL_F_N}) slots.")
    # Pad up to PANEL_F_N with None - those slots are left empty below.
    source_keys += [None] * (PANEL_F_N - len(source_keys))

    table = pd.read_csv(SELECTION_TABLE_CSV).set_index("source_key")
    missing = [k for k in source_keys if k is not None and k not in table.index]
    if missing:
        raise ValueError(f"{PANEL_F_CPD_IDS_CSV} lists source_key(s) not found in "
                          f"{SELECTION_TABLE_CSV}'s source_key column: {missing}")

    filled = [k for k in source_keys if k is not None]
    print(f"Panel f - {len(filled)}/{PANEL_F_N} slot(s) filled: "
          + ", ".join(table.loc[filled, "cpd_id"]))

    # Slot's own physical aspect (see CORNER_CELL_FRAC) - only 1.0 when PANEL_F_N_ROWS ==
    # PANEL_F_N_COLS; corner_h_frac below corrects for the general (non-square-slot) case.
    slot_w_in, slot_h_in = _slot_size_in(size, PANEL_F_N_ROWS, PANEL_F_N_COLS)
    corner_h_frac = CORNER_CELL_FRAC * slot_w_in / slot_h_in

    # Each raster is autocropped to its own content bounding box inside _mol_image (see that
    # function's own docstring) - mol_zoom is then rebased on the LARGEST cropped raster among
    # this panel's own filled compounds (not _mol_image()'s fixed nominal canvas, as before).
    mol_imgs = {k: _mol_image(table.loc[k, "smiles"]) for k in filled}
    max_raster_h = max(img.shape[0] for img in mol_imgs.values())
    max_raster_w = max(img.shape[1] for img in mol_imgs.values())
    mol_zoom = _mol_image_zoom(slot_w_in, slot_h_in, max_raster_w, max_raster_h)

    fig = plt.figure(figsize=size)
    fig.patch.set_facecolor("white")
    outer = fig.add_gridspec(PANEL_F_N_ROWS, PANEL_F_N_COLS, wspace=GRIDSPEC_SPACE, hspace=GRIDSPEC_SPACE)

    for i, source_key in enumerate(source_keys):
        ax = stylize(fig.add_subplot(outer[i // PANEL_F_N_COLS, i % PANEL_F_N_COLS]))
        ax.axis("off")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

        # Square (not rounded) border - drawn for every one of the PANEL_F_N slots, filled or not
        # (same always-show-all-outlines convention as this panel's earlier card version).
        ax.add_patch(FancyBboxPatch(
            (0, 0), 1, 1, transform=ax.transAxes, clip_on=False,
            boxstyle="square,pad=0", facecolor="none",
            edgecolor="black", linewidth=CARD_BORDER_LINEWIDTH))

        if source_key is None:
            continue

        row = table.loc[source_key]

        # OffsetImage/AnnotationBbox instead of imshow - places the RDKit raster at a uniform,
        # dynamically-fit mol_zoom (see _mol_image_zoom; never stretched/distorted to fill an
        # arbitrary box, unlike imshow's extent+aspect="auto"), centered in the slot (was
        # top-anchored - user request). zorder=1 (below the corner grids' zorder=2, added below)
        # so a molecule that overlaps a corner never covers it (user-flagged: "some grids are
        # not being seen now" - centering makes overlap with the taller-than-wide molecules more
        # likely, not less, so this needs an explicit z-order rather than relying on draw order).
        imagebox = OffsetImage(mol_imgs[source_key], zoom=mol_zoom)
        ab = AnnotationBbox(imagebox, (0.5, 0.5), xycoords="axes fraction",
                             frameon=False, box_alignment=(0.5, 0.5), pad=0, annotation_clip=False,
                             zorder=1)
        ax.add_artist(ab)

        # One bordered CORNER_GRID_N-wide horizontal row of cells in each of the 4 corners,
        # inset from the card edge by CARD_MARGIN - replaces the earlier version's pair of wide
        # center grids. Real data (see CORNER_GRID_SPECS/_gradient_hex_rgb, user-confirmed): each
        # cell is one of CORNER_GRID_PROTEINS, colored by that compound's own value in this
        # corner's metric.
        corner_w = CORNER_GRID_N * CORNER_CELL_FRAC
        corner_h = corner_h_frac
        corner_xy = {
            "top_left": (CARD_MARGIN, 1 - CARD_MARGIN - corner_h),
            "top_right": (1 - CARD_MARGIN - corner_w, 1 - CARD_MARGIN - corner_h),
            "bottom_left": (CARD_MARGIN, CARD_MARGIN),
            "bottom_right": (1 - CARD_MARGIN - corner_w, CARD_MARGIN),
        }
        for corner, (x0, y0) in corner_xy.items():
            col_suffix, white_at, red_at, hex_color = CORNER_GRID_SPECS[corner]
            # Individual vector Rectangle per cell (was a single imshow raster) - an
            # embedded raster and a vector path don't necessarily land on the exact same
            # output pixel when a PDF viewer/rasterizer resamples the page, however closely
            # their extent/xy matches on paper (confirmed: matplotlib's own window_extent for
            # the old imshow and its border Rectangle were bit-identical, yet the rendered PDF
            # still showed the raster fill bleeding past the border - user-flagged). Drawing
            # every cell as its own Rectangle makes fill and border both vector paths through
            # the same renderer, so they're pixel-exact at any zoom, in any viewer.
            for j, protein in enumerate(CORNER_GRID_PROTEINS):
                color = _gradient_hex_rgb(row[f"{protein} ({col_suffix})"], white_at, red_at, hex_color)
                cell_x0 = x0 + j * CORNER_CELL_FRAC  # j=0 (first of CORNER_GRID_PROTEINS) is
                # the leftmost cell.
                ax.add_patch(Rectangle(
                    (cell_x0, y0), CORNER_CELL_FRAC, corner_h, transform=ax.transData, clip_on=False,
                    facecolor=color, edgecolor="black", linewidth=CARD_BORDER_LINEWIDTH, zorder=2))
            ax.add_patch(Rectangle(
                (x0, y0), corner_w, corner_h, transform=ax.transData, clip_on=False,
                facecolor="none", edgecolor="black", linewidth=CARD_BORDER_LINEWIDTH, zorder=2))

        # Compound number (1-6, in PANEL_F_CPD_IDS_CSV's own row order - user request), placed
        # by _label_xy_for_mol in whichever part of THIS molecule's own structure has the most
        # free space, clear of the 4 corner grids. zorder=3 (above both the molecule's zorder=1
        # and the corner grids' zorder=2) as a legibility backstop, though the placement itself
        # is already chosen to avoid overlapping either.
        label_text = str(filled.index(source_key) + 1)
        label_x, label_y = _label_xy_for_mol(mol_imgs[source_key], mol_zoom, slot_w_in, slot_h_in,
                                              corner_xy, corner_w, corner_h, label_text)
        ax.text(label_x, label_y, label_text, transform=ax.transAxes,
                ha="center", va="center", zorder=3, **LABEL_FONT)

    # top=0.88/bottom=0.12 (was the 0.99/0.01 default) - gives this panel the same ~0.48cm
    # absolute top/bottom inset that panel b currently has at its own padding=0.12 (see panel e's
    # own save_panel call for the full rationale), so its top row clears the bold panel-letter
    # (which collided with the top-left compound card's own corner badges at the old near-zero
    # top margin) and both its top and bottom edges line up with panel b's, which spans the same
    # page height as panel c + the d/e row + this panel combined. left stays at the original
    # near-zero default; right=0.95 (was 0.99) adds a small right-hand margin, same as panels d/e
    # (user request). letter_margin as in panels c/d/e (user request: "c,d,e labels... a bit to
    # the left").
    save_panel(fig, letter, use_tight_layout=False, padding=padding,
               subplots_adjust=dict(left=0.01, right=0.95, top=0.88, bottom=0.12),
               letter_margin=0.01)


def main(rerun=False, subpanels=None):
    if subpanels is None:
        subpanels = PANEL_LETTERS
    sizes = load_panel_sizes(subpanels)
    paddings = load_panel_padding(subpanels)

    if "a" in subpanels:
        plot_structure_panel("a", sizes["a"], rerun=rerun, padding=paddings["a"])
    if "b" in subpanels:
        plot_pocket_grid_panel("b", sizes["b"], rerun=rerun, padding=paddings["b"])
    if "c" in subpanels:
        plot_upset_panel("c", sizes["c"], padding=paddings["c"])
    if "d" in subpanels:
        plot_sibling_robustness_panel("d", sizes["d"], padding=paddings["d"])
    if "e" in subpanels:
        plot_affinity_panel("e", sizes["e"], padding=paddings["e"])
    if "f" in subpanels:
        plot_compound_cards_panel("f", sizes["f"], padding=paddings["f"])

    merge_panels()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rerun", action="store_true", default=False,
                         help="Force PyMOL re-rendering of every panel's PNG(s), even if they "
                              "already exist (default: False, reusing whatever's already rendered)")
    parser.add_argument("--subpanels", type=str, default=None,
                         help=f"Comma-separated subset of panels to generate, from "
                              f"{{{','.join(PANEL_LETTERS)}}}. Default: all.")
    args = parser.parse_args()

    if args.subpanels is None:
        subpanels = PANEL_LETTERS
    else:
        subpanels = [p.strip() for p in args.subpanels.split(",")]
        invalid = [p for p in subpanels if p not in PANEL_LETTERS]
        if invalid:
            parser.error(f"Unknown panel(s) {invalid}, must be from {PANEL_LETTERS}")
    main(rerun=args.rerun, subpanels=subpanels)
