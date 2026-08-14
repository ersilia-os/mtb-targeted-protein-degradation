"""
Figure 3: multi-target docking hits from figure_3_calculations.py.

Each named panel (a-f) is rendered into its OWN standalone figure and saved as its own PDF
(Fig_3a.pdf ... Fig_3f.pdf, see main()/save_panel()), rather than assembled into one composite
figure - panels are meant to be laid out externally (e.g. in a manuscript figure tool). Each
panel's figsize comes from output/plots/figure_3/panel_layout.csv, a user-authored CSV with columns
panel, x, delta_x, y, delta_y, padding. delta_x/delta_y (cm) are each panel's exact saved page size
(converted to inches for matplotlib - see load_panel_sizes()); x/y are reserved for future
external-layout positioning; padding (0-1, see apply_padding()) is a per-dimension linear shrink
factor controlling deliberate blank margin around a panel's content WITHOUT affecting any font/
tick/label size (those are always absolute point values in matplotlib, independent of an axes' box
size) - 0 fills the page (every panel's current default), 1 shrinks the content to a zero-size
point at the center.

Panel f is still ONE file containing all 5 of its internal sub-columns (see its own paragraph
below), sized as one block; panel d and e (also internally multi-column) are otherwise independent
figures - see main() for how each is built.

Panel a follows the nested-overlaid-bars convention from
notebooks/46_docking_exploration_deliverables.ipynb: one bar per protein per cutoff, drawn loosest
cutoff (most compounds) first and strictest (fewest compounds) last, so each stricter cutoff's bar
sits fully inside the previous one instead of a standard grouped/stacked layout. The population is
figure_3_calculations.py's SELECTED_SET_CUTOFF multi-target hit set (compounds already selected as
hitting >= MIN_TARGETS genes at that cutoff), not the full screened library - so the y-axis ceiling
is that set's own size, and every bar height is a count out of that fixed total. Proteins sorted by
their compound count at cutoff -11, descending. Bars are colored with stylia's FadingColormap
("crimson" preset), sampled evenly across CUTOFFS - near-white for the loosest/outermost cutoff
(-8, most compounds) fading to vivid red for the strictest/innermost cutoff (-12, fewest but
best-scoring compounds), per request over the previous arbitrary categorical palette.

Panel b is a chord (Circos) diagram of pairwise protein hit overlap, at CIRCOS_CUTOFF, within that
same selected set - following notebooks/46_docking_exploration_IIIa_colors.ipynb's precedent
(pycirclize's Circos.chord_diagram, rasterized via FigureCanvasAgg into a stylia axes), colored with
figure_1's own gene->color mapping so protein colors stay consistent across every figure. Genes with
zero hits at the cutoff are dropped outright (a chord diagram sector needs a nonzero span); among
the rest, low-degree genes get a blank (padded-space) label instead of their real name - same
declutter heuristic as that notebook - since their sector is too small for real text without
overlapping its neighbors.

Panel c (plot_tier_grid_placeholder, name kept despite no longer being a pure placeholder) is a
TIER_GRID_ROW_LABELS rows (P2Rank prob., Docking HL, Docking REAL, Cross-pharm.) x 21 columns (all
genes) grid, colored mostly from figure_3_calculations.py's figure_3_gene_summary_stats.json - see
TIER_ROW_FIELDS for the row->field mapping. All 4 rows are continuous and graded by equal terciles
across the 21 genes (green/amber/red for best/middle/worst third, _tercile_colors) - a data-driven
split, no invented absolute cutoff, per request over hand-picked thresholds. Cross-pharm. is the
odd one out: it has no field in figure_3_gene_summary_stats.json, instead reusing
plot_circos_overlap's own node_strength metric (total pairwise multi-target-hit compound overlap
with every other gene, at CIRCOS_CUTOFF - see _compute_cross_pharmacology_scores) so this column's
ranking agrees with panel b's chord-diagram ordering (e.g. pheS > glyS > alaS > ...), per request.
The earlier Exp. tractability/Novelty boolean rows were dropped per request (mostly-NaN data, see
figure_3_calculations.py's NOVELTY_OVERRIDES/EXPERIMENTAL_TRACTABILITY_OVERRIDES). Gene rows are
ordered by TIER_GRID_SCORE (each gene's summed +1/0/-1 across its own green/amber/red cells), best
at top, worst at bottom -
per request over the previous alphabetical order.

Panel d (previously split across two analogous panels, d and e, merged into one per request) is
two stacked full-width rows, one per showcase compound -
figure_3_calculations.py's TOP_AVG_SCORE_COMPOUND_IDS, the 2 compounds (of the cutoff-12
multi-target hit tier, each hitting exactly 4 genes there) with the best average docking score
across their hit genes, a criterion the user picked explicitly (AskUserQuestion) over
best-single-score or broadest-gene-coverage alternatives. Each row is a 2D RDKit structure
depiction (render_2d_structure, no PyMOL/stylia, per project convention for plain structure images,
labeled with its own compound_id since two different compounds are on screen) followed by that
compound's 4 docking-pose renders, one per hit gene, at each gene's ACTUAL best-scoring pocket
(figure_3_calculations.py's compute_top_avg_score_compounds_pockets - the same per-gene min-score
pocket compute_gene_min_scores uses to call a hit, not merely the highest-P2Rank-probability pocket
notebooks/46_pocket_visualisation.ipynb simplifies to, and falling back to a gene's next-best pocket
- flagged with a trailing "*" on that gene's badge - when the true best pocket has no retained
pose; both compounds' hit genes turned out to already have a verified retained pose, no omissions
needed here unlike the single-compound panel this replaced). Rendering follows
scripts/47b_reference_pocket_visualization.py's established PyMOL recipe (reference structure
translucent cartoon, docked ligand pose extracted from docking.tar.gz and colored with the project's
standard orange-carbon convention), simplified to one ligand/no pocket sphere/no PDB or AlphaFill
overlays, camera oriented + zoomed tightly on the ligand alone (cmd.orient + a fixed-size zoom box,
no per-image autocrop) so every panel is framed at an identical, directly comparable physical scale.
Also shown: pocket atoms within SHOWCASE_NEAR_LIGAND_CUTOFF A of the ligand as thin context sticks
(hydrogens excluded), and detected protein-ligand H-bonds as yellow dashes (_add_hbond_dashes,
adapted from 47b's own unused draft of the same logic). Each pose render carries a gene-name badge
(top-left, colored via figure_1's gene_to_color) and the docking score (bottom, plain text in a
white-boxed frame) - the InterPro CAT/Other/NA classification and each gene's UniProt AC are still
computed/saved by compute_top_avg_score_compounds_pockets but not drawn here.

Panel f is ONE merged panel (previously two separate panels, e and f, merged per request), split
internally into 5 equally thin columns with a small gap between them (RIGHT_BLOCK_WSPACE) - the
P2Rank probability curve, then the 4 domain bands. All 5 are oriented vertically to fit that
tall/narrow shape: pocket_rank runs down the y-axis (rank 1, highest P2Rank probability, at the
top), not along x, with each domain band's title rotated 90 degrees at the top of its column
(plain horizontal text would overlap the next column over, this thin).

The 1st column is a rank-ordered area plot of P2Rank pocket probabilities across all 276 detected
pockets (output/pocket_detection_data.csv, via figure_3_calculations.py's compute_pocket_scores),
sorted probability descending down the y-axis - that sort rank (1-276) is an arbitrary pocket
identity, not a biological quantity.

The remaining 4 columns (DOMAIN_STRIP_COLUMNS/plot_domain_strip, all real data) are domain-strip
bands: each is one colored cell per pocket (all 276, same pocket_rank y-order as the probability
column), in that band's own color (red/blue/amber/green) where the condition holds, white
otherwise. Catalytic uses catalytic_confidence >= CATALYTIC_CONFIDENCE_MIN (strong ligand evidence
for the Catalytic Domain (ATP/ligase) label); the other 3 (tRNA binding, Editing, Anticodon
binding) have no analogous confidence score, so "present" means "InterPro label present, AND
catalytic_confidence < CATALYTIC_CONFIDENCE_MIN" - the latter clause so a pocket already flagged
in the Catalytic band (3 pockets carry both a catalytic and a non-catalytic label) isn't also
flagged in these. All from figure_3_calculations.py's compute_pocket_scores, sourced from
output/77_pocket_annotation/pocket_detection_interpro_updated.csv.

Usage:
    python figure_3_plot.py [--rerun]
"""
import argparse
import json
import os
import sys
import tarfile
import tempfile

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

os.environ["QT_QPA_PLATFORM"] = "offscreen"

import matplotlib as mpl
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymol
import stylia
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from pycirclize import Circos
from pymol import cmd
from rdkit import Chem
from rdkit.Chem import rdCoordGen
from rdkit.Chem.Draw import rdMolDraw2D
from stylia.config import get_fg_color
from stylia.figure.figure import stylize

from docking_utils import LIBRARIES, lookup_smiles

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

# pycirclize sets this to "tight" as an import-time side effect, which would silently override
# every plt.savefig() call below (even ones that don't pass bbox_inches="tight" themselves) and
# defeat save_panel()'s exact-figsize saving - reset to matplotlib's own default.
mpl.rcParams["savefig.bbox"] = None

plots_dir = os.path.join(root, "..", "..", "output", "plots", "figure_3")
data_dir = os.path.join(plots_dir, "data")
renders_dir = os.path.join(plots_dir, "renders")
os.makedirs(plots_dir, exist_ok=True)
os.makedirs(data_dir, exist_ok=True)
os.makedirs(renders_dir, exist_ok=True)

PANEL_LETTERS = ["a", "b", "c", "d", "e"]  # former "f" (tier grid) relabeled "e", per request
panel_layout_path = os.path.join(plots_dir, "panel_layout.csv")


def load_panel_sizes(panels):
    """panel -> (width_in, height_in), from output/plots/figure_3/panel_layout.csv (columns:
    panel, x, delta_x, y, delta_y, all in cm - x/y unused for now, reserved for future
    external-layout positioning). Only requires/returns rows for `panels` (e.g. a --subpanels
    subset), not every row in PANEL_LETTERS. Raises a clear error if the file or a needed row is
    missing."""
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
    already-tuned zero-margin behavior, unchanged)."""
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
    of an axes' box size, so this only moves/shrinks the plotted geometry, not the text."""
    scale = 1 - padding
    for ax in fig.axes:
        pos = ax.get_position()
        x0 = 0.5 + (pos.x0 - 0.5) * scale
        x1 = 0.5 + (pos.x1 - 0.5) * scale
        y0 = 0.5 + (pos.y0 - 0.5) * scale
        y1 = 0.5 + (pos.y1 - 0.5) * scale
        ax.set_position([x0, y0, x1 - x0, y1 - y0])


PANEL_LABEL_MARGIN = 0.02  # figure-fraction inset from the page corner - small breathing room


def add_panel_label(fig, letter):
    """Bold panel letter at the top-left of the FIGURE (page) itself, not the axes - stays at a
    fixed page corner regardless of panel_layout.csv's padding, unlike stylia.label()'s own abc
    mechanism (ax.set_title(loc="left"), anchored to the axes' own top edge, which would drift
    inward as padding grows). Same style values as that mechanism (bold, FONTSIZE_BIG, stylia's
    own foreground color for the active style)."""
    fig.text(PANEL_LABEL_MARGIN, 1 - PANEL_LABEL_MARGIN, letter, fontweight="bold",
             fontsize=stylia.FONTSIZE_BIG, color=get_fg_color(), ha="left", va="top",
             transform=fig.transFigure)


# Max page width a saved panel PDF may reach - stylia.SIZE is documented as Nature's two-column
# width (7.09in / ~18.0cm, dynamically resolved for the active format). No max height for now.
MAX_WIDTH_IN = stylia.SIZE


def save_panel(fig, letter, padding=0.0):
    # No bbox_inches="tight" (unlike stylia.save_figure) - that recomputes the saved page to fit
    # the actual rendered content, ignoring figsize. Panels must save at EXACTLY their
    # panel_layout.csv delta_x/delta_y, per request - plt.tight_layout() (also on via stylia's
    # figure.autolayout rcParam) still reflows labels/legends to fit within that fixed canvas, it
    # just doesn't grow the canvas to match them.
    plt.tight_layout()
    apply_padding(fig, padding)  # after tight_layout, so it isn't undone by it for a/b/c
    add_panel_label(fig, letter)
    output_path = os.path.join(plots_dir, f"Fig_3{letter}.pdf")
    plt.savefig(output_path, dpi=600, transparent=False)
    plt.close(fig)
    print(f"Saved Fig_3{letter}.pdf")


# Panel f's 5 sub-columns (probability curve + 4 domain bands) are packed with a gap between them
# (RIGHT_BLOCK_WSPACE) - independent of MOSAIC's own (default) row/column gaps. Widths are 4:1:1:1:1
# (RIGHT_BLOCK_WIDTH_RATIOS), not equal - the probability curve carries a real value per pocket, the
# 4 domain bands are just present/absent, so the curve gets more room to read, per request.
RIGHT_BLOCK_WSPACE = 0.45
RIGHT_BLOCK_WIDTH_RATIOS = [4, 1, 1, 1, 1]

# Fraction of each of panel f's 5 sub-axes reserved, at the top, for their rotated column-header
# xlabel (see reserve_top_header). tight_layout() can't budget space for this label itself -
# ax.xaxis.set_label_coords() (needed to keep the rotated label's alignment stable across
# redraws) disables matplotlib's automatic layout awareness of it, so without a manual
# reservation the label gets clipped flush against the saved page's top edge (confirmed
# empirically both by removing the reservation entirely, and by panel f's own delta_y shrinking
# from 11.7 to 8.78cm - "P2Rank prob." and "Anticodon binding" both visibly truncated at the old
# 0.07-0.12 range). This is a FRACTION of panel f's own height, calibrated against the longest
# header ("Anticodon binding") - re-verify (render + inspect, don't just eyeball) whenever
# delta_y changes, since the same fraction doesn't necessarily still clear the header at a
# different absolute height.
RIGHT_BLOCK_HEADER_FRAC = 0.13

# Fixed absolute left margin for panel c's 5 stacked rows, wide enough to fit the longest row
# label ("Anticodon binding") at its own fixed point size - see main()'s "c" block for why this
# is an absolute width (converted to a fraction there), not a fixed fraction. Measured directly
# off a render at panel c's delta_x=8cm (300dpi, left fraction 0.317/1.0in margin): the axis spine
# lands at pixel column 300 but "Anticodon binding" only starts at column 84, i.e. ~0.28in of that
# 1.0in was dead space - 0.75in leaves it flush with a small buffer instead.
ROW_LABEL_LEFT_MARGIN_IN = 0.75

# Vertical gap between panel c's P2Rank prob. row and the 4-domain-row block below it -
# independent of MOSAIC's own (default) row/column gaps, same role as RIGHT_BLOCK_WSPACE above
# but for rows.
ROW_BLOCK_HSPACE = 0.4

# Vertical gap among just the 4 domain rows themselves, via their own inner subgridspec (see
# main()'s "c" block) - tighter than ROW_BLOCK_HSPACE but not as tight as an earlier 0.3, per
# request.
DOMAIN_ROWS_HSPACE = 0.45

# Row heights for panel c's 5 stacked rows: P2Rank prob. gets 4 parts, each of the 4 domain bands
# gets 1 part (8 parts total) - so the probability curve occupies exactly 4/8 = 1/2 of the total
# stack height, per request, versus an equal 1/5 split.
ROW_BLOCK_HEIGHT_RATIOS = [4, 1, 1, 1, 1]

# Extra headroom (axes-fraction) added on top of each label's own anchor y, so the label sits
# with real empty space above the plot box instead of flush against it (or, for plot_pocket_scores,
# flush against its "0"/"1" tick labels). Tune by eye.
LABEL_GAP = 0.02


def reserve_top_header(ax, frac):
    """Shrink `ax` from the top by `frac` of its own height, leaving its bottom edge fixed - so a
    later ax.set_xlabel(..., top-positioned) has real, allocated room to render into instead of
    overflowing past the axes into whatever sits above it (see RIGHT_BLOCK_HEADER_FRAC)."""
    pos = ax.get_position()
    ax.set_position([pos.x0, pos.y0, pos.width, pos.height * (1 - frac)])
    return ax


# Loosest -> strictest, matching the draw order (background bar first, tightest cutoff on top) and
# figure_3_calculations.py's PANEL_A_CUTOFFS/SELECTED_SET_CUTOFF.
CUTOFFS = [-8, -9, -10, -11, -12]
SELECTED_SET_CUTOFF = -11
SORT_BY_CUTOFF = -11
YLIM_MAX = 475


def plot_protein_hit_counts(ax):
    counts = pd.read_csv(os.path.join(data_dir, "figure_3_selected_set_protein_hit_counts.csv"))
    counts = counts.sort_values(f"count_cutoff{abs(SORT_BY_CUTOFF)}", ascending=False).reset_index(drop=True)

    palette = stylia.FadingColormap("crimson").sample(len(CUTOFFS))
    for x, row in counts.iterrows():
        for c, cutoff in enumerate(CUTOFFS):
            label = str(cutoff) if x == 0 else None
            ax.bar(x, row[f"count_cutoff{abs(cutoff)}"], color=palette[c], edgecolor="black",
                   linewidth=stylia.LINEWIDTH, zorder=2, label=label)

    # Headroom above the tallest bar (362) so the bars themselves only fill ~80% of the axis
    # height, leaving blank space between the bars and the top spine/legend/"a" label - without
    # moving the legend or title from their original positions.
    ax.set_ylim(0, YLIM_MAX)
    ax.set_xlim(-1, 21)
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts["gene"], rotation=90, fontsize=stylia.FONTSIZE_SMALL)
    # Outside the axes to the right, vertically centered - no xlim padding trick needed since it
    # doesn't compete with the plot area at all.
    ax.legend(title="Docking score", loc="upper center", bbox_to_anchor=(0.5, 1), fontsize=stylia.FONTSIZE_SMALL, framealpha=0.8, ncol=5)
    stylia.label(ax, xlabel="", ylabel="Number of compounds")


CIRCOS_CUTOFF = -11

# Genes contributing less than this fraction of the diagram's total pairwise overlap would get a
# blank (padded-space, still correctly colored) sector label instead of their real name - per
# request, every gene keeps its real label now except CIRCOS_HIDE_LABELS below, so this is 0 (never
# triggers). Was 1/25, matching notebooks/46_docking_exploration_IIIa_colors.ipynb's threshold.
CIRCOS_DECLUTTER_FRACTION = 0

# Blanked regardless of the threshold above - every gene here has a tiny sector (count_cutoff11 of
# 1-15, vs. 45+ for the smallest of the labeled genes) and they all cluster together at the same
# point on the circle, so their real-text labels stack on top of each other into an unreadable
# smear rather than reading as separate sectors.
CIRCOS_HIDE_LABELS = {"trpS", "serS", "gltS", "leuS", "metS", "pheT", "proS", "thrS"}

CIRCOS_DPI = 200
# stylia's medium tier (FONTSIZE, not FONTSIZE_SMALL), per request. pycirclize's own label_kws
# size lives outside stylia's point system entirely (it's a creation-time value in a separate
# 8x8in/200dpi figure that then gets rasterized and rescaled into panel b's actual axes box), so
# hitting this exactly requires computing the real scale factor at runtime (see
# plot_circos_overlap's calibration pass) rather than a fixed literal - the correct label_kws size
# depends on panel b's own declared size in panel_layout.csv and on the matrix/crop content, both
# of which can change.
CIRCOS_LABEL_TARGET_FONTSIZE = stylia.FONTSIZE


def _build_cropped_circos_image(matrix, cmap, label_size):
    """Builds the circos chord diagram at `label_size` (pycirclize's own creation-time label size,
    not a final on-page point size) and crops to its non-white content bounding box - the render
    part of plot_circos_overlap's two-pass calibration, called once as a probe and once at the
    calibrated size."""
    circos = Circos.chord_diagram(
        matrix, space=2, cmap=cmap,
        label_kws=dict(size=label_size),
        link_kws=dict(ec="white", lw=0, alpha=1),
    )
    # chord_diagram() hardcodes each sector's outer arc border (ec="black", lw=0.5) with no kwarg
    # to override it - drop straight to the underlying Patch objects (2 per sector's outer track:
    # a borderless facecolor patch, then a fill-less edgecolor patch drawn on top - see pycirclize's
    # Track.axis()) and zero their linewidth before rendering, to match the borderless links above.
    for sector in circos.sectors:
        for track in sector.tracks:
            for patch in track.patches:
                patch.set_linewidth(0)
    fig_circos = circos.plotfig(figsize=(8, 8), dpi=CIRCOS_DPI)
    canvas = FigureCanvas(fig_circos)
    canvas.draw()
    w, h = fig_circos.canvas.get_width_height()
    img = np.frombuffer(canvas.tostring_argb(), dtype=np.uint8).reshape(h, w, 4)[:, :, [1, 2, 3, 0]]
    plt.close(fig_circos)

    # plotfig()'s own canvas leaves a lot of blank margin around the actual diagram - crop to the
    # non-white content bounding box (as figure_1_plot.py's autocrop_to_content does for the PyMOL
    # renders) so the diagram fills panel b instead of sitting small in a mostly-empty box.
    nonwhite = np.where((img[..., :3] < 250).any(axis=2))
    y0, y1 = nonwhite[0].min(), nonwhite[0].max()
    x0, x1 = nonwhite[1].min(), nonwhite[1].max()
    return img[y0:y1 + 1, x0:x1 + 1]


def plot_circos_overlap(ax):
    hits_path = os.path.join(data_dir, f"figure_3_multi_target_hits_cutoff{abs(CIRCOS_CUTOFF)}.csv")
    selected = pd.read_csv(hits_path)
    gene_cols = [c.removeprefix("score_") for c in selected.columns if c.startswith("score_")]

    with open(os.path.join(root, "..", "..", "output", "plots", "figure_1", "color_mapping.json")) as f:
        gene_to_color = json.load(f)["gene_to_color"]

    hit_sets = {g: set(selected.loc[selected[f"score_{g}"] <= CIRCOS_CUTOFF, "compound_id"]) for g in gene_cols}
    genes = [g for g in gene_cols if hit_sets[g]]

    matrix = pd.DataFrame(0, index=genes, columns=genes, dtype=int)
    for g1 in genes:
        for g2 in genes:
            if g1 != g2:
                matrix.loc[g1, g2] = len(hit_sets[g1] & hit_sets[g2])

    node_strength = matrix.sum(axis=0) + matrix.sum(axis=1)
    order = node_strength.sort_values(ascending=False).index
    matrix = matrix.loc[order, order]
    node_strength = node_strength.loc[order]

    # Declutter: low-degree genes' sectors are too small for real text without the labels
    # overlapping each other, so give them a unique blank label instead (still correctly colored
    # via cmap, just not named) - same heuristic and threshold as the notebook precedent.
    cmap = {g: gene_to_color[g] for g in genes}
    threshold = (node_strength.sum() / 2) * CIRCOS_DECLUTTER_FRACTION
    rename_map = {}
    for i, g in enumerate(matrix.index, start=1):
        if node_strength[g] < threshold or g in CIRCOS_HIDE_LABELS:
            blank = " " * i
            rename_map[g] = blank
            cmap[blank] = gene_to_color[g]
    matrix = matrix.rename(index=rename_map, columns=rename_map)

    # Chord diagrams double-count a symmetric matrix's upper+lower triangles - zero the lower
    # triangle (diagonal is already 0, self-overlap is meaningless) so each pair's link is drawn once.
    matrix.values[np.tril_indices_from(matrix.values, k=0)] = 0

    # Two-pass calibration: the scale factor between pycirclize's own label_kws size and a real
    # on-page point size can only be known AFTER a render is cropped and placed into `ax` (whose
    # own final box size only resolves once the layout engine actually runs) - so probe once at an
    # arbitrary size, measure the true points-per-raster-pixel via ax.transData directly (not an
    # assumed aspect-ratio match - dpi-independent of whatever plt.savefig(dpi=...) is used later,
    # since points are tied to fig.get_size_inches(), not raster resolution), then re-render at the
    # size that actually lands on CIRCOS_LABEL_TARGET_FONTSIZE for THIS panel's own current
    # delta_x/delta_y and THIS matrix's own crop shape.
    PROBE_LABEL_SIZE = 30
    probe_img = _build_cropped_circos_image(matrix, cmap, PROBE_LABEL_SIZE)
    ax.imshow(probe_img)
    ax.axis("off")
    ax.figure.tight_layout()  # matches save_panel()'s own call, so the measured box is the real one
    ax.figure.canvas.draw()

    p0x = ax.transData.transform((0, 0))[0]
    p1x = ax.transData.transform((1, 0))[0]
    pt_per_raster_px = abs(p1x - p0x) / ax.figure.dpi * 72
    probe_em_height_px = PROBE_LABEL_SIZE * CIRCOS_DPI / 72
    probe_effective_pt = probe_em_height_px * pt_per_raster_px
    calibrated_size = PROBE_LABEL_SIZE * (CIRCOS_LABEL_TARGET_FONTSIZE / probe_effective_pt)
    print(f"  Circos labels: size={PROBE_LABEL_SIZE} probe measured {probe_effective_pt:.2f}pt "
          f"-> calibrated size={calibrated_size:.2f} for target {CIRCOS_LABEL_TARGET_FONTSIZE}pt")

    final_img = _build_cropped_circos_image(matrix, cmap, calibrated_size)
    ax.images[0].remove()
    ax.imshow(final_img)
    stylia.label(ax, xlabel="", ylabel="")


# Typical P2Rank pocket-probability cutoff, per request - source/cite before this value appears
# in a manuscript legend (e.g. the P2Rank paper or its own documentation). Only used to stop the
# fill in plot_pocket_scores_row below the cutoff; no other visual change (no reference line, no
# tick relabeling), per request.
P2RANK_CUTOFF = 0.2


def plot_pocket_scores_row(ax):
    """Panel c's 1st row (of 5, alongside the 4 plot_domain_strip_row bands, per request) -
    pocket_rank on the x-axis (not y), probability on y, to match the other stacked rows'
    horizontal orientation. Transposed from plot_pocket_scores' original narrow/tall column
    layout, used when this content lived under panel f. Only fills above P2RANK_CUTOFF, per
    request - sub-cutoff pockets are left unfilled."""
    scores = pd.read_csv(os.path.join(data_dir, "figure_3_pocket_scores.csv"))
    nc = stylia.NamedColors()
    ax.plot(scores["pocket_rank"], scores["pocket_probability"], color=nc.orchid, linewidth=stylia.LINEWIDTH, zorder=2)
    above_cutoff = scores["pocket_probability"].where(scores["pocket_probability"] >= P2RANK_CUTOFF)
    ax.fill_between(scores["pocket_rank"], above_cutoff, color=nc.orchid, alpha=0.3, zorder=1)
    ax.set_xlim(scores["pocket_rank"].min(), scores["pocket_rank"].max())
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    # Explicit tick at 0.5 (unlabeled) in addition to the 0/1 endpoints, per request - a visible
    # midpoint reference without cluttering the axis with a "0.5" label.
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(["0", "", "1"])
    stylia.label(ax, xlabel="", ylabel="P2Rank prob.")
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.label.set_ha("right")
    ax.yaxis.label.set_va("center")
    # ha alone doesn't control distance from the axis - only how the text aligns around its own
    # anchor point, and this row's real "0"/"1" tick labels throw off matplotlib's automatic
    # anchor placement (Axis._autolabelpos) when combined with rotation=0, pushing the anchor (and
    # so the whole label) too far left. set_label_coords pins the anchor explicitly, close to the
    # tick labels, and disables that automatic repositioning so it stays put on every redraw.
    ax.yaxis.set_label_coords(-0.022, 0.5, transform=ax.transAxes)


def plot_pocket_scores(ax):
    # Vertical orientation (pocket_rank on the y-axis, probability on x) - this panel sits in a
    # tall/narrow column, so pockets run top-to-bottom instead of left-to-right. y-axis inverted
    # so rank 1 (highest probability) is at the top, matching the original left-to-right reading
    # order.
    scores = pd.read_csv(os.path.join(data_dir, "figure_3_pocket_scores.csv"))
    nc = stylia.NamedColors()
    ax.plot(scores["pocket_probability"], scores["pocket_rank"], color=nc.orchid, linewidth=stylia.LINEWIDTH, zorder=2)
    ax.fill_betweenx(scores["pocket_rank"], scores["pocket_probability"], color=nc.orchid, alpha=0.3, zorder=1)
    ax.set_ylim(scores["pocket_rank"].max(), scores["pocket_rank"].min())
    ax.set_xlim(0, 1)
    # Explicit tick at 0.5 (unlabeled) in addition to the 0/1 endpoints matplotlib already picks,
    # per request - a visible midpoint reference without cluttering the axis with a "0.5" label.
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(["0", "", "1"])
    ax.set_yticks([])
    stylia.label(ax, xlabel="P2Rank prob.", ylabel="")
    # Moved to the top and rotated 90 degrees, per request, to match the domain-strip columns'
    # own xlabels (plot_domain_strip) - keeps all 5 of panel f's column headers on the same edge,
    # in the same orientation, as real xlabels rather than free-floating ax.text.
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    ax.xaxis.label.set_rotation(90)
    # matplotlib rotates first, then aligns the ALREADY-ROTATED box - so ha (not va) centers a
    # 90-degree label horizontally here, and va="bottom" is what anchors its rotated-bottom edge
    # at the axis and lets it grow upward.
    ax.xaxis.label.set_ha("center")
    ax.xaxis.label.set_va("bottom")
    # set_label_coords disables matplotlib's automatic label repositioning (Axis._autolabelpos),
    # which otherwise resets this label's alignment/position on every draw - without this, the
    # ha/va above get silently overwritten and the label ends up off-center. y=0.08 clears
    # tick_top()'s "0"/"1" tick labels, plus LABEL_GAP on top of that so the label text itself
    # doesn't sit flush against the tick labels either (va="bottom" anchors its own bottom edge
    # here, so this gap is real empty space, not just tick clearance).
    ax.xaxis.set_label_coords(0.5, 1.03 + LABEL_GAP, transform=ax.transAxes)


# Panel f's 4 bands: column in figure_3_pocket_scores.csv -> band title -> its own "present" color
# (one per domain type, per request), against a shared white "absent" background. Catalytic uses a
# ligand-evidence confidence threshold (CATALYTIC_CONFIDENCE_MIN); the other 3 use plain InterPro
# label presence, mutually exclusive with Catalytic - see figure_3_calculations.py's
# compute_pocket_scores/CURATED_LABEL_COLUMNS. Each band is its own thin column within panel f's
# cell (see main()), not stacked inside one shared column - per request.
DOMAIN_STRIP_COLUMNS = [
    ("is_catalytic", "Catalytic", "crimson"),
    ("is_trna_binding", "tRNA binding", "cobalt"),
    ("is_editing", "Editing", "amber"),
    ("is_anticodon_binding", "Anticodon binding", "lime"),
]

# Each domain-row "present" mark's width, in pocket_rank units (1 unit = 1 pocket's own rank
# spacing) - 2x a pocket's natural 1-unit spacing, per request, for better legibility.
DOMAIN_ROW_BAR_WIDTH = 2.0


def plot_domain_strip_row(ax, column, title, color_name):
    """One row of panel c (4 stacked rows total, 1 per DOMAIN_STRIP_COLUMNS entry, per request) -
    a thin colored cell per pocket (all 276, same pocket_rank x-order - sorted by P2Rank
    probability descending, rank 1 leftmost), in this band's own color where `column` is True,
    white otherwise - figure_3_calculations.py's compute_pocket_scores, itself from
    output/77_pocket_annotation/pocket_detection_interpro_updated.csv. Horizontal orientation
    (pocket_rank on x, not y) to fit panel c's wide/short box - transposed from this band's
    original narrow/tall side-by-side-column layout, used when this content lived under panel f.
    Row label sits at the standard left-of-axis y-position (horizontal, not rotated) since a
    stacked row has room for it, unlike the packed narrow columns the top-rotated label convention
    was designed for."""
    scores = pd.read_csv(os.path.join(data_dir, "figure_3_pocket_scores.csv"))
    nc = stylia.NamedColors()
    present_color = getattr(nc, color_name)
    colors = [present_color if v else nc.white for v in scores[column]]
    # DOMAIN_ROW_BAR_WIDTH=2 (double each pocket's own 1-unit rank spacing), per request, so a
    # "present" mark is visually thicker/more legible - bars now overlap half their width into
    # each neighboring pocket's own column rather than sitting flush edge-to-edge.
    ax.bar(scores["pocket_rank"], 1, width=DOMAIN_ROW_BAR_WIDTH, color=colors, edgecolor="none", zorder=2)
    ax.set_xlim(scores["pocket_rank"].min(), scores["pocket_rank"].max())
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    stylia.label(ax, xlabel="", ylabel=title)
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.label.set_ha("right")
    ax.yaxis.label.set_va("center")


# Row label -> figure_3_gene_summary_stats.json field -> how to grade it. "continuous_high_better"/
# "continuous_low_better" rows are graded by equal terciles across the 21 genes (7 green/7 amber/7
# red, ranked by that field) - a data-driven split with no invented absolute cutoff, per request
# (AskUserQuestion) over hand-picked thresholds. "cross_pharmacology" has no field of its own (it
# isn't in figure_3_gene_summary_stats.json) - see _compute_cross_pharmacology_scores below.
TIER_ROW_FIELDS = [
    ("P2Rank prob.", "max_p2rank_prob", "continuous_high_better"),
    ("Docking HL", "best_hl_docking_score", "continuous_low_better"),
    ("Docking REAL", "best_real10b_docking_score", "continuous_low_better"),
    ("Cross-pharm.", None, "cross_pharmacology"),
]
TIER_GRID_ROW_LABELS = [label for label, _, _ in TIER_ROW_FIELDS]
TIER_GRID_GREEN = "lime"
TIER_GRID_AMBER = "amber"
TIER_GRID_RED = "crimson"


def _tercile_colors(values_by_gene, higher_is_better):
    """Split genes into equal thirds (7/7/7 of 21) by `values_by_gene`, best third green, middle
    third amber, worst third red - see TIER_ROW_FIELDS' comment for why terciles, not fixed
    cutoffs."""
    ordered_genes = sorted(values_by_gene, key=lambda g: values_by_gene[g], reverse=higher_is_better)
    third = len(ordered_genes) // 3
    colors = {}
    for i, gene in enumerate(ordered_genes):
        if i < third:
            colors[gene] = TIER_GRID_GREEN
        elif i < 2 * third:
            colors[gene] = TIER_GRID_AMBER
        else:
            colors[gene] = TIER_GRID_RED
    return colors


def _compute_cross_pharmacology_scores(cutoff):
    """Per-gene node strength (summed pairwise multi-target-hit compound overlap with every other
    gene) at `cutoff` - the same node_strength metric plot_circos_overlap (panel b) sorts its chord
    diagram by, so panel c's Cross-pharm. column agrees with panel b on gene order (e.g.
    pheS > glyS > alaS > ...). Computed independently of plot_circos_overlap (small duplicated
    block, same per-file-helper convention as autocrop_to_content) over ALL 21 genes, not just the
    ones with a nonzero chord sector there - a gene with zero multi-target hits at this cutoff
    legitimately gets a node strength of 0 (worst tercile), rather than being dropped."""
    hits_path = os.path.join(data_dir, f"figure_3_multi_target_hits_cutoff{abs(cutoff)}.csv")
    selected = pd.read_csv(hits_path)
    gene_cols = [c.removeprefix("score_") for c in selected.columns if c.startswith("score_")]
    hit_sets = {g: set(selected.loc[selected[f"score_{g}"] <= cutoff, "compound_id"]) for g in gene_cols}
    matrix = pd.DataFrame(0, index=gene_cols, columns=gene_cols, dtype=int)
    for g1 in gene_cols:
        for g2 in gene_cols:
            if g1 != g2:
                matrix.loc[g1, g2] = len(hit_sets[g1] & hit_sets[g2])
    return (matrix.sum(axis=0) + matrix.sum(axis=1)).to_dict()



# Per-cell-color point value for plot_tier_grid_placeholder's row ordering (green best/+1, amber
# neutral/0, red worst/-1) - summed across a gene's TIER_ROW_FIELDS columns to rank genes best to
# worst top to bottom, per request over the previous alphabetical row order.
TIER_GRID_SCORE = {TIER_GRID_GREEN: 1, TIER_GRID_AMBER: 0, TIER_GRID_RED: -1}


def plot_tier_grid_placeholder(ax):
    """21 genes (rows) x TIER_GRID_ROW_LABELS columns grid of rectangles, colored from
    figure_3_calculations.py's figure_3_gene_summary_stats.json (plus Cross-pharm., see
    _compute_cross_pharmacology_scores) - see TIER_ROW_FIELDS for the column->field mapping and
    grading rule. Rows are ordered by each gene's summed TIER_GRID_SCORE across all columns,
    best (greenest) at the top, worst (reddest) at the bottom - ties (there will be several, since
    the score is a small integer) are broken by raw cross_pharm_scores node_strength (descending),
    then alphabetically as a final tiebreak, per request. Genes-as-rows (transposed from this
    panel's original genes-as-columns layout) to fit panel f's narrow/tall box after the c/f
    content switch - 21 gene rows read fine down a tall axis, whereas 21 gene columns didn't fit
    this box's narrow width. A second per-row label column sits to the right of the grid (a twinx
    axes, so it's a real tick label rather than plotted text) reading "Sel." for genes with at
    least one row in output/selected_pockets.csv (its gene_name column), "Nov." for a gene whose
    own figure_3_gene_summary_stats.json novelty field is explicitly False (currently only ileS -
    see figure_3_calculations.py's NOVELTY_OVERRIDES), "Exp." for a gene whose experimental_
    tractability field is explicitly False (currently only glyS - see that same script's
    EXPERIMENTAL_TRACTABILITY_OVERRIDES; earlier labels win if a gene were ever more than one of
    these), blank otherwise, per request. Only labeled rows get a tick mark on this side."""
    with open(os.path.join(root, "..", "..", "output", "plots", "figure_1", "color_mapping.json")) as f:
        genes = sorted(json.load(f)["gene_to_color"].keys())
    with open(os.path.join(data_dir, "figure_3_gene_summary_stats.json")) as f:
        stats = json.load(f)
    selected_genes = set(pd.read_csv(os.path.join(root, "..", "..", "output", "selected_pockets.csv"))["gene_name"])

    nc = stylia.NamedColors()
    cross_pharm_scores = None
    col_colors = []
    for _, field, kind in TIER_ROW_FIELDS:
        if kind == "cross_pharmacology":
            cross_pharm_scores = _compute_cross_pharmacology_scores(CIRCOS_CUTOFF)
            values_by_gene = {gene: cross_pharm_scores[gene] for gene in genes}
            col_colors.append(_tercile_colors(values_by_gene, higher_is_better=True))
        else:
            values_by_gene = {gene: stats[gene][field] for gene in genes}
            col_colors.append(_tercile_colors(values_by_gene, higher_is_better=(kind == "continuous_high_better")))

    # Row order: primarily each gene's summed TIER_GRID_SCORE (best/greenest first); ties broken
    # by its raw (pre-tercile) cross_pharm_scores node_strength, descending (still no invented
    # cutoff - reuses data already computed above), then alphabetically as a final tiebreak for
    # true 0-vs-0 node-strength ties (e.g. gatA/gatB, neither of which shares any multi-target hit
    # with any other gene) - reproduces both explicitly requested tie orders (pheS > glyS > lysS >
    # valS; pheT > gatA > gatB) without hardcoding either as a literal gene-name list.
    gene_scores = {gene: sum(TIER_GRID_SCORE[colors[gene]] for colors in col_colors) for gene in genes}
    genes = sorted(genes, key=lambda g: (-gene_scores[g], -cross_pharm_scores[g], g))

    for col, colors in enumerate(col_colors):
        for row, gene in enumerate(genes):
            facecolor = getattr(nc, colors[gene])
            ax.add_patch(plt.Rectangle((col, row), 1, 1, facecolor=facecolor,
                                        edgecolor="white", linewidth=stylia.LINEWIDTH))

    ax.set_xlim(0, len(TIER_GRID_ROW_LABELS))
    ax.set_ylim(0, len(genes))
    ax.invert_yaxis()  # first gene (highest TIER_GRID_SCORE) at the top

    ax.set_yticks([i + 0.5 for i in range(len(genes))])
    ax.set_yticklabels(genes, fontsize=stylia.FONTSIZE_SMALL)
    ax.set_xticks([i + 0.5 for i in range(len(TIER_GRID_ROW_LABELS))])
    ax.set_xticklabels(TIER_GRID_ROW_LABELS, rotation=90, fontsize=stylia.FONTSIZE_SMALL)
    stylia.label(ax, xlabel="", ylabel="")
    # No outer axes border, per request - only the individual cells' own white edges (set above)
    # remain visible.
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Second per-row label column, to the right of the grid (twinx shares ax's x-axis/transform,
    # gets its own independent y-axis - ylim copied from ax, post-invert_yaxis, so rows line up) -
    # "Sel." for genes in output/selected_pockets.csv, "Nov." for a gene with an explicit
    # novelty=False, "Exp." for a gene with an explicit experimental_tractability=False (currently
    # only glyS - see figure_3_calculations.py's EXPERIMENTAL_TRACTABILITY_OVERRIDES), blank
    # otherwise, per request. Only rows with a non-blank label get a tick mark (also per request) -
    # toggled per-tick below since set_yticks/tick_params only offer one visibility setting for
    # every tick at once.
    row_labels = []
    for gene in genes:
        if gene in selected_genes:
            row_labels.append("Sel.")
        elif stats[gene]["novelty"] is False:
            row_labels.append("Nov.")
        elif stats[gene]["experimental_tractability"] is False:
            row_labels.append("Exp.")
        else:
            row_labels.append("")

    right_ax = ax.twinx()
    right_ax.set_ylim(ax.get_ylim())
    right_ax.grid(False)
    right_ax.set_yticks([i + 0.5 for i in range(len(genes))])
    right_ax.set_yticklabels(row_labels, fontsize=stylia.FONTSIZE_SMALL)
    # twinx() doesn't go through stylize(ax) above, so its ticks default to matplotlib's own
    # ytick.major.width (0.8) instead of stylia.LINEWIDTH (0.5, what ax's own left-side ticks use)
    # - set explicitly so the two sides match.
    right_ax.yaxis.set_tick_params(width=stylia.LINEWIDTH)
    for tick, label in zip(right_ax.yaxis.get_major_ticks(), row_labels):
        tick.tick2line.set_visible(bool(label))
    for spine in right_ax.spines.values():
        spine.set_visible(False)


# Matches scripts/47b_reference_pocket_visualization.py's COLOR_LIGAND_DOCKED / notebooks/
# 46_pocket_visualisation.ipynb's convention - orange carbons for a docked (as opposed to
# experimental/AlphaFill/homolog) ligand pose, kept consistent across every PyMOL render in the repo.
COLOR_LIGAND_DOCKED = [0xF5 / 255, 0xA6 / 255, 0x3A / 255]
COLOR_STRUCTURE = [0.7804, 0.8275, 0.8667]

# Angstrom cube centered on the ligand alone (not the whole structure, unlike figure_1_plot.py's/
# 47b's renders) so every panel-g sub-image is framed at a comparable, tightly-zoomed scale
# regardless of each protein's own size - a framing choice, tune by eye if too tight or too sparse.
SHOWCASE_ZOOM_BOX = 16.0

# Distance cutoff for the near-ligand context sticks (not the H-bond geometry below, which has its
# own fixed 3.5 A / 45 degree criteria).
SHOWCASE_NEAR_LIGAND_CUTOFF = 5.0


def _color_ligand(selection, carbon_color, color_name):
    cmd.util.cbag(selection)
    cmd.set_color(color_name, carbon_color)
    cmd.color(color_name, f"{selection} and elem C")


def _add_hbond_dashes(receptor_sel, ligand_sel, dash_name, cutoff=3.5, angle=45):
    """Detect donor/acceptor H-bond contacts between ligand_sel and receptor_sel, drawn as a
    yellow dashed-line object. Adapted from scripts/47b_reference_pocket_visualization.py's own
    (commented-out, unused there) add_hbond_dashes - same criteria (cutoff/angle) and
    find_pairs/distance mechanics, just actually wired in here and colored yellow per request
    instead of that draft's black. Returns the number of contacts found."""
    nearby = f"({receptor_sel} within 6 of {ligand_sel}) and not solvent"
    pairs = []
    pairs += cmd.find_pairs(f"{ligand_sel} and donor", f"{nearby} and acceptor", mode=1, cutoff=cutoff, angle=angle)
    pairs += cmd.find_pairs(f"{ligand_sel} and acceptor", f"{nearby} and donor", mode=1, cutoff=cutoff, angle=angle)
    pairs = list(dict.fromkeys(pairs))
    if not pairs:
        return 0
    cmd.delete(dash_name)
    for a1, a2 in pairs:
        cmd.distance(dash_name, a1, a2)
    cmd.hide("labels", dash_name)
    cmd.color("yellow", dash_name)
    cmd.set("dash_width", 3, dash_name)
    return len(pairs)


def _zoom_to_fixed_box(selection, box_size, box_name="_zoom_box"):
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


def render_showcase_pocket(row, rerun=False):
    # Keyed by compound_id AND gene, not gene alone - panels c and d dock different compounds into
    # some of the same genes (glyS, lysS, pheS, alaS all recur between the two), which would
    # otherwise silently collide on a gene-only cache key.
    png_path = os.path.join(renders_dir, f"figure_3_showcase_{row['compound_id']}_{row['gene']}.png")
    if os.path.exists(png_path) and not rerun:
        print(f"Reusing existing render for {row['compound_id']}/{row['gene']}: {png_path}")
        return png_path

    pymol.finish_launching(["pymol", "-cq"])
    cmd.reinitialize()

    structure_path = os.path.join(root, "..", "..", "output", "aligned_relaxed_structures",
                                   row["uniprot_ac"], f"{row['structure_name']}.pdb")
    cmd.load(structure_path, "structure")
    cmd.set_color("structure_color", COLOR_STRUCTURE)
    cmd.color("structure_color", "structure")
    cmd.show("cartoon", "structure")
    cmd.hide("lines", "structure")
    cmd.hide("sticks", "structure")
    # No pocket-sphere overlay (per request) - gene identity/color is already carried by each
    # render's own corner badge, so the cartoon's only job here is to stay out of the ligand's way.
    # Higher transparency than the original 0.3 for exactly that reason.
    # cartoon_transparency, not the generic transparency setting - the latter governs surface
    # representation and is a no-op on a cartoon (confirmed by the first render still showing a
    # fully-opaque loop occluding the ligand despite "transparency" being set).
    cmd.set("cartoon_transparency", 0.6, "structure")

    tar_path = os.path.join(LIBRARIES["REAL"], row["pocket_name"], "docking.tar.gz")
    with tarfile.open(tar_path, "r|gz") as tf:
        for member in tf:
            if member.name != f"docking/{row['compound_id']}_out.sdf":
                continue
            data = tf.extractfile(member).read()
            with tempfile.NamedTemporaryFile(suffix=".sdf", mode="wb", delete=False) as tmp:
                tmp.write(data)
                tmp_path = tmp.name
            cmd.load(tmp_path, "ligand")
            os.remove(tmp_path)
            break
    cmd.show("sticks", "ligand")
    cmd.hide("lines", "ligand")
    _color_ligand("ligand", COLOR_LIGAND_DOCKED, "ligC_docked")

    # Atoms within SHOWCASE_NEAR_LIGAND_CUTOFF A of the ligand, as sticks (element-colored, carbons
    # tinted to match the cartoon so they read as "protein" rather than a second ligand) - literal
    # per-atom distance, not byres-expanded to whole residues (tried that first - full side chains
    # reaching in from outside the cutoff radius buried the ligand in an unreadable thicket).
    # Explicit hydrogens (present in these relaxed structures) are excluded - at this zoom every C-H
    # bond became its own visible stick, still overwhelming the ligand even without byres. Thinner
    # than the ligand's own default stick_radius so the pocket reads as context, not competing with
    # the ligand for attention. cartoon_transparency only affects the cartoon representation, so
    # these sticks stay fully opaque regardless of that setting.
    cmd.select("near_ligand", f"(structure within {SHOWCASE_NEAR_LIGAND_CUTOFF} of ligand) and not hydro")
    cmd.show("sticks", "near_ligand")
    cmd.hide("lines", "near_ligand")
    cmd.set("stick_radius", 0.1, "near_ligand")
    _color_ligand("near_ligand", COLOR_STRUCTURE, "near_ligand_color")
    cmd.delete("near_ligand")

    # Protein-ligand H-bond interactions, as yellow dashes (see _add_hbond_dashes) - the identified
    # "interactions" requested, on top of (not instead of) the near-ligand context sticks above.
    _add_hbond_dashes("structure", "ligand", "hbonds")

    cmd.bg_color("white")
    # orthoscopic=1 -> true orthographic (parallel-projection) camera, not perspective - required at
    # this render's tight zoom, since a perspective camera would visibly distort/skew the ligand's
    # apparent geometry as the view gets this close.
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_trace_fog", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_trace_gain", 0.02)
    cmd.set("spec_reflect", 0)
    cmd.set("specular", 0)
    # PyMOL's default transparency_mode doesn't alpha-blend cartoon geometry correctly during ray
    # tracing at this render's tight, ligand-centered zoom - a nearby loop rendered edge-on showed up
    # as a solid opaque black blob over the ligand despite cartoon_transparency being set. Mode 2
    # fixes the blending (confirmed by isolated testing) and also makes the transparency actually
    # visible instead of still reading as near-opaque.
    cmd.set("transparency_mode", 2)
    cmd.set("antialias", 2)
    # orient (not just zoom) to the ligand - without this, the camera keeps whatever default view
    # angle the loaded structure's own coordinate frame happens to have, which for some proteins
    # (e.g. lysS) looks down the ligand's long axis almost end-on instead of showing it spread out.
    # cmd.orient() rotates to align the selection's own principal axes with the screen, which reads
    # as "flat"/extended for a small molecule regardless of the structure's native orientation.
    cmd.orient("ligand")
    _zoom_to_fixed_box("ligand", box_size=SHOWCASE_ZOOM_BOX)
    cmd.ray(1200, 1200)
    cmd.png(png_path, dpi=600)

    return png_path


def render_2d_structure(compound_id, smiles, rerun=False):
    # No PyMOL, no stylia (RDKit structure images are drawn plain, per project convention) - just a
    # cached RDKit depiction, same cache-by-file-existence idiom as render_showcase_pocket. Drawing
    # style matches the molecule-auditing skill's svg_for() (CoordGen coordinates), via the Cairo
    # PNG backend instead of SVG to embed like every other raster image in this file -
    # padding/bondLineWidth deviate from that skill's exact values per request (smaller molecule
    # within the same canvas, thicker bonds).
    png_path = os.path.join(renders_dir, f"figure_3_structure2d_{compound_id}.png")
    if os.path.exists(png_path) and not rerun:
        print(f"Reusing existing 2D structure render for {compound_id}: {png_path}")
        return png_path
    mol = Chem.MolFromSmiles(smiles)
    rdCoordGen.AddCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(500, 500)
    o = d.drawOptions()
    # True (opaque white), not the skill's own transparent canvas - _draw_compound_pose_row now
    # autocrops this raster and picks a number-label position by detecting non-white "ink"
    # (same convention as figure_4_plot.py's _mol_image/autocrop_to_content); a transparent
    # margin reads as black rather than white to that detection, breaking both. Same reasoning as
    # figure_4_plot.py's own _mol_image, which found no visible/measurable difference between a
    # transparent and opaque canvas on a white page.
    o.clearBackground = True
    # Back to RDKit's own default padding (0.05 per side, molecule fills 1-2*0.05=0.9 of the
    # canvas) - the earlier 0.1625 (shrunk to 75% size) was reverted per request ("my bad").
    o.padding = 0.05
    o.bondLineWidth = 3
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    d.WriteDrawingText(png_path)
    return png_path


def autocrop_to_content(img, padding_frac=0.05, background_frac=0.98):
    """Crop a rendered raster to its non-white bounding box (plus a small margin) - same helper
    as figure_1_plot.py's/figure_4_plot.py's own autocrop_to_content (duplicated per-file, same
    convention as this project's other small per-figure helpers)."""
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


# Panel d/e's structure column - shrunk to leave room for a compound-number label beside it
# (user request), styled the same way as figure_4 panel e's own numbered compound cards
# (figure_4_plot.py's MOL_IMAGE_FILL_FRAC/_mol_image_zoom/_label_xy_for_mol), ported here rather
# than imported since figure_N_plot.py files don't cross-import (see e.g. autocrop_to_content's
# own per-file duplication above).
STRUCT_IMAGE_FILL_FRAC = 0.9
# _draw_compound_pose_row's structure column is column 0 of plot_compound_pose_panel's own
# fig.add_gridspec(1, n_cols, wspace=POSE_ROW_WSPACE) - kept as a named constant here (rather than
# reading it back off the GridSpec) so the slot-size math below stays a plain, auditable formula.
POSE_ROW_WSPACE = 0.05
# Matches plot_compound_pose_panel's own fig.subplots_adjust(left=0.002, right=0.998, top=0.998,
# bottom=0.002) - symmetric margin on all 4 sides.
POSE_ROW_MARGIN = 0.002


def _pose_row_slot_size_in(fig_size, n_cols):
    """Physical (width, height) in inches of one column-slot within plot_compound_pose_panel's
    1 x n_cols row - same analytical approach as figure_4_plot.py's _slot_size_in, adapted for a
    single row (no hspace term needed)."""
    usable_w = fig_size[0] * (1 - 2 * POSE_ROW_MARGIN)
    usable_h = fig_size[1] * (1 - 2 * POSE_ROW_MARGIN)
    slot_w_in = usable_w / (n_cols + (n_cols - 1) * POSE_ROW_WSPACE)
    return slot_w_in, usable_h


def _struct_image_zoom(slot_w_in, slot_h_in, raster_w, raster_h):
    """OffsetImage zoom filling STRUCT_IMAGE_FILL_FRAC of the tighter dimension of a
    slot_w_in x slot_h_in slot - verbatim port of figure_4_plot.py's _mol_image_zoom."""
    return STRUCT_IMAGE_FILL_FRAC * min(slot_w_in * 72 / raster_w, slot_h_in * 72 / raster_h)


# Compound-number label style - ported from figure_4_plot.py's own LABEL_FONT/LABEL_REACH_*/
# LABEL_FOOTPRINT_MARGIN (user request: "bold Arial", same as figure_4 panel e's own compound
# numbering).
STRUCT_LABEL_FONT = {"family": "Arial", "fontweight": "bold", "fontsize": 7}
STRUCT_LABEL_REACH_MIN = 0.15
STRUCT_LABEL_REACH_MAX = 0.95
STRUCT_LABEL_REACH_STEP = 0.05
STRUCT_LABEL_FOOTPRINT_MARGIN = 2.0


def _struct_text_half_size_in(s):
    """Rendered (half_width, half_height) of `s` in STRUCT_LABEL_FONT, in inches - verbatim port
    of figure_4_plot.py's _text_half_size_in."""
    fig = plt.figure()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    t = fig.text(0, 0, s, **STRUCT_LABEL_FONT)
    fig.canvas.draw()
    bbox = t.get_window_extent(renderer=renderer)
    plt.close(fig)
    return (STRUCT_LABEL_FOOTPRINT_MARGIN * bbox.width / fig.dpi / 2,
            STRUCT_LABEL_FOOTPRINT_MARGIN * bbox.height / fig.dpi / 2)


def _label_xy_for_struct(img, mol_zoom, slot_w_in, slot_h_in, label_text):
    """Axes-fraction (x, y) placement for the compound-number label next to (not on top of) the
    centered structure - ported from figure_4_plot.py's _label_xy_for_mol, dropping that
    function's corner-grid avoidance (panel d/e's structure cell has no corner grids, only the
    molecule's own ink to avoid). Splits the raster into 2x2 quadrants, picks the emptiest one as
    the label's general direction, then walks outward from STRUCT_LABEL_REACH_MIN to
    STRUCT_LABEL_REACH_MAX until the label's own measured footprint is ink-free."""
    h, w = img.shape[:2]
    max_val = 1.0 if np.issubdtype(img.dtype, np.floating) else 255
    gray = img[..., :3].mean(axis=2)
    ink = gray < 0.98 * max_val
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

    half_w_in, half_h_in = _struct_text_half_size_in(label_text)
    half_x_frac = half_w_in / slot_w_in
    half_y_frac = half_h_in / slot_h_in
    half_px_x = max(1, round(half_x_frac / disp_w_frac * w))
    half_px_y = max(1, round(half_y_frac / disp_h_frac * h))

    def _xy_at(reach):
        return (0.5 + x_sign * reach * disp_w_frac / 2,
                0.5 + y_sign * reach * disp_h_frac / 2)

    def _clear_of_ink(reach):
        col = round(mid_x + x_sign * reach * mid_x)
        row = round(mid_y - y_sign * reach * mid_y)
        r0, r1 = max(0, row - half_px_y), min(h, row + half_px_y)
        c0, c1 = max(0, col - half_px_x), min(w, col + half_px_x)
        return not ink[r0:r1, c0:c1].any()

    reach = STRUCT_LABEL_REACH_MIN
    while reach <= STRUCT_LABEL_REACH_MAX:
        if _clear_of_ink(reach):
            return _xy_at(reach)
        reach += STRUCT_LABEL_REACH_STEP
    return _xy_at(STRUCT_LABEL_REACH_MAX)


# Hand-tuned (axes-fraction) nudge on top of _label_xy_for_struct's own auto-placed position, per
# compound_rank - user-requested fine adjustment for these 2 specific showcase compounds (now both
# rows within panel d: "1" further right/down, "2" slightly further up), same per-compound
# hand-tuning convention as figure_1's MY_VIEWS / figure_4's VIEWS_DIR hand-tuned camera views.
STRUCT_LABEL_NUDGE = {
    1: (0.06, -0.06),
    2: (-0.04, 0.12),
}


def _corner_badge(sub_ax, dot_color, label, loc, anchor_point, markersize=stylia.MARKERSIZE):
    # Reuses figure_1_plot.py's colored-circle Line2D-handle-as-legend trick (already used
    # project-wide for gene color coding), called multiple times per axes via add_artist() so each
    # corner badge survives the next call instead of replacing it (ax.legend() alone only keeps the
    # most recent one). anchor_point is in DATA (pixel) coordinates, not axes-fraction - with
    # adjustable="datalim" the axes BOX can be letterboxed around the actual image (blank padding on
    # the sides), so an axes-fraction anchor drifts off the visible image; pixel coordinates always
    # track the image's own corners regardless of that padding.
    handle = [Line2D([0], [0], marker="o", color="w", label=label,
                      markerfacecolor=dot_color, markeredgecolor="black", markeredgewidth=0.5,
                      markersize=markersize)]
    legend = sub_ax.legend(handles=handle, loc=loc, bbox_to_anchor=anchor_point,
                            bbox_transform=sub_ax.transData, frameon=True,
                            framealpha=0.85, edgecolor="black", fontsize=stylia.FONTSIZE_SMALL,
                            handletextpad=0.3, borderpad=0.3, labelspacing=0)
    sub_ax.add_artist(legend)


def _draw_compound_pose_row(fig, row_subgs, compound_rank, rerun=False):
    """Draws one compound's 2D structure + one docking-pose render per hit gene into the given
    1xN row_subgs slice (already positioned by the caller - this function only fills it in)."""
    # keep_default_na=False - otherwise pandas silently reads the literal string "NA" in
    # interpro_categories (unused here but present in the shared CSV schema) back as a real NaN.
    pockets = pd.read_csv(os.path.join(data_dir, "figure_3_top_avg_score_compounds_pockets.csv"),
                           keep_default_na=False)
    pockets = pockets[pockets["compound_rank"] == compound_rank].sort_values("gene").reset_index(drop=True)
    compound_id = pockets["compound_id"].iloc[0]

    with open(os.path.join(root, "..", "..", "output", "plots", "figure_1", "color_mapping.json")) as f:
        gene_to_color = json.load(f)["gene_to_color"]

    # Column 0: 2D structure - no PyMOL badges/score (they don't apply to the whole-compound cell),
    # no border, and no compound_id title, per request. Drawn smaller than the cell (via
    # OffsetImage/AnnotationBbox at STRUCT_IMAGE_FILL_FRAC, not a plain cell-filling imshow), with
    # a "1"/"2" (compound_rank) label placed next to it clear of the molecule's own ink - same
    # OffsetImage+ink-avoiding-label mechanism as figure_4 panel e's own numbered compound cards
    # (user request).
    smiles = lookup_smiles([compound_id], "REAL")[compound_id]
    struct_ax = stylize(fig.add_subplot(row_subgs[0, 0]))
    struct_ax.axis("off")
    struct_ax.set_xlim(0, 1)
    struct_ax.set_ylim(0, 1)
    stylia.label(struct_ax, xlabel="", ylabel="")

    struct_img = autocrop_to_content(mpimg.imread(render_2d_structure(compound_id, smiles, rerun=rerun)))
    slot_w_in, slot_h_in = _pose_row_slot_size_in(fig.get_size_inches(), row_subgs.ncols)
    raster_h, raster_w = struct_img.shape[:2]
    struct_zoom = _struct_image_zoom(slot_w_in, slot_h_in, raster_w, raster_h)

    struct_imagebox = OffsetImage(struct_img, zoom=struct_zoom)
    struct_ab = AnnotationBbox(struct_imagebox, (0.5, 0.5), xycoords="axes fraction",
                                frameon=False, box_alignment=(0.5, 0.5), pad=0,
                                annotation_clip=False, zorder=1)
    struct_ax.add_artist(struct_ab)

    label_text = str(compound_rank)
    label_x, label_y = _label_xy_for_struct(struct_img, struct_zoom, slot_w_in, slot_h_in, label_text)
    nudge_x, nudge_y = STRUCT_LABEL_NUDGE.get(compound_rank, (0.0, 0.0))
    label_x, label_y = label_x + nudge_x, label_y + nudge_y
    struct_ax.text(label_x, label_y, label_text, transform=struct_ax.transAxes,
                    ha="center", va="center", zorder=2, **STRUCT_LABEL_FONT)

    # Columns 1-N: one docking-pose render per hit gene (border, gene-name badge, docking-score
    # text) for visual consistency across the whole figure.
    for i, row in pockets.iterrows():
        sub_ax = stylize(fig.add_subplot(row_subgs[0, i + 1]))
        img = mpimg.imread(render_showcase_pocket(row, rerun=rerun))

        sub_ax.imshow(img)
        sub_ax.axis("off")
        sub_ax.set_aspect("equal", adjustable="datalim")

        img_h, img_w = img.shape[:2]
        # clip_on=False - each cell is very slightly taller than wide, so
        # set_aspect("equal", adjustable="datalim") expands ylim (blank margin above/below the
        # image) but leaves xlim locked exactly to the image's own edge; with the default
        # clip_on=True the left/right edges of this border sit exactly ON that (unmargined) x
        # clip boundary and lose their outer half to clipping, while the top/bottom edges (safely
        # inside the y margin) render at full width - confirmed by rasterizing the saved PDF and
        # measuring: left/right rendered at half the linewidth of top/bottom. clip_on=False
        # renders the full stroke on all 4 sides, confirmed by the same measurement.
        sub_ax.add_patch(plt.Rectangle((-0.5, -0.5), img_w, img_h, fill=False,
                                        edgecolor="black", linewidth=stylia.LINEWIDTH, clip_on=False))

        gene_label = row["gene"] + ("" if row["is_true_best_pocket"] else "*")
        _corner_badge(sub_ax, gene_to_color[row["gene"]], gene_label,
                      loc="upper left", anchor_point=(0, 0), markersize=stylia.MARKERSIZE * 0.5)

        sub_ax.text(img_w / 2, img_h * 0.97, f"Docking score: {row['score']:.2f}",
                    ha="center", va="bottom", fontsize=stylia.FONTSIZE, color="black",
                    bbox=dict(facecolor="white", edgecolor="black", alpha=0.85, boxstyle="square,pad=0.25"))

        stylia.label(sub_ax, xlabel="", ylabel="")


def plot_compound_pose_panel(letter, size, compound_ranks, padding=0.0, rerun=False):
    """Panel d (previously split across d and e, merged into one panel per request): one
    standalone figure with one stacked row per entry in `compound_ranks` (2D structure + one
    docking-pose render per hit gene), each row for TOP_AVG_SCORE_COMPOUND_IDS[compound_rank - 1]
    - see _draw_compound_pose_row."""
    fig = plt.figure(figsize=size)
    fig.patch.set_facecolor("white")
    # Every axes here uses set_aspect("equal", adjustable="datalim") to keep its square image
    # undistorted, which (like reserve_top_header's manual ax.set_position() in panel f) makes
    # matplotlib back off tight_layout() for this whole figure - so the stock default margins
    # (left=0.125, right=0.9, top=0.88, bottom=0.11) would otherwise govern the canvas, even
    # though nothing here (image axes with axis("off"), badges/score text drawn INSIDE each
    # image's own bounds) needs any margin at all.
    fig.subplots_adjust(left=0.002, right=0.998, top=0.998, bottom=0.002)
    pockets = pd.read_csv(os.path.join(data_dir, "figure_3_top_avg_score_compounds_pockets.csv"),
                           keep_default_na=False)
    outer_gs = fig.add_gridspec(len(compound_ranks), 1, hspace=0.01)
    for i, compound_rank in enumerate(compound_ranks):
        n_cols = 1 + (pockets["compound_rank"] == compound_rank).sum()  # 2D structure + one pose per hit gene
        row_gs = outer_gs[i, 0].subgridspec(1, n_cols, wspace=0.05)
        _draw_compound_pose_row(fig, row_gs, compound_rank, rerun=rerun)
    save_panel(fig, letter, padding)


def main(rerun=False, subpanels=None):
    # Each panel is its own standalone figure/PDF now (see module docstring) - sizes come from
    # the user-authored output/plots/figure_3/panel_layout.csv, not a shared mosaic. subpanels
    # restricts which of PANEL_LETTERS actually get (re)generated - see --subpanels.
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
        fig, ax = plt.subplots(figsize=sizes["a"])
        fig.patch.set_facecolor("white")
        stylize(ax)
        plot_protein_hit_counts(ax)
        save_panel(fig, "a", paddings["a"])

    if "b" in subpanels:
        fig, ax = plt.subplots(figsize=sizes["b"])
        fig.patch.set_facecolor("white")
        stylize(ax)
        plot_circos_overlap(ax)
        save_panel(fig, "b", paddings["b"])

    if "e" in subpanels:
        # Tier grid, labeled "e" (was "f", was "c") - relabeled to "e" (freed up when d/e merged
        # into one panel "d") so the panel letters read contiguously a-e again.
        fig, ax = plt.subplots(figsize=sizes["e"])
        fig.patch.set_facecolor("white")
        stylize(ax)
        plot_tier_grid_placeholder(ax)
        save_panel(fig, "e", paddings["e"])

    if "d" in subpanels:
        # Merged panel d/e (per request) - both showcase compounds now stack as 2 rows within one
        # panel, instead of 2 separate analogous panels.
        plot_compound_pose_panel("d", sizes["d"], compound_ranks=[1, 2], padding=paddings["d"], rerun=rerun)

    if "c" in subpanels:
        # P2Rank probability curve + 4 domain bands, labeled "c" (was "f") - see the swap note
        # above. Outer gridspec splits P2Rank prob. from the 4-domain-row block (ROW_BLOCK_HSPACE
        # gap, ROW_BLOCK_HEIGHT_RATIOS split); an inner subgridspec then packs the 4 domain rows
        # with their own tighter DOMAIN_ROWS_HSPACE gap, per request - replaces the narrow-box
        # side-by-side-column layout this used under panel f.
        fig = plt.figure(figsize=sizes["c"])
        fig.patch.set_facecolor("white")
        outer_gs = fig.add_gridspec(2, 1, height_ratios=[ROW_BLOCK_HEIGHT_RATIOS[0], sum(ROW_BLOCK_HEIGHT_RATIOS[1:])],
                                     hspace=ROW_BLOCK_HSPACE)
        plot_pocket_scores_row(stylize(fig.add_subplot(outer_gs[0, 0])))
        domain_gs = outer_gs[1, 0].subgridspec(len(DOMAIN_STRIP_COLUMNS), 1, hspace=DOMAIN_ROWS_HSPACE)
        for i, (column, title, color_name) in enumerate(DOMAIN_STRIP_COLUMNS):
            ax = stylize(fig.add_subplot(domain_gs[i, 0]))
            plot_domain_strip_row(ax, column=column, title=title, color_name=color_name)
        # save_panel()'s plt.tight_layout() warns "not compatible with tight_layout" for this
        # figure (the nested subgridspec) and falls back to matplotlib's stock default margins
        # (figure.subplot.left/right, normally 0.125/0.9) instead of shrinking to fit - leaving a
        # real blank strip on the right (no right-side decoration needs it), and clipping the
        # longest row labels ("P2Rank prob.", "tRNA binding", "Anticodon binding") on the left.
        # Force both margins explicitly. The left margin is a FIXED ABSOLUTE width
        # (ROW_LABEL_LEFT_MARGIN_IN), converted to the fraction-of-figure-width subplots_adjust
        # wants, rather than a fixed fraction - a fixed fraction's absolute clearance shrinks
        # whenever this panel's own delta_x does (panel_layout.csv is user-edited and has been
        # changing), reintroducing the clipping; an absolute width stays correct regardless.
        fig.subplots_adjust(left=ROW_LABEL_LEFT_MARGIN_IN / sizes["c"][0], right=0.995)
        save_panel(fig, "c", paddings["c"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rerun", action="store_true", default=False,
                         help="Force PyMOL/RDKit re-rendering of panels d/e's showcase-compound "
                              "PNGs, even if they already exist")
    parser.add_argument("--subpanels", type=str, default=None,
                         help="Comma-separated subset of panels to generate (e.g. 'a,c'), from "
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
