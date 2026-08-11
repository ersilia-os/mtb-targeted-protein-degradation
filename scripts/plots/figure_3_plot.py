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
their compound count at cutoff -11, descending.

Panel b is a chord (Circos) diagram of pairwise protein hit overlap, at CIRCOS_CUTOFF, within that
same selected set - following notebooks/46_docking_exploration_IIIa_colors.ipynb's precedent
(pycirclize's Circos.chord_diagram, rasterized via FigureCanvasAgg into a stylia axes), colored with
figure_1's own gene->color mapping so protein colors stay consistent across every figure. Genes with
zero hits at the cutoff are dropped outright (a chord diagram sector needs a nonzero span); among
the rest, low-degree genes get a blank (padded-space) label instead of their real name - same
declutter heuristic as that notebook - since their sector is too small for real text without
overlapping its neighbors.

Panel c (plot_tier_grid_placeholder) is a gray rectangle grid, TIER_GRID_ROW_LABELS rows
(Druggability, Docking scores HL, Docking Scores REAL, Exp. tractability, Novelty) x 21 columns
(all genes alphabetically) - the row/column structure and row labels are real (per request), but
every cell is still an unfilled gray placeholder - the actual per-gene/per-row values (e.g. from
the 5-tier vulnerability ranking in output/protein_prioritization/final_results.tsv, or wherever
each row's values end up coming from) haven't been decided yet.

Panels d and e are each one full-width row for a different showcase compound -
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
import math
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

PANEL_LETTERS = ["a", "b", "c", "d", "e", "f"]
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
YLIM_MAX = 425


def plot_protein_hit_counts(ax):
    counts = pd.read_csv(os.path.join(data_dir, "figure_3_selected_set_protein_hit_counts.csv"))
    counts = counts.sort_values(f"count_cutoff{abs(SORT_BY_CUTOFF)}", ascending=False).reset_index(drop=True)

    palette = stylia.CategoricalPalette("npg").get(len(CUTOFFS))
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
    ax.legend(title="Docking score", loc="upper center", bbox_to_anchor=(0.5, 1), fontsize=stylia.FONTSIZE_SMALL, framealpha=1, ncol=5)
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
        link_kws=dict(ec="black", lw=1, alpha=1),
    )
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


def plot_domain_strip(ax, column, title, color_name):
    """One band of panel f: a thin colored cell per pocket (all 276, same pocket_rank y-order as
    panel e - sorted by P2Rank probability descending, rank 1 at the top), in this band's own
    color where `column` is True, white otherwise - figure_3_calculations.py's
    compute_pocket_scores, itself from output/77_pocket_annotation/
    pocket_detection_interpro_updated.csv. Vertical orientation (pocket_rank on y, not x) to match
    the probability curve's own tall/narrow column, in the same panel f cell."""
    scores = pd.read_csv(os.path.join(data_dir, "figure_3_pocket_scores.csv"))
    nc = stylia.NamedColors()
    present_color = getattr(nc, color_name)
    colors = [present_color if v else nc.white for v in scores[column]]
    ax.barh(scores["pocket_rank"], 1, height=1.0, color=colors, edgecolor="none", zorder=2)
    ax.set_ylim(scores["pocket_rank"].max(), scores["pocket_rank"].min())
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    # Real xlabel (not free-floating ax.text), moved to the top and rotated 90 degrees - these
    # columns are packed edge-to-edge (RIGHT_BLOCK_WSPACE), too narrow for even one word of
    # horizontal text without overlapping the neighboring column's title, and top placement
    # matches plot_pocket_scores' own xlabel in the same panel f cell.
    stylia.label(ax, xlabel=title, ylabel="")
    ax.xaxis.set_label_position("top")
    ax.xaxis.label.set_rotation(90)
    # See plot_pocket_scores' identical block for why ha="center"/va="bottom" (not the reverse).
    ax.xaxis.label.set_ha("center")
    ax.xaxis.label.set_va("bottom")
    # Stops matplotlib's auto label positioning from overwriting this alignment on every draw.
    # LABEL_GAP leaves real empty space between the label's own bottom edge (va="bottom" anchors
    # it here) and the axes' top edge, instead of the two sitting flush against each other.
    ax.xaxis.set_label_coords(0.5, 1.0 + LABEL_GAP, transform=ax.transAxes)


# Row label -> figure_3_gene_summary_stats.json field -> how to grade it. "continuous_high_better"/
# "continuous_low_better" rows are graded by equal terciles across the 21 genes (7 green/7 amber/7
# red, ranked by that field) - a data-driven split with no invented absolute cutoff, per request
# (AskUserQuestion) over hand-picked thresholds. "boolean" rows (Novelty, Exp. tractability) are
# almost entirely NaN right now (only ileS/glyS have a real value - see figure_3_calculations.py's
# NOVELTY_OVERRIDES/EXPERIMENTAL_TRACTABILITY_OVERRIDES) - NaN stays the gray placeholder color,
# and True/False map to green/red respectively, per request.
TIER_ROW_FIELDS = [
    ("Druggability", "max_p2rank_prob", "continuous_high_better"),
    ("Docking scores HL", "best_hl_docking_score", "continuous_low_better"),
    ("Docking Scores REAL", "best_real10b_docking_score", "continuous_low_better"),
    ("Exp. tractability", "experimental_tractability", "boolean"),
    ("Novelty", "novelty", "boolean"),
]
TIER_GRID_ROW_LABELS = [label for label, _, _ in TIER_ROW_FIELDS]
TIER_GRID_FACECOLOR = "lightgray"
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


def _boolean_color(value):
    if isinstance(value, float) and math.isnan(value):
        return None  # no data yet - stays TIER_GRID_FACECOLOR
    return TIER_GRID_GREEN if value else TIER_GRID_RED


def plot_tier_grid_placeholder(ax):
    """TIER_GRID_ROW_LABELS rows x 21 (genes, alphabetical) columns grid of rectangles, colored
    from figure_3_calculations.py's figure_3_gene_summary_stats.json - see TIER_ROW_FIELDS for the
    row->field mapping and grading rule. Cells with no data yet (all of Novelty/Exp. tractability
    except ileS/glyS) stay the gray placeholder color rather than being guessed at."""
    with open(os.path.join(root, "..", "..", "output", "plots", "figure_1", "color_mapping.json")) as f:
        genes = sorted(json.load(f)["gene_to_color"].keys())
    with open(os.path.join(data_dir, "figure_3_gene_summary_stats.json")) as f:
        stats = json.load(f)

    nc = stylia.NamedColors()
    row_colors = []
    for _, field, kind in TIER_ROW_FIELDS:
        if kind == "boolean":
            row_colors.append({gene: _boolean_color(stats[gene][field]) for gene in genes})
        else:
            values_by_gene = {gene: stats[gene][field] for gene in genes}
            row_colors.append(_tercile_colors(values_by_gene, higher_is_better=(kind == "continuous_high_better")))

    for row, colors in enumerate(row_colors):
        for col, gene in enumerate(genes):
            color_name = colors[gene]
            facecolor = TIER_GRID_FACECOLOR if color_name is None else getattr(nc, color_name)
            ax.add_patch(plt.Rectangle((col, row), 1, 1, facecolor=facecolor,
                                        edgecolor="white", linewidth=stylia.LINEWIDTH))

    ax.set_xlim(0, len(genes))
    ax.set_ylim(0, len(TIER_GRID_ROW_LABELS))
    ax.invert_yaxis()  # row A at the top, matching standard matrix reading order

    ax.set_xticks([i + 0.5 for i in range(len(genes))])
    ax.set_xticklabels(genes, rotation=90, fontsize=stylia.FONTSIZE_SMALL)
    ax.set_yticks([i + 0.5 for i in range(len(TIER_GRID_ROW_LABELS))])
    ax.set_yticklabels(TIER_GRID_ROW_LABELS)
    stylia.label(ax, xlabel="", ylabel="")
    # No outer axes border, per request - only the individual cells' own white edges (set above)
    # remain visible.
    for spine in ax.spines.values():
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
    # style matches the molecule-auditing skill's svg_for() (CoordGen coordinates, transparent
    # background), via the Cairo PNG backend instead of SVG to embed like every other raster image
    # in this file - padding/bondLineWidth deviate from that skill's exact values per request
    # (smaller molecule within the same canvas, thicker bonds).
    png_path = os.path.join(renders_dir, f"figure_3_structure2d_{compound_id}.png")
    if os.path.exists(png_path) and not rerun:
        print(f"Reusing existing 2D structure render for {compound_id}: {png_path}")
        return png_path
    mol = Chem.MolFromSmiles(smiles)
    rdCoordGen.AddCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(500, 500)
    o = d.drawOptions()
    o.clearBackground = False
    # Back to RDKit's own default padding (0.05 per side, molecule fills 1-2*0.05=0.9 of the
    # canvas) - the earlier 0.1625 (shrunk to 75% size) was reverted per request ("my bad").
    o.padding = 0.05
    o.bondLineWidth = 3
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    d.WriteDrawingText(png_path)
    return png_path


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
    # no border, and no compound_id title, per request.
    smiles = lookup_smiles([compound_id], "REAL")[compound_id]
    struct_ax = stylize(fig.add_subplot(row_subgs[0, 0]))
    struct_img = mpimg.imread(render_2d_structure(compound_id, smiles, rerun=rerun))
    struct_ax.imshow(struct_img)
    struct_ax.axis("off")
    struct_ax.set_aspect("equal", adjustable="datalim")
    stylia.label(struct_ax, xlabel="", ylabel="")

    # Columns 1-N: one docking-pose render per hit gene (border, gene-name badge, docking-score
    # text) for visual consistency across the whole figure.
    for i, row in pockets.iterrows():
        sub_ax = stylize(fig.add_subplot(row_subgs[0, i + 1]))
        img = mpimg.imread(render_showcase_pocket(row, rerun=rerun))

        sub_ax.imshow(img)
        sub_ax.axis("off")
        sub_ax.set_aspect("equal", adjustable="datalim")

        img_h, img_w = img.shape[:2]
        sub_ax.add_patch(plt.Rectangle((-0.5, -0.5), img_w, img_h, fill=False,
                                        edgecolor="black", linewidth=stylia.LINEWIDTH))

        gene_label = row["gene"] + ("" if row["is_true_best_pocket"] else "*")
        _corner_badge(sub_ax, gene_to_color[row["gene"]], gene_label,
                      loc="upper left", anchor_point=(0, 0), markersize=stylia.MARKERSIZE * 0.5)

        sub_ax.text(img_w / 2, img_h * 0.97, f"Docking score: {row['score']:.2f}",
                    ha="center", va="bottom", fontsize=stylia.FONTSIZE, color="black",
                    bbox=dict(facecolor="white", edgecolor="black", alpha=0.85, boxstyle="square,pad=0.25"))

        stylia.label(sub_ax, xlabel="", ylabel="")


def plot_compound_pose_panel(letter, size, compound_rank, padding=0.0, rerun=False):
    """Panel d or e: one standalone figure (2D structure + one docking-pose render per hit gene)
    for TOP_AVG_SCORE_COMPOUND_IDS[compound_rank - 1] - see _draw_compound_pose_row."""
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
    n_cols = 1 + (pockets["compound_rank"] == compound_rank).sum()  # 2D structure + one pose per hit gene
    gs = fig.add_gridspec(1, n_cols, wspace=0.05)
    _draw_compound_pose_row(fig, gs, compound_rank, rerun=rerun)
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

    if "c" in subpanels:
        fig, ax = plt.subplots(figsize=sizes["c"])
        fig.patch.set_facecolor("white")
        stylize(ax)
        plot_tier_grid_placeholder(ax)
        save_panel(fig, "c", paddings["c"])

    if "d" in subpanels:
        plot_compound_pose_panel("d", sizes["d"], compound_rank=1, padding=paddings["d"], rerun=rerun)
    if "e" in subpanels:
        plot_compound_pose_panel("e", sizes["e"], compound_rank=2, padding=paddings["e"], rerun=rerun)

    if "f" in subpanels:
        # Panel f: probability curve + 4 domain bands, packed into one figure as 5 tight columns
        # (RIGHT_BLOCK_WSPACE). RIGHT_BLOCK_WIDTH_RATIOS gives the probability curve (a real
        # value, not a present/absent band like the other 4) more room to read, per request.
        fig = plt.figure(figsize=sizes["f"])
        fig.patch.set_facecolor("white")
        # Nothing is anchored at the bottom here (plot_pocket_scores moves ticks to the top,
        # plot_domain_strip has no ticks at all) - matplotlib's stock default bottom margin
        # (rcParams["figure.subplot.bottom"], normally 0.11) would otherwise go unused, since
        # reserve_top_header's manual ax.set_position() calls make tight_layout() back off for
        # this whole figure (see the "not compatible with tight_layout" warning) and leave stock
        # defaults in charge of every edge it doesn't itself touch.
        fig.subplots_adjust(bottom=0.02)
        right_block_gs = fig.add_gridspec(1, 1 + len(DOMAIN_STRIP_COLUMNS),
                                           width_ratios=RIGHT_BLOCK_WIDTH_RATIOS, wspace=RIGHT_BLOCK_WSPACE)
        axes = [stylize(fig.add_subplot(right_block_gs[0, i])) for i in range(1 + len(DOMAIN_STRIP_COLUMNS))]
        fig.canvas.draw()  # positions must be final before reserve_top_header reads them via get_position()
        plot_pocket_scores(reserve_top_header(axes[0], RIGHT_BLOCK_HEADER_FRAC))
        for i, (column, title, color_name) in enumerate(DOMAIN_STRIP_COLUMNS):
            ax = reserve_top_header(axes[1 + i], RIGHT_BLOCK_HEADER_FRAC)
            plot_domain_strip(ax, column=column, title=title, color_name=color_name)
        save_panel(fig, "f", paddings["f"])


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
