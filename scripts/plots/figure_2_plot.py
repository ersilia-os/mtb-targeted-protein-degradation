import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

import json
import re
import subprocess
import tempfile
from collections import defaultdict

import numpy as np
import pandas as pd
import stylia
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
from PIL import Image, ImageDraw, ImageFont
from scipy.stats import gaussian_kde
from stylia.config import get_fontsize_big, get_size

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

plots_dir = os.path.join(root, "..", "..", "plots", "figure_2")

PROPERTIES = ["MW", "cLogP", "TPSA", "RotBonds", "HBA", "AromaticRings", "QED"]
DISCRETE_PROPS = {"RotBonds", "HBA", "AromaticRings"}
PROP_XLIMS = {
    "MW":            (200, 500),
    "cLogP":         (-3, 6),
    "TPSA":          (0, 160),
    "RotBonds":      (0, 11),
    "HBA":           (0, 11),
    "AromaticRings": (0, 5),
    "QED":           (0.2, 1),
}
PROP_XLABELS = {"AromaticRings": "Aromatic rings", "RotBonds": "Rot. bonds", "HBA": "HBA"}
PROP_XTICKS = {
    "MW":            [300, 400],
    "cLogP":         [-2, 0, 2, 4],
    "TPSA":          [50, 100],
    "AromaticRings": [0, 2, 4],
    "QED":           [0.4, 0.6, 0.8],
}
PROP_BINSIZE = {"RotBonds": 2, "HBA": 2}

# (label, csv filename in plots/figure_2/, NamedColors name)
LIBRARIES = [
    ("HL",       "figure_2a_HL.csv",       "crimson"),
    ("REAL 10M", "figure_2a_REAL10M.csv",  "turquoise"),
    ("REAL 10B", "figure_2b_REAL10B.csv",  "amber"),
]

# Relative vertical sizes of panel A (scheme.svg), B (property row), C (gene grid rows), D (new,
# still-empty row). B's share was halved from the original 3:1:5 (A:B:C) design to shrink just
# B's plot area - per request, this should shrink the KDE/bar drawing area only, not the tick
# label text below it (fontsize is independent of axes-box height, so that already holds without
# extra work) - while keeping C's rows at their original absolute size (see
# UNIT_HEIGHT_FRACTION in main()). D's ratio just matches a single C row for now - a placeholder,
# easy to retune once D has real content.
ROW_HEIGHT_RATIOS = [1.2, 5, 5, 5, 5]

# Panel A (scheme.pdf) is a hand-authored Inkscape schematic, exported to PDF and rasterized via
# Poppler (see pdf_page_size_mm()/render_pdf_to_png()) rather than embedded as vector SVG - it
# can't use stylia.label(abc=...) directly, so its "A" label is drawn by hand onto the composited
# raster in merge_with_scheme(). Panels B/C are embedded natively via stylia.label(..., abc=...).
# Label A's own x is NOT fixed here - it's set dynamically to match wherever matplotlib actually
# put B/C's title (see main()/merge_with_scheme()), so all three stay aligned regardless of the
# distributions block's own left margin (y-axis label/tick width). Only its y (top) margin is a
# fixed constant, since there's no equivalent matplotlib position to match it against.
LABEL_FONT_PATH = "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"
LABEL_SIZE_PT = 8
# One shared margin: label A's own offset from the schematic's top edge, AND the breathing room
# left on both sides of the B+C block (it otherwise fills the full scheme width - see panel_bc_x
# in merge_with_scheme() - which put its content flush against both edges with zero margin).
# Keeping both uses on the same constant is what makes the top/left spacing read as one
# consistent margin throughout the figure, rather than two independently-tuned numbers.
MARGIN_MM = 5

# scheme.pdf/distributions are rasterized (via Poppler and matplotlib's own Agg backend,
# respectively) and composited with plain PIL pasting - both engines rasterize their own text
# themselves, so neither goes through the svgutils-compose text-anchor bug that a raw-SVG-embed
# approach hit (see docs/ for the diagnosis: every <text> in scheme.svg shifted horizontally once
# embedded via svgutils.compose and re-rendered, while paths/shapes stayed pixel-exact).
MERGED_DPI = 600


def pdf_page_size_mm(pdf_path):
    """Parses `pdfinfo`'s "Page size: W x H pts" line and converts pt -> mm."""
    result = subprocess.run(["pdfinfo", pdf_path], check=True, capture_output=True, text=True)
    match = re.search(r"Page size:\s+([\d.]+)\s*x\s*([\d.]+)\s*pts", result.stdout)
    width_pt, height_pt = float(match.group(1)), float(match.group(2))
    return width_pt * 25.4 / 72, height_pt * 25.4 / 72


def render_pdf_to_png(pdf_path, dpi, out_path):
    """pdftoppm (Poppler) - the same renderer backing Inkscape's own PDF export/preview, so this
    reproduces exactly what Inkscape shows, unlike embedding scheme's SVG text via svgutils."""
    out_prefix = out_path[:-4] if out_path.endswith(".png") else out_path
    subprocess.run(
        ["pdftoppm", "-png", "-singlefile", "-r", str(dpi), pdf_path, out_prefix],
        check=True, capture_output=True,
    )


def content_bottom_mm(png_path, dpi):
    """Row of the bottom-most non-white pixel, in mm - used to find how much blank margin a
    hand-authored SVG has below its actual visual content, so panel spacing can be computed from
    real content edges instead of trusting each file's own (variable, editor-dependent) declared
    canvas size."""
    arr = np.array(Image.open(png_path).convert("L"))
    nonwhite_rows = np.where((arr < 250).any(axis=1))[0]
    if len(nonwhite_rows) == 0:
        return 0.0
    return nonwhite_rows[-1] / (dpi / 25.4)


def boost_abc_fontsize(ax, shrink_factor):
    """Panels B/C's titles get shrunk by shrink_factor when composited under the schematic in
    merge_with_scheme() (the whole B+C block is scaled down to fit the scheme's width minus
    margins) - boost the source fontsize here so they end up the same physical size as panel A's
    own (unscaled) label after that compositing shrink, instead of visibly smaller."""
    ax._left_title.set_fontsize(get_fontsize_big() / shrink_factor)


def expected_shrink_factor(scheme_width_mm):
    """How much merge_with_scheme() will shrink the distributions block's width (and everything
    in it, fonts included) to fit scheme.pdf's width minus MARGIN_MM on each side - computed here,
    before the matplotlib figure is even built, purely from scheme.pdf's page width and stylia's
    own configured figure width, so plot_property()/plot_gene_panel() can pre-compensate the B/C
    title fontsize (see boost_abc_fontsize()) instead of needing a second render pass."""
    target_width_mm = scheme_width_mm - 2 * MARGIN_MM
    natural_width_mm = get_size() * 25.4  # stylia figure width in inches, before any compositing scale
    return target_width_mm / natural_width_mm


def load_libraries():
    loaded = []
    for label, filename, color_name in LIBRARIES:
        df = pd.read_csv(os.path.join(plots_dir, filename))
        loaded.append((label, color_name, df))
    return loaded


def load_gene_colors():
    """Gene->color mapping from figure 1 - reused here so protein colors stay consistent across
    figures instead of a second palette."""
    with open(os.path.join(root, "..", "..", "plots", "figure_1", "color_mapping.json")) as f:
        return json.load(f)["gene_to_color"]


def plot_property(ax, prop, libraries, abc=None, shrink_factor=None):
    nc = stylia.NamedColors()
    lo, hi = PROP_XLIMS[prop]

    # Matplotlib draws axes left-to-right, and each axes' own opaque white background patch is
    # painted before its own artists - which silently covers anything from a PRECEDING axes that
    # overlaps into this one's screen region.
    ax.patch.set_visible(False)

    if prop in DISCRETE_PROPS:
        bin_size = PROP_BINSIZE.get(prop, 1)

        def binned(df):
            return (df[prop].dropna().astype(int) // bin_size) * bin_size

        all_vals = set()
        for _, _, df in libraries:
            all_vals |= set(binned(df))
        vals = sorted(all_vals)
        edges = vals + [vals[-1] + bin_size]
        for label, color_name, df in libraries:
            data = binned(df)
            freq = pd.Series(data).value_counts(normalize=True).reindex(vals, fill_value=0)
            color = nc.get(color_name)
            # stairs() draws one continuous step-outline per series, unfilled - so a line only
            # ever appears at a bar's outer edge (where it meets zero) or at a height change
            # between adjacent bars, never as a divider between two same-height neighbors.
            ax.stairs(freq.values, edges, baseline=0, color=color, linewidth=stylia.LINEWIDTH * 1.5,
                      label=f"{label} (n={len(df) // 1000}k)")
        ax.set_xlim(lo - 0.5, hi + 0.5)
        stylia.label(ax, xlabel="", ylabel="", title=PROP_XLABELS.get(prop, prop), abc=abc)
    else:
        for label, color_name, df in libraries:
            data = df[prop].dropna().values
            kde = gaussian_kde(data)
            x = np.linspace(lo, hi, 300)
            y = kde(x)
            color = nc.get(color_name)
            # Plain line, no fill_between() below the curve, per request.
            ax.plot(x, y, color=color, linewidth=stylia.LINEWIDTH * 1.5, label=f"{label} (n={len(df) // 1000}k)")
        ax.set_xlim(lo, hi)
        stylia.label(ax, xlabel="", ylabel="", title=PROP_XLABELS.get(prop, prop), abc=abc)
        if prop == "MW":
            # Density values here are naturally small (peak ~0.0125) - fewer, rounder ticks
            # read better than the default 4-decimal-heavy tick set.
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))

    if prop in PROP_XTICKS:
        ax.set_xticks(PROP_XTICKS[prop])

    # y-axis border shown (left spine) alongside the x-axis border (bottom spine), matching
    # panel C's box - but still no y tick marks/numbers, since the cross-panel y-scales aren't
    # meant to be compared.
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if abc is not None:
        boost_abc_fontsize(ax, shrink_factor)

    ax.grid(False)


def load_docking_percentiles():
    return pd.read_csv(os.path.join(plots_dir, "figure_2c_docking_percentiles.csv"))


def gene_docking_stats(docking_percentiles, gene):
    """Per-pocket median and p1 docking score for one gene - columns are a (stat, library)
    MultiIndex, pockets sorted by REAL 10B's p1 ascending (most negative/best first) so both
    stats share one x-axis order across all three libraries instead of each sorting
    independently. Empty DataFrame if the gene has no docking data yet (all genes but alaS, for
    now) - plot_gene_panel() falls back to its dummy-dot placeholder in that case."""
    gene_df = docking_percentiles[docking_percentiles["gene"] == gene]
    if gene_df.empty:
        return gene_df
    pivoted = gene_df.pivot(index="pocket", columns="library", values=["median", "p1"])
    return pivoted.sort_values(("p1", "REAL 10B"))


def load_docking_snapshot_index():
    """Ordered rows (rank, gene, library, pocket, compound, score, filename) for panel D's 7
    snapshot cells - figure_2_calculations.py's compute_docking_snapshots() writes this index
    alongside the PNGs since panel D shows a fixed top-7 by score (repeats allowed), not one slot
    per gene, so the same gene's snapshots need distinguishing beyond just its name."""
    index_path = os.path.join(plots_dir, "docking_snapshots", "index.csv")
    if not os.path.isfile(index_path):
        return pd.DataFrame()
    return pd.read_csv(index_path).sort_values("rank")


LIBRARY_DISPLAY_NAMES = {"HL": "Hit Locator", "REAL 10M": "REAL 10M", "REAL 10B": "REAL 10B"}


def plot_docking_snapshot(ax, rank, image_path, library, score):
    """One panel D cell: a best-compound PyMOL docking snapshot, titled with the panel's own
    plain-text number (outside the frame via loc="lower center", no color/marker), plus two more
    plain-text legends - which library (top center) and its docking score (bottom center)."""
    ax.imshow(np.array(Image.open(image_path)))
    ax.set_xticks([])
    ax.set_yticks([])

    number_handle = [Line2D([0], [0], linestyle="None", marker="None", label=str(rank))]
    number_legend = ax.legend(handles=number_handle, loc="lower center", frameon=False,
                               bbox_to_anchor=(0.5, 1.0), fontsize=stylia.FONTSIZE_BIG, handlelength=0,
                               handletextpad=0, borderaxespad=0)
    # Each subsequent ax.legend() call below would otherwise silently replace the previous one -
    # add_artist() keeps all but the last one on the axes at once.
    ax.add_artist(number_legend)

    # handlelength=0 collapses the (invisible) marker's reserved space so only the text shows.
    library_handle = [Line2D([0], [0], linestyle="None", marker="None",
                              label=LIBRARY_DISPLAY_NAMES.get(library, library))]
    library_legend = ax.legend(handles=library_handle, loc="upper center", frameon=True, framealpha=0.8,
                                edgecolor="black", fontsize=stylia.FONTSIZE, handlelength=0, handletextpad=0)
    ax.add_artist(library_legend)

    score_handle = [Line2D([0], [0], linestyle="None", marker="None", label=f"Dock. score: {score:.2f}")]
    ax.legend(handles=score_handle, loc="lower center", frameon=True, framealpha=0.8, edgecolor="black",
              fontsize=stylia.FONTSIZE, handlelength=0, handletextpad=0)


def plot_gene_panel(ax, gene, color, docking_df=None, abc=None, shrink_factor=None, show_yaxis=True,
                     snapshot_rows=None):
    """One cell of panel C's 3x7 gene grid, titled with a colored circle + gene name analogously
    to figure 1 panel A's per-structure gene legend (render_structure_panel()'s legend_handles
    pattern), reusing the same gene->color mapping. Plots each library's per-pocket median and p1
    docking score (x = pocket identity, unlabeled - only their relative order matters here) when
    docking_df is given; otherwise a dummy single-dot scatter placeholder (for a gene with no
    docking data at all) in the gene's own color. show_yaxis=False drops the "Docking score"
    ylabel and y tick numbers - the 3x7 grid only needs them on each row's first column."""
    # Matplotlib draws axes in creation order, and each axes' own opaque white background patch
    # is painted before its own artists - silently covering anything from a PRECEDING axes that
    # overlaps into this one's screen region (the same bug plot_property() hit with panel B's
    # legend). With C's rows pulled close together, a cross-reference label near the top of one
    # row's panel can visually reach into the row below's rectangle, drawn right after it.
    ax.patch.set_visible(False)
    if docking_df is not None and not docking_df.empty:
        nc = stylia.NamedColors()
        n = len(docking_df)
        x = list(range(n))
        # Reference structure is this gene's AlphaFold2 model_0 (AF2 has exactly one model per
        # protein, confirmed across all 21 genes - a single unambiguous structure per gene) -
        # give a thin black outline (on both dots and dashes, drawn slightly larger and behind
        # the normal library-colored markers) to every pocket from that structure, wherever they
        # land in the REAL 10B p1 sort order. Uniprot AC is parsed the same way
        # figure_2_calculations.py's pocket_percentiles() does.
        uniprot_ac = docking_df.index[0].split("_model_")[0].split("_")[-1]
        best_structure = f"alphafold2_{uniprot_ac}_model_0"
        highlight = [i for i, pocket in enumerate(docking_df.index)
                     if pocket.rsplit("_pocket_", 1)[0] == best_structure]
        non_highlight = [i for i in range(n) if i not in highlight]
        highlight_x = [x[i] for i in highlight]
        non_highlight_x = [x[i] for i in non_highlight]
        for label, _, color_name in LIBRARIES:
            lib_color = nc.get(color_name)
            # No connecting line between pockets - a plain scatter of per-pocket stats. Median:
            # horizontal dash. p1: filled dot, smaller for the non-highlighted majority and
            # bigger for the highlighted (winning-structure) pockets, so the highlight reads as
            # emphasis rather than just an outline.
            p1_size = stylia.MARKERSIZE * 0.25
            p1_size_highlight = stylia.MARKERSIZE * 1.6
            outline_lw = 0.4
            ax.scatter(highlight_x, docking_df[("median", label)].iloc[highlight], marker="_",
                       s=p1_size, linewidths=stylia.LINEWIDTH * 1.5 + outline_lw * 2, color="black",
                       zorder=2)
            ax.scatter(highlight_x, docking_df[("p1", label)].iloc[highlight],
                       s=p1_size_highlight * 1.3, color="black", linewidth=0, zorder=2)
            # Dash length matches the (non-highlighted) p1 dot's own diameter (same s -> same
            # linear scale) instead of a wider, independently-sized dash.
            ax.scatter(x, docking_df[("median", label)], marker="_", s=p1_size,
                       linewidths=stylia.LINEWIDTH * 1.5, color=lib_color, zorder=3)
            ax.scatter(non_highlight_x, docking_df[("p1", label)].iloc[non_highlight], s=p1_size,
                       color=lib_color, linewidth=0, zorder=3)
            ax.scatter(highlight_x, docking_df[("p1", label)].iloc[highlight],
                       s=p1_size_highlight, color=lib_color, linewidth=0, zorder=3)
        # Margin between the axis border and the first/last dot scales with the pocket count, so
        # it reads as the same visual proportion of the panel across genes with very different
        # pocket counts (5 to 40) instead of a fixed absolute margin shrinking in relative terms
        # as more dots are packed into the same width.
        margin = 0.08 * max(n - 1, 1)
        ax.set_xlim(-margin, (n - 1) + margin)
        # -15 (not -12) at the bottom so every panel D cross-reference star fits on-scale - the
        # most extreme one (pheS's rank-1/2 snapshots) sits at -14.5/-14.3, both more negative
        # than a single pocket's own 1st-percentile cutoff (the previous -12 floor was tuned to
        # that, not to single-best-compound scores).
        ax.set_ylim(-15, -6)

        # Cross-reference to panel D: a black star at the exact best-compound score (sharper
        # than even the p1 dot, since it's the single best compound, not a 1st-percentile cutoff)
        # for every panel D snapshot that came from this gene, labeled with that snapshot's own
        # panel D number - so the two panels can be read together. Every label sits to the right
        # of its own star. Stars for the same gene can land very close in score (e.g. pheS's
        # rank-1/2, 0.27 apart, same pocket/x) - MIN_LABEL_SEP pushes a stacked label straight up
        # (same x as the one below it) just far enough that the numbers don't touch.
        points = []
        for row in snapshot_rows or []:
            if row["pocket"] not in docking_df.index:
                continue
            points.append((docking_df.index.get_loc(row["pocket"]), row["score"], int(row["rank"])))
        points.sort(key=lambda p: p[1])  # ascending score - bottom to top

        # Same idea as the margin above: a fixed offset would look cramped on a 5-pocket panel
        # and lost on a 40-pocket one, so scale it with pocket count instead.
        LABEL_BASE_DX = 0.05 * max(n - 1, 1)
        MIN_LABEL_SEP = 1.1
        prev_label_y = None
        for x_pos, score, rank in points:
            label_y = score if prev_label_y is None else max(score, prev_label_y + MIN_LABEL_SEP)
            prev_label_y = label_y
            label_x = x_pos + LABEL_BASE_DX
            ax.scatter([x_pos], [score], marker="*", s=stylia.MARKERSIZE * 1.5, color="black",
                       linewidth=0, zorder=6)
            ax.text(label_x, label_y, str(rank), fontsize=stylia.FONTSIZE, ha="left", va="center",
                    zorder=6)

        ax.set_xticks([])
        if not show_yaxis:
            # tick_params(left=False, labelleft=False), not set_yticks([]) - clearing the ticks
            # outright also removes the tick *locations* gridlines are drawn at, which silently
            # dropped the y-grid on every column but the first.
            ax.tick_params(axis="y", left=False, labelleft=False)
        stylia.label(ax, xlabel="", ylabel="Docking score" if show_yaxis else "", abc=abc)
        # Horizontal gridlines only, to help read off docking-score values - panel C only.
        # visible=True is required, not just axis="y": with visible left as the matplotlib
        # default (None), grid() toggles off the axis it's given whenever rcParams already
        # started it enabled (stylia's default), which is exactly this case.
        ax.grid(True, axis="y")
    else:
        ax.scatter([0.5], [0.5], color=color, edgecolor="black", linewidth=0.5,
                   s=stylia.MARKERSIZE * 3, zorder=2)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        stylia.label(ax, xlabel="", ylabel="", abc=abc)
        ax.grid(False)
    # Full rectangle border (all four spines - default visible, so nothing to disable) instead of
    # just the bottom+left "L" panel B uses.
    if abc is not None:
        boost_abc_fontsize(ax, shrink_factor)

    legend_handles = [Line2D([0], [0], marker="o", color="w", label=gene, markerfacecolor=color,
                              markeredgecolor="black", markeredgewidth=0.5, markersize=stylia.MARKERSIZE)]
    # loc="lower center" anchors the legend's own BOTTOM edge at bbox_to_anchor's y, so the whole
    # legend sits above that line (fully outside the axes rectangle) instead of straddling it -
    # "upper center" was anchoring the legend's top edge there, letting its body hang down into
    # the plot.
    ax.legend(handles=legend_handles, loc="lower center", frameon=False, bbox_to_anchor=(0.5, 1.05),
              fontsize=stylia.FONTSIZE, handletextpad=0.3)


def main():
    libraries = load_libraries()
    for label, _, df in libraries:
        print(f"{label}: {len(df):,} compounds")

    scheme_pdf_path = os.path.join(plots_dir, "scheme.pdf")
    if not os.path.isfile(scheme_pdf_path):
        raise FileNotFoundError(
            f"{scheme_pdf_path} not found - export scheme.svg from Inkscape as scheme.pdf "
            "(File > Save As > PDF); the merged figure reads that PDF rather than the SVG "
            "directly (embedding the SVG's text via svgutils.compose shifts label positions)."
        )
    scheme_width_mm, _ = pdf_page_size_mm(scheme_pdf_path)

    shrink_factor = expected_shrink_factor(scheme_width_mm)
    gene_colors = load_gene_colors()
    genes = sorted(gene_colors.keys())
    n_gene_rows = 3
    assert len(genes) == n_gene_rows * len(PROPERTIES), (
        f"Expected {n_gene_rows * len(PROPERTIES)} genes for a {n_gene_rows}x{len(PROPERTIES)} "
        f"grid, got {len(genes)}."
    )

    # UNIT_HEIGHT_FRACTION is the per-ratio-unit height that made a C row come out at its
    # original absolute size (5 units -> 0.12375*6/18*5 = 0.20625 of total height). Keeping this
    # unit fixed while shrinking just ROW_HEIGHT_RATIOS[0] (B's ratio) means B's plot area shrinks
    # without touching C's rows at all - create_figure()'s height is the *total* figure height,
    # not a per-row height, so it's recomputed from the (now smaller) ratio sum.
    UNIT_HEIGHT_FRACTION = 0.12375 * 6 / 18
    n_d_panels = len(PROPERTIES)
    fig, axs = stylia.create_figure(1 + n_gene_rows + 1, len(PROPERTIES),
                                     height=UNIT_HEIGHT_FRACTION * sum(ROW_HEIGHT_RATIOS),
                                     height_ratios=ROW_HEIGHT_RATIOS)
    row1_first_ax = None
    for i, prop in enumerate(PROPERTIES):
        ax = axs.next()
        # Panel B label sits on the first (MW) panel only - titles draw outside the axes
        # bounding box, so this doesn't shift its alignment with the other 4 panels in the row.
        plot_property(ax, prop, libraries, abc="B" if i == 0 else None, shrink_factor=shrink_factor)
        if i == 0:
            row1_first_ax = ax

    # Panel C: 3 rows x 7 columns, one cell per gene (alphabetical), each titled with its
    # gene-colored circle - content deferred (empty), per request. gene_row_axes groups the 7
    # axes per row so all of a row's axes can be repositioned together below (a single row's
    # axes all share the same y0/y1, coming from one gridspec row).
    # Docking percentiles cover all 21 genes (compute_docking_percentiles() runs over every
    # pocket) - now plotted for every gene, not just alaS.
    docking_percentiles = load_docking_percentiles()
    # snapshot_index rows are also used here (not just in panel D below) - each panel D
    # rendered snapshot's exact pocket/score gets a matching marker in that gene's own panel C
    # plot, so a reader can trace one directly from the other.
    snapshot_index = load_docking_snapshot_index()
    snapshot_rows_by_gene = defaultdict(list)
    for _, row in snapshot_index.iterrows():
        snapshot_rows_by_gene[row["gene"]].append(row)

    gene_row_axes = [[] for _ in range(n_gene_rows)]
    for i, gene in enumerate(genes):
        ax = axs.next()
        docking_df = gene_docking_stats(docking_percentiles, gene)
        plot_gene_panel(ax, gene, gene_colors[gene], docking_df=docking_df,
                         abc="C" if i == 0 else None, shrink_factor=shrink_factor,
                         show_yaxis=(i % len(PROPERTIES) == 0),
                         snapshot_rows=snapshot_rows_by_gene.get(gene))
        gene_row_axes[i // len(PROPERTIES)].append(ax)
    ax_bottom = gene_row_axes[0][0]

    # Panel D: 1 row x 7 columns, the top-7 best-compound PyMOL docking snapshots overall (see
    # figure_2_calculations.py's compute_docking_snapshots() - a fixed count by global score,
    # repeats allowed, not one slot per gene) - any leftover columns (if the index has fewer than
    # 7 rows) stay empty, visible-frame placeholders. d_row_axes is kept (mirroring
    # gene_row_axes) so this row can be repositioned below, independent of whatever hspace gives
    # the C3-D gap. snapshot_index itself was already loaded above, for panel C's markers.
    d_row_axes = []
    for i in range(n_d_panels):
        ax = axs.next()
        ax.set_xticks([])
        ax.set_yticks([])
        if i < len(snapshot_index):
            row = snapshot_index.iloc[i]
            image_path = os.path.join(plots_dir, "docking_snapshots", row["filename"])
            plot_docking_snapshot(ax, i + 1, image_path, row["library"], row["score"])
        stylia.label(ax, xlabel="", ylabel="", abc="D" if i == 0 else None)
        if i == 0:
            boost_abc_fontsize(ax, shrink_factor)
        d_row_axes.append(ax)

    # stylia.save_figure() always calls plt.tight_layout() with no args right before saving,
    # which would silently undo any fig.subplots_adjust() call made beforehand - so set the
    # gap between the two rows AFTER tight_layout() runs, then save directly (matching
    # stylia.save_figure()'s own dpi/transparent/bbox_inches settings) instead of going through it.
    plt.tight_layout()
    # wspace pulled in from matplotlib's ~0.2 default so more of the row's fixed total width goes
    # to each panel's actual plot area rather than the gaps between the 5 columns (panel B's
    # subplots were reading as too narrow).
    fig.subplots_adjust(hspace=0.6, wspace=0.08)

    # Measure, directly from matplotlib's own layout, (1) how far panel B's title sits from the
    # left edge of what will become the saved SVG's own viewBox, and (2) the real gap between
    # row 1 (B) and row 2 (C) - both expressed in points relative to the tight bbox that
    # bbox_inches="tight" crops to (pad_inches=0 below makes that crop's origin exactly the SVG's
    # own (0, 0), so no extra offset needs accounting for). merge_with_scheme() then only has to
    # multiply these by its own mm-per-point scale factor - this is what keeps the A/B/C labels
    # aligned and the A-B/B-C gaps equal without hardcoding any panel-specific measurements that
    # would go stale the moment fonts, panel count, or layout tweaks change.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    tight_bbox = fig.get_tightbbox(renderer)  # inches, uncropped frame
    # stylia.label(abc=...) sets the title via ax.set_title(loc="left", ...), which matplotlib
    # stores in ax._left_title - ax.title itself is the (here, empty/unused) centered title and
    # has a zero-size window extent, which silently produced a degenerate offset before this fix.
    title_bbox_px = row1_first_ax._left_title.get_window_extent(renderer=renderer)  # display pixels
    title_x0_in = title_bbox_px.x0 / fig.dpi
    b_title_offset_pt = (title_x0_in - tight_bbox.x0) * 72

    fig_width_in, fig_height_in = fig.get_size_inches()
    row1_pos = row1_first_ax.get_position()
    row2_pos = ax_bottom.get_position()

    # By default, a loc="left" title sits flush with the axes frame's own left edge - but row 2's
    # y-axis (tick numbers + "Docking score") needs more horizontal room than that frame edge
    # allows, so it extends b_title_offset_pt further left than the title itself, "surpassing" it.
    # Shift both titles left by that same offset (converted into each axes' own fraction-of-width
    # units, since row 2's merged axes is 5x wider than row 1's) so they land at the tight bbox's
    # actual left edge - i.e. exactly where the widest y-axis content reaches - instead of at the
    # frame edge. After this, label A can simply match panel_bc_x with no extra offset needed.
    row1_first_ax._left_title.set_x(-(b_title_offset_pt / 72) / (row1_pos.width * fig_width_in))
    ax_bottom._left_title.set_x(-(b_title_offset_pt / 72) / (row2_pos.width * fig_width_in))
    # Panel D's first cell never had this shift applied, so its "D" label sat flush with its own
    # frame edge (default loc="left" behavior) instead of lining up with A/B/C's shifted x - apply
    # the same correction here.
    row_d_pos = d_row_axes[0].get_position()
    d_row_axes[0]._left_title.set_x(-(b_title_offset_pt / 72) / (row_d_pos.width * fig_width_in))

    gap_fraction = row1_pos.y0 - row2_pos.y1  # figure-fraction; unaffected by the later tight crop
    gap_b_c_pt = gap_fraction * fig_height_in * 72

    # distributions_svg_path's own top edge (y=0) is the top of the tight bbox, i.e. the "B"
    # title's top - NOT row 1's own axes/frame top, which sits below the title. merge_with_scheme()
    # needs row 1's frame-top offset from that same SVG origin (not the title's) so the A-B gap
    # (schematic content -> row 1's frame) is measured the same way as the B-C gap (row 1's frame
    # -> row 2's frame), instead of accidentally comparing a label-top-referenced gap to a
    # frame-referenced one.
    row1_frame_top_in = row1_pos.y1 * fig_height_in
    row1_frame_top_offset_pt = (tight_bbox.y1 - row1_frame_top_in) * 72

    # natural_width_mm is the distributions block's own (unscaled) tight-cropped physical width;
    # scale is how much merge_with_scheme() needs to shrink it to fit scheme's width minus
    # MARGIN_MM on each side - the same ratio expected_shrink_factor() predicted earlier, just
    # recomputed from the actual rendered layout instead of the pre-build estimate.
    natural_width_mm = tight_bbox.width * 25.4
    scale = (scheme_width_mm - 2 * MARGIN_MM) / natural_width_mm

    # Panel C's 3 rows inherited whatever gap fig.subplots_adjust(hspace=...) produced above -
    # since that same hspace also sets the B-to-C1 gap (row1_pos/row2_pos, already measured), and
    # C's rows are taller than B's (height_ratios), the same hspace fraction yields a bigger
    # absolute C-internal gap than the B-C1 one. Tighten just the C-internal gaps by explicitly
    # repositioning rows 2/3 upward - independent of the B-C gap, which keeps whatever hspace gave
    # row 1 (already captured in gap_b_c_pt above).
    C_ROW_GAP_SHRINK = 0.65  # fraction of the current C-row gap to keep - now that the gene-name
    # legend sits fully above each row's frame (loc="lower center", see plot_gene_panel()), the
    # gap needs enough room for that legend, not just breathing room between two bare frames.
    c_row_positions = [row[0].get_position() for row in gene_row_axes]
    gap_c1_c2 = c_row_positions[0].y0 - c_row_positions[1].y1
    gap_c2_c3 = c_row_positions[1].y0 - c_row_positions[2].y1
    target_gap = gap_c1_c2 * C_ROW_GAP_SHRINK
    shift_row2 = gap_c1_c2 - target_gap
    shift_row3 = shift_row2 + (gap_c2_c3 - target_gap)
    for ax in gene_row_axes[1]:
        pos = ax.get_position()
        ax.set_position([pos.x0, pos.y0 + shift_row2, pos.width, pos.height])
    for ax in gene_row_axes[2]:
        pos = ax.get_position()
        ax.set_position([pos.x0, pos.y0 + shift_row3, pos.width, pos.height])

    # Panel D inherited the same raw hspace gap as C's rows (both ratio 5), so it needs the same
    # tightening, cumulative on top of C3's own shift_row3 move above.
    d_row_pos = d_row_axes[0].get_position()
    gap_c3_d = c_row_positions[2].y0 - d_row_pos.y1
    shift_d = shift_row3 + (gap_c3_d - target_gap)
    for ax in d_row_axes:
        pos = ax.get_position()
        ax.set_position([pos.x0, pos.y0 + shift_d, pos.width, pos.height])

    # Panels B+C alone aren't a deliverable - only the schematic+B+C combination (figure_2.png/
    # .pdf) is. Saving directly to PNG at dpi=scale*MERGED_DPI makes matplotlib's own tight-bbox
    # crop land exactly at the target (post-scale) physical size in one step, so
    # merge_with_scheme() can paste it into the composite with no further resizing.
    tmp_png_fd, tmp_png_path = tempfile.mkstemp(suffix=".png")
    os.close(tmp_png_fd)
    plt.savefig(tmp_png_path, dpi=scale * MERGED_DPI, transparent=False, bbox_inches="tight", pad_inches=0)
    try:
        merge_with_scheme(tmp_png_path, scale, gap_b_c_pt, row1_frame_top_offset_pt,
                           scheme_pdf_path, scheme_width_mm)
    finally:
        os.remove(tmp_png_path)


def merge_with_scheme(distributions_png_path, scale, gap_b_c_pt, row1_frame_top_offset_pt,
                       scheme_pdf_path, scheme_width_mm):
    """Combines panel A (scheme.pdf, exported from the hand-authored Inkscape schematic) with
    panels B+C (the property-distributions and gene-grid rows just rendered, already carrying
    their own "B"/"C" labels via stylia.label(abc=...)) into the final figure_2.png/.pdf.

    Both source images are already rasters at this point - scheme.pdf via Poppler
    (render_pdf_to_png(), the same renderer behind Inkscape's own PDF preview/export) and the
    distributions block via matplotlib's own Agg backend, already sized to its final physical
    scale (see the dpi=scale*MERGED_DPI savefig in main()). Compositing is therefore plain PIL
    pasting, with no vector-SVG-embedding step - deliberately, since embedding scheme.svg's raw
    text via svgutils.compose was found to shift every text label's position on re-render (shapes
    stayed pixel-exact; text did not), while rasterizing each panel with its own native renderer
    and pasting reproduces Inkscape's own PDF export exactly.

    gap_b_c_pt (from main(), in points relative to distributions_png_path's own origin) places the
    A-B gap to match B/C's actual rendered row spacing, instead of an independent hardcoded
    constant that would drift out of sync whenever the matplotlib layout changes. Label A itself
    just matches panel_bc_x directly - main() already shifted B/C's own titles to sit at that same
    origin (see the row1_first_ax._left_title.set_x(...) call), so no separate offset is needed
    here.
    """
    scheme_png_fd, scheme_png_path = tempfile.mkstemp(suffix=".png")
    os.close(scheme_png_fd)
    try:
        render_pdf_to_png(scheme_pdf_path, MERGED_DPI, scheme_png_path)
        scheme_img = Image.open(scheme_png_path).convert("RGB")
        # scheme.pdf's declared page height often has incidental blank space below its actual
        # drawing (varies as the file gets re-edited in Inkscape) - measure the real content edge
        # from the same raster instead of trusting the page size for spacing.
        scheme_content_bottom_mm = content_bottom_mm(scheme_png_path, MERGED_DPI)
    finally:
        os.remove(scheme_png_path)

    distributions_img = Image.open(distributions_png_path).convert("RGB")
    distributions_width_mm = distributions_img.width / MERGED_DPI * 25.4
    distributions_height_mm = distributions_img.height / MERGED_DPI * 25.4

    # panel_bc_y positions distributions_png_path's own origin (the "B" title's top, not row 1's
    # frame - see the row1_frame_top_offset_pt docstring note above), so it's offset backwards by
    # row 1's own frame-top offset - that's what makes (scheme content -> row 1 frame) end up the
    # same distance apart as (row 1 frame -> row 2 frame), instead of comparing a label-top gap to
    # a frame-top gap. gap_b_c_pt/row1_frame_top_offset_pt are real points (72/inch) from
    # matplotlib's renderer, so they need the standard pt->mm factor before the (dimensionless)
    # compositing shrink is applied on top.
    points_to_mm = 25.4 / 72
    gap_b_c_mm = gap_b_c_pt * points_to_mm * scale
    row1_frame_top_offset_mm = row1_frame_top_offset_pt * points_to_mm * scale
    panel_bc_x = MARGIN_MM  # flush against the same margin as label A, instead of centered (which
                            # pushed B/C inward by half of the old fixed-fraction leftover width)
    panel_bc_y = scheme_content_bottom_mm + gap_b_c_mm - row1_frame_top_offset_mm

    total_width_mm = scheme_width_mm
    # MARGIN_MM of blank space below panel D too, matching the breathing room already used on
    # the top/left/right (panel A's own top margin, panel_bc_x) instead of the content ending
    # flush with the canvas edge.
    total_height_mm = panel_bc_y + distributions_height_mm + MARGIN_MM
    mm_to_px = lambda v: round(v * MERGED_DPI / 25.4)  # noqa: E731

    canvas = Image.new("RGB", (mm_to_px(total_width_mm), mm_to_px(total_height_mm)), "white")
    canvas.paste(scheme_img, (0, 0))
    canvas.paste(distributions_img, (mm_to_px(panel_bc_x), mm_to_px(panel_bc_y)))

    draw = ImageDraw.Draw(canvas)
    label_font = ImageFont.truetype(LABEL_FONT_PATH, round(LABEL_SIZE_PT * MERGED_DPI / 72))
    draw.text((mm_to_px(panel_bc_x), mm_to_px(MARGIN_MM)), "A", font=label_font, fill="black")

    png_path = os.path.join(plots_dir, "figure_2.png")
    pdf_path = os.path.join(plots_dir, "figure_2.pdf")
    canvas.save(png_path)
    canvas.save(pdf_path, "PDF", resolution=float(MERGED_DPI))

    print(f"Panel A (scheme): {scheme_width_mm:.1f} mm wide (content to {scheme_content_bottom_mm:.1f} mm)")
    print(f"Panels B+C (distributions): {distributions_width_mm:.1f} x {distributions_height_mm:.1f} mm (scaled by {scale:.4f})")
    print(f"A-B gap: {gap_b_c_mm:.2f} mm, B-C gap: {gap_b_c_mm:.2f} mm")
    print(f"Total: {total_width_mm:.1f} x {total_height_mm:.1f} mm")
    print(f"Saved to {png_path}")
    print(f"Saved to {pdf_path}")


if __name__ == "__main__":
    main()
