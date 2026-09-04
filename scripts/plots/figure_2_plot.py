"""
Figure 2: compound-property distributions (panel b), per-gene docking-percentile grid
one column per canonical (deduplicated) pocket (panel c), and top docking-pose snapshots
(panel d), from figure_2_calculations.py. Panel a is a hand-authored Inkscape schematic
(scheme.svg, exported to scheme.pdf) - not built by this script at all, just pasted in.

Built on figure_1_plot.py's panel-per-file architecture instead of the original, single
shared b+c+d matplotlib block PIL-composited onto a rasterized scheme.pdf: b/c/d each
save as their own standalone PDF at an EXACT physical size read from
output/plots/figure_2/panel_layout.csv (columns: panel, x, delta_x, y, delta_y, padding,
in cm), with a baked-in panel letter (add_panel_label). Panel a (scheme.pdf) is pasted
the same way via pypdf, scaled to fit its declared box in case the Inkscape export
doesn't match panel_layout.csv exactly - no rasterization, no Poppler dependency, and its
own "a" label is baked on separately (_render_panel_a_label(), a small transparent
overlay PDF merged on top of it) since scheme.pdf itself carries no label. Merging into
one positioned master PDF (Fig_2_full.pdf) happens in this same file (merge_panels(),
called at the end of main()), matching figure_1_plot.py's choice to fold plot+merge into
one script. merge_panels() also exports a flattened Fig_2_full.png of that same merged
page (user request) via pymupdf - a real Python dependency, but not Poppler, so this
doesn't reintroduce the dependency panel a's own vector merge above deliberately avoids.

panel_layout.csv is a first-draft starting point (a's box matches scheme.pdf's own
measured size; b/c/d are rough guesses), to be tuned iteratively by rendering and
looking, not solved analytically up front.

Usage:
    python figure_2_plot.py [--subpanels b,c,d]
"""
import argparse
import json
import os
import sys
from collections import defaultdict

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
from PIL import Image
from pypdf import PdfReader, PdfWriter, Transformation
import pymupdf
from scipy.stats import gaussian_kde
import stylia
from stylia.config import get_fg_color
from stylia.figure.figure import stylize

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

plots_dir = os.path.join(root, "..", "..", "output", "plots", "figure_2")
os.makedirs(plots_dir, exist_ok=True)

PANEL_LETTERS = ["a", "b", "c", "d"]
panel_layout_path = os.path.join(plots_dir, "panel_layout.csv")

# ===========================================================================
# Panel b data + helper: compound-property distributions
# (ported verbatim from the original shared-figure version of this script, abc/shrink_factor dropped - panel-lettering
# is now centralized in save_panel()/add_panel_label() instead of per-axes)
# ===========================================================================

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


def load_libraries():
    loaded = []
    for label, filename, color_name in LIBRARIES:
        df = pd.read_csv(os.path.join(plots_dir, filename))
        loaded.append((label, color_name, df))
    return loaded


def load_gene_colors():
    """Gene->color mapping from figure 1 - reused here so protein colors stay consistent
    across figures instead of a second palette."""
    with open(os.path.join(root, "..", "..", "output", "plots", "figure_1", "color_mapping.json")) as f:
        return json.load(f)["gene_to_color"]


def plot_property(ax, prop, libraries):
    nc = stylia.NamedColors()
    lo, hi = PROP_XLIMS[prop]

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
            ax.stairs(freq.values, edges, baseline=0, color=color, linewidth=stylia.LINEWIDTH * 1.5,
                      label=f"{label} (n={len(df) // 1000}k)")
        ax.set_xlim(lo - 0.5, hi + 0.5)
        stylia.label(ax, xlabel="", ylabel="", title=PROP_XLABELS.get(prop, prop))
    else:
        for label, color_name, df in libraries:
            data = df[prop].dropna().values
            kde = gaussian_kde(data)
            x = np.linspace(lo, hi, 300)
            y = kde(x)
            color = nc.get(color_name)
            ax.plot(x, y, color=color, linewidth=stylia.LINEWIDTH * 1.5, label=f"{label} (n={len(df) // 1000}k)")
        ax.set_xlim(lo, hi)
        stylia.label(ax, xlabel="", ylabel="", title=PROP_XLABELS.get(prop, prop))
        if prop == "MW":
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))

    if prop in PROP_XTICKS:
        ax.set_xticks(PROP_XTICKS[prop])

    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)


# ===========================================================================
# Panel c data + helper: per-gene docking-percentile grid
# (ported verbatim from the original shared-figure version of this script, abc/shrink_factor dropped)
# ===========================================================================

def load_docking_percentiles():
    return pd.read_csv(os.path.join(plots_dir, "figure_2c_docking_percentiles.csv"))


def load_pocket_clusters():
    """structure-level pocket id -> spatial_cluster_id, from the same 6.14 A greedy
    centroid dedup figure_1_calculations.py uses for gene_to_unique_pocket_count
    (scripts/77_pocket_annotation/09_cluster_pockets.py's persisted assignments) -
    collapses panel c's per-structure pockets down to one column per canonical site."""
    path = os.path.join(root, "..", "..", "output", "77_pocket_annotation", "pocket_clusters.csv")
    df = pd.read_csv(path)
    df["pocket_id"] = df["File name"].str.replace(".pdb", "", regex=False) + "_pocket_" + df["Pocket number"].astype(str)
    return df.set_index("pocket_id")["spatial_cluster_id"].to_dict()


def gene_docking_stats(docking_percentiles, gene, pocket_to_cluster):
    """Per-canonical-pocket median and p1 docking score for one gene - columns are a
    (stat, library) MultiIndex, canonical pockets sorted by REAL 10B's p1 ascending (most
    negative/best first) so both stats share one x-axis order across all three libraries
    instead of each sorting independently. A canonical pocket can have several
    structure-level candidates (same physical site detected on different structures) -
    only the one with the best (most negative) REAL 10B p1 represents it here, so each
    cluster contributes exactly one column instead of one per redundant detection. Empty
    DataFrame if the gene has no docking data - plot_gene_panel() falls back to its
    dummy-dot placeholder in that case."""
    gene_df = docking_percentiles[docking_percentiles["gene"] == gene].copy()
    if gene_df.empty:
        return gene_df
    gene_df["spatial_cluster_id"] = gene_df["pocket"].map(pocket_to_cluster)
    real10b = gene_df[gene_df["library"] == "REAL 10B"]
    winners = real10b.sort_values("p1").drop_duplicates("spatial_cluster_id")["pocket"]
    gene_df = gene_df[gene_df["pocket"].isin(winners)]
    pivoted = gene_df.pivot(index="pocket", columns="library", values=["median", "p1"])
    return pivoted.sort_values(("p1", "REAL 10B"))


def load_docking_snapshot_index():
    """Ordered rows (rank, gene, library, pocket, compound, score, filename) for panel d's 7
    snapshot cells - figure_2_calculations.py's compute_docking_snapshots() writes this index
    alongside the PNGs since panel d shows a fixed top-7 by score (repeats allowed), not one slot
    per gene, so the same gene's snapshots need distinguishing beyond just its name."""
    index_path = os.path.join(plots_dir, "docking_snapshots", "index.csv")
    if not os.path.isfile(index_path):
        return pd.DataFrame()
    return pd.read_csv(index_path).sort_values("rank")


def plot_gene_panel(ax, gene, color, docking_df=None, show_ylabel=False, show_yticks=False,
                     snapshot_rows=None):
    """One cell of panel c's 3x7 gene grid, titled with a colored circle + gene name.
    Plots each library's per-pocket median and p1 docking score (x = pocket identity,
    unlabeled - only their relative order matters here) when docking_df is given;
    otherwise a dummy single-dot scatter placeholder (for a gene with no docking data at
    all) in the gene's own color. show_ylabel (row's first column) draws the "Docking
    score" text with no tick marks/numbers; show_yticks (row's LAST column) draws the
    tick marks/numbers on that column's right edge with no label text - splitting the two
    across opposite ends of the row (rather than stacking both on the first column) is
    what lets the row's own left margin shrink to just the label's width, per request.
    Every column is one canonical pocket (gene_docking_stats() already dedups to that),
    so all points render at the same size."""
    ax.patch.set_visible(False)
    if docking_df is not None and not docking_df.empty:
        nc = stylia.NamedColors()
        n = len(docking_df)
        x = list(range(n))
        for label, _, color_name in LIBRARIES:
            lib_color = nc.get(color_name)
            marker_size = stylia.MARKERSIZE * 1.6
            ax.scatter(x, docking_df[("median", label)], marker="_", s=marker_size,
                       linewidths=stylia.LINEWIDTH * 1.5, color=lib_color, zorder=3)
            ax.scatter(x, docking_df[("p1", label)], s=marker_size, color=lib_color, linewidth=0, zorder=3)
        # A margin that's a pure multiple of (n - 1) shrinks to near-zero for genes with only
        # 2-3 canonical pockets, so the outermost big points sit right on the frame - the
        # 0.35 floor keeps a readable gap regardless of n, only kicking in below ~n=5 (0.08 *
        # (n-1) already exceeds it from n=6 up, so larger grids are unaffected).
        margin = max(0.08 * (n - 1), 0.35)
        ax.set_xlim(-margin, (n - 1) + margin)
        ax.set_ylim(-15, -6)

        points = []
        for row in snapshot_rows or []:
            if row["pocket"] not in docking_df.index:
                continue
            points.append((docking_df.index.get_loc(row["pocket"]), row["score"], int(row["rank"])))
        points.sort(key=lambda p: p[1])

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
        ax.tick_params(axis="y", left=False, labelleft=False, right=False, labelright=False)
        if show_yticks:
            # Numbers on the RIGHT edge of the row's last column, not attached to the labeled
            # first column - lets the first column's own left margin shrink to just the
            # "Docking score" text's width instead of text+numbers stacked together (same idea as
            # figure_1's boxplot tick_right(), which keeps a label on its own separate side).
            ax.yaxis.tick_right()
            ax.tick_params(axis="y", right=True, labelright=True)
        stylia.label(ax, xlabel="", ylabel="Docking score" if show_ylabel else "")
        ax.grid(True, axis="y")
    else:
        ax.scatter([0.5], [0.5], color=color, edgecolor="black", linewidth=0.5,
                   s=stylia.MARKERSIZE * 3, zorder=2)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        stylia.label(ax, xlabel="", ylabel="")
        ax.grid(False)

    legend_handles = [Line2D([0], [0], marker="o", color="w", label=gene, markerfacecolor=color,
                              markeredgecolor="black", markeredgewidth=0.5, markersize=stylia.MARKERSIZE)]
    ax.legend(handles=legend_handles, loc="lower center", frameon=False, bbox_to_anchor=(0.5, 1.05),
              fontsize=stylia.FONTSIZE, handletextpad=0.3)


# ===========================================================================
# Panel d data + helper: top docking-pose snapshots
# (ported verbatim from the original shared-figure version of this script, abc dropped)
# ===========================================================================

LIBRARY_DISPLAY_NAMES = {"HL": "Hit Locator", "REAL 10M": "REAL 10M", "REAL 10B": "REAL 10B"}


def plot_docking_snapshot(ax, rank, image_path, library, score):
    """One panel d cell: a best-compound PyMOL docking snapshot, titled with the panel's own
    plain-text number (outside the frame via loc="lower center", no color/marker), plus two more
    plain-text legends - which library (top center) and its docking score (bottom center)."""
    ax.imshow(np.array(Image.open(image_path)))
    ax.set_xticks([])
    ax.set_yticks([])

    number_handle = [Line2D([0], [0], linestyle="None", marker="None", label=str(rank))]
    number_legend = ax.legend(handles=number_handle, loc="lower center", frameon=False,
                               bbox_to_anchor=(0.5, 1.0), fontsize=stylia.FONTSIZE_BIG, handlelength=0,
                               handletextpad=0, borderaxespad=0)
    ax.add_artist(number_legend)

    library_handle = [Line2D([0], [0], linestyle="None", marker="None",
                              label=LIBRARY_DISPLAY_NAMES.get(library, library))]
    library_legend = ax.legend(handles=library_handle, loc="upper center", frameon=True, framealpha=0.8,
                                edgecolor="black", fontsize=stylia.FONTSIZE, handlelength=0, handletextpad=0)
    ax.add_artist(library_legend)

    score_handle = [Line2D([0], [0], linestyle="None", marker="None", label=f"Score: {score:.2f}")]
    ax.legend(handles=score_handle, loc="lower center", frameon=True, framealpha=0.8, edgecolor="black",
              fontsize=stylia.FONTSIZE, handlelength=0, handletextpad=0)


# ===========================================================================
# Panel-saving scaffolding (ported near-verbatim from figure_1_plot.py)
# ===========================================================================

PANEL_LABEL_MARGIN = 0.02


def add_panel_label(fig, letter):
    """Bold panel letter at the top-left of the FIGURE (page), fixed regardless of
    padding - from figure_2_plot.py. Used for b/c/d only - panel a's label is baked on
    separately by _render_panel_a_label(), since scheme.pdf isn't a matplotlib figure."""
    fig.text(PANEL_LABEL_MARGIN, 1 - PANEL_LABEL_MARGIN, letter, fontweight="bold",
              fontsize=stylia.FONTSIZE_BIG, color=get_fg_color(), ha="left", va="top",
              transform=fig.transFigure)


MAX_WIDTH_IN = stylia.SIZE  # Nature two-column guideline (~7.09in/18cm) - from figure_2_plot.py


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
    """Saves at exactly panel_layout.csv's delta_x/delta_y (no bbox_inches="tight").
    use_tight_layout=False + subplots_adjust is the escape hatch for panels where
    tight_layout() warns/fails - from figure_2_plot.py. fig.set_tight_layout(False)
    defeats stylia's own auto-relayout rcParam, which would otherwise silently override
    this layout at savefig time."""
    fig.set_tight_layout(False)
    if use_tight_layout:
        plt.tight_layout(pad=tight_pad, w_pad=tight_w_pad)
    else:
        fig.subplots_adjust(**(subplots_adjust or dict(left=0.01, right=0.99, top=0.99, bottom=0.01)))
    apply_padding(fig, padding)
    add_panel_label(fig, letter)
    output_path = os.path.join(plots_dir, f"Fig_2{letter}.pdf")
    plt.savefig(output_path, dpi=600, transparent=False)
    plt.close(fig)
    print(f"Saved Fig_2{letter}.pdf")


CM_TO_PT = 72 / 2.54


def _render_panel_a_label(delta_x_cm, delta_y_cm):
    """Bakes a bold "a" label (same style/corner-inset as add_panel_label()) onto its own
    transparent-background PDF, sized to match panel a's declared box - scheme.pdf itself
    has no label drawn into it, so this overlay supplies one instead of requiring a manual
    edit in Inkscape."""
    fig = plt.figure(figsize=(delta_x_cm / 2.54, delta_y_cm / 2.54))
    fig.patch.set_alpha(0)
    add_panel_label(fig, "a")
    label_path = os.path.join(plots_dir, "_panel_a_label.pdf")
    fig.savefig(label_path, transparent=True)
    plt.close(fig)
    return label_path


def merge_panels():
    """Pastes each panel onto one positioned master PDF (Fig_2_full.pdf) per
    panel_layout.csv's x/y. b/c/d are pure translation (already saved at their exact
    declared size via save_panel()). Panel a (scheme.pdf, hand-authored in Inkscape, not
    script-generated) is additionally scaled to fit its declared box, in case the
    Inkscape export doesn't match panel_layout.csv pixel-for-pixel. Ported from
    figure_1_plot.py's own merge_panels(), kept in this same file per the same
    plot+merge-in-one-script convention now shared by all 4 figures."""
    df = pd.read_csv(panel_layout_path).set_index("panel")
    missing = [p for p in PANEL_LETTERS if p not in df.index]
    if missing:
        raise ValueError(f"{panel_layout_path} is missing row(s) for panel(s): {missing}")

    total_width_cm = max(df.loc[p, "x"] + df.loc[p, "delta_x"] for p in PANEL_LETTERS)
    total_height_cm = max(df.loc[p, "y"] + df.loc[p, "delta_y"] for p in PANEL_LETTERS)

    writer = PdfWriter()
    master_page = writer.add_blank_page(width=total_width_cm * CM_TO_PT, height=total_height_cm * CM_TO_PT)

    for p in PANEL_LETTERS:
        panel_path = os.path.join(plots_dir, "scheme.pdf" if p == "a" else f"Fig_2{p}.pdf")
        if not os.path.exists(panel_path):
            print(f"  Skipping {p}: {panel_path} not found - Fig_2_full.pdf will be missing it.")
            continue
        panel_page = PdfReader(panel_path).pages[0]

        x_top_cm, y_top_cm = df.loc[p, "x"], df.loc[p, "y"]
        delta_x_cm, delta_y_cm = df.loc[p, "delta_x"], df.loc[p, "delta_y"]
        x_pt = x_top_cm * CM_TO_PT
        y_bottom_pt = (total_height_cm - y_top_cm - delta_y_cm) * CM_TO_PT

        transform = Transformation()
        if p == "a":
            actual_w_pt = float(panel_page.mediabox.width)
            actual_h_pt = float(panel_page.mediabox.height)
            target_w_pt = delta_x_cm * CM_TO_PT
            target_h_pt = delta_y_cm * CM_TO_PT
            transform = transform.scale(target_w_pt / actual_w_pt, target_h_pt / actual_h_pt)
        transform = transform.translate(tx=x_pt, ty=y_bottom_pt)

        master_page.merge_transformed_page(panel_page, transform)
        print(f"  Placed {os.path.basename(panel_path)} at x={x_top_cm}cm, y={y_top_cm}cm (top-left)")

        if p == "a":
            label_path = _render_panel_a_label(delta_x_cm, delta_y_cm)
            label_page = PdfReader(label_path).pages[0]
            master_page.merge_transformed_page(label_page, Transformation().translate(tx=x_pt, ty=y_bottom_pt))
            os.remove(label_path)

    output_path = os.path.join(plots_dir, "Fig_2_full.pdf")
    with open(output_path, "wb") as f:
        writer.write(f)
    print(f"Saved merged master figure ({total_width_cm:.2f} x {total_height_cm:.2f} cm) to {output_path}")

    # Flattened PNG alongside the vector PDF (user request), rendered at the same dpi=600
    # save_panel's own PDFs already use for their embedded raster content - pymupdf renders the
    # already-merged page directly (no external Poppler binary, unlike pdftoppm/pdf2image, so
    # this doesn't reintroduce the Poppler dependency panel a's own merge deliberately avoids -
    # see module docstring).
    png_path = os.path.join(plots_dir, "Fig_2_full.png")
    pdf_doc = pymupdf.open(output_path)
    zoom = 600 / 72  # pymupdf's Pixmap render defaults to 72dpi
    pix = pdf_doc[0].get_pixmap(matrix=pymupdf.Matrix(zoom, zoom))
    pix.save(png_path)
    pdf_doc.close()
    print(f"Saved {png_path}")


# ===========================================================================
# Panel builders
# ===========================================================================

# Panels c and d are both 7-column grids at the same panel_layout.csv delta_x, and are meant to
# read as one continuous set of columns (panel d's numbered snapshots cross-reference panel c's
# per-gene plots via the same rank stars) - so their 7 columns must land at IDENTICAL x-positions/
# widths. left/right/wspace here are explicit and shared between both build_panel_c/build_panel_d
# (each with save_panel(..., use_tight_layout=False, ...)) instead of each panel's own independent
# plt.tight_layout() call, which would otherwise size each panel's left margin to just its own
# content (panel c's leftmost "Docking score" y-axis label needs real room; panel d's cells have no
# axis decoration at all) and misalign the two grids. left is set wide enough for panel c's y-axis
# label + tick numbers; panel d just carries that same margin as blank space on its own left edge.
POCKET_ROW_LEFT = 0.035
POCKET_ROW_RIGHT = 0.95
POCKET_ROW_WSPACE = 0.08


def build_panel_b(size, padding):
    libraries = load_libraries()
    fig, axs = plt.subplots(1, len(PROPERTIES), figsize=size)
    fig.patch.set_facecolor("white")
    for ax, prop in zip(axs, PROPERTIES):
        stylize(ax)
        plot_property(ax, prop, libraries)
    save_panel(fig, "b", use_tight_layout=False,
               subplots_adjust=dict(left=POCKET_ROW_LEFT, right=POCKET_ROW_RIGHT, top=0.62, bottom=0.32,
                                     wspace=POCKET_ROW_WSPACE),
               padding=padding)


def build_panel_c(size, padding):
    gene_colors = load_gene_colors()
    genes = sorted(gene_colors.keys())
    n_gene_rows = 3
    assert len(genes) == n_gene_rows * len(PROPERTIES), (
        f"Expected {n_gene_rows * len(PROPERTIES)} genes for a {n_gene_rows}x{len(PROPERTIES)} "
        f"grid, got {len(genes)}.")

    docking_percentiles = load_docking_percentiles()
    pocket_to_cluster = load_pocket_clusters()
    snapshot_index = load_docking_snapshot_index()
    snapshot_rows_by_gene = defaultdict(list)
    for _, row in snapshot_index.iterrows():
        snapshot_rows_by_gene[row["gene"]].append(row)

    fig, axs = plt.subplots(n_gene_rows, len(PROPERTIES), figsize=size)
    fig.patch.set_facecolor("white")
    for i, gene in enumerate(genes):
        row, col = divmod(i, len(PROPERTIES))
        ax = stylize(axs[row, col])
        docking_df = gene_docking_stats(docking_percentiles, gene, pocket_to_cluster)
        plot_gene_panel(ax, gene, gene_colors[gene], docking_df=docking_df,
                         show_ylabel=(col == 0), show_yticks=(col == len(PROPERTIES) - 1),
                         snapshot_rows=snapshot_rows_by_gene.get(gene))
    save_panel(fig, "c", use_tight_layout=False,
               subplots_adjust=dict(left=POCKET_ROW_LEFT, right=POCKET_ROW_RIGHT, top=0.88, bottom=0.02,
                                     wspace=POCKET_ROW_WSPACE, hspace=0.35),
               padding=padding)


def build_panel_d(size, padding):
    snapshot_index = load_docking_snapshot_index()
    fig, axs = plt.subplots(1, len(PROPERTIES), figsize=size)
    fig.patch.set_facecolor("white")
    for i, ax in enumerate(axs):
        stylize(ax)
        stylia.label(ax, xlabel="", ylabel="")
        ax.set_xticks([])
        ax.set_yticks([])
        if i < len(snapshot_index):
            row = snapshot_index.iloc[i]
            image_path = os.path.join(plots_dir, "docking_snapshots", row["filename"])
            plot_docking_snapshot(ax, i + 1, image_path, row["library"], row["score"])
    # top=0.76 (not 0.93 like panels b/c) - imshow's default equal-aspect shrinks each square
    # snapshot to the column WIDTH (fixed, shared with panel c), so a nominal box much taller
    # than it is wide (as panel_layout.csv's old delta_y=4 gave at top=0.93) just pads blank
    # space above+below the visible image instead of enlarging it. top=0.76 with delta_y=3.1
    # (panel_layout.csv) sizes the nominal box to match the column width, while still leaving
    # enough absolute top margin for the "d" label plus each cell's rank-number legend.
    save_panel(fig, "d", use_tight_layout=False,
               subplots_adjust=dict(left=POCKET_ROW_LEFT, right=POCKET_ROW_RIGHT, top=0.76, bottom=0.05,
                                     wspace=POCKET_ROW_WSPACE),
               padding=padding)


def main(subpanels=None):
    if subpanels is None:
        subpanels = PANEL_LETTERS
    sizes = load_panel_sizes(subpanels)
    paddings = load_panel_padding(subpanels)

    for letter in subpanels:
        if sizes[letter][0] > MAX_WIDTH_IN:
            print(f"WARNING: panel '{letter}' delta_x is {sizes[letter][0] * 2.54:.2f}cm, exceeding "
                  f"the {MAX_WIDTH_IN * 2.54:.2f}cm Nature two-column guideline by "
                  f"{(sizes[letter][0] - MAX_WIDTH_IN) * 2.54:.2f}cm.")

    if "b" in subpanels:
        build_panel_b(sizes["b"], paddings["b"])

    if "c" in subpanels:
        build_panel_c(sizes["c"], paddings["c"])

    if "d" in subpanels:
        build_panel_d(sizes["d"], paddings["d"])

    merge_panels()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--subpanels", type=str, default=None,
                         help="Comma-separated subset of panels to (re)generate (e.g. 'b,c'), from "
                              f"{{{','.join(PANEL_LETTERS)}}}. Panel 'a' has no build step (it's "
                              "scheme.pdf, pasted as-is) but is still validated/merged. Default: all. "
                              "Fig_2_full.pdf is always re-merged from whatever panel files "
                              "currently exist on disk.")
    args = parser.parse_args()
    if args.subpanels is None:
        subpanels = PANEL_LETTERS
    else:
        subpanels = [p.strip() for p in args.subpanels.split(",")]
        invalid = [p for p in subpanels if p not in PANEL_LETTERS]
        if invalid:
            parser.error(f"Unknown panel(s) {invalid}, must be from {PANEL_LETTERS}")
    main(subpanels=subpanels)
