import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

import json
import subprocess
import tempfile

import numpy as np
import pandas as pd
import stylia
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from lxml import etree
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
from PIL import Image
from scipy.stats import gaussian_kde
from svgutils.compose import Element, Figure, SVG, Text, Unit

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

plots_dir = os.path.join(root, "..", "..", "plots", "figure_2")

PROPERTIES = ["MW", "cLogP", "TPSA", "AromaticRings", "QED"]
DISCRETE_PROPS = {"AromaticRings"}
PROP_XLIMS = {
    "MW":            (200, 500),
    "cLogP":         (-3, 6),
    "TPSA":          (0, 160),
    "AromaticRings": (0, 5),
    "QED":           (0, 1),
}
PROP_XLABELS = {"AromaticRings": "Aromatic rings"}

# (label, csv filename in plots/figure_2/, NamedColors name)
LIBRARIES = [
    ("HL",       "figure_2a_HL.csv",       "crimson"),
    ("REAL 10M", "figure_2a_REAL10M.csv",  "turquoise"),
    ("REAL 10B", "figure_2b_REAL10B.csv",  "amber"),
]

# Panel A (scheme.svg) is a hand-authored Inkscape schematic, not a matplotlib axes, so it can't
# use stylia.label(abc=...) directly - these mirror that same style (stylia/figure/figure.py's
# label(): bold, Arial, print-format font size) for the one SVG-text label this script still has
# to inject by hand. Panels B/C are embedded natively via stylia.label(..., abc=...) instead.
SVG_NS = "http://www.w3.org/2000/svg"
LABEL_SIZE_MM = 8 * 25.4 / 72  # 8pt converted to mm, since the merged canvas's user units are mm
LABEL_MARGIN_MM = 5  # from the schematic's top-left corner
LABEL_KWARGS = dict(size=LABEL_SIZE_MM, weight="bold", font="Arial", color="black")
PANEL_GAP_MM = 6  # breathing room between the schematic and the distributions block below it
DISTRIBUTIONS_SIZE_FACTOR = 0.85  # render the distributions block a bit smaller than the scheme

# Rendering the composed merge SVG to PNG/PDF: cairosvg was tried first, but it mispositions
# scheme.svg's embedded base64 images relative to a spec-compliant renderer (confirmed by
# comparing against a Chrome headless screenshot of the identical file) - so Chrome headless
# renders the PNG instead, and the PDF is just that already-correct PNG repackaged as a
# single-page PDF via Pillow (avoids Chrome's print-to-pdf page-size/margin quirks entirely).
MERGED_DPI = 600
CSS_DPI = 96  # browsers render physical (mm/in) units at a fixed 96 CSS px per inch


def load_libraries():
    loaded = []
    for label, filename, color_name in LIBRARIES:
        df = pd.read_csv(os.path.join(plots_dir, filename))
        loaded.append((label, color_name, df))
    return loaded


def load_pocket_stats():
    """Per-pocket docking-score percentiles (figure_2c_calculations.py) for all three libraries,
    pocket order fixed by HL's 1st percentile (p1) ascending so REAL 10M/REAL 10B share the same
    x-axis instead of each sorting independently. Also loads figure 1's gene->color mapping -
    reused here so protein colors stay consistent across figures instead of a second palette."""
    df = pd.read_csv(os.path.join(plots_dir, "figure_2c_docking_percentiles.csv"))
    hl = df[df["library"] == "HL"].sort_values("p1").reset_index(drop=True)
    pocket_order = hl["pocket"].tolist()

    by_library = {}
    for label, _, _ in LIBRARIES:
        lib_df = df[df["library"] == label].set_index("pocket").reindex(pocket_order)
        by_library[label] = lib_df

    with open(os.path.join(root, "..", "..", "plots", "figure_1", "color_mapping.json")) as f:
        gene_to_color = json.load(f)["gene_to_color"]
    return hl, by_library, gene_to_color


def plot_property(ax, prop, libraries, show_legend, show_ylabel, abc=None):
    nc = stylia.NamedColors()
    lo, hi = PROP_XLIMS[prop]

    if prop in DISCRETE_PROPS:
        all_vals = set()
        for _, _, df in libraries:
            all_vals |= set(df[prop].dropna().astype(int))
        vals = sorted(all_vals)
        for label, color_name, df in libraries:
            data = df[prop].dropna().astype(int)
            freq = pd.Series(data).value_counts(normalize=True).reindex(vals, fill_value=0)
            color = nc.get(color_name)
            ax.bar(vals, freq.values, facecolor=mcolors.to_rgba(color, alpha=0.5), edgecolor=color,
                   linewidth=stylia.LINEWIDTH * 1.5, label=f"{label} (n={len(df) // 1000}k)")
        ax.set_xlim(lo - 0.5, hi + 0.5)
        stylia.label(ax, xlabel=PROP_XLABELS.get(prop, prop), ylabel="Frequency" if show_ylabel else "", abc=abc)
    else:
        for label, color_name, df in libraries:
            data = df[prop].dropna().values
            kde = gaussian_kde(data)
            x = np.linspace(lo, hi, 300)
            y = kde(x)
            color = nc.get(color_name)
            ax.plot(x, y, color=color, linewidth=stylia.LINEWIDTH * 1.5, label=f"{label} (n={len(df) // 1000}k)")
            ax.fill_between(x, y, color=color, alpha=0.5)
        ax.set_xlim(lo, hi)
        stylia.label(ax, xlabel=PROP_XLABELS.get(prop, prop), ylabel="Density" if show_ylabel else "", abc=abc)
        if prop == "MW":
            # Density values here are naturally small (peak ~0.0125) - fewer, rounder ticks
            # read better than the default 4-decimal-heavy tick set.
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))

    # Only the first panel keeps a y-axis label ("Density") - ticks/ticklabels are
    # dropped everywhere since the cross-panel y-scales aren't meant to be compared.
    ax.tick_params(axis="y", left=False, labelleft=False)

    if show_legend:
        ax.legend(loc="upper left", fontsize="small")


def plot_medians(ax):
    """Row 2: docking-score percentiles per pocket for all three libraries - top-1,000 (p1, solid)
    and top-100 (p0.1, dashed), sharing one x-axis ordered by HL's p1 ascending. Each library is
    colored the same as its row-1 panel (HL=crimson, REAL 10M=turquoise, REAL 10B=amber). A
    protein-identity color strip sits embedded at the very bottom of the same axes (not a separate
    axes below it), colored via figure 1's gene->color mapping so protein groupings stay visible
    without a 276-pocket x-tick label."""
    hl, by_library, gene_to_color = load_pocket_stats()
    n = len(hl)
    colors = [gene_to_color.get(gene, "#cccccc") for gene in hl["gene"]]

    nc = stylia.NamedColors()
    all_values = np.concatenate([
        by_library[label][col].values
        for label, _, _ in LIBRARIES for col in ("p1", "p0_1")
    ])
    all_values = all_values[~np.isnan(all_values)]
    lo, hi = all_values.min(), all_values.max()
    data_range = hi - lo
    strip_height = 0.08 * data_range
    strip_gap = 0.04 * data_range
    y_bottom = lo - strip_gap - strip_height
    y_top = hi + 0.1 * data_range

    ax.bar(range(n), [strip_height] * n, bottom=y_bottom, width=1, color=colors, edgecolor="none")

    lw = stylia.LINEWIDTH * 1.5
    for label, _, color_name in LIBRARIES:
        color = nc.get(color_name)
        lib_df = by_library[label]
        ax.plot(range(n), lib_df["p1"].values, color=color, linewidth=lw, linestyle="-",
                label=f"{label} - top 1,000")
        ax.plot(range(n), lib_df["p0_1"].values, color=color, linewidth=lw, linestyle="--",
                label=f"{label} - top 100")

    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(y_bottom, y_top)
    ax.set_xticks([])
    stylia.label(ax, xlabel="Pocket structures", ylabel="Docking score", abc="C")
    ax.legend(loc="upper left", fontsize="small", ncols=3)


def main():
    libraries = load_libraries()
    for label, _, df in libraries:
        print(f"{label}: {len(df):,} compounds")

    # height fraction chosen so 5 panels across the fixed Nature two-column width come out
    # roughly square, instead of stylia's multi-panel default (0.5 * width, calibrated for
    # fewer columns) which stretches each panel into a tall, narrow rectangle. Doubled here
    # (0.22 -> 0.44) so the median-score row below matches the property row's per-panel height —
    # create_figure()'s height is the *total* figure height, not a per-row height.
    fig, axs = stylia.create_figure(2, len(PROPERTIES), height=0.44)
    for i, prop in enumerate(PROPERTIES):
        ax = axs.next()
        # Panel B label sits on the first (MW) panel only - titles draw outside the axes
        # bounding box, so this doesn't shift its alignment with the other 4 panels in the row.
        plot_property(ax, prop, libraries, show_legend=(prop == "QED"), show_ylabel=(i == 0),
                      abc="B" if i == 0 else None)

    # Second row (panel C): merge the 5 individual cells into one axes spanning the full row
    # width (same x-space as the 5 property panels above), instead of 5 separate empty panels.
    row2_axes = [axs.next() for _ in range(len(PROPERTIES))]
    gridspec = row2_axes[0].get_gridspec()
    ax_bottom = row2_axes[0]
    for ax in row2_axes[1:]:
        ax.remove()
    ax_bottom.set_subplotspec(gridspec[1, :])
    plot_medians(ax_bottom)

    # stylia.save_figure() always calls plt.tight_layout() with no args right before saving,
    # which would silently undo any fig.subplots_adjust() call made beforehand - so tighten the
    # gap between the two rows AFTER tight_layout() runs, then save directly (matching
    # stylia.save_figure()'s own dpi/transparent/bbox_inches settings) instead of going through it.
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.3)

    output_path = os.path.join(plots_dir, "figure_2.png")
    plt.savefig(output_path, dpi=600, transparent=False, bbox_inches="tight")
    print(f"Saved to {output_path}")

    svg_path = os.path.join(plots_dir, "figure_2.svg")
    plt.savefig(svg_path, dpi=600, transparent=False, bbox_inches="tight")
    print(f"Saved to {svg_path}")

    # merge_with_scheme() needs its own SVG file on disk (svgutils.compose.SVG() can't embed a
    # matplotlib Figure object directly) - kept as a throwaway temp copy, separate from the
    # persisted figure_2.svg artifact above, and deleted once the merge is done.
    tmp_svg_fd, tmp_svg_path = tempfile.mkstemp(suffix=".svg")
    os.close(tmp_svg_fd)
    plt.savefig(tmp_svg_path, dpi=600, transparent=False, bbox_inches="tight")
    try:
        merge_with_scheme(tmp_svg_path)
    finally:
        os.remove(tmp_svg_path)


def white_background(width_mm, height_mm):
    """Opaque white backing rect - avoids relying on viewer/consumer background color."""
    rect = etree.Element(f"{{{SVG_NS}}}rect", {
        "x": "0", "y": "0", "width": str(width_mm), "height": str(height_mm), "fill": "white",
    })
    return Element(rect)


def merge_with_scheme(distributions_svg_path):
    """Combines panel A (scheme.svg, a hand-authored Inkscape schematic) with panels B+C (the
    property-distributions and docking-score rows just saved, already carrying their own "B"/"C"
    labels via stylia.label(abc=...)) into one figure, saved as PNG + PDF (not SVG - the
    composition itself is still built as vector SVG internally via svgutils, then rasterized/
    converted with cairosvg since that's the simplest way to get PNG/PDF out of an svgutils
    composition, which has no native non-SVG export).

    svgutils.compose.SVG() strips each source file's own <svg> wrapper and inlines its raw content,
    so embedding is always 1 (child viewBox unit) = 1 (parent user unit) - no implicit rescaling.
    The parent Figure below is declared directly in mm, so its user-unit space is literally
    millimeters; scheme.svg's own viewBox is already mm-equivalent (its declared width/height
    matches its viewBox numbers 1:1), so it needs no scaling. distributions_svg_path is matplotlib
    output (fix_mpl=True treats its declared pt-labelled width as bare viewBox units), so it's
    explicitly scaled to match the parent's mm-based width.
    """
    scheme_path = os.path.join(plots_dir, "scheme.svg")
    if not os.path.isfile(scheme_path):
        raise FileNotFoundError(
            f"{scheme_path} not found - scheme.svg is a hard requirement for the merged figure "
            "(hand-authored in Inkscape, not generated by this script)."
        )
    scheme = SVG(scheme_path)
    distributions = SVG(distributions_svg_path, fix_mpl=True)

    # scheme's own viewBox is already mm-equivalent, so this round-trips to the same numbers -
    # scheme itself needs no scaling.
    scheme_width_mm = Unit(scheme.width).to("mm").value
    scheme_height_mm = Unit(scheme.height).to("mm").value

    # distributions.width/.height are figure_2.svg's raw (post-fix_mpl) viewBox units, not mm -
    # scale so its width is a fraction of scheme's width in the shared mm-based parent canvas.
    scale = (scheme_width_mm * DISTRIBUTIONS_SIZE_FACTOR) / distributions.width
    distributions.scale(scale)
    distributions_width_mm = distributions.width * scale
    distributions_height_mm = distributions.height * scale

    panel_bc_x = (scheme_width_mm - distributions_width_mm) / 2  # center the smaller panel
    panel_bc_y = scheme_height_mm + PANEL_GAP_MM
    distributions.move(panel_bc_x, panel_bc_y)

    label_a = Text("A", LABEL_MARGIN_MM, LABEL_MARGIN_MM + LABEL_SIZE_MM, **LABEL_KWARGS)

    total_width_mm = scheme_width_mm
    total_height_mm = panel_bc_y + distributions_height_mm

    fig = Figure(
        f"{total_width_mm}mm",
        f"{total_height_mm}mm",
        white_background(total_width_mm, total_height_mm),
        scheme,
        distributions,
        label_a,
    )

    merged_svg_fd, merged_svg_path = tempfile.mkstemp(suffix=".svg")
    os.close(merged_svg_fd)
    png_path = os.path.join(plots_dir, "figure_2_merged.png")
    pdf_path = os.path.join(plots_dir, "figure_2_merged.pdf")
    try:
        fig.save(merged_svg_path)

        # Chrome renders mm/in-based documents at a fixed 96 CSS px/inch regardless of
        # --window-size, so --window-size must match that natural size (or content renders
        # shrunk into a corner of a mostly-blank canvas) - --force-device-scale-factor then
        # scales the whole render up to the real target resolution without leaving any padding.
        window_width = round(total_width_mm * CSS_DPI / 25.4)
        window_height = round(total_height_mm * CSS_DPI / 25.4)
        device_scale_factor = MERGED_DPI / CSS_DPI
        subprocess.run([
            "google-chrome", "--headless", "--disable-gpu", "--hide-scrollbars",
            f"--force-device-scale-factor={device_scale_factor}",
            f"--window-size={window_width},{window_height}",
            f"--screenshot={png_path}",
            "--default-background-color=FFFFFFFF",
            f"file://{merged_svg_path}",
        ], check=True, capture_output=True)

        Image.open(png_path).convert("RGB").save(pdf_path, "PDF", resolution=float(MERGED_DPI))
    finally:
        os.remove(merged_svg_path)

    print(f"Panel A (scheme): {scheme_width_mm:.1f} x {scheme_height_mm:.1f} mm")
    print(f"Panels B+C (distributions): {scheme_width_mm:.1f} x {distributions_height_mm:.1f} mm (scaled by {scale:.4f})")
    print(f"Total: {total_width_mm:.1f} x {total_height_mm:.1f} mm")
    print(f"Saved to {png_path}")
    print(f"Saved to {pdf_path}")


if __name__ == "__main__":
    main()
