"""
Extended Data Figure 4: overall shape of the per-pocket docking-score distributions across
the 3 screening libraries (HLL, REAL 9.56M, REAL 9.92B), stepping back from Figure 2c's
per-gene grid to show all 276 pockets pooled together. Reuses figure_2_calculations.py's
already-computed per-pocket percentiles (output/plots/figure_2/figure_2c_docking_percentiles.csv)
rather than recomputing anything.

(a) Median docking score, 3 distributions (HLL / REAL 9.56M / REAL 9.92B), 276 pockets each.
(b) P1 (1st-percentile) docking score, same 3 distributions.

Each distribution: a boxplot (styling ported from figure_1_plot.py's tax.boxplot(), library
colors matching figure_2_plot.py's LIBRARIES) plus a jittered strip of all 276 raw
per-pocket points on top.

Usage:
    python ExtendedDataFigure4.py
"""
import os
import sys

os.environ["QT_QPA_PLATFORM"] = "offscreen"

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "..", "src"))

import numpy as np
import pandas as pd
import stylia
from stylia.config import get_fg_color

from default import RANDOM_SEED

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

output_dir = os.path.join(root, "..", "..", "..", "output")
plots_dir = os.path.join(output_dir, "plots", "FigSupp", "ExtendedDataFigure4")
os.makedirs(plots_dir, exist_ok=True)

DOCKING_PERCENTILES_CSV = os.path.join(output_dir, "plots", "figure_2", "figure_2c_docking_percentiles.csv")

# (CSV library value, display label, NamedColors name) - same 3 libraries/colors as
# figure_2_plot.py's LIBRARIES, so colors stay consistent with Figure 2.
LIBRARY_DISPLAY = [
    ("HL", "HLL", "crimson"),
    ("REAL 10M", "REAL 9.56M", "turquoise"),
    ("REAL 10B", "REAL 9.92B", "amber"),
]

JITTER_WIDTH = 0.12
BOX_WIDTH = 0.5

# Panel letter sits above each panel's own full rendered extent (axes + tick labels +
# title), same tightbbox-based approach as ExtendedDataFigure1's add_panel_label.
PANEL_LABEL_Y_PAD = 0.015


def add_panel_label(fig, ax, letter):
    renderer = fig.canvas.get_renderer()
    bbox = ax.get_tightbbox(renderer)
    x0 = fig.transFigure.inverted().transform((bbox.x0, 0))[0]
    top_y = fig.transFigure.inverted().transform((0, bbox.y1))[1]
    fig.text(x0, top_y + PANEL_LABEL_Y_PAD, letter, fontweight="bold",
              fontsize=stylia.FONTSIZE_BIG, color=get_fg_color(), ha="left", va="bottom",
              transform=fig.transFigure)


def plot_metric_panel(ax, df, metric, ylabel, rng):
    nc = stylia.NamedColors()
    positions = range(len(LIBRARY_DISPLAY))
    box_data = [df.loc[df["library"] == lib, metric].values for lib, _, _ in LIBRARY_DISPLAY]

    # Dots drawn first (zorder=1), boxplot on top (zorder=3) with a transparent face so the
    # dots underneath a box remain visible.
    for i, (_, _, color_name) in enumerate(LIBRARY_DISPLAY):
        values = box_data[i]
        jitter = rng.uniform(-JITTER_WIDTH, JITTER_WIDTH, size=len(values))
        ax.scatter(np.full(len(values), i) + jitter, values, color=nc.get(color_name),
                   s=stylia.MARKERSIZE * 0.5, alpha=1, linewidth=0, zorder=1)

    bp = ax.boxplot(
        box_data, positions=positions, widths=BOX_WIDTH, showfliers=False, patch_artist=True,
        zorder=3,
        boxprops=dict(edgecolor="black", linewidth=stylia.LINEWIDTH),
        whiskerprops=dict(color="black", linewidth=stylia.LINEWIDTH),
        capprops=dict(color="black", linewidth=stylia.LINEWIDTH),
        medianprops=dict(color="black", linewidth=stylia.LINEWIDTH),
    )
    for box in bp["boxes"]:
        box.set_facecolor("none")

    ax.set_xticks(list(positions))
    ax.set_xticklabels([label for _, label, _ in LIBRARY_DISPLAY])
    stylia.label(ax, xlabel="", ylabel=ylabel)
    ax.set_box_aspect(1)


def main():
    df = pd.read_csv(DOCKING_PERCENTILES_CSV)
    rng = np.random.default_rng(RANDOM_SEED)

    data_path = os.path.join(plots_dir, "ExtendedDataFigure4_data.csv")
    df[["pocket", "gene", "library", "n", "median", "p1"]].to_csv(data_path, index=False)
    print(f"Saved {len(df)} row(s) to {data_path}")

    fig, axs = stylia.create_figure(1, 2, width=0.6, height=0.3)
    panel_axes = []

    ax_a = axs.next()
    plot_metric_panel(ax_a, df, "median", "Median docking score", rng)
    panel_axes.append(ax_a)

    ax_b = axs.next()
    plot_metric_panel(ax_b, df, "p1", "P1 docking score", rng)
    panel_axes.append(ax_b)

    # Needs a real renderer for add_panel_label's tightbbox measurements.
    fig.canvas.draw()
    for ax, letter in zip(panel_axes, ["a", "b"]):
        add_panel_label(fig, ax, letter)

    pdf_path = os.path.join(plots_dir, "ExtendedDataFigure4.pdf")
    png_path = os.path.join(plots_dir, "ExtendedDataFigure4.png")
    stylia.save_figure(pdf_path)
    stylia.save_figure(png_path)
    print(f"Saved {pdf_path}")
    print(f"Saved {png_path}")


if __name__ == "__main__":
    main()
