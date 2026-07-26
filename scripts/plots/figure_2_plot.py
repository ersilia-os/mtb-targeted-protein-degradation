import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

import numpy as np
import pandas as pd
import stylia
import matplotlib.colors as mcolors
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
from scipy.stats import gaussian_kde

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


def load_libraries():
    loaded = []
    for label, filename, color_name in LIBRARIES:
        df = pd.read_csv(os.path.join(plots_dir, filename))
        loaded.append((label, color_name, df))
    return loaded


def plot_property(ax, prop, libraries, show_legend, show_ylabel):
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
        stylia.label(ax, xlabel=PROP_XLABELS.get(prop, prop), ylabel="Frequency" if show_ylabel else "")
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
        stylia.label(ax, xlabel=PROP_XLABELS.get(prop, prop), ylabel="Density" if show_ylabel else "")
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


def main():
    libraries = load_libraries()
    for label, _, df in libraries:
        print(f"{label}: {len(df):,} compounds")

    # height fraction chosen so 5 panels across the fixed Nature two-column width come out
    # roughly square, instead of stylia's multi-panel default (0.5 * width, calibrated for
    # fewer columns) which stretches each panel into a tall, narrow rectangle. Doubled here
    # (0.22 -> 0.44) so the new empty row below matches the property row's per-panel height —
    # create_figure()'s height is the *total* figure height, not a per-row height.
    fig, axs = stylia.create_figure(2, len(PROPERTIES), height=0.44)
    for i, prop in enumerate(PROPERTIES):
        ax = axs.next()
        plot_property(ax, prop, libraries, show_legend=(prop == "QED"), show_ylabel=(i == 0))

    # Second row: merge the 5 individual cells into one axes spanning the full row width
    # (same x-space as the 5 property panels above), instead of 5 separate empty panels.
    row2_axes = [axs.next() for _ in range(len(PROPERTIES))]
    gridspec = row2_axes[0].get_gridspec()
    ax_bottom = row2_axes[0]
    for ax in row2_axes[1:]:
        ax.remove()
    ax_bottom.set_subplotspec(gridspec[1, :])
    stylia.label(ax_bottom, xlabel="", ylabel="")

    output_path = os.path.join(plots_dir, "figure_2.png")
    stylia.save_figure(output_path)
    print(f"Saved to {output_path}")

    svg_path = os.path.join(plots_dir, "figure_2.svg")
    stylia.save_figure(svg_path)
    print(f"Saved to {svg_path}")


if __name__ == "__main__":
    main()
