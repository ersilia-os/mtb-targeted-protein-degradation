"""
Supplementary figure: the 2 Enamine-library-docking t-SNE embeddings (HLL docking, REAL
9.92B docking) split out of figure_2_plot.py's panel a, which now shows only the
PocketVec embedding. Same 3 reference canonical pockets (gltS_cluster1, tyrS_cluster1,
ileS_cluster1) highlighted by their gene's canonical color, with bold-lettered callout
badges pinned to fixed corners of the HLL panel only (as before - both panels show the
same 3 pockets, so a 2nd legend on REAL 9.92B would be redundant).

Usage:
    python tSNE_Enamine.py
"""
import json
import os
import sys

os.environ["QT_QPA_PLATFORM"] = "offscreen"

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "..", "src"))

import numpy as np
import stylia

from default import RANDOM_SEED
from docking_tsne import compute_docking_tsne_embedding
from pocket_group_plots import calculate_density, density_to_sizes, load_canonical_pocket_labels

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

output_dir = os.path.join(root, "..", "..", "..", "output")
plots_dir = os.path.join(output_dir, "plots", "FigSupp", "tSNE_Enamine")
os.makedirs(plots_dir, exist_ok=True)

TARGET_CLUSTERS = ["gltS_cluster1", "tyrS_cluster1", "ileS_cluster1"]
# Mathtext ($\bf{...}$) bolds only the letter, not "Pocket" or the gene name.
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

EMBEDDINGS = [
    ("HLL docking", lambda: compute_docking_tsne_embedding(
        os.path.join(output_dir, "unidock_docking", "docking_results"), RANDOM_SEED)),
    ("REAL 9.92B docking", lambda: compute_docking_tsne_embedding(
        os.path.join(output_dir, "unidock_REAL_docking_2", "docking_results"), RANDOM_SEED)),
]


def plot_tsne_panel(ax, coords, canonical, title, annotate):
    nc = stylia.NamedColors()
    canonical = np.array(canonical)
    is_target = np.isin(canonical, TARGET_CLUSTERS)

    # "Other" points: light gray, sized by local 2D density; the 3 reference clusters keep
    # their own fixed color on top so they stay identifiable regardless of local density.
    sizes = density_to_sizes(calculate_density(coords), min_size=POINT_MIN_SIZE, max_size=POINT_MAX_SIZE)
    light_gray = nc.get("silver", lighten=0.5)

    # Explicit z-order (background < connector line < each cluster's own dots < caption
    # badge): the line must sit above the background but below the dots of the cluster it
    # points to.
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
            # (zorder=1) but below this cluster's own dots (zorder=3).
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


def main():
    fig, axs = stylia.create_figure(1, 2, width=0.6, height=0.3)
    for title, embed_fn in EMBEDDINGS:
        keys, coords = embed_fn()
        canonical = load_canonical_pocket_labels(output_dir, keys)
        plot_tsne_panel(axs.next(), coords, canonical, title, annotate=(title == "HLL docking"))
    stylia.save_figure(os.path.join(plots_dir, "tSNE_Enamine.png"))


if __name__ == "__main__":
    main()
