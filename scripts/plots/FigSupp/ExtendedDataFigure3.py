"""
Extended Data Figure 3: PocketVec 128-dim rank-fingerprint heatmap gallery, one thin row per
individual pocket structure belonging to the 3 reference canonical pockets highlighted in
Figure 1's panel e (TARGET_CLUSTERS: ileS_cluster1, gltS_cluster1, tyrS_cluster1 - grouped
into 3 blocks in that order: Pocket X/Y/Z). Same probe column order and fixed vmin/vmax=1/131
(the rank range including dummy/penalty ranks) as pocketvec_molmap.py's own
gallery_by_canonical_pocket - reuses that script's cached probe_order.pkl directly for exact
between-figure consistency, rather than recomputing the correlation-clustering ordering here.
POCKETVEC_CMAP (pheS-anchored, reversed) colors it - low rank/strong hit renders as the
salient deep color.

Split out of the former Figure 2's panel b (which paired this gallery with a 7-column
physicochemical-descriptor heatmap alongside it) - the physicochemical columns are dropped
here per request, keeping only the PocketVec descriptors.

FIGSIZE is a first-draft starting point, meant to be tuned iteratively by rendering and
looking, same convention as every other panel in this project.

Usage:
    python ExtendedDataFigure3.py
"""
import json
import os
import pickle
import sys

os.environ["QT_QPA_PLATFORM"] = "offscreen"

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "..", "src"))

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import stylia
from stylia.figure.figure import stylize

from pocket_group_plots import load_canonical_pocket_labels

stylia.set_format("print")
stylia.set_style("article")

output_dir = os.path.join(root, "..", "..", "..", "output")
plots_dir = os.path.join(root, "..", "..", "..", "output", "plots", "FigSupp", "ExtendedDataFigure3")
os.makedirs(plots_dir, exist_ok=True)

# Order is X, Y, Z (ileS, gltS, tyrS) - matches Figure 1 panel e's own TARGET_CLUSTERS order.
TARGET_CLUSTERS = ["ileS_cluster1", "gltS_cluster1", "tyrS_cluster1"]
CLUSTER_LABELS = {
    "ileS_cluster1": "Pocket $\\bf{X}$\n(ileS)",
    "gltS_cluster1": "Pocket $\\bf{Y}$\n(gltS)",
    "tyrS_cluster1": "Pocket $\\bf{Z}$\n(tyrS)",
}

# Canonical gene ordering/coloring, reused from figure_1_calculations.py so genes line up
# (and share colors) with every other figure in this project.
with open(os.path.join(output_dir, "plots", "figure_1", "color_mapping.json")) as f:
    figure_1_mappings = json.load(f)
gene_to_color = figure_1_mappings["gene_to_color"]


def cluster_color(cluster):
    gene = cluster.split("_cluster")[0]
    return gene_to_color[gene]


# PocketVec heatmap colormap - pale tint -> pheS's own canonical gene color, then reversed so
# low rank/strong hit renders as the salient deep color (rank 1 near vmin -> deep blue; a
# penalty rank near vmax -> pale/white).
_pheS_blue = np.array(mcolors.to_rgb(gene_to_color["pheS"]))
_pheS_blue_pale = 0.92 * np.array([1.0, 1.0, 1.0]) + 0.08 * _pheS_blue
POCKETVEC_CMAP = mcolors.LinearSegmentedColormap.from_list("pheS_blue", [_pheS_blue_pale, _pheS_blue])

FIGSIZE = (7.0, 2.3)


def build_gallery(fig):
    fps = pickle.load(open(os.path.join(output_dir, "pocketvec_RUN", "fps_rank.pkl"), "rb"))
    keys = sorted(fps)
    key_to_canonical = dict(zip(keys, load_canonical_pocket_labels(output_dir, keys)))

    order = pickle.load(open(
        os.path.join(output_dir, "plots", "FigSupp", "pocketvec_molmap", "probe_order.pkl"), "rb"))

    cmap = POCKETVEC_CMAP.reversed()
    vmin, vmax = 1, 131

    blocks = [(cluster, sorted(k for k in keys if key_to_canonical[k] == cluster)) for cluster in TARGET_CLUSTERS]
    block_sizes = [len(members) for _, members in blocks]

    outer_gs = fig.add_gridspec(len(blocks), 1, height_ratios=block_sizes, hspace=0.15)

    # No physchem column here (dropped per request) - just label, swatch, heatmap.
    LABEL_W, GAP_W, SWATCH_W, HEATMAP_W = 0.065, 0.02, 0.03, 0.865
    BLOCK_WIDTH_RATIOS = [LABEL_W, GAP_W, SWATCH_W, GAP_W, HEATMAP_W]

    for i, (cluster, members) in enumerate(blocks):
        color = cluster_color(cluster)
        block_gs = outer_gs[i, 0].subgridspec(1, 5, width_ratios=BLOCK_WIDTH_RATIOS, wspace=0)

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
        for j, key in enumerate(members):
            ax = stylize(fig.add_subplot(heatmap_gs[j, 0]))
            row = np.asarray(fps[key])[order].reshape(1, -1)
            ax.imshow(row, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
            ax.set_xticks([])
            ax.set_yticks([])
            stylia.label(ax, xlabel="", ylabel="", title="")


def main():
    fig = plt.figure(figsize=FIGSIZE)
    fig.patch.set_facecolor("white")
    build_gallery(fig)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.03)

    pdf_path = os.path.join(plots_dir, "ExtendedDataFigure3.pdf")
    png_path = os.path.join(plots_dir, "ExtendedDataFigure3.png")
    fig.savefig(pdf_path, dpi=600, transparent=False, bbox_inches="tight")
    fig.savefig(png_path, dpi=600, transparent=False, bbox_inches="tight")
    print(f"Saved {pdf_path}")
    print(f"Saved {png_path}")


if __name__ == "__main__":
    main()
