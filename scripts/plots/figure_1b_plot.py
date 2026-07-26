import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

import json
import numpy as np
import pandas as pd
import stylia
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import gaussian_kde

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

output_dir = os.path.join(root, "..", "..", "output")
plots_dir = os.path.join(root, "..", "..", "plots", "figure_1")
os.makedirs(plots_dir, exist_ok=True)

# uniprot_to_gene / gene ordering computed by figure_1_calculations.py — reused here so
# gene order matches figure_1's composite grid.
with open(os.path.join(plots_dir, "color_mapping.json")) as f:
    mappings = json.load(f)
uniprot_to_gene = mappings["uniprot_to_gene"]
gene_to_color = mappings["gene_to_color"]
genes_sorted = sorted(uniprot_to_gene.values())

# scripts/14_calculate_SeqId.py's output: pairwise Needleman-Wunsch (BLOSUM62) sequence
# identity, computed over the 21 tRNA synthetases plus 25 background Mtb proteins. Only
# the 21 tRNA synthetases are kept for this heatmap.
seqid_df = pd.read_csv(os.path.join(output_dir, "sequences", "NW_SeqAlign", "SeqId_matrix.tsv"), sep="\t", index_col=0)
seqid_df.columns = seqid_df.columns.str.strip()
seqid_df.index = seqid_df.index.str.strip()

trna_labels = [label for label in seqid_df.index if "(tRNA)" in label]
uniprot_order = [label.replace(" (tRNA)", "") for label in trna_labels]
alphabetical_order = sorted(uniprot_order, key=lambda uid: uniprot_to_gene[uid])
alphabetical_labels = [f"{uid} (tRNA)" for uid in alphabetical_order]

matrix = seqid_df.loc[alphabetical_labels, alphabetical_labels].values
gene_labels = [uniprot_to_gene[uid] for uid in alphabetical_order]

# scripts/plots/figure_1c_calculations.py's output: pairwise structural comparison —
# upper triangle = mean RMSD, lower triangle = min RMSD (both pooled across all
# structure-file combinations for that protein pair), diagonal = placeholder 0 (same-
# protein comparison not yet computed). Same alphabetical order as the SeqId matrix
# above, so rows/columns line up between the two panels.
struct_df = pd.read_csv(os.path.join(plots_dir, "structural_similarity_matrix.tsv"), sep="\t", index_col=0)
struct_matrix = struct_df.loc[alphabetical_order, alphabetical_order].values

# scripts/plots/figure_1d_calculations.py's output: pocket-level PocketVec comparison,
# aggregated to protein level as min cosine distance between every deduplicated pocket
# of protein i and protein j (symmetric; diagonal = min distance between a protein's
# own distinct pockets). Same alphabetical order as the other panels.
pocketvec_df = pd.read_csv(os.path.join(plots_dir, "pocketvec_similarity_matrix.tsv"), sep="\t", index_col=0)
pocketvec_matrix = pocketvec_df.loc[alphabetical_order, alphabetical_order].values

CBAR_VMIN = 20
CBAR_VMAX = 40
STRUCT_VMIN = 0
STRUCT_VMAX = 25
POCKETVEC_VMIN = 0.12
POCKETVEC_VMAX = 0.20


def plot_comparison_heatmap(ax, matrix, labels, cmap, vmin, vmax, cbar_label, box_yticks, line_color, cbar_ticks, abc=None):
    im = ax.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax)

    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=stylia.FONTSIZE_SMALL)
    ax.set_yticklabels(labels, fontsize=stylia.FONTSIZE_SMALL)
    # Row labels move to the right side of the heatmap — the left side is now taken up
    # by the distribution + colorbar panels appended below.
    ax.yaxis.tick_right()

    divider = make_axes_locatable(ax)
    n = matrix.shape[0]

    # Top panel: one boxplot per protein, column-aligned with the heatmap, each showing
    # that protein's pairwise comparison values against the other proteins (self-
    # comparison excluded). whis=(0, 100) stretches whiskers to the actual min/max (no
    # outlier detection/fliers). No sharex here — boxplot()'s own tick-labeling clobbers
    # a shared axis's labels, so alignment is matched explicitly via xlim instead
    # (imshow's default extent for n categories is exactly [-0.5, n-0.5]).
    tax = divider.append_axes("top", size="20%", pad=0.05)
    box_data = [np.delete(matrix[:, i], i) for i in range(n)]
    bp = tax.boxplot(
        box_data, positions=range(n), widths=0.6, whis=(0, 100), showfliers=False, patch_artist=True,
        boxprops=dict(edgecolor="black", linewidth=stylia.LINEWIDTH),
        whiskerprops=dict(color="black", linewidth=stylia.LINEWIDTH),
        capprops=dict(color="black", linewidth=stylia.LINEWIDTH),
        medianprops=dict(color="black", linewidth=stylia.LINEWIDTH),
    )
    for i, box in enumerate(bp["boxes"]):
        box.set_facecolor(gene_to_color[labels[i]])
        box.set_alpha(0.8)
    tax.set_xlim(-0.5, n - 0.5)
    tax.set_ylim(vmin, vmax)
    tax.set_yticks(box_yticks)
    tax.tick_params(axis="x", labelbottom=False, bottom=False)
    tax.tick_params(axis="y", labelsize=stylia.FONTSIZE_SMALL, length=2, width=stylia.LINEWIDTH)
    for spine in ("top", "right"):
        tax.spines[spine].set_visible(False)
    # stylia.create_figure() applies this same thin/light grid to ax automatically, but
    # axes appended manually via make_axes_locatable don't get it — set it explicitly so
    # tax/lax match the heatmap's own grid instead of falling back to matplotlib's bolder default.
    tax.set_axisbelow(True)
    tax.grid(visible=True, linewidth=stylia.LINEWIDTH * 0.5, color="#DDDDDD", alpha=0.6)
    stylia.label(tax, xlabel="", ylabel="")

    # Colorbar axis carved out of the same divider as the heatmap axis (rather than
    # fig.colorbar's default shrink=<1 sizing), so its height matches the heatmap's
    # height/ylim exactly instead of floating shorter in the middle of the panel.
    # Appended first (before the distribution panel) so it's the one adjacent to the
    # heatmap; ticks face right, into the open gap left by ax's labels moving right.
    cax = divider.append_axes("left", size="5%", pad=0.20)
    cbar = ax.figure.colorbar(im, cax=cax, ticklocation="right")
    cbar.set_ticks(cbar_ticks)
    # Label lives on dax (further out) instead of here, so it reads as a single shared
    # label for the colorbar+distribution pair rather than being squeezed next to cax.
    # Explicit length/width — the default colorbar tick size renders disproportionately
    # thick/blocky at 600 dpi compared to the rest of the figure's line weight.
    cbar.ax.tick_params(labelsize=stylia.FONTSIZE_SMALL, length=2, width=stylia.LINEWIDTH)
    cbar.outline.set_linewidth(stylia.LINEWIDTH)

    # Distribution panel: single KDE of pairwise comparison values across all proteins
    # (upper-triangle values, i.e. one value per protein pair — excludes the diagonal
    # self-comparison). Appended after the colorbar, so it sits further left, outermost.
    # sharey=cax ties its value axis directly to the colorbar's own axis — no second,
    # redundant set of tick labels for the same quantity.
    dax = divider.append_axes("left", size="20%", pad=0.03, sharey=cax)
    pairwise_values = matrix[np.triu_indices(n, k=1)]
    kde = gaussian_kde(pairwise_values)
    y = np.linspace(vmin, vmax, 200)
    density = kde(y)
    # Baseline (zero density) sits against the colorbar side (dax's right edge); the
    # curve bulges away from it, toward the outer margin (dax's left edge).
    dax.plot(-density, y, color=line_color, linewidth=stylia.LINEWIDTH)
    # Fill under the curve with the same colormap/scale as the heatmap's colorbar
    # instead of a flat color, so the fill's color at a given height directly mirrors
    # what that value looks like on the heatmap/colorbar.
    gradient = np.linspace(vmin, vmax, 256).reshape(-1, 1)
    grad_im = dax.imshow(gradient, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto",
                          origin="lower", extent=[-density.max(), 0, vmin, vmax])
    clip_vertices = np.column_stack([np.concatenate([-density, [0, 0]]), np.concatenate([y, [y[-1], y[0]]])])
    grad_im.set_clip_path(Polygon(clip_vertices, closed=True, transform=dax.transData))
    dax.set_xlim(-density.max() * 1.02, 0)
    dax.set_xticks([])
    dax.tick_params(axis="y", labelleft=False, left=False)
    for spine in ("top", "right", "left"):
        dax.spines[spine].set_visible(False)
    dax.set_axisbelow(True)
    dax.grid(visible=True, linewidth=stylia.LINEWIDTH * 0.5, color="#DDDDDD", alpha=0.6)
    stylia.label(dax, xlabel="", ylabel=cbar_label)

    stylia.label(ax, xlabel="", ylabel="", abc=abc)


fig, axs = stylia.create_figure(3, 1, width=0.5, height=1.5)
nc = stylia.NamedColors()

seqid_cmap = stylia.FadingColormap("crimson", transformation=None).cmap
plot_comparison_heatmap(axs.next(), matrix, gene_labels, seqid_cmap, CBAR_VMIN, CBAR_VMAX,
                         "Sequence identity (%)", box_yticks=[30, 40], line_color=nc.crimson,
                         cbar_ticks=np.arange(CBAR_VMIN, CBAR_VMAX + 1, 5), abc="A")

# Different color family from panel A (cobalt vs crimson) so the two panels are
# visually distinct at a glance; still reversed so dark reads as "more similar" here
# too, even though for RMSD that means a LOW value.
struct_cmap = stylia.FadingColormap("cobalt", transformation=None).cmap.reversed()
plot_comparison_heatmap(axs.next(), struct_matrix, gene_labels, struct_cmap, STRUCT_VMIN, STRUCT_VMAX,
                         "Structural RMSD (Å)", box_yticks=[10, 20], line_color=nc.cobalt,
                         cbar_ticks=np.arange(STRUCT_VMIN, STRUCT_VMAX + 1, 5), abc="B")

# Third distinct color family (lime). Reversed for the same reason as panel B — for
# cosine distance, LOW value = more similar. box_yticks mark the method's own published
# lenient/stringent "similar pocket" thresholds (scripts/README.md: 0.17 / 0.14).
pocketvec_cmap = stylia.FadingColormap("lime", transformation=None).cmap.reversed()
plot_comparison_heatmap(axs.next(), pocketvec_matrix, gene_labels, pocketvec_cmap, POCKETVEC_VMIN, POCKETVEC_VMAX,
                         "PocketVec cosine distance", box_yticks=[0.14, 0.17], line_color=nc.lime,
                         cbar_ticks=np.arange(POCKETVEC_VMIN, POCKETVEC_VMAX + 0.001, 0.02), abc="C")

output_path = os.path.join(plots_dir, "figure_1b.png")
stylia.save_figure(output_path)
print(f"Saved to {output_path}")
