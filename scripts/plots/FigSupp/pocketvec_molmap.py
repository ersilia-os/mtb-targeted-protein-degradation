"""
Fixed 1D ordering of the 128 PocketVec probe ligands (correlation distance -> hierarchical
clustering with optimal leaf ordering), used to render every pocket's 128-dim rank fingerprint
(output/pocketvec_RUN/fps_rank.pkl) as a single-row heatmap strip on a shared column order -
inspired by Fig. 5G of Offensperger et al., Science 2024 (eadk5864) and the MolMap method (Shen
et al., Nat Mach Intell 2021; github.com/shenwanxiang/bidd-molmap), adapted to 1D rather than
MolMap's 2D grid so that many pockets' strips can be stacked into one directly comparable
clustermap-style gallery. Complements, rather than replaces, pocketvec_tsne.py (which embeds
pockets, not probes): here every pocket shares the same probe column order, so two pockets'
strips are directly comparable column-by-column, and the probe legend shows which reference
ligand sits at any given column.

Same outlier definition as pocketvec_tsne.py / notebooks/16_PocketVec_analyses.ipynb: a pocket
is an outlier if >=80 of its 128 rank values are dummy/penalty ranks (raw value > 128). Outliers
are excluded only when building the shared column order / PCA / color range (so failed-docking
noise doesn't distort them) - they are still rendered in the gallery/tiles, marked with a red
border, since a "mostly failed docking" row is itself an informative result.

In addition to the full 128-probe view, also renders a "compressed" 32-feature view: PCA (fit
with StandardScaler on the 264 valid pockets, applied to all 276) of the same rank fingerprint.
Unlike the raw 128 probes, PCA components are already ordered by decreasing explained variance
and are mutually uncorrelated by construction, so no correlation-based reordering step applies
here - the natural column order (PC1..PC32, left to right) already is the informative one. Signed
component values use a diverging colormap (plum = negative, mint = positive) rather than the
sequential one used for ranks.

Usage:
    python pocketvec_molmap.py
"""
import os
import sys

os.environ["QT_QPA_PLATFORM"] = "offscreen"

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "..", "src"))

import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import stylia
from PIL import Image, ImageDraw
from rdkit import Chem
from rdkit.Chem import Draw
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from default import RANDOM_SEED
from pocket_group_plots import (
    CANONICAL_POCKET_MIN_COUNT,
    load_canonical_pocket_labels,
    load_domain_labels,
    load_gene_labels,
    sanitize,
)

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

output_dir = os.path.join(root, "..", "..", "..", "output")
plots_dir = os.path.join(output_dir, "plots", "FigSupp", "pocketvec_molmap")

fps_path = os.path.join(output_dir, "pocketvec_RUN", "fps_rank.pkl")
probe_sdf_path = os.path.join(output_dir, "pocketvec_RUN", "TOP128_rDock_LLM", "ALL", "all.sdf")
pockets_csv_path = os.path.join(output_dir, "pocket_detection_data.csv")

OUTLIER_MAX_PENALTY_RANKS = 80  # same threshold as pocketvec_tsne.py / notebooks/16
RAW_VMIN, RAW_VMAX = 1, 131  # fixed rank range incl. penalty ranks, shared across every row
PCA_N_COMPONENTS = 32
LEGEND_CELL_PX = 120  # legend cell size in pixels
ROW_HEIGHT_INCHES = 0.28  # per-row height for stacked galleries, so rows/labels stay legible
                          # regardless of group size (stylia's own height default is fixed,
                          # independent of nrows, and squeezes stacked rows illegible otherwise)


def gallery_height(n_rows):
    return max(0.05, n_rows * ROW_HEIGHT_INCHES / stylia.figure.get_size())


def load_fps():
    fps = pickle.load(open(fps_path, "rb"))
    n_penalty = {p: int(np.sum(fps[p] > len(fps[p]))) for p in fps}
    is_outlier = {p: n_penalty[p] >= OUTLIER_MAX_PENALTY_RANKS for p in fps}
    return fps, is_outlier


def build_probe_order(fps, is_outlier):
    valid_keys = sorted(p for p in fps if not is_outlier[p])
    X = np.array([fps[p] for p in valid_keys])  # (264, 128)
    corr = np.corrcoef(X.T)  # probe-probe Pearson correlation across the 264 valid pockets
    dist = 1 - corr
    np.fill_diagonal(dist, 0)
    dist = (dist + dist.T) / 2  # guard against asymmetric float noise
    np.fill_diagonal(dist, 0)

    Z = linkage(squareform(dist, checks=False), method="average", optimal_ordering=True)
    return leaves_list(Z)  # (128,) probe indices in an order where correlated probes are adjacent


def build_pca(fps, is_outlier):
    all_keys = sorted(fps)
    valid_keys = [k for k in all_keys if not is_outlier[k]]
    X_valid = np.array([fps[k] for k in valid_keys])
    X_all = np.array([fps[k] for k in all_keys])

    scaler = StandardScaler().fit(X_valid)
    pca = PCA(n_components=PCA_N_COMPONENTS, random_state=RANDOM_SEED).fit(scaler.transform(X_valid))
    pcs_all = pca.transform(scaler.transform(X_all))
    pcs = {key: pcs_all[i] for i, key in enumerate(all_keys)}

    # color range from the 264 valid pockets only, so outlier PCA projections (fit on noise)
    # don't wash out the color scale for everyone else - outliers simply clip to the extremes
    pcs_valid = pca.transform(scaler.transform(X_valid))
    max_abs = float(np.max(np.abs(pcs_valid)))
    return pcs, pca, max_abs


def make_row(values, order):
    return np.asarray(values)[order].reshape(1, -1)


def get_cmap_sequential():
    # stylia's FadingColormap.cmap always maps low->pale, high->deep regardless of the
    # `ascending` constructor arg (that arg only affects .transform(), not .cmap itself) -
    # reverse explicitly so low rank (strong predicted hit) renders as the salient deep color.
    cm = stylia.FadingColormap("plum")
    return cm.cmap.reversed()


def get_cmap_diverging():
    return stylia.DivergingColormap("plum_mint").cmap  # negative=plum, 0=white, positive=mint


def plot_row(ax, row, cmap, vmin, vmax, outlier=False):
    ax.imshow(row, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks([])
    ax.set_yticks([])
    if outlier:
        for spine in ax.spines.values():
            spine.set_edgecolor("red")
            spine.set_linewidth(2)
            spine.set_visible(True)


def save_single_tiles(values_dict, order, cmap, vmin, vmax, is_outlier, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for key in sorted(values_dict):
        fig, axs = stylia.create_figure(1, 1)  # wide strip, not square - omit width/height
        ax = axs.next()
        plot_row(ax, make_row(values_dict[key], order), cmap, vmin, vmax, outlier=is_outlier.get(key, False))
        title = key + (" (outlier)" if is_outlier.get(key, False) else "")
        stylia.label(ax, xlabel="", ylabel="", title=title)
        stylia.save_figure(os.path.join(out_dir, f"{key}.png"))
        plt.close(fig)


def key_to_short_label(pockets_df):
    labels = {}
    for _, row in pockets_df.iterrows():
        key = row["File name"].replace(".pdb", "") + "_pocket_" + str(row["Pocket number"])
        labels[key] = f"{row['Prediction type']}_p{row['Pocket number']}"
    return labels


def save_gallery(keys, values_dict, order, cmap, vmin, vmax, is_outlier, key_labels, out_path, min_count=1):
    keys = sorted(set(keys) & set(values_dict))
    if len(keys) < min_count:
        return
    fig, axs = stylia.create_figure(len(keys), 1, height=gallery_height(len(keys)))
    for key in keys:
        ax = axs.next()
        plot_row(ax, make_row(values_dict[key], order), cmap, vmin, vmax, outlier=is_outlier.get(key, False))
        label = key_labels.get(key, key) + (" *" if is_outlier.get(key, False) else "")
        # horizontal ax.text, not stylia.label's rotated ylabel - a rotated label per row
        # doesn't fit legibly once rows are this thin (many rows stacked in one gallery)
        ax.text(-0.01, 0.5, label, transform=ax.transAxes, ha="right", va="center",
                fontsize=stylia.FONTSIZE_SMALL)
        stylia.label(ax, xlabel="", ylabel="")
    stylia.save_figure(out_path)
    plt.close(fig)


def save_grouped_galleries(values_dict, order, cmap, vmin, vmax, is_outlier, key_labels, out_dir):
    keys = sorted(values_dict)

    genes, gene_to_color = load_gene_labels(output_dir, keys)
    gene_dir = os.path.join(out_dir, "gallery_by_gene")
    os.makedirs(gene_dir, exist_ok=True)
    for gene in sorted(set(genes)):
        gene_keys = [k for k, g in zip(keys, genes) if g == gene]
        save_gallery(gene_keys, values_dict, order, cmap, vmin, vmax, is_outlier, key_labels,
                     os.path.join(gene_dir, f"pocketvec_molmap_{sanitize(gene)}.png"))

    domains = load_domain_labels(output_dir, keys)
    domain_dir = os.path.join(out_dir, "gallery_by_domain")
    os.makedirs(domain_dir, exist_ok=True)
    for domain in sorted(set(domains)):
        domain_keys = [k for k, d in zip(keys, domains) if d == domain]
        save_gallery(domain_keys, values_dict, order, cmap, vmin, vmax, is_outlier, key_labels,
                     os.path.join(domain_dir, f"pocketvec_molmap_{sanitize(domain)}.png"))

    canonical = load_canonical_pocket_labels(output_dir, keys)
    canon_dir = os.path.join(out_dir, "gallery_by_canonical_pocket")
    os.makedirs(canon_dir, exist_ok=True)
    for group in sorted(set(canonical)):
        group_keys = [k for k, c in zip(keys, canonical) if c == group]
        save_gallery(group_keys, values_dict, order, cmap, vmin, vmax, is_outlier, key_labels,
                     os.path.join(canon_dir, f"pocketvec_molmap_{sanitize(group)}.png"),
                     min_count=CANONICAL_POCKET_MIN_COUNT)


def save_probe_legend(order, out_path):
    suppl = Chem.SDMolSupplier(probe_sdf_path)
    mols = list(suppl)
    assert len(mols) == 128, f"expected 128 probes in {probe_sdf_path}, got {len(mols)}"

    canvas = Image.new("RGB", (len(order) * LEGEND_CELL_PX, LEGEND_CELL_PX), "white")
    for col, probe_idx in enumerate(order):
        mol = mols[probe_idx]
        if mol is None:
            continue
        img = Draw.MolToImage(mol, size=(LEGEND_CELL_PX, LEGEND_CELL_PX))
        canvas.paste(img, (col * LEGEND_CELL_PX, 0))
    canvas.save(out_path)


def save_text_row_legend(labels, out_path, cell_px=LEGEND_CELL_PX):
    canvas = Image.new("RGB", (len(labels) * cell_px, cell_px), "white")
    draw = ImageDraw.Draw(canvas)
    for col, text in enumerate(labels):
        x0 = col * cell_px
        draw.rectangle([x0, 0, x0 + cell_px, cell_px], outline="black")
        draw.text((x0 + 10, 10), text, fill="black")
    canvas.save(out_path)


def main():
    fps, is_outlier = load_fps()
    n_excluded = sum(is_outlier.values())
    print(f"{len(fps)} pockets total, {n_excluded} excluded from order/PCA construction (outliers)")

    pockets_df = pd.read_csv(pockets_csv_path)
    key_labels = key_to_short_label(pockets_df)
    seq_cmap = get_cmap_sequential()

    # --- full 128-probe view ---
    raw_dir = os.path.join(plots_dir, "raw128")
    order = build_probe_order(fps, is_outlier)
    pickle.dump(order, open(os.path.join(plots_dir, "probe_order.pkl"), "wb"))
    save_single_tiles(fps, order, seq_cmap, RAW_VMIN, RAW_VMAX, is_outlier, os.path.join(raw_dir, "tiles"))
    save_grouped_galleries(fps, order, seq_cmap, RAW_VMIN, RAW_VMAX, is_outlier, key_labels, raw_dir)
    save_probe_legend(order, os.path.join(raw_dir, "probe_row_legend.png"))

    # --- compressed 32-feature PCA view ---
    pca_dir = os.path.join(plots_dir, "pca32")
    pcs, pca, max_abs = build_pca(fps, is_outlier)
    pca_order = np.arange(PCA_N_COMPONENTS)  # already meaningfully ordered (PC1..PC32), no reordering
    pickle.dump({"pca": pca, "max_abs": max_abs}, open(os.path.join(plots_dir, "pca32_model.pkl"), "wb"))

    var_ratio = pca.explained_variance_ratio_
    pd.DataFrame({"component": [f"PC{i + 1}" for i in range(PCA_N_COMPONENTS)],
                  "explained_variance_ratio": var_ratio,
                  "cumulative_variance_ratio": np.cumsum(var_ratio)}
                 ).to_csv(os.path.join(plots_dir, "pca32_explained_variance.csv"), index=False)
    print(f"PCA-32 cumulative explained variance: {np.cumsum(var_ratio)[-1]:.1%} "
          f"(of the 264 valid pockets' 128-dim rank fingerprint)")

    div_cmap = get_cmap_diverging()
    pcs_df = pd.DataFrame({k: v for k, v in pcs.items()}).T
    pcs_df.columns = [f"PC{i + 1}" for i in range(PCA_N_COMPONENTS)]
    pcs_df.to_csv(os.path.join(plots_dir, "pocketvec_pca32.csv"))

    save_single_tiles(pcs, pca_order, div_cmap, -max_abs, max_abs, is_outlier, os.path.join(pca_dir, "tiles"))
    save_grouped_galleries(pcs, pca_order, div_cmap, -max_abs, max_abs, is_outlier, key_labels, pca_dir)
    save_text_row_legend([f"PC{i + 1} ({v:.1%})" for i, v in enumerate(var_ratio)],
                          os.path.join(pca_dir, "pca32_row_legend.png"))


if __name__ == "__main__":
    main()
