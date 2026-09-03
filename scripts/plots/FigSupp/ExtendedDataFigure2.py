"""
Extended Data Figure 2: relationship between local and global structural comparisons of tRNA
synthetases, classified by Class pairs. All structure pairs are used for comparison, and the
best (lowest) result (RMSD) is kept per protein pair. Structural pairs having a coverage lower
than 10% are discarded. A whole protein pair is discarded if all their associated structure
pairs are discarded for this reason.

Reconstructs a scatter plot originally produced ad hoc (never saved as a script) while tuning
figure_1_calculations.py's GLOBAL_COVERAGE_MIN=0.10 threshold for its structural RMSD matrix -
this figure shows the per-pair values behind that matrix's upper/lower triangles directly, split
by aaRS class combination.

Usage:
    python ExtendedDataFigure2.py
"""
import itertools
import json
import os
import sys

os.environ["QT_QPA_PLATFORM"] = "offscreen"

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "..", "src"))

import pandas as pd
import stylia

# Format: print | Style: article - matches ExtendedDataFigure1.py (publication supplementary figure)
stylia.set_format("print")
stylia.set_style("article")

output_dir = os.path.join(root, "..", "..", "..", "output")
plots_dir = os.path.join(root, "..", "..", "..", "output", "plots", "FigSupp", "ExtendedDataFigure2")
os.makedirs(plots_dir, exist_ok=True)

COMPARISONS_DIR = os.path.join(output_dir, "structural_comparisons_full")

# Canonical gene ordering/coloring, reused from figure_1_calculations.py so genes line up with
# every other figure in this project.
with open(os.path.join(output_dir, "plots", "figure_1", "color_mapping.json")) as f:
    figure_1_mappings = json.load(f)
uniprot_to_gene = figure_1_mappings["uniprot_to_gene"]

# Same aaRS Class I/II mapping as figure_1_plot.py (AARS_CLASS_LABELS) - gatA/gatB have no aaRS
# class (transamidases, not ligases) and are excluded from this figure entirely.
AARS_CLASS_LABELS = {
    "alaS": "Class II", "argS": "Class I", "aspS": "Class II", "cysS1": "Class I",
    "gltS": "Class I", "glyS": "Class II", "hisS": "Class II", "ileS": "Class I",
    "leuS": "Class I", "lysS": "Class II", "metS": "Class I", "pheS": "Class II",
    "pheT": "Class II", "proS": "Class II", "serS": "Class II", "thrS": "Class II",
    "trpS": "Class I", "tyrS": "Class I", "valS": "Class I",
}
classified_uniprots = sorted(uid for uid, gene in uniprot_to_gene.items() if gene in AARS_CLASS_LABELS)

# Same coverage cutoff and rationale as figure_1_calculations.py's structural RMSD matrix: a
# plain unfiltered min is misleading (low-coverage alignments can produce spuriously low RMSD),
# so the global (cmd.super) min is restricted to combinations clearing this coverage - a pair
# with no combination clearing it is dropped (NaN), not silently backfilled.
GLOBAL_COVERAGE_MIN = 0.10

CATEGORY_ORDER = ["Class I vs Class I", "Class II vs Class II", "Class I vs Class II"]


def pair_category(gene_1, gene_2):
    classes = sorted([AARS_CLASS_LABELS[gene_1], AARS_CLASS_LABELS[gene_2]])
    if classes[0] == classes[1]:
        return f"{classes[0]} vs {classes[0]}"
    return f"{classes[0]} vs {classes[1]}"


rows = []
dropped_pairs = []
for uid1, uid2 in itertools.combinations(classified_uniprots, 2):
    ac_lo, ac_hi = sorted((uid1, uid2))
    path = os.path.join(COMPARISONS_DIR, f"{ac_lo}_{ac_hi}_global_local.csv")
    df = pd.read_csv(path)

    gene_1, gene_2 = uniprot_to_gene[uid1], uniprot_to_gene[uid2]
    filtered = df[df["global_coverage"] >= GLOBAL_COVERAGE_MIN]
    if not len(filtered):
        dropped_pairs.append((gene_1, gene_2))
        continue

    global_rmsd = filtered["global_rmsd"].min()
    local_rmsd = df["local_rmsd"].min()
    category = pair_category(gene_1, gene_2)
    rows.append((gene_1, gene_2, category, global_rmsd, local_rmsd))

print(f"Pairs dropped for no combination clearing {GLOBAL_COVERAGE_MIN:.0%} global coverage: {dropped_pairs}")

pair_data = pd.DataFrame(rows, columns=["gene_1", "gene_2", "category", "global_rmsd", "local_rmsd"])
print(f"Category counts: {pair_data['category'].value_counts().to_dict()}")

data_path = os.path.join(plots_dir, "ExtendedDataFigure2_data.csv")
pair_data.to_csv(data_path, index=False)
print(f"Saved {len(pair_data)} row(s) to {data_path}")

# ===========================================================================
# Plot
# ===========================================================================
pal = stylia.CategoricalPalette("npg")
category_colors = dict(zip(CATEGORY_ORDER, pal.get(len(CATEGORY_ORDER))))


def plot_global_vs_local_rmsd(ax, pair_data):
    for category in CATEGORY_ORDER:
        subset = pair_data[pair_data["category"] == category]
        ax.scatter(subset["global_rmsd"], subset["local_rmsd"], color=category_colors[category],
                   label=f"{category} (n={len(subset)})")
    # Two-line title (fits over the narrower square box now that the legend sits outside it)
    stylia.label(ax, xlabel="Global RMSD, 10%-coverage-filtered min (Å)", ylabel="Local RMSD, min (Å)",
                 title="Global vs. local structural RMSD\nper protein pair, by Class I/II category")
    # box_aspect(1) forces a physically square plot box even with the legend occupying space
    # outside the axes to the right - without it, the layout engine steals width from the axes
    # to make room for the legend, stretching the box into a tall rectangle.
    ax.set_box_aspect(1)
    # Single-column legend, outside the axes to the right (not overlapping any data point).
    ax.legend(loc="upper left", bbox_to_anchor=(1.03, 1.0), frameon=False)


def main():
    fig, axs = stylia.create_figure(1, 1, width=0.5, height=0.5)
    ax = axs.next()
    plot_global_vs_local_rmsd(ax, pair_data)

    png_path = os.path.join(plots_dir, "ExtendedDataFigure2.png")
    stylia.save_figure(png_path)
    print(f"Saved {png_path}")


if __name__ == "__main__":
    main()
