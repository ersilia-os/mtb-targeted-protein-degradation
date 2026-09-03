"""
Extended Data Figure 1: protein-level domain presence for Catalytic/tRNA binding/Editing/
Anticodon binding domains (panel a), per-residue AlphaFold2 pLDDT class breakdown as a
stacked bar (panel b), and structure inventory by source (panel c), across all 21 Mtb
tRNA synthetases.
All 3 panels share the same x-axis: genes sorted alphabetically, same order/colors as
figure_1_calculations.py's canonical gene ordering (output/plots/figure_1/color_mapping.json).

Unlike FigSupp/figure_supp_pocket.py's domain strips (pocket-level, mutually exclusive
with a catalytic-confidence threshold), panel a here is a protein-level annotation: a
domain counts as present if it's annotated anywhere in the protein at all, regardless of
whether a P2Rank pocket was ever detected there.

FIGSIZE and the row height ratios are a first-draft starting point, meant to be tuned
iteratively by rendering and looking, same convention as figure_supp_pocket.py.

Usage:
    python ExtendedDataFigure1.py
"""
import os
import sys

os.environ["QT_QPA_PLATFORM"] = "offscreen"

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "..", "src"))

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import stylia
from Bio.PDB import PDBParser
from matplotlib.colors import ListedColormap, to_hex
from matplotlib.transforms import Bbox
from stylia.figure.figure import stylize

stylia.set_format("print")
stylia.set_style("article")

data_dir = os.path.join(root, "..", "..", "..", "data")
output_dir = os.path.join(root, "..", "..", "..", "output")
plots_dir = os.path.join(root, "..", "..", "..", "output", "plots", "FigSupp", "ExtendedDataFigure1")
os.makedirs(plots_dir, exist_ok=True)

# Canonical gene ordering/coloring, reused from figure_1_calculations.py so genes line up
# (and share colors) with every other figure in this project.
with open(os.path.join(output_dir, "plots", "figure_1", "color_mapping.json")) as f:
    figure_1_mappings = json.load(f)
uniprot_to_gene = figure_1_mappings["uniprot_to_gene"]
gene_to_color = figure_1_mappings["gene_to_color"]
gene_to_uniprot = {gene: uid for uid, gene in uniprot_to_gene.items()}
GENES = sorted(gene_to_color.keys())

# ===========================================================================
# Panel c data: structure count per protein, split by source
# (output/trna_synthetases_data.csv - one row per structure; source is the file_name's
# leading prefix, e.g. "alphafold2_P9WFS9_model_0.pdb")
# ===========================================================================
SOURCE_LABELS = {
    "alphafold2": "AlphaFold2",
    "alphafold3": "AlphaFold3",
    "chai1": "Chai-1",
    "swissmodel": "SwissModel",
}
# stylia's npg categorical palette, greedy-farthest-point selection for max perceptual
# distinctness between the 4 sources (vs. SpectralColormap's continuous-hue sampling).
SOURCE_COLORS = {k: to_hex(v) for k, v in zip(SOURCE_LABELS.keys(), stylia.CategoricalPalette("npg").get(len(SOURCE_LABELS)))}

structures = pd.read_csv(os.path.join(output_dir, "trna_synthetases_data.csv"))
structures["source"] = structures["file_name"].str.extract(r"^([a-z0-9]+)_")
structure_counts = structures.groupby(["uniprot_ac", "source"]).size().unstack(fill_value=0)
structure_counts.index = structure_counts.index.map(uniprot_to_gene)
structure_counts = structure_counts.reindex(GENES).fillna(0).astype(int)[list(SOURCE_LABELS.keys())]

structure_counts_path = os.path.join(plots_dir, "ExtendedDataFigure1_structure_counts.csv")
structure_counts.to_csv(structure_counts_path)
print(f"Saved structure counts to {structure_counts_path}")

# ===========================================================================
# Panel a data: protein-level domain presence
# A domain is "present" for a protein if any InterPro accession in that protein's own
# scripts/77_pocket_annotation/<uid>_annotation_table_categorized.csv carries that
# category - independent of pocket detection (per request), unlike figure_supp_pocket.py's
# pocket-level strips. Uses the 77_pocket_annotation categorization (steps 02-03), NOT the
# older output/sequences/interpro_summary_curated.tsv: that older file's coverage<0.60
# filter (scripts/10_organize_sequence_info.py) drops large-footprint Class II aaRS
# catalytic domains before they can even be keyword-classified, silently showing no
# Catalytic Domain for alaS/aspS/pheS/glyS/hisS/lysS/serS - see
# scripts/77_pocket_annotation/README.md's "Bugs fixed vs. the old pipeline".
# ===========================================================================
DOMAIN_LABELS = {
    "Catalytic Domain (ATP/ligase)": "Catalytic",
    "tRNA Binding Domain": "tRNA binding",
    "Editing Domain": "Editing",
    "Anticodon Binding Domain": "Anticodon binding",
}
CATEGORIZED_TABLE_TEMPLATE = os.path.join(output_dir, "77_pocket_annotation", "{uid}_annotation_table_categorized.csv")

domain_presence = pd.DataFrame(index=list(DOMAIN_LABELS.values()), columns=GENES, dtype=int)
for gene in GENES:
    protein_categories = set(pd.read_csv(CATEGORIZED_TABLE_TEMPLATE.format(uid=gene_to_uniprot[gene]))["category"])
    for category, row_label in DOMAIN_LABELS.items():
        domain_presence.loc[row_label, gene] = int(category in protein_categories)

domain_presence_path = os.path.join(plots_dir, "ExtendedDataFigure1_domain_presence.csv")
domain_presence.to_csv(domain_presence_path)
print(f"Saved domain presence to {domain_presence_path}")

# ===========================================================================
# Panel b data: per-residue AlphaFold2 pLDDT
# One AF2 model per protein - pLDDT is read off the CA atom's B-factor column, where
# AlphaFold stores per-residue pLDDT (replicated across every atom in that residue).
# ===========================================================================
AF2_PDB_TEMPLATE = os.path.join(data_dir, "structures", "alphafold2_database", "{uid}", "AF-{uid}-F1-model_v4.pdb")

parser = PDBParser(QUIET=True)
plddt_by_gene = {}
for gene in GENES:
    uid = gene_to_uniprot[gene]
    structure = parser.get_structure(uid, AF2_PDB_TEMPLATE.format(uid=uid))
    plddt_by_gene[gene] = [res["CA"].get_bfactor() for res in structure[0].get_residues() if "CA" in res]

print(f"Extracted per-residue pLDDT for {len(plddt_by_gene)} proteins "
      f"(residue counts: { {g: len(v) for g, v in plddt_by_gene.items()} })")

# AlphaFold DB's own per-residue confidence bands and colors (Jumper et al. 2021, Nature;
# AlphaFold DB FAQ). pd.cut's own `labels` must list ascending bin edge order
# (PLDDT_BAND_EDGES_ASCENDING); PLDDT_BAND_LABELS is the separate display/stacking order
# (bottom-to-top, most to least confident) used by the plot and the fractions table below.
PLDDT_BAND_EDGES = [0, 50, 70, 90, 100]
PLDDT_BAND_LABELS_ASCENDING = ["Very low (<50)", "Low (50-70)", "Confident (70-90)", "Very high (>90)"]
PLDDT_BAND_LABELS = list(reversed(PLDDT_BAND_LABELS_ASCENDING))
PLDDT_BAND_COLORS = dict(zip(PLDDT_BAND_LABELS, ["#0053D6", "#65CBF3", "#FFDB13", "#FF7D45"]))

plddt_fractions = pd.DataFrame(index=GENES, columns=PLDDT_BAND_LABELS, dtype=float)
for gene in GENES:
    binned = pd.cut(plddt_by_gene[gene], bins=PLDDT_BAND_EDGES, labels=PLDDT_BAND_LABELS_ASCENDING, include_lowest=True)
    plddt_fractions.loc[gene] = pd.Series(binned).value_counts(normalize=True).reindex(PLDDT_BAND_LABELS)

plddt_fractions_path = os.path.join(plots_dir, "ExtendedDataFigure1_plddt_fractions.csv")
plddt_fractions.to_csv(plddt_fractions_path)
print(f"Saved pLDDT class fractions to {plddt_fractions_path}")

# ===========================================================================
# Plot
# ===========================================================================
FIGSIZE = (stylia.SIZE, 6.0)
# Panel order (top to bottom): a) domain presence, b) pLDDT stacked bar, c) structure
# counts - heights follow suit ([1, 2, 2]). A single gridspec hspace applies the same gap
# to every row, but both gaps here need extra room for the legend sitting atop panels b
# and c - so panels sit in rows 0, 2, 4 of a 5-row gridspec with hspace=0, and rows 1/3
# are blank spacer rows sized by GAP_AB/GAP_BC instead.
HEIGHT_RATIOS = [1, 2, 2]
GAP_AB = 0.7
GAP_BC = 0.7
BAR_WIDTH = 0.7

# Panel letters sit at a fixed figure-fraction x (near the page's own left edge, clearing
# every panel's y-axis label/ticks, whatever their width) and at the top of that panel's
# full rendered extent - axes + tick labels + y-label + its own legend, whichever sits
# highest - via get_tightbbox()/legend.get_window_extent() (needs a real renderer, hence
# the fig.canvas.draw() in main() before this runs), not just ax.get_position()'s bare
# axes box, since the top-of-panel legends below would otherwise sit above the letter.
PANEL_LABEL_X = 0.01
PANEL_LABEL_Y_PAD = 0.015
# Legend row, shared by every panel that has one: one horizontal row just above the axes
# (outside the plot, not inside/right), per request.
LEGEND_KWARGS = dict(frameon=False, fontsize=stylia.FONTSIZE_SMALL, loc="lower left",
                      bbox_to_anchor=(0.0, 1.0), borderaxespad=0.2, handletextpad=0.3, columnspacing=1.2)


def add_panel_label(fig, ax, letter):
    renderer = fig.canvas.get_renderer()
    bbox = ax.get_tightbbox(renderer)
    legend = ax.get_legend()
    if legend is not None:
        bbox = Bbox.union([bbox, legend.get_window_extent(renderer)])
    top_y = fig.transFigure.inverted().transform((0, bbox.y1))[1]
    fig.text(PANEL_LABEL_X, top_y + PANEL_LABEL_Y_PAD, letter,
              fontweight="bold", fontsize=stylia.FONTSIZE_BIG, ha="left", va="bottom",
              transform=fig.transFigure)


def plot_structure_counts(ax):
    x = np.arange(len(GENES))
    bottom = np.zeros(len(GENES))
    for source, color in SOURCE_COLORS.items():
        values = structure_counts[source].values
        ax.bar(x, values, bottom=bottom, width=BAR_WIDTH, color=color,
               label=SOURCE_LABELS[source], edgecolor="none", zorder=2)
        bottom += values
    ax.set_xlim(-0.5, len(GENES) - 0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(GENES, rotation=90)
    stylia.label(ax, xlabel="", ylabel="Number of structures")
    ax.legend(ncol=len(SOURCE_LABELS), **LEGEND_KWARGS)


# The "*" glyph's own ink sits high in its font line-box, so va="center" alone leaves it
# looking above the cell's true center - nudged down by eye (data units, half-integer cell
# grid) to compensate.
ASTERISK_Y_OFFSET = 0.08


# Matches stylize()'s own default axes.linewidth (0.5), so the inner grid and outer
# spines - both set to this width below - read as one consistent frame.
BORDER_LINEWIDTH = 0.5


def plot_domain_presence(ax):
    matrix = domain_presence.values.astype(float)
    nc = stylia.NamedColors()
    cmap = ListedColormap([nc.white, nc.cobalt])
    ax.imshow(matrix, cmap=cmap, vmin=0, vmax=1, aspect="auto", zorder=1)
    # Cell borders in black - visible against both the cobalt "present" fill and the
    # white "absent" background. Plain axhline/axvline rather than ax.grid(), since
    # stylia's "article" rcParams draw grid lines below patches (axisbelow) regardless of
    # an explicit zorder - invisible against imshow's own fill. The outer spines are set
    # to the same color/width so the surrounding frame reads as part of the same grid.
    for i in range(len(domain_presence.index) + 1):
        ax.axhline(i - 0.5, color="black", linewidth=BORDER_LINEWIDTH, zorder=3)
    for j in range(len(GENES) + 1):
        ax.axvline(j - 0.5, color="black", linewidth=BORDER_LINEWIDTH, zorder=3)
    for spine in ax.spines.values():
        spine.set_linewidth(BORDER_LINEWIDTH)
        spine.set_color("black")
    ax.set_xticks(np.arange(len(GENES)))
    ax.set_xticklabels(GENES, rotation=90)
    # gatA/gatB flagged in the Catalytic row, per request - GatA/GatB's catalytic mechanism
    # is transamidation (Glu-tRNA(Gln)/Asp-tRNA(Asn) + glutamine), not the aminoacyl-adenylate
    # ligase mechanism the other 19 targets share (see scripts/77_pocket_annotation/03_categorize.py).
    catalytic_row = list(domain_presence.index).index("Catalytic")
    for gene in ["gatA", "gatB"]:
        col = GENES.index(gene)
        text_color = nc.white if domain_presence.loc["Catalytic", gene] == 1 else nc.black
        ax.text(col, catalytic_row + ASTERISK_Y_OFFSET, "*", ha="center", va="center",
                 color=text_color, fontsize=stylia.FONTSIZE_BIG, fontweight="bold", zorder=3)
    ax.set_yticks(np.arange(len(domain_presence.index)))
    ax.set_yticklabels(domain_presence.index)
    stylia.label(ax, xlabel="", ylabel="")


def plot_plddt_stacked_bar(ax):
    x = np.arange(len(GENES))
    bottom = np.zeros(len(GENES))
    for label in PLDDT_BAND_LABELS:
        values = plddt_fractions[label].values.astype(float)
        ax.bar(x, values, bottom=bottom, width=BAR_WIDTH, color=PLDDT_BAND_COLORS[label],
               label=label, edgecolor="none", zorder=2)
        bottom += values
    ax.set_xlim(-0.5, len(GENES) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_xticks(x)
    ax.set_xticklabels(GENES, rotation=90)
    stylia.label(ax, xlabel="", ylabel="Fraction of residues\n(AlphaFold2 pLDDT)")
    ax.legend(ncol=len(PLDDT_BAND_LABELS), **LEGEND_KWARGS)


def main():
    fig = plt.figure(figsize=FIGSIZE)
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(5, 1, height_ratios=[HEIGHT_RATIOS[0], GAP_AB, HEIGHT_RATIOS[1], GAP_BC, HEIGHT_RATIOS[2]], hspace=0)
    ax_a = stylize(fig.add_subplot(gs[0, 0]))
    plot_domain_presence(ax_a)
    ax_b = stylize(fig.add_subplot(gs[2, 0]))
    plot_plddt_stacked_bar(ax_b)
    ax_c = stylize(fig.add_subplot(gs[4, 0]))
    plot_structure_counts(ax_c)

    # Needs a real renderer for add_panel_label's tightbbox/legend-extent measurements.
    fig.canvas.draw()
    for ax, letter in [(ax_a, "a"), (ax_b, "b"), (ax_c, "c")]:
        add_panel_label(fig, ax, letter)

    pdf_path = os.path.join(plots_dir, "ExtendedDataFigure1.pdf")
    png_path = os.path.join(plots_dir, "ExtendedDataFigure1.png")
    fig.savefig(pdf_path, dpi=600, transparent=False, bbox_inches="tight")
    fig.savefig(png_path, dpi=600, transparent=False, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {pdf_path}")
    print(f"Saved {png_path}")


if __name__ == "__main__":
    main()
