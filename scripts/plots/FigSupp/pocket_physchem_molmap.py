"""
Separate, independent counterpart to pocketvec_molmap.py: instead of the 128-dim PocketVec
docking fingerprint, this builds a fixed 1D ordering (correlation distance -> hierarchical
clustering with optimal leaf ordering, same recipe as pocketvec_molmap.py) over 9 human-readable
Tier-1 physicochemical pocket descriptors derived from residue composition - the pocket-level
analogue of the interpretable physicochemical descriptors used for fragment interpretability in
Offensperger et al., Science 2024 (eadk5864, Fig. 5G), adapted because molecule-style descriptors
(LogP, HBD/HBA, rotatable bonds) don't apply to a protein cavity, and rendered as a single-row
strip (rather than MolMap's 2D grid) so many pockets stack into one directly comparable gallery.

Descriptors (9, in a single row - no padding needed):
    pocket size (# residues), fraction hydrophobic (AVLIMFWC), fraction aromatic (FYW),
    fraction positively charged (KR), fraction negatively charged (DE), average hydropathy
    (Kyte & Doolittle, J Mol Biol 1982), average per-residue B-factor/pLDDT, P2Rank score,
    P2Rank probability.

Residue *identity* (not just chain_resnum) is not stored in pocket_detection_data.csv, so this
script parses each pocket's source structure once (Biopython, cached per file - 178 unique
structures for 276 pockets) to resolve chain_resnum -> three-letter amino acid code. Structure
paths are reconstructed as output/aligned_relaxed_structures/<Uniprot AC>/<File name> rather than
trusting the CSV's stale "Full path" column, which still has the pre-rename "processed/" prefix.

No outlier exclusion here (unlike pocketvec_molmap.py): the PocketVec docking-failure mode does
not apply to residue-composition descriptors, so all 276 pockets are used throughout.

Usage:
    python pocket_physchem_molmap.py
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
from Bio.PDB import PDBParser
from PIL import Image, ImageDraw
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

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
plots_dir = os.path.join(output_dir, "plots", "FigSupp", "pocket_physchem_molmap")
tiles_dir = os.path.join(plots_dir, "tiles")
os.makedirs(tiles_dir, exist_ok=True)

pockets_csv_path = os.path.join(output_dir, "pocket_detection_data.csv")
structures_dir = os.path.join(output_dir, "aligned_relaxed_structures")

ROW_HEIGHT_INCHES = 0.28  # per-row height for stacked galleries, see pocketvec_molmap.py

HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "CYS"}
AROMATIC = {"PHE", "TYR", "TRP"}
POSITIVE = {"LYS", "ARG"}
NEGATIVE = {"ASP", "GLU"}
# Kyte & Doolittle, J Mol Biol 157:105-132 (1982)
KYTE_DOOLITTLE = {
    "ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5, "GLN": -3.5, "GLU": -3.5,
    "GLY": -0.4, "HIS": -3.2, "ILE": 4.5, "LEU": 3.8, "LYS": -3.9, "MET": 1.9, "PHE": 2.8,
    "PRO": -1.6, "SER": -0.8, "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2,
}
DESCRIPTOR_NAMES = [
    "pocket size", "frac. hydrophobic", "frac. aromatic", "frac. positive", "frac. negative",
    "avg. hydropathy", "avg. B-factor/pLDDT", "P2Rank score", "P2Rank probability",
]


def gallery_height(n_rows):
    return max(0.05, n_rows * ROW_HEIGHT_INCHES / stylia.figure.get_size())


def make_key(file_name, pocket_number):
    return file_name.replace(".pdb", "") + "_pocket_" + str(pocket_number)


def parse_residue_map(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_path)
    resmap = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":  # skip heteroatoms/waters
                    continue
                resmap[(chain.id, residue.id[1])] = residue.resname
        break  # predicted structures here are single-model
    return resmap


def resolve_residue_types(chain_resn_str, resmap, pocket_key):
    resnames = []
    for token in chain_resn_str.split():
        chain, resnum = token.split("_", 1)
        try:
            resname = resmap[(chain, int(resnum))]
        except (KeyError, ValueError):
            print(f"WARNING: could not resolve residue {token} for {pocket_key}, skipping it")
            continue
        if resname not in KYTE_DOOLITTLE:
            print(f"WARNING: non-standard residue {resname} at {token} for {pocket_key}, skipping it")
            continue
        resnames.append(resname)
    return resnames


def compute_descriptors(resnames, bfactors, score, probability):
    n = len(resnames)
    frac = lambda s: sum(r in s for r in resnames) / n
    return np.array([
        n,
        frac(HYDROPHOBIC),
        frac(AROMATIC),
        frac(POSITIVE),
        frac(NEGATIVE),
        np.mean([KYTE_DOOLITTLE[r] for r in resnames]),
        np.mean(bfactors),
        score,
        probability,
    ])


def build_descriptor_matrix():
    pockets_df = pd.read_csv(pockets_csv_path)
    resmap_cache = {}
    keys, rows = [], []

    for _, row in pockets_df.iterrows():
        key = make_key(row["File name"], row["Pocket number"])
        pdb_path = os.path.join(structures_dir, row["Uniprot AC"], row["File name"])
        if pdb_path not in resmap_cache:
            resmap_cache[pdb_path] = parse_residue_map(pdb_path)

        chain_resn = row["Pocket residues (chain_resn)"]
        bfactors = [float(b) for b in row["B-factors"].split()]
        resnames = resolve_residue_types(chain_resn, resmap_cache[pdb_path], key)
        if not resnames:
            print(f"WARNING: no resolvable residues for {key}, dropping it from the physchem analysis")
            continue

        keys.append(key)
        rows.append(compute_descriptors(resnames, bfactors, row["Pocket score"], row["Pocket probability"]))

    X = np.array(rows)
    return keys, X


def build_descriptor_order(X):
    corr = np.corrcoef(X.T)  # descriptor-descriptor Pearson correlation across all pockets
    dist = 1 - corr
    np.fill_diagonal(dist, 0)
    dist = (dist + dist.T) / 2
    np.fill_diagonal(dist, 0)

    Z = linkage(squareform(dist, checks=False), method="average", optimal_ordering=True)
    return leaves_list(Z)  # (9,) descriptor indices, correlated ones adjacent


def normalize_columns(X):
    return (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))


def row_from_values(row_norm, order):
    return row_norm[order].reshape(1, -1)


def get_cmap():
    cm = stylia.FadingColormap("plum")  # here high (normalized) value = deep, no reversal needed
    return cm.cmap


def plot_row(ax, row, cmap):
    ax.imshow(row, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks([])
    ax.set_yticks([])


def save_single_tiles(keys, X_norm, order, cmap):
    for key, row_norm in zip(keys, X_norm):
        fig, axs = stylia.create_figure(1, 1)  # wide strip, not square - omit width/height
        ax = axs.next()
        plot_row(ax, row_from_values(row_norm, order), cmap)
        stylia.label(ax, xlabel="", ylabel="", title=key)
        stylia.save_figure(os.path.join(tiles_dir, f"{key}.png"))
        plt.close(fig)


def key_to_short_label(pockets_df):
    labels = {}
    for _, row in pockets_df.iterrows():
        key = make_key(row["File name"], row["Pocket number"])
        labels[key] = f"{row['Prediction type']}_p{row['Pocket number']}"
    return labels


def save_gallery(group_keys, key_to_row, cmap, key_labels, out_path, min_count=1):
    group_keys = sorted(set(group_keys) & set(key_to_row))
    if len(group_keys) < min_count:
        return
    fig, axs = stylia.create_figure(len(group_keys), 1, height=gallery_height(len(group_keys)))
    for key in group_keys:
        ax = axs.next()
        plot_row(ax, key_to_row[key], cmap)
        label = key_labels.get(key, key)
        # horizontal ax.text, not stylia.label's rotated ylabel - see pocketvec_molmap.py
        ax.text(-0.01, 0.5, label, transform=ax.transAxes, ha="right", va="center",
                fontsize=stylia.FONTSIZE_SMALL)
        stylia.label(ax, xlabel="", ylabel="")
    stylia.save_figure(out_path)
    plt.close(fig)


def save_grouped_galleries(keys, key_to_row, cmap, key_labels):
    genes, gene_to_color = load_gene_labels(output_dir, keys)
    gene_dir = os.path.join(plots_dir, "gallery_by_gene")
    os.makedirs(gene_dir, exist_ok=True)
    for gene in sorted(set(genes)):
        group_keys = [k for k, g in zip(keys, genes) if g == gene]
        save_gallery(group_keys, key_to_row, cmap, key_labels,
                     os.path.join(gene_dir, f"pocket_physchem_molmap_{sanitize(gene)}.png"))

    domains = load_domain_labels(output_dir, keys)
    domain_dir = os.path.join(plots_dir, "gallery_by_domain")
    os.makedirs(domain_dir, exist_ok=True)
    for domain in sorted(set(domains)):
        group_keys = [k for k, d in zip(keys, domains) if d == domain]
        save_gallery(group_keys, key_to_row, cmap, key_labels,
                     os.path.join(domain_dir, f"pocket_physchem_molmap_{sanitize(domain)}.png"))

    canonical = load_canonical_pocket_labels(output_dir, keys)
    canon_dir = os.path.join(plots_dir, "gallery_by_canonical_pocket")
    os.makedirs(canon_dir, exist_ok=True)
    for group in sorted(set(canonical)):
        group_keys = [k for k, c in zip(keys, canonical) if c == group]
        save_gallery(group_keys, key_to_row, cmap, key_labels,
                     os.path.join(canon_dir, f"pocket_physchem_molmap_{sanitize(group)}.png"),
                     min_count=CANONICAL_POCKET_MIN_COUNT)


def save_descriptor_legend(order):
    cell_px = 200
    canvas = Image.new("RGB", (len(order) * cell_px, cell_px), "white")
    draw = ImageDraw.Draw(canvas)
    for col, descriptor_idx in enumerate(order):
        text = DESCRIPTOR_NAMES[descriptor_idx]
        x0 = col * cell_px
        draw.rectangle([x0, 0, x0 + cell_px, cell_px], outline="black")
        # simple manual word-wrap so long descriptor names don't overflow one cell
        words, lines, line = text.split(), [], ""
        for w in words:
            trial = (line + " " + w).strip()
            if len(trial) > 14:
                lines.append(line)
                line = w
            else:
                line = trial
        lines.append(line)
        for i, l in enumerate(lines):
            draw.text((x0 + 10, 10 + i * 18), l, fill="black")
    canvas.save(os.path.join(plots_dir, "descriptor_row_legend.png"))


def main():
    keys, X = build_descriptor_matrix()
    pd.DataFrame(X, index=keys, columns=DESCRIPTOR_NAMES).to_csv(
        os.path.join(plots_dir, "pocket_physchem_descriptors.csv"))
    print(f"{len(keys)}/276 pockets resolved with usable residue identities")

    order = build_descriptor_order(X)
    pickle.dump(order, open(os.path.join(plots_dir, "descriptor_order.pkl"), "wb"))

    X_norm = normalize_columns(X)
    cmap = get_cmap()
    key_to_row = {key: row_from_values(row_norm, order) for key, row_norm in zip(keys, X_norm)}

    pockets_df = pd.read_csv(pockets_csv_path)
    key_labels = key_to_short_label(pockets_df)

    save_single_tiles(keys, X_norm, order, cmap)
    save_grouped_galleries(keys, key_to_row, cmap, key_labels)
    save_descriptor_legend(order)


if __name__ == "__main__":
    main()
