"""
Shared gene / canonical-pocket / domain labeling for pocket-level (keys, coords) embeddings,
plus the density-based point-sizing helpers - used by scripts/plots/figure_1_5_plot.py and
by the FigSupp molmap scripts (pocket_physchem_molmap.py, pocketvec_molmap.py).

Domain color choices (crimson/cobalt/amber/lime) match figure_1_5_plot.py's own
DOMAIN_STRIP_COLUMNS; the gene color mapping is output/plots/figure_1's canonical
gene_to_color.
"""
import json
import os

import pandas as pd
from scipy.stats import gaussian_kde

DENSITY_BANDWIDTH = 0.15  # same as notebooks/16_PocketVec_analyses.ipynb's calculate_density
MIN_SIZE, MAX_SIZE = 10, 70

CATALYTIC_CONFIDENCE_MIN = 3
DOMAIN_COLUMNS = [
    ("is_catalytic", "Catalytic", "crimson"),
    ("is_trna_binding", "tRNA binding", "cobalt"),
    ("is_editing", "Editing", "amber"),
    ("is_anticodon_binding", "Anticodon binding", "lime"),
]
CURATED_LABEL_NAMES = {
    "is_trna_binding": "tRNA Binding Domain",
    "is_editing": "Editing Domain",
    "is_anticodon_binding": "Anticodon Binding Domain",
}
CANONICAL_POCKET_MIN_COUNT = 5  # skip per-canonical-pocket highlight plots for smaller groups


def make_key(file_name, pocket_number):
    return file_name.replace(".pdb", "") + "_pocket_" + str(pocket_number)


def load_gene_labels(output_dir, keys):
    pockets = pd.read_csv(os.path.join(output_dir, "pocket_detection_data.csv"))
    pockets["key"] = [make_key(fn, pn) for fn, pn in zip(pockets["File name"], pockets["Pocket number"])]
    key_to_uniprot = dict(zip(pockets["key"], pockets["Uniprot AC"]))

    with open(os.path.join(output_dir, "plots", "figure_1", "color_mapping.json")) as f:
        mappings = json.load(f)
    uniprot_to_gene = mappings["uniprot_to_gene"]
    gene_to_color = mappings["gene_to_color"]

    genes = [uniprot_to_gene[key_to_uniprot[k]] for k in keys]
    return genes, gene_to_color


def load_canonical_pocket_labels(output_dir, keys):
    interpro = pd.read_csv(os.path.join(output_dir, "77_pocket_annotation", "pocket_detection_interpro_updated.csv"))
    interpro["key"] = [make_key(fn, pn) for fn, pn in zip(interpro["File name"], interpro["Pocket number"])]
    interpro["canonical_pocket"] = interpro["Gene"] + "_cluster" + interpro["spatial_cluster_id"].astype(str)
    key_to_canonical = dict(zip(interpro["key"], interpro["canonical_pocket"]))
    return [key_to_canonical[k] for k in keys]


def load_domain_labels(output_dir, keys):
    interpro = pd.read_csv(os.path.join(output_dir, "77_pocket_annotation", "pocket_detection_interpro_updated.csv"),
                            keep_default_na=False)
    interpro["key"] = [make_key(fn, pn) for fn, pn in zip(interpro["File name"], interpro["Pocket number"])]
    key_to_row = interpro.set_index("key")

    labels = []
    for k in keys:
        row = key_to_row.loc[k]
        if row["catalytic_confidence"] >= CATALYTIC_CONFIDENCE_MIN:
            labels.append("Catalytic")
            continue
        curated = row["curated_labels"].split("|")
        assigned = "None"
        for col, title, _ in DOMAIN_COLUMNS:
            if col == "is_catalytic":
                continue
            if CURATED_LABEL_NAMES[col] in curated:
                assigned = title
                break
        labels.append(assigned)
    return labels


def sanitize(label):
    return label.replace(" ", "_").replace("/", "-")


def calculate_density(coords, bandwidth=DENSITY_BANDWIDTH):
    xy = coords.T
    kde = gaussian_kde(xy, bw_method=bandwidth)
    return kde(xy)


def density_to_sizes(density, min_size=MIN_SIZE, max_size=MAX_SIZE):
    return min_size + (density - density.min()) / (density.max() - density.min()) * (max_size - min_size)

