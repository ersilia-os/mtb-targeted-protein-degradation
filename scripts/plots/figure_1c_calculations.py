import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

import itertools
import json

import numpy as np
import pandas as pd

output_dir = os.path.join(root, "..", "..", "output")
plots_dir = os.path.join(root, "..", "..", "plots", "figure_1")
os.makedirs(plots_dir, exist_ok=True)

COMPARISONS_DIR = os.path.join(output_dir, "structural_comparisons")


def banner(title):
    line = "=" * (len(title) + 10)
    print(line)
    print(f"==== {title} ====")
    print(line)


# uniprot_to_gene / gene ordering computed by figure_1_calculations.py — reused here so
# this matrix's rows/columns line up exactly with figure 1b's SeqId heatmap.
with open(os.path.join(plots_dir, "color_mapping.json")) as f:
    mappings = json.load(f)
uniprot_to_gene = mappings["uniprot_to_gene"]
uniprot_order = sorted(uniprot_to_gene.keys(), key=lambda uid: uniprot_to_gene[uid])

banner("DIAGONAL: placeholder 0 RMSD")
# Same-protein min-RMSD (comparing each protein's own structure files pairwise via
# PyMOL) isn't computed anywhere in the pipeline yet (scripts/15_calculate_StSim.py only
# does cross-protein pairs) and is slow (~1600 superpositions) — placeholder for now.
diag_min = {uid: 0.0 for uid in uniprot_order}

banner("OFF-DIAGONAL: pooling existing output/structural_comparisons/ CSVs (both directions)")
pair_rmsds = {}
for uid1, uid2 in itertools.combinations(uniprot_order, 2):
    values = []
    for a, b in [(uid1, uid2), (uid2, uid1)]:
        path = os.path.join(COMPARISONS_DIR, f"{a}_{b}_rmsd.csv")
        if os.path.exists(path):
            values.extend(pd.read_csv(path)["rmsd"].tolist())
    if not values:
        raise FileNotFoundError(f"No structural comparison file found for {uid1}/{uid2}")
    pair_rmsds[frozenset((uid1, uid2))] = values
    print(f"{uniprot_to_gene[uid1]}/{uniprot_to_gene[uid2]}: {len(values)} pooled comparisons "
          f"(mean={np.mean(values):.2f}, min={np.min(values):.2f})")

banner("ASSEMBLING MATRIX")
# Upper triangle: mean RMSD between the two proteins (pooled, both directions).
# Diagonal: min RMSD between structures of the same protein.
# Lower triangle: min RMSD between the two proteins (pooled, both directions).
n = len(uniprot_order)
matrix = np.zeros((n, n))
for i, uid_i in enumerate(uniprot_order):
    for j, uid_j in enumerate(uniprot_order):
        if i == j:
            matrix[i, j] = diag_min[uid_i]
        elif i < j:
            matrix[i, j] = np.mean(pair_rmsds[frozenset((uid_i, uid_j))])
        else:
            matrix[i, j] = np.min(pair_rmsds[frozenset((uid_i, uid_j))])

gene_labels = [uniprot_to_gene[uid] for uid in uniprot_order]
matrix_df = pd.DataFrame(matrix, index=uniprot_order, columns=uniprot_order)
output_path = os.path.join(plots_dir, "structural_similarity_matrix.tsv")
matrix_df.to_csv(output_path, sep="\t")
print(f"Saved to {output_path}")
