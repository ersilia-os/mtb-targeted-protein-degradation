import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

import itertools
import json
import pickle

import numpy as np
import pandas as pd
from scipy.spatial.distance import cosine

output_dir = os.path.join(root, "..", "..", "output")
plots_dir = os.path.join(root, "..", "..", "plots", "figure_1")
os.makedirs(plots_dir, exist_ok=True)


def banner(title):
    line = "=" * (len(title) + 10)
    print(line)
    print(f"==== {title} ====")
    print(line)


# uniprot_to_gene / gene ordering computed by figure_1_calculations.py — reused here so
# this matrix's rows/columns line up with figure 1b's other panels.
with open(os.path.join(plots_dir, "color_mapping.json")) as f:
    mappings = json.load(f)
uniprot_to_gene = mappings["uniprot_to_gene"]
gene_to_unique_pocket_count = mappings["gene_to_unique_pocket_count"]
uniprot_order = sorted(uniprot_to_gene.keys(), key=lambda uid: uniprot_to_gene[uid])

# Same greedy pocket dedup as figure_1_calculations.py (by pocket-centroid distance,
# sorted by Pocket score descending) — reused here so the 93 pockets feeding this
# matrix are exactly the ones already reported as gene_to_unique_pocket_count.
POCKET_DEDUP_DISTANCE_THRESHOLD = 6.14
pocket_detection_data = pd.read_csv(os.path.join(output_dir, "pocket_detection_data.csv"))

with open(os.path.join(output_dir, "pocketvec_RUN", "fps_rank.pkl"), "rb") as f:
    fps_rank = pickle.load(f)


def pocketvec_key(file_name, pocket_number):
    return f"{file_name.replace('.pdb', '')}_pocket_{pocket_number}"


banner("DEDUPLICATING POCKETS PER PROTEIN AND MAPPING TO POCKETVEC DESCRIPTORS")
uid_to_vectors = {}
for uid in uniprot_order:
    gene = uniprot_to_gene[uid]
    df = pocket_detection_data[pocket_detection_data["Uniprot AC"] == uid].sort_values("Pocket score", ascending=False)
    accepted_centroids = []
    vectors = []
    for _, row in df.iterrows():
        centroid = np.array([float(v) for v in row["Pocket centroid coordinate (x y z)"].split()])
        if all(np.linalg.norm(centroid - c) > POCKET_DEDUP_DISTANCE_THRESHOLD for c in accepted_centroids):
            key = pocketvec_key(row["File name"], row["Pocket number"])
            if key not in fps_rank:
                print(f"  WARNING: {key} not found in fps_rank.pkl — skipping this pocket")
                continue
            accepted_centroids.append(centroid)
            vectors.append(fps_rank[key])
    uid_to_vectors[uid] = vectors
    expected = gene_to_unique_pocket_count[gene]
    match = "OK" if len(vectors) == expected else "MISMATCH"
    print(f"{gene} ({uid}): {len(vectors)} pockets with descriptors (expected {expected} from color_mapping.json) [{match}]")

banner("COMPUTING MIN COSINE DISTANCE — PROTEIN-LEVEL MATRIX")
# Mirrors notebooks/16_PocketVec_analyses.ipynb's own convention: off-diagonal = min
# cosine distance between every pocket of protein i and every pocket of protein j;
# diagonal = min cosine distance between a protein's own distinct pockets (every
# protein here has >=2 deduplicated pockets, so this is always well-defined).
n = len(uniprot_order)
matrix = np.zeros((n, n))
for i, uid_i in enumerate(uniprot_order):
    for j, uid_j in enumerate(uniprot_order):
        if j < i:
            continue
        if i == j:
            pairs = itertools.combinations(uid_to_vectors[uid_i], 2)
        else:
            pairs = itertools.product(uid_to_vectors[uid_i], uid_to_vectors[uid_j])
        dist = min(cosine(va, vb) for va, vb in pairs)
        matrix[i, j] = dist
        matrix[j, i] = dist

gene_labels = [uniprot_to_gene[uid] for uid in uniprot_order]
matrix_df = pd.DataFrame(matrix, index=uniprot_order, columns=uniprot_order)
print(f"Matrix value range: min={matrix_df.values.min():.4f}, max={matrix_df.values[~np.eye(n, dtype=bool)].max():.4f} (excluding diagonal)")

output_path = os.path.join(plots_dir, "pocketvec_similarity_matrix.tsv")
matrix_df.to_csv(output_path, sep="\t")
print(f"Saved to {output_path}")
