import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

import itertools
import json
import pickle

import numpy as np
import pandas as pd
import stylia
from matplotlib.colors import to_hex
from scipy.spatial.distance import cosine

# Format: print — change with stylia.set_format()
stylia.set_format("print")

data_dir = os.path.join(root, "..", "..", "data")
output_dir = os.path.join(root, "..", "..", "output")
plots_dir = os.path.join(root, "..", "..", "plots", "figure_1")
os.makedirs(plots_dir, exist_ok=True)


def banner(title):
    line = "=" * (len(title) + 10)
    print(line)
    print(f"==== {title} ====")
    print(line)


uniprot_to_gene_df = pd.read_csv(os.path.join(data_dir, "mtb_trna_synthetases_bosch_2021_fig5.csv"))
uniprot_to_gene = {i: j for i, j in zip(uniprot_to_gene_df['uniprot_ac'], uniprot_to_gene_df['gene_name_in_bosch_2021'])}

uniprot_ids = ["P9WFW5", "P9WFW7", "P9WFW3", "P9WQA1", "P9WN61", "P9WFT5", "P9WFV3", "P9WFS9", "P9WFV1", "P9WFV9", "P9WFT9", "P9WFV7", "P9WFT7",
               "P9WFW1", "P9WFU5", "P9WFU9", "P9WFV5", "P9WFT3", "P9WFU3", "P9WFU1", "P9WFT1"]

# Single canonical ordering (alphabetical by gene name), reused for every section below —
# also what figure_1b_plot.py's SeqId/structural/PocketVec panels sort by, so all outputs
# from this script line up row-for-row.
uniprot_order = sorted(uniprot_ids, key=lambda uid: uniprot_to_gene[uid])
genes_sorted = [uniprot_to_gene[uid] for uid in uniprot_order]

# stylia SpectralColormap ("npg" preset), one color per protein, assigned in
# alphabetical gene-name order so the gradient lines up with PROTEINS in figure_1_plot.py.
palette = stylia.SpectralColormap("npg").sample(len(genes_sorted))
gene_to_color = {gene: to_hex(palette[i]) for i, gene in enumerate(genes_sorted)}

banner("LOADING BOSCH 2021 VULNERABILITY INDEX")
# Bosch et al. 2021, Supplementary Data S2, sheet "(1) Mtb H37Rv" — CRISPRi-derived
# Vulnerability Index (VI) per H37Rv gene, matched here on gene common name ("name").
vi_df = pd.read_excel(os.path.join(data_dir, "bosch_2021_DataS2.xlsx"), sheet_name="(1) Mtb H37Rv")
gene_to_vi_all = {i: j for i, j in zip(vi_df['name'], vi_df['Vulnerability Index'])}
gene_to_vi = {gene: gene_to_vi_all[gene] for gene in genes_sorted}
vi_genome_min = float(vi_df['Vulnerability Index'].min())
vi_genome_max = float(vi_df['Vulnerability Index'].max())
print(f"Matched Vulnerability Index for {len(gene_to_vi)}/{len(genes_sorted)} genes")
print(f"Genome-wide VI range (n={len(vi_df)}): [{vi_genome_min}, {vi_genome_max}]")

banner("COUNTING EXPERIMENTAL PDB STRUCTURES")
# data/structures/pdbe_database/<uniprot_ac>/<uniprot_ac>_updated/*_updated.cif — one file per
# experimental PDB entry (the sibling *_archive-PDB folder holds the same entries as legacy .ent
# files, so counting the "_updated" folder alone avoids double-counting).
pdbe_dir = os.path.join(data_dir, "structures", "pdbe_database")
gene_to_pdb_count = {}
for uid in uniprot_order:
    gene = uniprot_to_gene[uid]
    updated_dir = os.path.join(pdbe_dir, uid, f"{uid}_updated")
    if os.path.isdir(updated_dir):
        pdb_codes = {f.split("_updated")[0] for f in os.listdir(updated_dir) if f.endswith(".cif")}
        gene_to_pdb_count[gene] = len(pdb_codes)
    else:
        gene_to_pdb_count[gene] = 0
print(f"PDB counts by gene: {gene_to_pdb_count}")

banner("COUNTING CHEMBL BINDERS")
# data/ligands/chembl/<uniprot_ac>.json — IC50 activities fetched by script 07 for the
# 3 proteins with a ChEMBL target mapping (data/chembl_uniprot_mapping.txt). A "binder" is
# a unique molecule with an exact (standard_relation "=") IC50 <= 10 uM (10000 nM); genes
# with no ChEMBL target mapping at all have no entry here (rendered as "No data").
CHEMBL_BINDER_CUTOFF_NM = 10000
chembl_dir = os.path.join(data_dir, "ligands", "chembl")
gene_to_chembl_binders = {}
gene_to_chembl_total = {}
for uid in uniprot_order:
    gene = uniprot_to_gene[uid]
    chembl_file = os.path.join(chembl_dir, f"{uid}.json")
    if not os.path.exists(chembl_file):
        continue
    with open(chembl_file) as f:
        activities = json.load(f)
    chembl_df = pd.DataFrame(activities)
    total_molecules = chembl_df['molecule_chembl_id'].nunique()
    exact = chembl_df[chembl_df['standard_relation'] == '='].copy()
    exact['standard_value'] = pd.to_numeric(exact['standard_value'], errors='coerce')
    n_binders = exact[exact['standard_value'] <= CHEMBL_BINDER_CUTOFF_NM]['molecule_chembl_id'].nunique()
    gene_to_chembl_binders[gene] = int(n_binders)
    gene_to_chembl_total[gene] = int(total_molecules)
print(f"ChEMBL binders (IC50 <= {CHEMBL_BINDER_CUTOFF_NM} nM) by gene: {gene_to_chembl_binders}")
print(f"ChEMBL total unique molecules by gene: {gene_to_chembl_total}")

banner("DEDUPLICATING POCKETS PER PROTEIN")
# pocket_detection_data.csv has one row per pocket per structure, so a protein with 13
# structures can log the same physical pocket 13 times over. Structures are aligned into
# a shared coordinate frame (output/aligned_relaxed_structures), so pocket centroids are
# directly comparable across structures. Greedy dedup: sort by Pocket score descending,
# accept a pocket as "new" only if its centroid is farther than 6.14 A from every
# already-accepted centroid. 6.14 A matches the empirical same-pocket distance cutoff
# from notebooks/08_coherence_detected_pockets.ipynb (pairwise pocket-centroid distance
# analysis across aligned structures for the same protein).
# Computed once here — both the count (gene_to_unique_pocket_count, below) and the
# accepted rows themselves (uid_to_accepted_rows, reused for the PocketVec matrix later)
# come from this single pass.
POCKET_DEDUP_DISTANCE_THRESHOLD = 6.14
pocket_detection_data = pd.read_csv(os.path.join(output_dir, "pocket_detection_data.csv"))
gene_to_unique_pocket_count = {}
uid_to_accepted_rows = {}
for uid in uniprot_order:
    gene = uniprot_to_gene[uid]
    df = pocket_detection_data[pocket_detection_data['Uniprot AC'] == uid].sort_values('Pocket score', ascending=False)
    accepted_centroids = []
    accepted_rows = []
    for _, row in df.iterrows():
        centroid = np.array([float(v) for v in row['Pocket centroid coordinate (x y z)'].split()])
        if all(np.linalg.norm(centroid - c) > POCKET_DEDUP_DISTANCE_THRESHOLD for c in accepted_centroids):
            accepted_centroids.append(centroid)
            accepted_rows.append(row)
    gene_to_unique_pocket_count[gene] = len(accepted_rows)
    uid_to_accepted_rows[uid] = accepted_rows
print(f"Distinct pocket counts by gene: {gene_to_unique_pocket_count}")

banner("SAVING MAPPINGS")
mappings = {
    "uniprot_to_gene": uniprot_to_gene,
    "gene_to_color": gene_to_color,
    "gene_to_vulnerability_index": gene_to_vi,
    "vulnerability_index_genome_min": vi_genome_min,
    "vulnerability_index_genome_max": vi_genome_max,
    "gene_to_pdb_count": gene_to_pdb_count,
    "gene_to_chembl_binders": gene_to_chembl_binders,
    "gene_to_chembl_total": gene_to_chembl_total,
    "gene_to_unique_pocket_count": gene_to_unique_pocket_count,
}
color_mapping_path = os.path.join(plots_dir, "color_mapping.json")
with open(color_mapping_path, "w") as f:
    json.dump(mappings, f, indent=2)
print(f"Saved to {color_mapping_path}")

banner("STRUCTURAL RMSD MATRIX: pooling output/structural_comparisons/ CSVs (both directions)")
# Upper triangle: mean RMSD between the two proteins (pooled, both directions).
# Diagonal: min RMSD between structures of the same protein — not yet computed anywhere
# in the pipeline (scripts/15_calculate_StSim.py only does cross-protein pairs) and slow
# (~1600 superpositions), so left as a placeholder 0 for now.
# Lower triangle: min RMSD between the two proteins (pooled, both directions).
COMPARISONS_DIR = os.path.join(output_dir, "structural_comparisons")
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

n = len(uniprot_order)
struct_matrix = np.zeros((n, n))
for i, uid_i in enumerate(uniprot_order):
    for j, uid_j in enumerate(uniprot_order):
        if i == j:
            struct_matrix[i, j] = 0.0
        elif i < j:
            struct_matrix[i, j] = np.mean(pair_rmsds[frozenset((uid_i, uid_j))])
        else:
            struct_matrix[i, j] = np.min(pair_rmsds[frozenset((uid_i, uid_j))])

struct_matrix_df = pd.DataFrame(struct_matrix, index=uniprot_order, columns=uniprot_order)
struct_output_path = os.path.join(plots_dir, "structural_similarity_matrix.tsv")
struct_matrix_df.to_csv(struct_output_path, sep="\t")
print(f"Saved to {struct_output_path}")

banner("POCKETVEC MATRIX: min cosine distance over deduplicated pockets")
# Mirrors notebooks/16_PocketVec_analyses.ipynb's own convention: off-diagonal = min
# cosine distance between every pocket of protein i and every pocket of protein j;
# diagonal = min cosine distance between a protein's own distinct pockets (every protein
# here has >=2 deduplicated pockets, so this is always well-defined). Reuses
# uid_to_accepted_rows from the dedup pass above — no second dedup pass needed.
with open(os.path.join(output_dir, "pocketvec_RUN", "fps_rank.pkl"), "rb") as f:
    fps_rank = pickle.load(f)


def pocketvec_key(file_name, pocket_number):
    return f"{file_name.replace('.pdb', '')}_pocket_{pocket_number}"


uid_to_vectors = {}
for uid in uniprot_order:
    vectors = []
    for row in uid_to_accepted_rows[uid]:
        key = pocketvec_key(row["File name"], row["Pocket number"])
        if key not in fps_rank:
            print(f"  WARNING: {key} not found in fps_rank.pkl — skipping this pocket")
            continue
        vectors.append(fps_rank[key])
    uid_to_vectors[uid] = vectors
    print(f"{uniprot_to_gene[uid]} ({uid}): {len(vectors)} pockets with descriptors")

pocketvec_matrix = np.zeros((n, n))
for i, uid_i in enumerate(uniprot_order):
    for j, uid_j in enumerate(uniprot_order):
        if j < i:
            continue
        if i == j:
            pairs = itertools.combinations(uid_to_vectors[uid_i], 2)
        else:
            pairs = itertools.product(uid_to_vectors[uid_i], uid_to_vectors[uid_j])
        dist = min(cosine(va, vb) for va, vb in pairs)
        pocketvec_matrix[i, j] = dist
        pocketvec_matrix[j, i] = dist

pocketvec_matrix_df = pd.DataFrame(pocketvec_matrix, index=uniprot_order, columns=uniprot_order)
print(f"Matrix value range: min={pocketvec_matrix_df.values.min():.4f}, "
      f"max={pocketvec_matrix_df.values[~np.eye(n, dtype=bool)].max():.4f} (excluding diagonal)")
pocketvec_output_path = os.path.join(plots_dir, "pocketvec_similarity_matrix.tsv")
pocketvec_matrix_df.to_csv(pocketvec_output_path, sep="\t")
print(f"Saved to {pocketvec_output_path}")
