import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

import json
import numpy as np
import pandas as pd
import stylia
from matplotlib.colors import to_hex

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

# stylia SpectralColormap ("npg" preset), one color per protein, assigned in
# alphabetical gene-name order so the gradient lines up with PROTEINS in figure_1_plot.py.
genes_sorted = sorted(uniprot_to_gene[uid] for uid in uniprot_ids)
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
for uid in uniprot_ids:
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
for uid in uniprot_ids:
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

banner("COUNTING DISTINCT POCKETS PER PROTEIN")
# pocket_detection_data.csv has one row per pocket per structure, so a protein with 13
# structures can log the same physical pocket 13 times over. Structures are aligned into
# a shared coordinate frame (output/aligned_relaxed_structures), so pocket centroids are
# directly comparable across structures. Greedy dedup: sort by Pocket score descending,
# accept a pocket as "new" only if its centroid is farther than 6.14 A from every
# already-accepted centroid. 6.14 A matches the empirical same-pocket distance cutoff
# from notebooks/08_coherence_detected_pockets.ipynb (pairwise pocket-centroid distance
# analysis across aligned structures for the same protein).
POCKET_DEDUP_DISTANCE_THRESHOLD = 6.14
pocket_detection_data = pd.read_csv(os.path.join(output_dir, "pocket_detection_data.csv"))
gene_to_unique_pocket_count = {}
for uid in uniprot_ids:
    gene = uniprot_to_gene[uid]
    df = pocket_detection_data[pocket_detection_data['Uniprot AC'] == uid].sort_values('Pocket score', ascending=False)
    centroids = [np.array([float(v) for v in c.split()]) for c in df['Pocket centroid coordinate (x y z)']]
    accepted = []
    for c in centroids:
        if all(np.linalg.norm(c - a) > POCKET_DEDUP_DISTANCE_THRESHOLD for a in accepted):
            accepted.append(c)
    gene_to_unique_pocket_count[gene] = len(accepted)
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
output_path = os.path.join(plots_dir, "color_mapping.json")
with open(output_path, "w") as f:
    json.dump(mappings, f, indent=2)
print(f"Saved to {output_path}")
