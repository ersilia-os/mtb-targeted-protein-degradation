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
plots_dir = os.path.join(root, "..", "..", "output", "plots", "figure_1")
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
# Live UniProt PDB cross-references, fetched by scripts/77_pocket_annotation/04_query_pdb_xrefs.py
# (output/77_pocket_annotation/pdb_xrefs.json). Supersedes the older local
# data/structures/pdbe_database/ snapshot, which undercounts: that snapshot predates several
# PDB entries deposited since it was downloaded (e.g. pheS/pheT: 7 in the stale local snapshot
# vs. 12 in a live UniProt query, confirmed 2026-08).
pdb_xrefs_path = os.path.join(output_dir, "77_pocket_annotation", "pdb_xrefs.json")
with open(pdb_xrefs_path) as f:
    pdb_xrefs = json.load(f)
gene_to_pdb_count = {uniprot_to_gene[uid]: len(pdb_xrefs[uid]["pdb_ids"]) for uid in uniprot_order}
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

banner("STRUCTURAL RMSD MATRIX: global (coverage-filtered min, upper) vs local (min, lower)")
# scripts/15c_calculate_StSim_global_local_coverage.py's output: every structure
# combination for each pair, both cmd.super (global) and cealign (local) RMSD, each with
# its own coverage (aligned_residues / min(res_1, res_2)). Files are named by UniProt AC
# STRING sort (ac_lo < ac_hi), NOT by uniprot_order's gene-name order, so that pair must
# be re-sorted independently for the file lookup below.
#
# Upper triangle: global (cmd.super) min RMSD, restricted to combinations with
# global_coverage >= 10% - a plain unfiltered min is misleading (e.g. cysS1/lysS hit
# 0.63 A from a ~5% coverage fluke, see conversation/notebook notes). 30% was tried
# first but discards 73/210 pairs entirely (every combination fails the cutoff); 10%
# only fails for 4 pairs (tyrS/glyS, proS/gatA, pheT/gatA, aspS/gatA - three involve
# gatA, which genuinely doesn't share fold architecture with the rest of the set) -
# those 4 cells are left as NaN (no reliable global comparison), not silently
# backfilled with the unfiltered value.
# Lower triangle: local (cealign) min RMSD, no coverage filter (cealign doesn't force a
# global correspondence, so it isn't prone to the same low-coverage/low-RMSD artifact).
# Diagonal: placeholder 0 (self-comparison data now exists in the same source files as
# self-pairs, but incorporating it here is a separate, not-yet-made decision).
GLOBAL_COVERAGE_MIN = 0.10
COMPARISONS_DIR = os.path.join(output_dir, "structural_comparisons_full")
pair_global_min = {}
pair_local_min = {}
no_reliable_global = []
for uid1, uid2 in itertools.combinations(uniprot_order, 2):
    ac_lo, ac_hi = sorted((uid1, uid2))
    path = os.path.join(COMPARISONS_DIR, f"{ac_lo}_{ac_hi}_global_local.csv")
    df = pd.read_csv(path)

    filtered = df[df["global_coverage"] >= GLOBAL_COVERAGE_MIN]
    if len(filtered):
        pair_global_min[frozenset((uid1, uid2))] = filtered["global_rmsd"].min()
    else:
        pair_global_min[frozenset((uid1, uid2))] = np.nan
        no_reliable_global.append((uniprot_to_gene[uid1], uniprot_to_gene[uid2]))

    pair_local_min[frozenset((uid1, uid2))] = df["local_rmsd"].min()

    print(f"{uniprot_to_gene[uid1]}/{uniprot_to_gene[uid2]}: {len(df)} combinations "
          f"(global min @>={GLOBAL_COVERAGE_MIN:.0%} cov={pair_global_min[frozenset((uid1, uid2))]:.2f}, "
          f"local min={pair_local_min[frozenset((uid1, uid2))]:.2f})")

print(f"Pairs with no combination clearing {GLOBAL_COVERAGE_MIN:.0%} global coverage "
      f"(left as NaN): {no_reliable_global}")

n = len(uniprot_order)
struct_matrix = np.full((n, n), np.nan)
for i, uid_i in enumerate(uniprot_order):
    for j, uid_j in enumerate(uniprot_order):
        if i == j:
            struct_matrix[i, j] = 0.0
        elif i < j:
            struct_matrix[i, j] = pair_global_min[frozenset((uid_i, uid_j))]
        else:
            struct_matrix[i, j] = pair_local_min[frozenset((uid_i, uid_j))]

struct_matrix_df = pd.DataFrame(struct_matrix, index=uniprot_order, columns=uniprot_order)
struct_output_path = os.path.join(plots_dir, "structural_similarity_matrix.tsv")
struct_matrix_df.to_csv(struct_output_path, sep="\t")
print(f"Saved to {struct_output_path}")

banner("SEQUENCE IDENTITY MATRIX: global (upper) vs local (lower)")
# Upper triangle: global (Needleman-Wunsch) identity, scripts/14_calculate_SeqId.py's
# SeqId_matrix.tsv. Lower triangle: local (Smith-Waterman) identity, the same script's
# LocalSeqId_matrix.tsv. Diagonal: self-identity (~100%, meaningfully computed by both
# aligners already - no placeholder needed, unlike the structural matrix's diagonal).
global_seqid_df = pd.read_csv(os.path.join(output_dir, "sequences", "NW_SeqAlign", "SeqId_matrix.tsv"), sep="\t", index_col=0)
global_seqid_df.columns = global_seqid_df.columns.str.strip()
global_seqid_df.index = global_seqid_df.index.str.strip()
local_seqid_df = pd.read_csv(os.path.join(output_dir, "sequences", "SW_LocalSeqAlign", "LocalSeqId_matrix.tsv"), sep="\t", index_col=0)
local_seqid_df.columns = local_seqid_df.columns.str.strip()
local_seqid_df.index = local_seqid_df.index.str.strip()

seqid_labels = [f"{uid} (tRNA)" for uid in uniprot_order]
global_seqid_matrix = global_seqid_df.loc[seqid_labels, seqid_labels].values
local_seqid_matrix = local_seqid_df.loc[seqid_labels, seqid_labels].values

n = len(uniprot_order)
seqid_matrix = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        if i == j:
            seqid_matrix[i, j] = global_seqid_matrix[i, j]
        elif i < j:
            seqid_matrix[i, j] = global_seqid_matrix[i, j]
        else:
            seqid_matrix[i, j] = local_seqid_matrix[i, j]

seqid_matrix_df = pd.DataFrame(seqid_matrix, index=uniprot_order, columns=uniprot_order)
seqid_output_path = os.path.join(plots_dir, "sequence_identity_matrix.tsv")
seqid_matrix_df.to_csv(seqid_output_path, sep="\t")
print(f"Saved to {seqid_output_path}")

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

banner("POCKET SCORES: P2Rank probability + domain-strip annotations (276 pockets)")
# Migrated from figure_3_calculations.py, along with figure_1_plot.py's own panel d (figure 3's
# panel c no longer needs this data). A pocket counts as "catalytic" at catalytic_confidence >= 3
# (either strong direct-PDB or strong AlphaFill ligand evidence for the Catalytic Domain
# (ATP/ligase) label - see scripts/77_pocket_annotation/10_assemble_final_table.py's confidence
# scale, 0-4), vs. the rest (weak/no evidence or no catalytic label at all, confidence 0-2) - a
# threshold the user specified directly, not chosen here.
CATALYTIC_CONFIDENCE_MIN = 3

# The other 3 domain rows (tRNA binding, Editing, Anticodon binding) have no ligand-evidence
# confidence score like catalytic_confidence - only a binary "this InterPro domain label is
# present for this pocket" (curated_labels). Red there means label present AND
# catalytic_confidence < CATALYTIC_CONFIDENCE_MIN - mutually exclusive with the Catalytic row,
# since a pocket can carry both a catalytic and a non-catalytic label at once.
CURATED_LABEL_COLUMNS = {
    "is_trna_binding": "tRNA Binding Domain",
    "is_editing": "Editing Domain",
    "is_anticodon_binding": "Anticodon Binding Domain",
}

pockets = pd.read_csv(os.path.join(output_dir, "pocket_detection_data.csv"))
pockets = pockets.sort_values("Pocket probability", ascending=False).reset_index(drop=True)

interpro = pd.read_csv(os.path.join(output_dir, "77_pocket_annotation", "pocket_detection_interpro_updated.csv"),
                        keep_default_na=False)
pockets = pockets.merge(
    interpro[["Uniprot AC", "File name", "Pocket number", "catalytic_confidence", "curated_labels"]],
    on=["Uniprot AC", "File name", "Pocket number"], how="left",
)

pocket_scores = pockets[["Uniprot AC", "File name", "Pocket number", "Pocket probability", "catalytic_confidence"]].copy()
pocket_scores.insert(0, "pocket_rank", range(1, len(pocket_scores) + 1))
pocket_scores = pocket_scores.rename(columns={"Pocket probability": "pocket_probability"})
pocket_scores["is_catalytic"] = pocket_scores["catalytic_confidence"] >= CATALYTIC_CONFIDENCE_MIN
has_label = {
    col: pockets["curated_labels"].apply(lambda labels, label=label: label in labels.split("|"))
    for col, label in CURATED_LABEL_COLUMNS.items()
}
for col, mask in has_label.items():
    pocket_scores[col] = mask & ~pocket_scores["is_catalytic"]

pocket_scores_path = os.path.join(plots_dir, "figure_1_pocket_scores.csv")
pocket_scores.to_csv(pocket_scores_path, index=False)
print(f"Saved {len(pocket_scores):,} row(s) to {pocket_scores_path}")
