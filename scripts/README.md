# Scripts

One entry per script, in pipeline execution order. Explicit thresholds and cutoffs are stated with their rationale where the original design decision recorded one. See the root [README.md](../README.md) for the high-level pipeline overview and [CLAUDE.md](../CLAUDE.md) for repository-wide conventions.

## Structure preparation

### `00_generate_fasta_files.py`
Generates a FASTA file for each of the 21 tRNA synthetases (from `data/mtb_trna_synthetases_bosch_2021_fig5_annotated.csv`), used to query external structure/sequence resources.

### `01_download_alphafill.py`
Programmatically downloads AF2 structures with co-crystallized ligands from [AlphaFill](https://alphafill.eu/) into `/data/structures/alphafill_database/`.

### `02_organize_structures.py`
Organizes structure files from all sources (PDBe, AlphaFold2, AlphaFold3, Chai-1, AlphaFill, Swiss-Model) into `processed/structures`, stored in both `.cif` and `.pdb` formats. Ensures only one chain is saved per file and that sequences are not chunked. PDBe files are omitted from this automated processing.

### `03_align_structures.py`
Aligns all structures to simplify visualization. Based on RMSD, structures that seemed far apart from the rest were removed.

### `04_relax_structures.py`
Prepares structures for docking: protonation with PDB2PQR and relaxation with PyRosetta. Computationally intensive.

### `05_align_relaxed_structures.py`
Re-aligns structures after relaxation, using their unrelaxed counterparts as the alignment reference. No structures were removed at this stage, even those with high RMSD against the unrelaxed structures.

**Aggregate output (scripts 02–05):** `processed/trna_synthetases_data.csv` has one row per processed structure:

| **Field** | **Description** |
|-----------|------------------|
| `file_name` | Name of the processed PDB structure file |
| `uniprot_ac` | Uniprot AC identifier |
| `n_residues` | Number of residues |
| `start_resid` | First residue number (first residue is 1) of the sequence available in the PDB file, with respect to the Uniprot full sequence |
| `end_resid` | Last residue number (first residue is 1) of the sequence available in the PDB file, with respect to the Uniprot full sequence |
| `coverage` | Percentage sequence coverage |
| `sequence_structure` | Sequence found in the PDB file |
| `full_sequence` | Sequence found in UniProt |

## Sequence & pocket analysis

### `06_sequence_annotation.py`
Processes protein family/domain annotations downloaded from [InterPro](https://www.ebi.ac.uk/interpro/) (files in `data/sequences/interpro`).

### `07_fetch_from_chembl.py`
Fetches bioactivity data from [ChEMBL](https://www.ebi.ac.uk/chembl/) using UniProt AC identifiers. Data was found for only 3 of the 21 tRNA synthetases.

### `08_detect_pockets.py`
Detects pockets in AF2, AF3, Chai-1 and SwissModel predicted models. Following the [authors' recommendations](https://github.com/rdk/p2rank/issues/76): pockets are kept if probability (K) > 0.2, with at least the Top-3 (N) per structure retained regardless. Pockets with at least one residue below pLDDT < 65 (AF2, AF3, Chai-1) or QSQE < 0.65 (SwissModel) are then filtered out, discarding about 25–30% of pockets. These cut-offs are arbitrary; usual recommendations are 70 & 0.7 — this pipeline is slightly less restrictive. See also script 48, which applies the same P2Rank approach to experimental multimeric structures.

Summary output `processed/pocket_detection_data.csv` (one row per detected pocket and structure):

| **Field** | **Description** |
|-----------|------------------|
| `Uniprot_AC` | Uniprot AC identifier |
| `File name` | PDB file in which pockets have been detected |
| `Prediction type` | The method used for protein structure prediction |
| `Full path` | The full directory path where the PDB file is stored |
| `Pocket number` | The identified pocket number within the structure |
| `Pocket score` | The score assigned to the detected pocket |
| `Pocket probability` | The probability value indicating confidence in pocket detection |
| `Pocket centroid coordinate (x y z)` | The (x, y, z) coordinates of the pocket's centroid |
| `Pocket residues (chain_resn)` | List of residues forming the pocket, with chain and residue number |
| `B-factors` | Confidence measures: pLDDT (AF2, AF3, Chai) or QSQE (SM) |

### `09_prepare_pymol_visualizations.py`
Generates one PyMOL session file (`.pse`) per protein to visualize detected pockets and their residues:

| **Element** | **Description** | **Displayed?** |
|-------------|------------------|--------------|
| Reference structure (AF2) | Wheat color, surface + cartoon representation | Yes |
| Pockets detected in reference structure (AF2) | Sky-blue spheres with arbitrary size (pocket detection provides a single 3D coordinate) | Yes |
| Residues defining each pocket in AF2 | Orange color, surface + cartoon representation | Yes |
| Aligned structures (all but AF2) | Gray color, cartoon representation | No |
| Pockets detected in aligned structures | Gray-colored points | No |
| InterPro annotations | Conserved sites, domains, families, homologous superfamilies etc. (red color, surface representation) | No |

### `10_organize_sequence_info.py`
Organizes sequence annotation information. Output `processed/sequences/interpro_summary.tsv`:

| **Field** | **Description** |
|-----------|------------------|
| `Accession` | InterPro Accession identifier |
| `Name` | Name of the InterPro entry |
| `Type` | Type of annotation (e.g., domain, conserved site, homologous superfamily) |
| `Number of proteins` | Number of proteins associated with the InterPro entry |
| `Average Coverage` | Average sequence coverage across the associated proteins |

The file also contains binary (1/0) columns indicating presence/absence of each annotation per protein. A manually curated companion file (Catalytic domain / ATP binding site, Anticodon Binding Domain, etc.) is at `processed/sequences/interpro_summary_curated.tsv`.

### `11_organize_pocket_info.py`
Organizes pocket information together with InterPro annotations. Output `processed/pocket_detection_data_interpro.csv`:

| **Field** | **Description** |
|-----------|------------------|
| `Uniprot_AC` | Uniprot AC identifier |
| `File name` | PDB file in which pockets have been detected |
| `Prediction type` | The method used for protein structure prediction (e.g., AlphaFold2, AlphaFold3) |
| `Full path` | The full directory path where the PDB file is stored |
| `Pocket number` | The identified pocket number within the structure |
| `Pocket score` | The score assigned to the detected pocket |
| `Pocket probability` | The probability value indicating confidence in pocket detection |
| `Pocket centroid coordinate (x y z)` | The (x, y, z) coordinates of the pocket's centroid |
| `Pocket residues (chain_resn)` | List of residues forming the pocket, with chain and residue number |
| `Confidence score` | Confidence measures: pLDDT (AF2, AF3, Chai) or QSQE (SM) |
| `Interpro ID` | Identifier of the matched InterPro domain |
| `Interpro name` | Name of the matched InterPro domain |
| `Interpro Matches` | Residue ranges corresponding to the InterPro domain |
| `Residues in pocket` | Number of residues forming the pocket |
| `Residues in Interpro` | Number of residues forming the InterPro domain |
| `Residues overlap` | Number of residues shared between the pocket and the InterPro domain |
| `Coverage pocket` | Fraction of pocket residues that overlap with an InterPro domain |
| `Coverage domain` | Fraction of InterPro domain residues overlapping with the pocket |

Pocket–InterPro pairs with no overlapping residues are omitted from the file.

### `12_align_PDB_structures.py`
Aligns experimental structures (e.g., from PDBe) to predicted models to evaluate spatial coherence between structure sources. Also checks for overlaps between detected pockets and known ligand binding sites. Results in `processed/pdbe_annotation_report.csv`.

### `13_prepare_PocketVec.py`
*Not yet documented — pending individual review.*

### `14_calculate_SeqId.py`
Computes pairwise sequence identities between the 21 tRNA synthetases using global alignment (Needleman–Wunsch algorithm). Outputs a sequence identity matrix (`processed/sequences/NW_SeqAlign/SeqId_matrix.tsv`) and corresponding aligned proportions (`processed/sequences/NW_SeqAlign/Prop_matrix.tsv`).

### `15_calculate_StSim.py`
Evaluates structural similarity between tRNA synthetases using Pymol. Outputs CSV files of RMSDs between all structure pairs across all 21 tRNA synthetases, in `processed/structural_comparisons`.

### `16_calculate_PocketVec.py`
Calculates PocketVec descriptors ([Comajuncosa-Creus et al., Nat Commun 2024](https://www.nature.com/articles/s41467-024-52146-3)) from raw docking scores (docking calculations were run on an HPC cluster with CPU parallelization). 12/276 PocketVec descriptors were filtered out due to excessive presence of outlier compounds (>80), following the authors' recommendations. A PocketVec distance threshold of 0.17 is used to classify a pocket pair as similar. Descriptors are stored in `/processed/pocketvec_RUN/fps_rank.pkl`.

## Protein prioritization

### `17_protein_prioritization.py`
This script stubs out the prioritization pipeline; the logic itself is implemented in notebooks `17a_Protein_prioritization.ipynb`, `17b_Protein_prioritization_pairs.ipynb`/`17b_Protein_prioritization_triplets.ipynb`, and `17c_Protein_prioritization_individual.ipynb`. In brief: protein pairs (210) and triplets (1,330) were exhaustively enumerated and classified by PocketVec distance, global similarity (structural and sequential) and the number of pockets mapped to the catalytic domain. Global similarity thresholds: 35% sequence identity and 10 Å RMSD for sequential and structural similarity, respectively. Comparisons were extended to the pocket level (32,561 pairs, 2,499,258 triplets), collecting lenient (PocketVec distance < 0.17) and stringent (PocketVec distance < 0.14) sets of pairs (1,481 / 76) and triplets (3,880 / 807). An intra-set normalization then derives a "Prioritization score" per protein based on how often it appears across the collected sets — indicating similarity to other tRNA synthetases and potential polypharmacology.

Final output `processed/protein_prioritization/final_results.tsv`:

| **Field** | **Description** |
|-----------|------------------|
| `Uniprot_AC` | Uniprot AC identifier |
| `Gene name` | Standard gene name associated with the Uniprot AC identifier |
| `Vulnerability score` | Vulnerability score derived from [Bosch et al, 2021](https://pubmed.ncbi.nlm.nih.gov/34297925/) |
| `Score` | Prioritization score |
| `Tier` | Protein-associated tier |
| `sim_Tier1-5` | Number of proteins in the tier N having a PocketVec distance < 0.14 |

## 100k compound screening (Enamine Diversity Library HLL-100)

After structural characterization of the 21 essential tRNA synthetases, potential small-molecule binders are searched for in an active-learning fashion: first docking a chemically diverse 100k-compound set against each pocket structure, then training a surrogate model with [LazyQSAR](https://github.com/ersilia-os/lazy-qsar).

Computational resources used throughout this and later screening stages: `herbert` (default host unless noted), `norrsken-gpu-wsl` (NVIDIA GeForce RTX 4090), `SBNB-IRB cluster` (HPC, `/aloy/home`).

### `18_enamine_mfps.py`
Calculates Morgan Count Fingerprints for the 100,157 compounds of the ENAMINE Diversity Library HLL-100 (sublibrary of HLL-460). Stored in `/processed/enamine_characterization`.

### `19_enamine_conformations.py`
Generates 3D conformations (ETKDGv3 + UFF energy minimization) for the 100k library, yielding 100,154 conformations (`/processed/enamine_characterization/conformations.tar.gz`). 90% of these compounds fall in the [270, 400] Da range.

### `20_unidock_ligandprep.py`
Prepares the 100,154 3D-conformer compounds for docking via `unidocktools`' `ligandprep` (run on `norrsken-gpu-wsl`, `unidock` conda env). Output: `processed/unidock_docking/conformations_prepared`.

### `21_unidock_proteinprep.py`
Prepares the 178 protein structures (276 pocket structures, see `processed/pocket_detection_data.csv`) via `unidocktools`' `proteinprep` (run on `norrsken-gpu-wsl`, `unidock` conda env with `ambertools` and `reduce`). 8 structures required manual intervention to fix PDB2PQR/PyRosetta protonation issues. Output: `processed/unidock_docking/structures_prepared.tar.gz`.

### `22_unidock_docking.py`
Docks the 100k library against all 276 pocket structures with [Uni-Dock](https://github.com/dptech-corp/Uni-Dock), search mode _fast_, _vina_ scoring, on `norrsken-gpu-wsl` (~1 week of compute). Output: `processed/unidock_docking/docking_results` (`docking.tar.gz` per pocket directory, 276 total).

### `23_unidock_checks.py`
Systematically scans docking logs (`logs.tar.gz` files) for errors and warnings once docking (script 22) completes.

### `24_enamine_chemeleon.py`
Calculates [CheMeleon](https://github.com/JacksonBurns/chemeleon) embeddings for the same 100,157 compounds, using Ersilia's NVIDIA GeForce RTX 4090. Stored alongside the fingerprints from script 18 in `/processed/enamine_characterization`.

### `25_enamine_binarization.py`
Binarizes docking scores at 3 percentiles (0.1%, 0.5%, 1%) per pocket structure, producing reports mapping compound IDs, docking scores and the 3 bin flags. Output: `processed/unidock_docking/binarized_reports`.

### `26_train_models/`
Trains a [LazyQSAR](https://github.com/ersilia-os/lazy-qsar) ML model per pocket structure (x276) using the binarized data (script 25) and CheMeleon descriptors (script 24), at all 3 binarization percentiles. Run on the `SBNB-IRB cluster` for parallelization (3 bins × 276 pockets = 828 jobs). Models saved as joblib files in `processed/unidock_docking/models`. After analysis, only the 0.1%-binarized models were used downstream.

## Enamine REAL 9.56M: characterization, screening & prioritization

### `27_enamine_REAL_chemeleon.py`
Calculates CheMeleon embeddings for the 9.56M-compound Enamine REAL library using Ersilia's NVIDIA GeForce RTX 4090 (`norrsken-gpu-wsl`). Stored at `/processed/enamine_REAL_characterization/embeddings` (96 files, ~100k compounds each, ~25GB total).

### `28_transfer_embeddings.sh`
Transfers the embeddings from script 27 to the `/aloy/home` partition on the `SBNB-IRB cluster` via `rsync` through `dante`.

### `29_inference_enamine_REAL/`
Uses the trained surrogate models (script 26) to screen Enamine REAL 9.56M compounds, leveraging `SBNB-IRB cluster` parallelization. Output: `processed/unidock_docking/inferences`.

### `30_enamine_REAL_conformations.py`
Generates 3D conformations (ETKDGv3 + UFF energy minimization) for all Enamine REAL 9.56M compounds. Stored in pkl format per chunk under `/processed/enamine_REAL_characterization/conformations`.

### `31_enamine_REAL_checks_conformations.py`
Identifies compounds for which conformation generation (script 30) failed; IDs stored in `processed/enamine_REAL_characterization/failed_REAL.csv`. Also samples 25k compounds for subsequent analyses, stored in `processed/enamine_REAL_characterization/random_sample_REAL.csv`.

### `32_enamine_REAL_checks_inferences.py`
Systematically scans inference logs (`processed/unidock_docking/inferences_output`) for errors/warnings once script 29 completes. Only the inference probabilities for compounds that survived conformation generation (scripts 30, 31) are kept, stored in `processed/unidock_docking/inference_probs` — this is why the pipeline moves from script 29 to script 32 at this step rather than directly to selection.

### `33_enamine_REAL_selection.py`
Prioritizes a curated set of 100k compounds per pocket structure (x276), ranked by that pocket's inference probability, while diversifying by synthons (unweighting compounds with systematically high inference probabilities). Also defines a global "inactive" set (~13k compounds) whose synthons never appear in any pocket's active set. Combined list (100k+13k per pocket) stored in `processed/unidock_REAL_docking/input_ligands`.

### `34_enamine_REAL_mfps.py`
Calculates Morgan Count Fingerprints for all Enamine REAL 9.56M compounds, stored in `processed/enamine_REAL_characterization/enamine_REAL_ECFP6.h5`. Must be run with RDKit version 2025.9.1 — the ECFP6 bit ordering must match the version used in the related [ready-to-screen-enamine-real](https://github.com/ersilia-os/ready-to-screen-enamine-real) repository, which characterizes the 10-billion-compound library.

## Enamine REAL 9.56M: docking (100k prioritized compounds per pocket)

### `35_unidock_REAL_ligandprep.py`
Prepares the successfully-processed Enamine REAL 9.56M compounds (3D conformations) for docking via `unidocktools`' `ligandprep` (`norrsken-gpu-wsl`, `unidock` conda env). Output: `processed/unidock_REAL_docking/conformations_prepared`.

### `36_unidock_REAL_docking.py`
Docks the 276 pocket structures against ~113k compounds per pocket with Uni-Dock (_fast_ search, _vina_ scoring, `norrsken-gpu-wsl`, ~1 week compute). Output: `processed/unidock_REAL_docking/docking_results` (`docking.tar.gz` per pocket directory, 276 total). Protein structures reuse the preparation from script 21.

### `37_surrogate_model.py`
Trains a Naive Bayes Classifier per pocket structure (x276) using Morgan Count Fingerprints (RDKit 2025.9.1) at a 1% binarization threshold on docking scores (1,130 actives per pocket). Models stored at `processed/unidock_REAL_docking/training_results/models`.

## Enamine REAL 10B: characterization, screening & prioritization

Characterization (ECFP6 fingerprints for the full 10B library) was performed in the related [ready-to-screen-enamine-real](https://github.com/ersilia-os/ready-to-screen-enamine-real/tree/main) repository. Screening (surrogate models from script 37 applied to the full library, split into 994 chunks) was performed in the related [gcadda4tb-enamine-real-screening](https://github.com/ersilia-os/gcadda4tb-enamine-real-screening) repository. Per-chunk, per-pocket top-1% inference results are the starting point for the scripts below.

### `38_summarize_screening_results.py`
Checks whether compounds shared between pockets' top-1% hits reflect genuine multi-target/multi-protein binding or are an artifact of overlapping pocket geometry. For each chunk, shared hits between pocket pairs are classified as `SAME_POCKET` (centroid distance < 6.14 Å, i.e. effectively the same physical site detected twice), `SAME_PROTEIN` (different pocket, same protein) or `DIFF_PROTEIN`. Output: `processed/unidock_REAL_docking/inference_10B/shared_compounds/` (per chunk).

### `39_reduce_n_hits_I.py`
Per-chunk selection step of the 10B hit-reduction pipeline. A compound "hits" a target (pocket, or protein when any of its pockets hits) if it falls in the top 1% of surrogate-model-predicted scores for that target, computed independently within its 10M-compound chunk (994 chunks make up the 10B library) — a per-chunk, per-pocket relative cutoff, not a global rank or fixed absolute score. See script 40 for the cross-chunk aggregation and the two selection conditions applied.

### `40_reduce_n_hits_II.py`
Cross-chunk aggregation step, computed both at the pocket level (276 pockets) and protein level (21 proteins, pockets collapsed), using two complementary selection strategies:

* **Condition A ("promiscuous")**: compounds ranked by total number of targets hit, keeping the top 10,000 per chunk and reducing to the top **250,000** globally (both pocket- and protein-level sets).
* **Condition B ("selective"\*)**: compounds must first clear a promiscuity floor — hitting at least 50 pockets (pocket-level) or 2 proteins (protein-level) among their top-100 hits. Within that already-promiscuous eligible pool, the compounds with the *fewest additional* target hits are kept: 1,000 per pocket (**276,000** total) and 13,000 per protein (**273,000** total). \*The "selective" label means "least promiscuous among an already-promiscuous pool," not "hits few targets" — it follows the naming used in the pipeline scripts themselves, which flag the term with an asterisk.

A fixed random seed (42) is used for tie-breaking throughout. Final sets: `A_pockets.csv`, `B_pockets.csv`, `A_proteins.csv`, `B_proteins.csv` in `processed/unidock_REAL_docking/inference_10B/`.

### `41_download_and_map.py`
Since the full 10B-compound SMILES mapping is too large to store locally, downloads per-chunk ID→SMILES files from a shared Google Drive folder (via a Google service account), merges the four selection sets from script 40 (deduplicating overlaps), and maps each selected compound to its SMILES string. Output: `processed/unidock_REAL_docking/inference_10B/selected_compounds/` (per chunk).

### `42_annotate_and_filter.py`
Computes physicochemical descriptors (MW, LogP, QED) with RDKit and merges in a natural-product/synthetic-likeness score (NSPS) from the Ersilia Model Hub model `eos12x7`. Filters to a drug-like chemical space: MW in [250, 450], LogP in [−1, 5], QED > 0.4, NSPS in [10, 40]. Filtered compounds (~965k) stored in `processed/unidock_REAL_docking/inference_10B/filtered_compounds.csv`.

### `43_clustering.ipynb`
Reduces the ~965k filtered compounds to a final validation set. First applies a synthon-diversity cap (max 3 occurrences per synthon, analogous to script 33), then clusters the remaining compounds' ECFP6 fingerprints with [BitBirch](https://github.com/mqcomplab/bitbirch) (`bblean` package; similarity threshold 0.3, branching factor 50, 5 recluster iterations), keeping one representative per cluster. Final set: **~99k** compounds, stored in `processed/unidock_REAL_docking/inference_10B/clustered_compounds.csv`.

## Final validation docking and hit selection

### `44_generate_conformations.py`
Generates 3D conformations (ETKDGv3 + UFF energy minimization) for the ~99k clustered compounds. Output: `processed/unidock_REAL_docking_2/conformations/`.

### `45_unidock_REAL_2_ligandprep.py`
Prepares the ~99k compound conformations for docking via `unidocktools`' `ligandprep` (`norrsken-gpu-wsl`, `unidock` conda env). Output: `processed/unidock_REAL_docking_2/conformations_prepared/`.

### `46_unidock_REAL_2_docking.py`
Docks the 276 pocket structures against the ~99k compound set with Uni-Dock (_fast_ search, _vina_ scoring, `norrsken-gpu-wsl`). Output: `processed/unidock_REAL_docking_2/docking_results`, one `report.csv` per pocket.

### `47_docking_summary.py`
Console reporting tool, per gene, summarizing for every candidate pocket of a given tRNA synthetase: P2Rank pocket probability, InterPro domain annotation, AlphaFill co-crystallized-ligand evidence, structural confidence (minimum pLDDT for AF2/AF3/Chai-1 models, or GMQE for SwissModel models), and docking score percentiles (0.01%, 0.1%, 1%) from both docking libraries (100k Enamine DL and the ~99k Enamine REAL set). Used to manually curate one reference pocket per tRNA synthetase, recorded in `output/reference_pocket.csv`, which is consumed by downstream hit-selection scripts.

### `48_detect_pocket_multimers.py`
Complements the UniProt-AC-keyed, single-chain pipeline (script 08) with a standalone script keyed by PDB code, characterizing pockets in real, potentially multi-chain experimental structures (e.g. to check whether a candidate pocket sits at a subunit interface). Invoked as `python 48_detect_pocket_multimers.py --pdb-codes 6XYZ,7ABC`. For each PDB code:

- Downloads the RCSB-annotated **biological assembly** (not just the as-deposited asymmetric unit) live from `files.rcsb.org`, falling back to the asymmetric unit if no assembly is annotated. Since RCSB entries and annotations can be revised, the aggregate report records the access date per row.
- Strips ligands, waters and other heteroatoms while keeping **all** protein chains (unlike script 02, which reduces structures to a single chosen chain).
- Detects pockets with P2Rank, as script 08 does, but using P2Rank's **default** config instead of `-c alphafold` (these are experimental, not AlphaFold-predicted, structures), applying only the probability/rank filter (probability ≥ 0.2 or Top-3, unchanged from script 08). Script 08's per-residue B-factor/pLDDT confidence gate is **not** applied here: the PDB B-factor column encodes crystallographic atomic displacement, not prediction confidence, so that gate's semantics don't transfer to experimental structures.

Unlike script 04, this script does **not** run PDB2PQR protonation or PyRosetta relaxation: experimental structures are already validated against experimental data by crystallographic/cryo-EM refinement, so an unconstrained relax risks moving pocket residues away from their observed conformation. Pocket detection runs directly on the ligand-stripped structure.

Outputs are stored under `output/48_detect_pocket_multimers/`. Aggregate `pocket_detection_data.csv` (one row per detected pocket):

| **Field** | **Description** |
|-----------|------------------|
| `PDB code` | RCSB identifier |
| `File name` | Stripped structure file name |
| `Full path` | Path to the structure used for pocket detection |
| `Chains` | Protein chain IDs kept after stripping (multimer composition) |
| `Pocket number` | P2Rank rank |
| `Pocket score` | P2Rank score |
| `Pocket probability` | P2Rank probability |
| `Pocket centroid coordinate (x y z)` | Centroid coordinates |
| `Pocket residues (chain_resn)` | Residues forming the pocket, with chain and residue number |
| `Experimental method` | From the RCSB entry API |
| `Resolution (A)` | From the RCSB entry API |
| `Biological assembly info` | Oligomeric state and author/software provenance, from the RCSB assembly API |
| `RCSB access date` | Date the structure/annotations were fetched |

Verified on three real PDB codes spanning a homodimer (1HXW), a heterotetramer (7K98) and a monomer (5W25): ligand/water stripping, multi-chain retention, cross-chain interface pocket detection and RCSB metadata all behaved as expected, and re-running the script is a no-op on already-completed stages.

---

## Pending review (no README text yet)

The following scripts exist on disk but haven't been individually reviewed against documentation yet: `13_prepare_PocketVec.py` (placeholder above), `47b_reference_pocket_visualization.py`, `49_unidock_proteinprep_multimers.py` through `76_boltz2_docking_supervisor.sh` (52 `CAT_promiscuous`/`CAT_selective`, 54–61 `NONCAT`/`REAL10M`/`REAL10B` selection and docking, 62–70 hit aggregation/merging/characterization/plotting/filtering, 70b compound GIFs, 71–76 Boltz-2 validation). These will get entries as we go through them one at a time.

## Stale entries (pending correctness pass)

The following five entries are carried over **verbatim** from the previous root README. They describe an **older numbering** of the pipeline — the scripts that exist under these numbers today are different (`49_unidock_proteinprep_multimers.py`, `50_unidock_docking_multimers.py`, `51_selected_pockets_visualization.py`, `54_NONCAT_promiscuous.py`; there is no second `48` file — script 48 is `48_detect_pocket_multimers.py`, documented above). Kept here for now rather than deleted, until each of the actual current scripts 48–76 is reviewed and given an accurate entry.

### `48_docking_hits_raw.py` (as originally described — does not match the current script 48)
Using each gene's reference pocket, quantifies raw overlap between genes' top-100 and top-1,000 docking hits (visualized as UpSet plots) and compares observed vs. expected-by-chance multi-target binder counts. Compounds hitting at least 2 targets within the top 1,000 are collected into a CSV with per-gene scores/ranks and physicochemical properties (MW, cLogP, TPSA, HBD, HBA, rotatable bonds, aromatic rings, QED, PAINS flag), benchmarked against a 25,000-compound random background. Outputs were stored in `output/48_docking_hits_raw/`.

### `49_docking_hits_selective.py` (as originally described — does not match the current script 49)
Builds on the same reference pockets to select a final set of up to **500** compounds balancing potency against selectivity, using five complementary ranking metrics (m1–m5, each contributing up to 50 compounds unless otherwise capped):

* **m1** — max potency, low selectivity: top target rank, excluding compounds whose non-target 50th-percentile rank falls below 20,000.
* **m2** — max potency, high selectivity: same as m1, but requiring a non-target 50th-percentile rank of at least 50,000.
* **m3** — high potency, selectivity gap: ranked by the gap between the non-target 10th-percentile rank and the top target rank (target rank ≤ 20,000, non-target 50th-percentile rank ≥ 20,000).
* **m4** — high potency, max selectivity: ranked by the non-target 1st-percentile rank (target rank ≤ 20,000).
* **m5** — diversity rescue: compounds binding at least 2 targets well (2nd-best target rank ≤ 20,000, non-target 50th-percentile rank ≥ 50,000) and not already selected by m1–m4, deduplicated by Murcko scaffold, topping the total selection up to 500.

Outputs (per-compound metric table, score/profiling plots) were stored in `output/49_docking_hits_selective/`.

### `50_ersilia_eos42ez.sh` (as originally described — does not match the current script 50)
SMILES from the (formerly-)48 and 49 scripts are merged and deduplicated, then run through the Ersilia Model Hub model `eos42ez` (HepG2, HSkMC and IMR90 cytotoxicity prediction) via the Ersilia CLI. Predictions were stored in `output/50_ersilia_eos42ez/eos42ez.csv`.

### `51_filter_hits.py` (as originally described — does not match the current script 51)
Applies a sequential filter to the combined raw and selective hit sets: QED > 0.5, then exclusion of PAINS structural alerts, then all three cytotoxicity scores < 0.3, printing the number of compounds surviving each stage. The final shortlist was stored in `output/51_filtered_hits/filtered_hits.csv`:

| **Field** | **Description** |
|-----------|------------------|
| `key` | Compound identifier |
| `smiles` | Compound SMILES string |
| `cytotoxicity_hepg2` | Predicted HepG2 cytotoxicity score (Ersilia model `eos42ez`) |
| `cytotoxicity_hskmc` | Predicted HSkMC cytotoxicity score (Ersilia model `eos42ez`) |
| `cytotoxicity_imr90` | Predicted IMR90 cytotoxicity score (Ersilia model `eos42ez`) |
| `QED` | Quantitative Estimate of Drug-likeness (RDKit) |
| `is_pains` | Whether the compound matches a PAINS structural alert (RDKit) |

### `54_merge_scores_select_hits.py` (as originally described — does not match the current script 54)
Consolidated round-1 (113k, `output/unidock_REAL_docking`) and round-2 (~99k, `output/unidock_REAL_docking_2`) docking scores for the `alphafold3_P9WFU3_model_2_pocket_2` reference pocket into a single ranked table (compound_id, smiles, source, docking_score), sorted by docking score ascending (112,953 `10M` + 99,100 `10B` + 5 `Both` = 212,058 rows; compounds present in both campaigns kept their best/lowest score and were labeled `Both`). It then greedily selected the top 10,000 compounds while capping chemical redundancy: compound IDs encode synthons (building blocks) as `{m|s}_{reaction_id}____{synthon_1}____{synthon_2}[____{synthon_3}]`, and no synthon could appear in more than `MAX_SYN = 3` selected compounds (same diversity-capping idea as script 33 and notebook 43). Two columns were added: `observed_synthons` (semicolon-separated per-synthon counts across ALL previously considered rows, selected or not — a diagnostic column, not itself capped) and `select` (1/0, driven by a separate selected-only counter enforcing the MAX_SYN cap). The output was truncated to the rows actually considered (not the full 212,058) — reaching 10,000 selections took 10,080 considered rows at this cap. Output: `output/54_merge_scores_select_hits/alphafold3_P9WFU3_model_2_pocket_2_merged_docking_scores.csv`.
