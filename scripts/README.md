# Scripts

One entry per script, in pipeline execution order. Explicit thresholds and cutoffs are stated with their rationale where the original design decision recorded one. See the root [README.md](../README.md) for the high-level pipeline overview and [CLAUDE.md](../CLAUDE.md) for repository-wide conventions.

## Structure preparation

### `00_generate_fasta_files.py`
Generates a FASTA file for each of the 21 tRNA synthetases (from `data/mtb_trna_synthetases_bosch_2021_fig5_annotated.csv`), used to query external structure/sequence resources.

### `01_download_alphafill.py`
Programmatically downloads AF2 structures with co-crystallized ligands from [AlphaFill](https://alphafill.eu/) into `/data/structures/alphafill_database/`.

### `02_organize_structures.py`
Organizes structure files from all sources (AlphaFold2, AlphaFold3, Chai-1, AlphaFill, SwissModel; PDBe is omitted from this automated processing) into `output/structures/<uniprot_ac>/`, converting each from `.cif` to `.pdb` via PyMOL and then deleting the source `.cif` — only `.pdb` files remain on disk. For multi-chain files, keeps the single chain with the highest sequence coverage against the UniProt reference sequence, requiring that maximum coverage to be ≥95%; structures whose selected-chain sequence is <80% of the reference length are dropped outright, and those in the 80–95% range are kept only if the sequence is exactly continuous with (a contiguous substring of) the UniProt reference, otherwise removed. Non-protein atoms (waters, ligands, other heteroatoms) are then stripped from the retained chain. Produces the per-structure lookup table `output/trna_synthetases_data.csv` (see aggregate table below).

### `03_align_structures.py`
Aligns each protein's structures (Cα superimposition over the residue range shared with the reference, via Biopython) to a single reference — the first-listed AlphaFold3 model (`model_0`) for that protein — skipping AlphaFill structures (identical to AlphaFold2, so redundant to align). Per-structure RMSD against the reference is saved to `output/alignment_rmsd_data.csv`; structures with RMSD > 10 Å are removed from `output/structures`, `output/aligned_structures`, and `output/trna_synthetases_data.csv`.

### `04_relax_structures.py`
Prepares structures for docking: protonation via PDB2PQR (AMBER force field, pH 7.0, run in the `adda4tb` conda env) followed by relaxation with PyRosetta's `FastRelax` (`ref2015` score function). Each structure is relaxed 3 times independently and the lowest-scoring (best) pose is kept. Reads the RMSD-filtered file list from `output/alignment_rmsd_data.csv`, inputs from `output/aligned_structures/`, outputs to `output/relaxed_structures/<uniprot_ac>/`. Computationally intensive; re-running skips structures whose relaxed output already exists.

### `05_align_relaxed_structures.py`
Re-aligns each relaxed structure to its own unrelaxed counterpart (`output/aligned_structures/`) as the alignment reference (same Cα-superimposition approach as script 03). RMSDs are saved to `output/alignment_relaxed_rmsd_data.csv`; no structures are removed at this stage — the removal code path is present but commented out, since post-relaxation RMSD against the unrelaxed reference was manually confirmed to stay under the 10 Å threshold used in script 03 for every structure.

**Aggregate output (scripts 02–05):** `output/trna_synthetases_data.csv` has one row per processed structure:

| **Field** | **Description** |
|-----------|------------------|
| `file_name` | Name of the processed PDB structure file |
| `chain_id` | Chain ID retained in the processed structure |
| `uniprot_ac` | Uniprot AC identifier |
| `n_residues` | Number of residues in the UniProt full sequence |
| `start_resid` | First residue number (first residue is 1) of the sequence available in the PDB file, with respect to the Uniprot full sequence |
| `end_resid` | Last residue number (first residue is 1) of the sequence available in the PDB file, with respect to the Uniprot full sequence |
| `coverage` | Percentage sequence coverage |
| `structure_sequence_length` | Length of the sequence found in the PDB file |
| `full_sequence_length` | Length of the sequence found in UniProt |
| `sequence_structure` | Sequence found in the PDB file |
| `full_sequence` | Sequence found in UniProt |

## Sequence & pocket analysis

### `06_sequence_annotation.py`
Aggregates protein family/domain annotations downloaded from [InterPro](https://www.ebi.ac.uk/interpro/) (per-protein files `data/sequences/interpro/entry-matching-<uniprot_ac>.tsv`) across all 21 tRNA synthetases into one combined table, `output/sequences/interpro_data.csv` (columns: `uniprot_ac`, `interpro_ac`, `name`, `source_database`, `type`, `integrated_into`, `integrated_signatures`, `go_terms`, `protein_length`, `matches`). Input for script 10's summary.

### `07_fetch_from_chembl.py`
Fetches ChEMBL bioactivity data (IC50-type activities only, via the `chembl_webresource_client` package) for the tRNA synthetases with a mapped ChEMBL target ID in `data/chembl_uniprot_mapping.txt`, saving one raw JSON dump per UniProt AC to `data/ligands/chembl/<uniprot_ac>.json`. Data was found for only 3 of the 21 tRNA synthetases.

### `08_detect_pockets.py`
Detects pockets in AF2, AF3, Chai-1 and SwissModel predicted models (`output/aligned_relaxed_structures/`, per `output/alignment_relaxed_rmsd_data.csv`; AlphaFill is excluded upstream since it's identical to AlphaFold2). Runs P2Rank with the `-c alphafold` config. Following the [authors' recommendations](https://github.com/rdk/p2rank/issues/76) (probability ≥ 0.2, at least top-K by rank) and keeping the total number of pockets tractable, we limited K to 2. Pockets with at least one residue below pLDDT < 65 (AF2, AF3, Chai-1) or QSQE < 0.65 (SwissModel) are then filtered out, discarding about 25–30% of pockets. These cut-offs are arbitrary; usual recommendations are 70 & 0.7 — this pipeline is slightly less restrictive. See also script 48, which applies the same P2Rank approach to experimental multimeric structures.

Summary output `output/pocket_detection_data.csv` (one row per detected pocket and structure):

| **Field** | **Description** |
|-----------|------------------|
| `Uniprot AC` | Uniprot AC identifier |
| `File name` | PDB file in which pockets have been detected |
| `Prediction type` | The method used for protein structure prediction |
| `Full path` | The full directory path where the PDB file is stored |
| `Pocket number` | The identified pocket number within the structure |
| `Pocket score` | The score assigned to the detected pocket |
| `Pocket probability` | The probability value indicating confidence in pocket detection |
| `Pocket centroid coordinate (x y z)` | The (x, y, z) coordinates of the pocket's centroid |
| `Pocket residues (chain_resn)` | List of residues forming the pocket, with chain and residue number |
| `B-factors` | Confidence measures: pLDDT (AF2, AF3, Chai) or QSQE (SM) |

Note: the `Full path` column in the currently-committed `output/pocket_detection_data.csv` still records the pre-rename `processed/aligned_relaxed_structures/...` prefix (from before this repo's `processed/` → `output/` rename, now fixed in the scripts themselves) — treat existing values in that column as informational, not necessarily the literal current path.

### `09_prepare_pymol_visualizations.py`
Generates one gzipped PyMOL session file (`output/pymol_sessions/<uniprot_ac>.pse.gz`) per protein, using the post-relax/aligned structures (`output/aligned_relaxed_structures/`) with the AlphaFold2 model as reference, to visualize detected pockets and their residues:

| **Element** | **Description** | **Displayed?** |
|-------------|------------------|--------------|
| Reference structure (AF2) | Wheat color, surface + cartoon representation | Yes |
| Pockets detected in reference structure (AF2) | Sky-blue spheres with arbitrary size (pocket detection provides a single 3D coordinate) | Yes |
| Residues defining each pocket in AF2 | Orange color, surface + cartoon representation | Yes |
| Aligned structures (all but AF2) | Gray color, cartoon representation | No (loaded but disabled) |
| Pockets detected in aligned structures | Gray-colored points | No (loaded but disabled) |
| InterPro annotations | One object per matched entry: white surface for the whole reference structure, with residues in the entry's matched range(s) colored red | No (loaded but disabled) |

### `10_organize_sequence_info.py`
Aggregates the per-protein InterPro data from script 06 into two files. `output/sequences/interpro_summary.tsv` — one row per InterPro accession, across all 21 tRNA synthetases:

| **Field** | **Description** |
|-----------|------------------|
| `Accession` | InterPro Accession identifier |
| `Name` | Name of the InterPro entry |
| `Type` | Type of annotation (e.g., domain, conserved site, homologous superfamily) |
| `Number of proteins` | Number of proteins carrying this InterPro entry |
| `Average Coverage` | Average, across those proteins, of (residues covered by the entry's matches) / (protein length) |
| `Proteins` | Comma-separated UniProt ACs carrying this InterPro entry |
| `Protein Coverages` | Per-protein coverage values, comma-separated in the same order as `Proteins` |

Plus one binary (1/0) column per protein (named by UniProt AC) indicating presence/absence of that InterPro entry. A manually curated companion file, `output/sequences/interpro_summary_curated.tsv`, restricts to entries with `Average Coverage` < 0.60 (excluding broad, low-specificity domains) and inserts a `Curated annotation` column classifying each entry's `Name` by keyword match into one of: Catalytic Domain (ATP Binding Site), Anticodon Binding Domain, Editing Domain, tRNA Binding Domain, RNA Binding Domain, Other too broad/unspecified functional entities, or Other.

The curation half of this script (from `### CURATE INTERPRO ANNOTATIONS ###` onward) previously hardcoded an absolute path (`/home/acomajuncosa/Documents/mtb-targeted-protein-degradation`) instead of reusing the script's own `root` variable; fixed to use `root` like the rest of the script.

### `11_organize_pocket_info.py`
Joins script 08's pocket data (`output/pocket_detection_data.csv`) with script 10's curated InterPro table (`output/sequences/interpro_summary_curated.tsv`) by UniProt AC, and computes residue overlap between each detected pocket and each of that protein's (already <60%-coverage-curated) InterPro domain matches. Output `output/pocket_detection_data_interpro.tsv` — one row per pocket × InterPro-entry pair with at least one overlapping residue, carrying all of script 08's columns (see above; note the field is `Uniprot AC` with a space, not an underscore) plus:

| **Field** | **Description** |
|-----------|------------------|
| `Interpro ID` | Identifier of the matched InterPro domain |
| `Interpro name` | Name of the matched InterPro domain |
| `Interpro curated annotation` | The functional category assigned in script 10 (Catalytic Domain (ATP Binding Site), Anticodon Binding Domain, etc.) |
| `Interpro Matches` | Residue ranges corresponding to the InterPro domain |
| `Residues in pocket` | Number of residues forming the pocket |
| `Residues in interpro` | Number of residues forming the InterPro domain |
| `Residues overlap` | Number of residues shared between the pocket and the InterPro domain |
| `Coverage pocket` | Fraction of pocket residues that overlap with the InterPro domain |
| `Coverage domain` | Fraction of InterPro domain residues overlapping with the pocket |
| `Overall coverage domain` | Fraction of the protein's full UniProt sequence length covered by the InterPro domain match |

Pocket–InterPro pairs with no overlapping residues are omitted from the file.

### `12_align_PDB_structures.py`
For every UniProt AC with PDBe structures under `data/structures/pdbe_database/`, loads that protein's PyMOL session from script 09 (`output/pymol_sessions/<uniprot_ac>.pse.gz`), aligns each PDBe structure onto the session's AlphaFold2 reference object (`cmd.align`), and for every non-water, non-standard-amino-acid heteroatom group (co-crystallized ligand) in the aligned structure, computes its minimum 3D distance to each of that protein's P2Rank pocket centroids. Saves a long-format report, `output/pdbe_annotation_report.csv` (tab-separated despite the `.csv` extension) with columns `Uniprot ID`, `PDB Structure`, `Pocket`, `HET RES`, `Min Distance`, plus a modified copy of each session (surfaces hidden, PDBe structures added) to `output/pymol_sessions_pdbe/<uniprot_ac>.pse.gz`.

### `13_prepare_PocketVec.py`
Prepares [PocketVec](https://github.com/sbnb-irb/pocketvec) inputs for every detected pocket: converts the full aligned/relaxed structure (`output/aligned_relaxed_structures/<uniprot_ac>/<file>.pdb`) to `.mol2`, and the single-point pocket-centroid PDB written by script 08 (`output/detected_pockets/<uniprot_ac>/<file>/pockets/pocket_<n>.pdb`) to `.sd`, both via OpenBabel, into one directory per pocket (`<file>_pocket_<n>/`) under `output/pocketvec_PRE/`. Unlike other scripts, its paths are relative to the current working directory rather than the script's own location, so it must be run with `scripts/` as the working directory. Writes into `output/pocketvec_RUN/pocketvec_PRE/` (script 22 reads from the same nested location) — the extra `pocketvec_RUN/` level, alongside the plain `processed/`→`output/` rename, matches where PocketVec's pre-computed inputs actually live on disk.

### `14_calculate_SeqId.py`
Computes global-alignment (Needleman–Wunsch, BLOSUM62, gap open −10/extend −1) pairwise sequence identity across the 21 tRNA synthetases **plus** 25 random Mtb proteome background proteins (`data/mtb_proteome_from_uniprot.tsv`, `random_state=42`, excluding the 21 tRNA synthetases) — a 46×46 matrix, not only the 21 tRNA synthetases against each other. Outputs `output/sequences/NW_SeqAlign/SeqId_matrix.tsv` (% identity) and `.../Prop_matrix.tsv` (fraction of aligned, non-gap positions relative to the shorter sequence), plus a heatmap plotted directly with matplotlib rather than stylia (this script predates that convention), saved to `output/plots/SeqId_matrix.png`.

### `15_calculate_StSim.py`
Evaluates structural similarity between the 21 tRNA synthetases at the structure level: for every ordered pair of proteins (computed in both directions, not just the upper triangle), aligns every combination of their `output/aligned_relaxed_structures/` files via PyMOL's sequence-independent `cmd.super` (unlike script 12's sequence-based `cmd.align`) and records the RMSD. Output: one CSV per ordered protein pair, `output/structural_comparisons/<uniprot_ac_1>_<uniprot_ac_2>_rmsd.csv` (columns `file_name_1`, `file_name_2`, `rmsd`).

### `15b_calculate_StSim_local.py`
Complements script 15's whole-chain comparison with a local, topology-independent one: for each of the 210 unordered pairs of the 21 tRNA synthetases, aligns only their AF2 model (`alphafold2_<uniprot_ac>_model_0.pdb`, one structure per protein, not every predicted structure like script 15) via PyMOL's `cmd.cealign` (Combinatorial Extension) and records the RMSD and aligned length. Unlike `cmd.super`/`cmd.align`, `cealign` doesn't attempt a global 1:1 correspondence across the whole chain, so it can find a well-matching substructure even when the full-chain comparison is poor — e.g. hisS/alaS have the highest pairwise sequence identity in the set (38%) but a poor whole-chain RMSD (~17 Å, script 15), while a residue-restricted comparison of just their shared Class II fold core superposes at ~1–1.6 Å; `cealign` independently finds the same region without needing any InterPro/Class I/II annotation, which is also why it needs no special-casing for gatA/gatB (the 2 of the 21 that aren't classical aaRS and have no Class I/II fold at all — it just reports a correspondingly weaker match for them). Output: a single table, `output/structural_comparisons_local/StSim_local_pairwise.csv` (columns `uniprot_ac_1`, `uniprot_ac_2`, `rmsd`, `alignment_length`), rather than per-pair files, since there's only one value per pair.

### `15c_calculate_StSim_global_local_coverage.py`
The exhaustive, coverage-aware version combining scripts 15 and 15b's methods: for every unordered pair of the 21 tRNA synthetases, **including self-pairs** (a protein's own structures against each other — AF2 vs. AF3 vs. Chai-1 vs. SwissModel, etc. — filling a gap `scripts/plots/figure_1_calculations.py` flags in a comment but never computes), and every combination of their `output/aligned_relaxed_structures/` files (231 pairs, ~33,700 structure-combinations), computes **both** `cmd.super` (global) and `cmd.cealign` (local) on the same loaded pair. Records not just RMSD but the alignment coverage for each method (`aligned_residues / min(res_1, res_2)`) — motivated by a real failure mode found during development: whole-chain `cmd.super` can report a deceptively low RMSD (e.g. 0.63 Å between cysS1 and lysS, two different Class I/II folds) driven by only a tiny fraction of residues (24/449, ~5%) actually superimposing. Output: one CSV per pair, `output/structural_comparisons_full/<uniprot_ac_1>_<uniprot_ac_2>_global_local.csv` (columns `file_name_1`, `file_name_2`, `res_1`, `res_2`, `global_rmsd`, `global_aligned_res`, `global_coverage`, `local_rmsd`, `local_aligned_res`, `local_coverage`), plus a shared `errors.log` for any structure combination that fails to align. Runtime is on the order of hours (~5, benchmarked), so it's meant to run as a background/overnight job; since a pair's CSV is only written once that pair's combinations are all done, restarting the script simply skips any pair whose output file already exists.

### `16_calculate_PocketVec.py`
Computes rank-transformed [PocketVec](https://github.com/sbnb-irb/pocketvec) descriptors ([Comajuncosa-Creus et al., Nat Commun 2024](https://www.nature.com/articles/s41467-024-52146-3)) from each pocket's raw rDock docking scores against a fixed 128-compound reference set (`output/pocketvec_RUN/pocketvec_POST/<pocket>/rDock_results_LLM/*scores*`; the docking itself ran upstream on an HPC cluster with CPU parallelization, not in this script). For each of the 276 pockets, builds a raw-score vector (ordered by `output/pocketvec_RUN/TOP128_rDock_LLM/ALL/all.pkl`) and rank-transforms it (`scipy.stats.rankdata`, `method='max'`); raw scores in (0, 50), (50, 100), and above 100 are additionally pushed to increasingly high dummy ranks (`len+1`, `len+2`, `len+3` respectively) — the script doesn't state why, but this is consistent with demoting docking runs that returned an invalid/penalty score rather than a genuine (negative) binding score. Only the rank-transformed fingerprints are saved, to `output/pocketvec_RUN/fps_rank.pkl` — the raw scores are not persisted.

The outlier-pocket filtering (12/276 descriptors excluded for having >80 outlier compounds) and the 0.17 PocketVec-distance similarity threshold used by script 17's prioritization are computed separately, in `notebooks/16_PocketVec_analyses.ipynb` — not by this script.

### `77_pocket_annotation/`
**Note on numbering**: despite its late number, this only depends on scripts 00-08 (structure prep + P2Rank pocket detection) and logically belongs here, alongside scripts 06-12. It's numbered 77 because it was added long after the rest of the pipeline was built, and no integer was free in the 06-17 range. It **coexists** with scripts 06/09/10/11/12 rather than replacing them — those and their outputs (`interpro_data.csv`, `interpro_summary_curated.tsv`, `pocket_detection_data_interpro.tsv`, `pymol_sessions/<uniprot_ac>.pse.gz`) are untouched, and script 47 / the notebooks that read the old file by its old column names keep working unchanged. This pipeline's own outputs live in separate, numbered folders (`output/77_pocket_annotation/`, `output/77_pymol_sessions/`) precisely to avoid any collision with scripts 06/09's outputs.

A corrected redo of the domain-annotation and pocket-domain-overlap analysis, fixing several real bugs found in scripts 06/09/10/11/12 (script 10's cross-protein coverage filter drops the correct catalytic annotations for compact proteins like pheS; its keyword curation has substring false-positives; gatA/gatB's different catalytic mechanism got zero domain coverage under the old GO-term map; script 11 never resolves for pheS at all; script 12's whole-chain rigid alignment breaks badly, 10-14 Å RMSD, on proteins with independently-flexible domains). Also properly parses AlphaFill's ligand-provenance metadata (`data/structures/alphafill_database/*.json`) for the first time anywhere in the pipeline — previously only raw coordinate proximity was used.

11 numbered steps (`01_fetch_interpro.py` through `11_build_pymol_sessions.py`, run via `run_all.sh`) produce:
- `output/77_pocket_annotation/pocket_detection_interpro_updated.csv` — one row per detected pocket (276 total), with curated domain label(s) (multi-valued), per-category entry-support counts, direct-PDB and AlphaFill ligand evidence (identity/source/distance/strength, plus AlphaFill's `local_rmsd`/homolog identity), whether *this pocket* has direct-PDB ligand evidence (`has_direct_ligand_evidence`), 3D spatial-cluster reproducibility (`spatial_cluster_id`/`n_models_in_cluster`, via the same fixed-6.14-Å greedy dedup as `scripts/plots/figure_1_calculations.py`), and a `catalytic_confidence` score (0-4: 0 whenever the pocket lacks the curated Catalytic Domain label — no ligand evidence can raise it without that label; otherwise 1 for the label alone, +1 for any weak ligand evidence, ceiling 4 for a strong ligand — substrate/ATP-cofactor/reaction-intermediate-mimic/curated-inhibitor — found directly in an experimental structure of that protein, ceiling 3 if the only strong evidence is AlphaFill-transplanted). A low `catalytic_confidence` does not mean a pocket is invalid — only that no catalytic-specific evidence was found for it.
- `output/77_pymol_sessions/<gene>_pocket_annotation.pse` — one comprehensive session per protein: category-colored domain objects, cluster-colored pocket-centroid spheres, aligned experimental PDB structures (ligands colored by weak/strong evidence), and an aligned AlphaFill object.

See `scripts/77_pocket_annotation/README.md` for the full per-step breakdown and the ligand weak/strong classification rules.

## Protein prioritization

### `17_protein_prioritization.py`
This script stubs out the prioritization pipeline; the logic itself is implemented in notebooks `17a_Protein_prioritization.ipynb`, `17b_Protein_prioritization_pairs.ipynb`/`17b_Protein_prioritization_triplets.ipynb`, and `17c_Protein_prioritization_individual.ipynb`. In brief: protein pairs (210) and triplets (1,330) were exhaustively enumerated and classified by PocketVec distance, global similarity (structural and sequential) and the number of pockets mapped to the catalytic domain. Global similarity thresholds: 35% sequence identity and 10 Å RMSD for sequential and structural similarity, respectively. Comparisons were extended to the pocket level (32,561 pairs, 2,499,258 triplets), collecting lenient (PocketVec distance < 0.17) and stringent (PocketVec distance < 0.14) sets of pairs (1,481 / 76) and triplets (3,880 / 807). An intra-set normalization then derives a "Prioritization score" per protein based on how often it appears across the collected sets — indicating similarity to other tRNA synthetases and potential polypharmacology.

Final output `output/protein_prioritization/final_results.tsv`:

| **Field** | **Description** |
|-----------|------------------|
| `Uniprot AC` | Uniprot AC identifier |
| `Gene name` | Standard gene name associated with the Uniprot AC identifier |
| `Vulnerability score` | Vulnerability score derived from [Bosch et al, 2021](https://pubmed.ncbi.nlm.nih.gov/34297925/) |
| `Score` | Prioritization score |
| `Tier` | Protein-associated tier (1–5) |
| `sim_Tier1` … `sim_Tier5` | Five separate columns (not one) — for the protein in this row, the number of proteins in tier *N* (1 through 5) having a PocketVec distance < 0.14 to it |

## 100k compound screening (Enamine Diversity Library HLL-100)

After structural characterization of the 21 essential tRNA synthetases, potential small-molecule binders are searched for in an active-learning fashion: first docking a chemically diverse 100k-compound set against each pocket structure, then training a surrogate model with [LazyQSAR](https://github.com/ersilia-os/lazy-qsar).

Computational resources used throughout this and later screening stages: `herbert` (default host unless noted), `nebula` (NVIDIA GeForce RTX 4090), `SBNB-IRB cluster` (HPC, `/aloy/home`).

### `18_enamine_mfps.py`
Calculates radius-3/2048-bit Morgan Count Fingerprints (RDKit `GetMorganGenerator(radius=3, fpSize=2048)`, i.e. ECFP6-equivalent) for the ENAMINE Diversity Library HLL-100 (`data/Enamine/Enamine_Hit_Locator_Library_100K_Set_plated_100160cmpds_20250629.smiles`, nominally 100,160 compounds; SMILES parsing/fingerprinting failures are silently skipped, yielding the 100,157 compounds actually saved to `X.npz`). Also saves InChIKey/SMILES lookup dictionaries (`IK_TO_SMI.pkl`, `ID_TO_SMI.pkl`, `ID_TO_IK.pkl`). Stored in `output/enamine_characterization/`.

### `19_enamine_conformations.py`
Generates 3D conformations (ETKDGv3, seed 42, + UFF energy minimization; 3 more compounds fail embedding, yielding 100,154 of the 100,157 fingerprinted compounds) for the 100k library, parallelized over 12 workers. Saves both a single gzipped multi-molecule SDF (`output/enamine_characterization/conformations.sdf.gz`) and a `.tar.gz` archive of one `.sdf` file per compound (`output/enamine_characterization/conformations.tar.gz`) — the uncompressed per-compound files are deleted after archiving. 90% of these compounds fall in the [270, 400] Da range.

### `20_unidock_ligandprep.py`
Extracts script 19's `conformations.tar.gz`, builds a ligand-index text file, and prepares the 100,154 3D-conformer compounds for docking via `unidocktools`' `ligandprep` (batches of 100; run on `nebula`, `unidock_env` conda env). The ligand-index file is then rewritten to list only the compounds that survived `ligandprep` (failures discarded). Output: `output/unidock_docking/conformations_prepared.tar.gz` plus `output/unidock_docking/input_ligands.txt` — the extracted/prepared directories themselves are deleted after tarring, so script 22 re-extracts the tarball before docking.

### `21_unidock_proteinprep.py`
Prepares the 178 protein structures (276 pocket structures, see `output/pocket_detection_data.csv`) via `unidocktools`' `proteinprep` (run on `nebula`, `unidock_tools` conda env with `ambertools` and `reduce`). 8 structures required manual intervention to fix PDB2PQR/PyRosetta protonation issues. Output: `output/unidock_docking/structures_prepared.tar.gz`.

### `22_unidock_docking.py`
Docks the 100k library against all 276 pocket structures with [Uni-Dock](https://github.com/dptech-corp/Uni-Dock) (`unidock_env` conda env), search mode _fast_, _vina_ scoring, box size 22.5 Å per side, seed 42, only the single top-scoring pose kept per compound (`num_modes=1`), on `nebula` (~1 week of compute). Pockets are processed in a fixed random order (`random_state=42` shuffle of `pocket_detection_data.csv`). Per pocket, also copies in that pocket's PocketVec `.sd` file (script 13) for reference alongside the docking results. Output: `output/unidock_docking/docking_results/<pocket_label>/` (`docking.tar.gz` of docked poses, `logs.tar.gz`, `report.csv` of per-compound Uni-Dock energies parsed from each pose's `<Uni-Dock RESULT>` block), 276 pocket directories total.

### `23_unidock_checks.py`
Console-only sanity check (writes nothing to disk): for each pocket's `logs.tar.gz` from script 22, prints the pocket name and the first log line containing "error" or "warning" (excluding the known-benign `in add_to_output_container` warning) — only the first matching line per pocket is shown, so later problems in the same log won't surface here.

### `24_enamine_chemeleon.py`
Calculates [CheMeleon](https://github.com/JacksonBurns/chemeleon) embeddings for the same 100,157 compounds, using Ersilia's NVIDIA GeForce RTX 4090 (hardcoded `device="cuda"` — needs manual editing to run CPU-only). Downloads the pretrained MPNN checkpoint from Zenodo ([record 15460715](https://zenodo.org/records/15460715/files/chemeleon_mp.pt)) to `data/CheMeleon/chemeleon_mp.pt` on first run, then extracts the mean-aggregated message-passing fingerprint (not a downstream property prediction — the regression head that ships with the checkpoint is instantiated but unused) in batches of 1,000, stored as `float16`. Stored alongside the fingerprints from script 18 in `output/enamine_characterization/` (`X_CheMeleon.npz`, `IDs_CheMeleon.txt`).

### `25_enamine_binarization.py`
Binarizes docking scores at 3 percentiles (0.1%, 0.5%, 1%) per pocket structure — `bin_01`/`bin_05`/`bin_1` = 1 if the compound's score is more negative (better) than that percentile of the pocket's own score distribution — producing one report per pocket, `output/unidock_docking/binarized_reports/report_bin_<pocket>.csv` (compound ID, docking score, 3 bin flags). Also writes `output/unidock_docking/pickle.pkl`, a flat list of 828 `(pocket, bin_name)` pairs (276 pockets × 3 bins) used to index script 26's training jobs.

### `26_train_models/`
`main.py` trains one [LazyQSAR](https://github.com/ersilia-os/lazy-qsar) random-forest classifier (`pca=False`, `num_trials=20`, `max_samples=10000`) per `(pocket, bin)` job — indexed by an integer `alpha` (0–827) into script 25's `pickle.pkl` — on CheMeleon embeddings (script 24) vs. that bin's flag. Reports 3-fold stratified CV AUROC (seed 42, `AUROCs.csv`) from a first pass, then refits on all data and saves the final model as `output/unidock_docking/models/<pocket>/<bin>/LQ_RF.joblib`. `run.sh` wraps `main.py` for the `SBNB-IRB cluster` (hardcoded to a specific user's conda env path, `/aloy/home/acomajuncosa/anaconda3/envs/lazyqsar`); `submission.sh` is the SLURM array job wrapper submitting all 3 bins × 276 pockets = 828 jobs (`#SBATCH --array=0-827%3`; a since-corrected, previously-checked-in version of this file only covered a partial rerun of 3 jobs), logging to `output/unidock_docking/training_outputs/` on the same SBNB-cluster checkout (`/aloy/home/acomajuncosa/Ersilia/mtb`) — its `processed/`→`output/` rename happened independently of this local checkout, since `data/`/`output/` are eosvc-synced, not git-tracked. After analysis, only the 0.1%-binarized (`bin_01`) models were used downstream.

## Enamine REAL 9.56M: characterization, screening & prioritization

### `27_enamine_REAL_chemeleon.py`
Loads the 9.56M-compound Enamine REAL library (`data/Enamine/2024.07_Enamine_REAL_DB_9.6M.cxsmiles`), sorts by InChIKey and saves `output/enamine_REAL_characterization/enamine_REAL.tsv` (InChIKey, id, smiles) — the printed "compounds in the original dataset" vs. "IDs in the final compound set" counts are a duplicate-ID diagnostic only; the saved table itself is not deduplicated. Calculates 2048-dim CheMeleon embeddings (same mechanism as script 24, including the hardcoded `device="cuda"`) using Ersilia's NVIDIA GeForce RTX 4090 (`nebula`), in chunks of 100,000 compounds each written to its own gzip-compressed HDF5 file (`X` and `ids` datasets). Stored at `output/enamine_REAL_characterization/embeddings/` (96 files, ~100k compounds each, ~25GB total).

### `28_transfer_embeddings.sh`
Run from `nebula` itself (not from this checkout): transfers script 27's embeddings from its local `output/enamine_REAL_characterization/embeddings/` to the same path on the `SBNB-IRB cluster` (`/aloy/home/.../Ersilia/mtb/output/...`) via `rsync` over ssh through `dante`.

### `29_inference_enamine_REAL/`
`main.py` applies one `bin_01` surrogate model (script 26) per pocket to all 9.56M Enamine REAL compounds: loads the pocket's `LQ_RF.joblib`, scores all 96 CheMeleon embedding chunks (script 27), and saves `output/unidock_docking/inferences/<pocket>_bin_01.csv.gz` (`ids`, `probs`). Indexed by `alpha` (0–275, one per pocket) into `output/unidock_docking/pickle_bins_01.pkl`, built in `notebooks/26_docking_models.ipynb` by filtering script 25's `pickle.pkl` down to its `bin_01` entries — the same notebook whose AUROC comparison across bins justified using only `bin_01` downstream (per script 26). `run.sh`/`submission.sh` mirror script 26's SBNB-cluster SLURM setup (276-job array, one per pocket).

### `30_enamine_REAL_conformations.py`
Generates 3D conformations (ETKDGv3, seed 42, + UFF energy minimization) for all Enamine REAL 9.56M compounds, in chunks of 100,000, parallelized over 32 workers. Each chunk's conformers are pickled together (gzip), `output/enamine_REAL_characterization/conformations/mols_<chunk>.pkl.gz`; a separate `id_to_split.tsv.gz` records which chunk each compound ID landed in.

### `31_enamine_REAL_checks_conformations.py`
Identifies compounds for which conformation generation (script 30) failed (`mols_<chunk>.pkl.gz` entries mapped to `None`); IDs stored in `output/enamine_REAL_characterization/failed_REAL.csv`. Also draws a seeded (42) random sample of 25k successfully-embedded compounds for subsequent analyses, stored in `output/enamine_REAL_characterization/random_sample_REAL.csv`. Reads script 30's conformations directly from its `nebula`/`Documents_GPU` output location rather than a local copy.

### `32_enamine_REAL_checks_inferences.py`
Scans inference logs (`output/unidock_docking/inferences_outputs` on the SBNB cluster) for errors/warnings once script 29 completes. Computes the surviving-compound set (all compounds minus script 31's `failed_REAL.csv`) and saves it (`success_mols.pkl`) as a fixed ordering; then, per pocket, reorders that pocket's raw inference probabilities (script 29) to match this fixed compound order and saves as a compressed `.npz` — stored in `output/unidock_docking/inference_probs`. This fixed ordering is what lets script 33 stack all 276 pockets' probabilities into one matrix directly.

### `33_enamine_REAL_selection.py`
Stacks all 276 pockets' `bin_01` inference probabilities (script 32) into one compound × pocket matrix, then z-score normalizes it twice — first per pocket (column-wise), then per compound (row-wise) — before ranking; the row-normalization is what "diversifies by synthons" in practice, by flattening the profile of any compound that scores well indiscriminately across most pockets rather than for a specific one. Per pocket, greedily selects the top `MAX_MOL=100,000` compounds by (doubly-normalized) score while capping any synthon (parsed from the compound ID's `____`-separated segments) at `MAX_SYN=5` occurrences. Also defines a global "inactive" set (~13k compounds) that are neither selected by any pocket nor share a synthon with any pocket's selection. Each pocket's `output/unidock_REAL_docking/input_ligands/input_ligands_<pocket>.txt` lists that pocket's 100k active molecules plus the full inactive set (same ~13k appended to every pocket's file).

### `34_enamine_REAL_mfps.py`
Calculates radius-3/2048-bit Morgan Count fingerprints (int8, clipped at 127) for all Enamine REAL 9.56M compounds, stored as one chunked, gzip-compressed HDF5 file (`output/enamine_REAL_characterization/enamine_REAL_ECFP6.h5`, `ids`/`fps` datasets, processed in 1M-row batches). Must be run with RDKit version 2025.9.1 — the ECFP6 bit ordering must match the version used in the related [ready-to-screen-enamine-real](https://github.com/ersilia-os/ready-to-screen-enamine-real) repository, which characterizes the 10-billion-compound library.

## Enamine REAL 9.56M: docking (100k prioritized compounds per pocket)

### `35_unidock_REAL_ligandprep.py`
Prepares script 30's 3D conformations for docking via `unidocktools`' `ligandprep` (`nebula`, `unidock` conda env), one of script 30's 96 chunks at a time rather than all at once (unlike script 20's single-batch approach for the 100k library): writes each chunk's non-failed (per script 31) compounds as individual SDFs to a temp folder, runs `ligandprep` (batch size 100) into that chunk's own numbered subdirectory, then deletes the temp folder before starting the next chunk, keeping disk usage bounded. Raises if a compound not listed in `failed_REAL.csv` has a missing conformer. Output: `output/unidock_REAL_docking/conformations_prepared/<chunk>/`.

### `36_unidock_REAL_docking.py`
Docks the 276 pocket structures against ~113k compounds per pocket with Uni-Dock (_fast_ search, _vina_ scoring, `nebula`, ~1 week compute). Output: `processed/unidock_REAL_docking/docking_results` (`docking.tar.gz` per pocket directory, 276 total). Protein structures reuse the preparation from script 21.

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
Prepares the ~99k compound conformations for docking via `unidocktools`' `ligandprep` (`nebula`, `unidock` conda env). Output: `processed/unidock_REAL_docking_2/conformations_prepared/`.

### `46_unidock_REAL_2_docking.py`
Docks the 276 pocket structures against the ~99k compound set with Uni-Dock (_fast_ search, _vina_ scoring, `nebula`). Output: `processed/unidock_REAL_docking_2/docking_results`, one `report.csv` per pocket.

### `47_docking_summary.py`
Console reporting tool, per gene, summarizing for every candidate pocket of a given tRNA synthetase: P2Rank pocket probability, InterPro domain annotation, AlphaFill co-crystallized-ligand evidence, structural confidence (minimum pLDDT for AF2/AF3/Chai-1 models, or GMQE for SwissModel models), and docking score percentiles (0.01%, 0.1%, 1%) from both docking libraries (100k Enamine DL and the ~99k Enamine REAL set). Used to manually curate one reference pocket per tRNA synthetase — the script prints instructions to record the choice in `output/reference_pocket_catalytic.csv` and `output/reference_pocket_noncatalytic.csv` (both `gene_name, pocket_name`; neither is written by this script itself), which are consumed by downstream hit-selection scripts.

### `47b_reference_pocket_visualization.py`
Builds one PyMOL session per manually-curated reference pocket (script 47's `output/reference_pocket_catalytic.csv`/`_noncatalytic.csv`, both checked automatically and merged into one session if a gene has both) for visual QC of that curation. For a given gene it loads the reference structure and pocket-centroid sphere, then overlays: the top-`TOP_N_LIGANDS` (5) best-scoring REAL round-2 docked poses (extracted from `docking.tar.gz`), every UniProt-cross-referenced experimental PDB structure (full biological assembly, aligned onto the reference via `cmd.align` restricted to one chain's pocket residues), the local AlphaFill structure with its transplanted ligands (same pocket-residue-restricted `cmd.align`), and optionally user-supplied homolog PDB structures (`--homologs` CSV) aligned via whole-chain `cmd.super` — chosen over `cmd.align`/`cmd.cealign` after empirically comparing RMSDs for this use case, per the script's own docstring.

**Inputs:** `output/reference_pocket_catalytic.csv`/`_noncatalytic.csv`, `output/pocket_detection_data.csv`, `output/aligned_relaxed_structures/`, `output/detected_pockets/`, `data/structures/alphafill_database/<uniprot_ac>/<uniprot_ac>.cif`, REAL round-2 docking results (`docking_utils.py`'s `LIBRARIES`/`load_gene_map`), plus live UniProt/RCSB REST API queries for cross-referenced PDB structures.

**Outputs:** `output/47b_reference_pocket_visualization/session_<uniprot_ac>_<gene>.pse`, with downloaded structures cached under `output/47b_reference_pocket_visualization/pdb_structures/<uniprot_ac>/` and `homolog_structures/<uniprot_ac>/` so repeat runs don't re-fetch them.

**Hardcoded parameters:** `TOP_N_LIGANDS = 5`; a fixed RGB per structure source (reference structure, docked/PDB/AlphaFill/homolog ligands); a 21-entry `UNIPROT_IDS_ORDERED` list giving each gene a consistent color from a tab20/tab20b palette; a fixed 100 Å zoom box.

Complete and documented (full `--help`); one H-bond-dashing helper is left commented out as an intentionally-disabled feature, not an unfinished TODO.

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

## Experimental multimer structures: receptor prep, docking & visualization

### `49_unidock_proteinprep_multimers.py`
Prepares script 48's stripped multimer structures for Uni-Dock, mirroring script 21's approach for the single-chain pipeline with one addition: script 48 skips PDB2PQR protonation/PyRosetta relaxation (pocket detection doesn't need them), so its stripped structures have no explicit hydrogens — unlike script 21's inputs, already protonated+relaxed by script 04. This script runs PDB2PQR (pH 7.0, AMBER force field, `adda4tb` conda env) first, then `unidocktools proteinprep` (`unidock_tools` conda env), so `proteinprep`'s RDKit-based charge/H-bond typing has explicit hydrogens to work with. Invoked as `python 49_unidock_proteinprep_multimers.py --pdb-codes 7K98,6XYZ`; resumable per PDB code (skips if the `.pdbqt` already exists). Depending on the installed `unidock_tools` build, `proteinprep` may shell out to AmberTools' `pdb4amber` — confirmed on `nebula`, resolved with `conda install -c conda-forge ambertools`.

**Inputs:** `output/48_detect_pocket_multimers/stripped_structures/<pdb_code>.pdb`.

**Outputs:** `output/49_unidock_proteinprep_multimers/<pdb_code>/<pdb_code>.pdb` (protonated) and `.pdbqt` (prepared receptor).

### `50_unidock_docking_multimers.py`
Docks one multimer receptor/pocket at a time against the same 99,105-compound Enamine REAL round-2 set used by scripts 44–46, so scores stay directly comparable — same box size (22.5 Å), seed (42), search mode (`fast`) and scoring function (`vina`) as script 46's `run_unidock` (also `max_gpu_memory=22000`, `num_modes=1`, `cpu=32`, `batch=500`). Resolves the receptor from script 49's output and the pocket centroid from `output/48_detect_pocket_multimers/pocket_detection_data.csv`; per-compound scores are parsed from the output SDFs' `ENERGY` field into `report.csv`, and the raw `docking/`/`logs.log` outputs are compressed to `.tar.gz` and deleted afterward. Resumable per pocket, via the same `REPORT_MIN_ROWS = 99105` completeness check script 46 uses to skip already-finished pockets. Must run on a GPU machine inside the `unidock_tools` conda env.

**Inputs:** `output/49_unidock_proteinprep_multimers/<pdb_code>/<pdb_code>.pdbqt`, `output/48_detect_pocket_multimers/pocket_detection_data.csv`, `output/unidock_REAL_docking_2/conformations_prepared/` (or its `.tar`).

**Outputs:** `output/50_unidock_docking_multimers/<pdb_code>_pocket_<n>/report.csv` (compound, score), plus `docking.tar.gz`, `logs.tar.gz`, `input_ligands.txt`.

### `51_selected_pockets_visualization.py`
Builds one combined PyMOL session covering every pocket curated in `output/selected_pockets.csv` (`gene_name, site_type, pocket_name, comment`) — the file that supersedes script 47's single-reference-pocket curation once multiple catalytic/non-catalytic sites per gene are being tracked (also consumed by `scripts/plots/figure_3_plot.py`/`figure_4_plot.py`). Each pocket becomes a merged structure + pocket-centroid-sphere object (`POCKET_SPHERE_COLOR = "blue"`, `POCKET_SPHERE_SCALE = 10`; surface representation by default, cartoon via `--no-surface`) plus its top-`--top-n` (default 5) best-scoring docked ligands, loaded but hidden (`cmd.disable`) so the session stays uncluttered by default. Branches per pocket on whether `pocket_name` contains `"_model_"` to resolve either the monomeric pipeline's structures/docking (scripts 04–08, 44–46) or the multimeric/experimental-PDB pipeline's (scripts 48–50).

**Inputs:** `output/selected_pockets.csv`; monomeric pockets: `output/aligned_relaxed_structures/`, `output/detected_pockets/`, `output/pocket_detection_data.csv`, `LIBRARIES["REAL"]`; multimeric pockets: `output/48_detect_pocket_multimers/{stripped_structures,detected_pockets,pocket_detection_data.csv}`, `output/50_unidock_docking_multimers/`.

**Outputs:** `output/51_selected_pockets_visualization/session_selected_pockets.pse`.

---

## Catalytic/non-catalytic hit selection (100k/9.56M/10M scale)

Scripts 52–56 share `src/docking_utils.py`'s helpers (`LIBRARIES`, `build_matrix`, `load_scores`, `load_real_negative_scores`/`load_real_positive_scores`, `lookup_smiles`, `sample_prescreened_smiles`, `compute_properties`, `save_upset`, `plot_score_boxplots`/`_multi`, `plot_profiling`, `plot_tsne`) and `src/default.py`'s `RANDOM_SEED = 42`. `compute_properties` computes MW/cLogP/TPSA/HBD/HBA/RotBonds/AromaticRings/QED (RDKit) and `is_pains` (RDKit PAINS A+B+C `FilterCatalog`).

### `52_CAT_promiscuous.py`
Multi-target hit-overlap ("promiscuity") analysis over the 4 catalytic (CAT) reference pockets (pheS, aspS, lysS, alaS — pheT has no CAT reference) against the 99,105-compound REAL round-2 screen. For each of `TOP_NS = [100, 1_000, 10_000]`: saves an UpSet plot of per-gene top-N overlap, and prints observed-vs-expected-by-chance counts for compounds hitting k=2,3,4 targets (expected = `top_n**k / M**(k-1)`, M = library size). Every compound hitting ≥2 targets within the top-1,000 (`PROFILING_TOP_N`/`PROFILING_MIN_TARGETS`) is saved with per-gene score/rank, SMILES and physicochemical properties, then profiled against a random `BG_SAMPLE_SIZE = 20,000`-compound REAL background (seed 42): a physchem plot, an ECFP4 t-SNE plot, and per-gene raw-score boxplots (fixed order `BOXPLOT_GENE_ORDER = ["pheS","aspS","lysS","alaS"]`) against the 100k DL library and the REAL round-1 negative set.

**Inputs:** `output/reference_pocket_catalytic.csv`; `LIBRARIES["REAL"]`/`LIBRARIES["DL"]` docking reports; REAL round-1 negative scores.

**Outputs**, `output/52_CAT_promiscuous/`:

| **File** | **Description** |
|-----------|------------------|
| `<genes>_CAT_upset_top{100,1000,10000}.png` | Per-top-N UpSet overlap plot |
| `<genes>_REAL_CAT_multi_target_top{100,1000,10000}.csv` | `compound, smiles, n_targets`, `score_<gene>`/`rank_<gene>` per gene, plus physicochemical properties |
| `<genes>_REAL_CAT_scores.png` / `_profiling.png` / `_tsne.png` | Score boxplots / physchem profile / t-SNE, for the ≥2-target top-1,000 subset |

### `53_CAT_selective.py`
Selectivity analysis for the same 4 CAT pockets, using rank-transformed docking scores (ascending, 0=best) against the 4 target pockets vs. the other 223 non-target pockets (also excluding pheT's own 13 pockets, since pheT is pheS's obligate PheRS partner). Selects ~500 compounds total via 5 complementary metrics, each capped at 50 except the last (`METRIC_CONFIGS`):

* **m1** — `max_target_rank` ascending, excluding `nontarget_p50 < 20,000` (max potency, low selectivity bar)
* **m2** — same, excluding `nontarget_p50 < 50,000` (max potency, high selectivity bar)
* **m3** — `nontarget_p10 − max_target_rank` descending, excluding `max_target_rank > 20,000`, a non-positive gap, or `nontarget_p50 < 20,000` (potency + selectivity gap)
* **m4** — `nontarget_p1` descending, excluding `max_target_rank > 20,000` (potency + max selectivity)
* **m5** — "diversity rescue": `top2_rank ≤ 20,000` and `nontarget_p50 ≥ 50,000`, not already selected by m1–m4, deduplicated to one compound per novel Murcko scaffold, capped at `500 − len(m1–m4 selection)`

Then builds UpSet plots (`UPSET_THRESHOLDS = [100, 1_000, 10_000]`, intersected with the selected set), score boxplots, physchem profiling and t-SNE for the full selection, plus a console report of overlap with script 52's top-1,000 multi-target list.

**Inputs:** `output/reference_pocket_catalytic.csv`, `LIBRARIES["REAL"]` (all 276 pockets, split target/non-target by UniProt AC), `LIBRARIES["DL"]`, REAL round-1 negative scores.

**Outputs**, `output/53_CAT_selective/`: `<genes>_REAL_CAT.csv` (`smiles, top2_rank, max_target_rank`, `rank_<gene>` per gene, `nontarget_p1/p5/p10/p50`, `m1`..`m5`, `selected` [which metric picked it], physicochemical properties), plus the same `_upset_top{100,1000,10000}.png`/`_scores.png`/`_profiling.png`/`_tsne.png` outputs as script 52.

### `54_NONCAT_promiscuous.py`
NON-CAT counterpart of script 52, over `output/selected_pockets.csv`'s 8 `site_type == "NON-CAT"` rows rather than a reference-pocket file — merging pheS+pheT into one "pheST" target (3 pockets: `alphafold3_P9WFU3_model_2_pocket_2`, `7K98_pocket_6`, `alphafold2_P9WFU1_model_0_pocket_1`) so there are 4 targets total, matching CAT's target count (pheST 3 pockets, aspS 2, lysS 1, alaS 2). A target's "top-N hit" is the union of its own pockets' top-N compound sets. For `TOP_NS = [100, 1_000]` (no top-10,000 here): per-target hit-set sizes, UpSet plots, observed-vs-expected-by-chance per exact target combination, and a final ≥2-target selection that additionally **excludes** compounds already CAT-promiscuous one level up — NON-CAT top-100 excludes script 52's CAT top-1,000 list, NON-CAT top-1,000 excludes CAT top-10,000 (`CAT_EXCLUSION_CSV`). Also reports (informational only) overlap with script 53's ~500 CAT-selective compounds, then profiling/boxplots/t-SNE for the top-1,000 selection.

**Inputs:** `output/selected_pockets.csv`, `output/52_CAT_promiscuous/alaS_aspS_lysS_pheS_REAL_CAT_multi_target_top{1000,10000}.csv`, `output/53_CAT_selective/alaS_aspS_lysS_pheS_REAL_CAT.csv`, `LIBRARIES["REAL"]` (monomeric pockets) and `output/50_unidock_docking_multimers/` (the `7K98_pocket_6` multimer pocket).

**Outputs**, `output/54_NONCAT_promiscuous/`: `<targets>_NONCAT_upset_top{100,1000}.png`, `<targets>_REAL_NONCAT_multi_target_top{100,1000}.csv` (per-pocket `score_<pocket>`/`rank_<pocket>` for all 8 pockets, `n_targets`, `targets_hit`, SMILES, physchem properties), plus `_scores.png`/`_profiling.png`/`_tsne.png` for the top-1,000 subset.

### `55_REAL10M_promiscuous.py`
Pocket-level (not protein-level) promiscuity screen over the `bin_01` surrogate model's predicted probabilities (scripts 29/37) for the full 9.56M-compound Enamine REAL library, across all 276 pockets independently. Per pocket, raw probabilities are rank-transformed to a [0,1] percentile (`scipy.stats.rankdata(method="min")`); per compound, percentiles-of-percentiles are computed across the 276 pockets at `PERCENTILE_LEVELS = [1,10,25,50,75,90,99]`. Filter: `p75 > 0.90` (`FILTER_LEVEL`/`FILTER_THRESHOLD`) — a compound ranking in the top 10% of predicted probability in ≥25% of the 276 pockets, deliberately more permissive than scripts 52–54's ≥2-distinct-gene requirement. Also reports (diagnostic only, doesn't affect the filter) `n_distinct_proteins_top10`: the number of distinct proteins (of 21) where the compound is top-10% in ≥1 of that protein's own pockets — UniProt AC parsed from each pocket filename via `PROTEIN_RE = r"^[^_]+_(P9W[A-Z0-9]+)_"`.

**Inputs:** per-pocket probability `.npz` files and `success_mols.pkl` (fixed compound order) under `output/unidock_docking/inference_probs/`; SMILES looked up from the `"REAL_ROUND1"` library.

**Outputs**, `output/55_identify_promiscuous_enamine_REAL/`: `promiscuous_hits.csv` (`id, smiles, p1, p10, p25, p50, p75, p90, p99, n_distinct_proteins_top10`) and `promiscuous_indices.pkl` (row indices into `success_mols.pkl`, reused by script 56).

### `56_NONCAT_top100_REAL10M.py`
Selects the top-100 non-promiscuous docking hits per NON-CAT pocket — 7 of the 8 `output/selected_pockets.csv` NON-CAT pockets, excluding the dimer-interface pocket `DIMER_POCKET = "7K98_pocket_6"`. Per pocket: loads its Enamine REAL **round-1** (10M) docking report, sorts ascending by score, drops any compound flagged by script 55's `promiscuous_hits.csv`, keeps the top `TOP_N = 100` and prints the 100th-place cutoff score. The 7 per-pocket tables are concatenated **without deduplication** — a compound landing in multiple pockets' top-100 appears once per pocket, reported explicitly rather than treated as an error. Also builds a per-pocket boxplot across 5 reference score distributions (Hit Locator DL, REAL round-1 negatives/positives, REAL round-2 pre-screened, and this script's own "REAL 1 – selected").

**Inputs:** `output/selected_pockets.csv`, `output/55_identify_promiscuous_enamine_REAL/promiscuous_hits.csv`, the REAL round-1 docking results directory, `LIBRARIES["DL"]`, `LIBRARIES["REAL"]`.

**Outputs**, `output/56_NONCAT_top100_selection/`: `top100_per_noncat_pocket.csv` (`gene_name, pocket_name, rank, compound, score, smiles`, physicochemical properties), plus `noncat_score_boxplots.png`, `noncat_top100_physchem_profile.png`, `noncat_top100_tsne.png`.

## Non-catalytic hit selection & docking at 10B scale

### `57_NONCAT_REAL10B_selective.py`
Identifies, for each of the 7 curated NON-CAT pockets (`output/selected_pockets.csv`, excluding the dimer pocket `DIMER_POCKET = "7K98_pocket_6"`), Enamine REAL 10B compounds selective for that pocket's target protein: top-1% predicted probability (`ind_1.npz`, the 99th-percentile threshold from the external screen) for that pocket, but *not* top-1% for any pocket belonging to a **different** protein — only set membership survives from the external screen, there's no continuous score to rank by. pheS/pheT are mutually exempted from counting as each other's "background" (`PARTNER_AC_OF = {"P9WFU3": "P9WFU1", "P9WFU1": "P9WFU3"}`, symmetric heterodimer exemption). Iterates all `N_CHUNKS = 994` pre-computed screening chunks from the external `gcadda4tb-enamine-real-screening` repo (read directly from `~/github/gcadda4tb-enamine-real-screening/output`, not copied in); for each chunk, opens `{chunk}.tar` and reads every pocket's `{pocket}_ind_1.npz` via `src/screening_10b_utils.py`'s `load_ind1`, computes `target_idx − union(background_idx)`, and writes one CSV per pocket per chunk. Resumable (skips chunks whose all-pocket output already exists); each chunk's ID→SMILES mapping (~90MB, downloaded on demand from a shared Google Drive folder via a service account) is deleted after use to keep `tmp/` bounded across all 994 chunks.

**Inputs:** `output/selected_pockets.csv`; `output/pocket_detection_data.csv` (via `get_pocket_to_ac()`); per-chunk `.tar` archives from the external screening repo's own output; per-chunk SMILES/ID mappings from Google Drive.

**Outputs:** `output/57_NONCAT_selective_10B/{gene}_{pocket}/{chunk}.csv` (`chunk, index, compound_id, smiles, ind1_threshold, gene_name, pocket_name`).

### `58_NONCAT_REAL10B_conformations.py`
Consolidates script 57's output: for each of the 7 NON-CAT pockets, concatenates all chunk CSVs, drops duplicate `compound_id`s (57's 994 chunks span 3 overlapping Enamine families), and randomly caps each pocket at `MAX_PER_POCKET = 100_000` (a random sample, `random_state=RANDOM_SEED`, since there's no score to rank by). Merges across pockets and drops any compound selective for more than one distinct pocket (kept for none). Generates a 3D conformer (RDKit ETKDGv3 + UFF minimization, seed 42) per surviving unique compound, using `N_WORKERS = 16` parallel workers.

**Inputs:** `output/selected_pockets.csv`; `output/57_NONCAT_selective_10B/{gene}_{pocket}/*.csv`.

**Outputs:** `output/58_generate_conformations_noncat_selective_10B/merged_selective_hits.csv`; per-compound SDFs at `output/58_generate_conformations_noncat_selective_10B/conformations/{compound_id}.sdf` (written atomically via a `.part` temp file + `os.replace`; resumable, skips SDFs already on disk including guarding against 0-byte truncated files from a killed prior run).

### `59_NONCAT_REAL10B_ligandprep.py`
Wraps script 58's conformers for docking via `unidocktools ligandprep` (`unidock_tools` conda env, `--batch_size 100`), mirroring script 45's approach: writes an index of all conformer SDF paths, runs `ligandprep` in batches, then rewrites the index restricted to files that survived preparation. Runs once over all conformers regardless of pocket — script 60 builds per-pocket subsets from this shared output.

**Inputs:** `output/58_generate_conformations_noncat_selective_10B/conformations/*.sdf`.

**Outputs:** `output/59_unidock_ligandprep_noncat_selective_10B/conformations_prepared/`, `output/59_unidock_ligandprep_noncat_selective_10B/input_ligands.txt`.

### `60_NONCAT_REAL10B_docking.py`
Docks script 59's prepared ligands against each of the 7 NON-CAT pockets as **per-compound subset docking**: each pocket only docks the compounds scripts 57/58 selected as selective for it, not a full cross-matrix. For each target, reuses an already-prepared `.pdbqt` receptor from the prior round-2 REAL docking (`output/unidock_REAL_docking_2/docking_results/`) if available, otherwise prepares one fresh via `unidocktools proteinprep` from `output/aligned_relaxed_structures/`. Box center comes from `output/pocket_detection_data.csv`; docking itself matches script 46's round-2 REAL docking settings — box size 22.5 Å, seed 42, `search_mode="fast"`, `scoring="vina"`, `cpu=32` — reusing script 46's own `run_unidock`/`extract_score_from_sdf`/`generate_report` helpers. Resumable per pocket (`report.csv` exists → skip) and per receptor (`.pdbqt` already prepared → skip). Must run inside the `unidock_tools` conda env on a GPU machine.

**Inputs:** `output/selected_pockets.csv`, `output/pocket_detection_data.csv`, `output/58_generate_conformations_noncat_selective_10B/merged_selective_hits.csv`, `output/59_unidock_ligandprep_noncat_selective_10B/conformations_prepared/`, `output/unidock_REAL_docking_2/docking_results/` (receptor reuse), `output/aligned_relaxed_structures/` (fallback receptor prep).

**Outputs:** `output/60_unidock_docking_noncat_selective_10B/docking_results/{pocket_name}/` — `report.csv` (compound, score), `docking.tar.gz`, `logs.tar.gz`.

### `61_NONCAT_top100_REAL10B.py`
Per NON-CAT pocket (again excluding the dimer pocket), sorts script 60's `report.csv` ascending by score and greedily selects the top `TARGET_N = 100` compounds under a strict "no synthon reused" constraint (`MAX_SYN = 1`, stricter than scripts 33/54's caps) — compound IDs' synthon segments parsed via `compound_id.split("____")[1:]`. Also builds a 6-group score-distribution boxplot (script 56's 5 groups — Hit Locator, REAL round-1 negatives/positives, REAL round-2 all, REAL round-1 selected — plus this script's own new "REAL 2 – selected" group), a physicochemical profiling plot, and an ECFP4 t-SNE plot of the selection vs. a background sampled evenly (`BG_SAMPLE_PER_POCKET = 2_000` per pocket) from script 57's *raw* per-pocket selective-hit pools (not script 58's capped pool).

**Inputs:** `output/selected_pockets.csv`; `output/57_NONCAT_selective_10B/`; `output/60_unidock_docking_noncat_selective_10B/docking_results/{pocket}/report.csv`; `output/58_generate_conformations_noncat_selective_10B/merged_selective_hits.csv` (SMILES lookup); `output/56_NONCAT_top100_selection/top100_per_noncat_pocket.csv`; `LIBRARIES["DL"]`, `LIBRARIES["REAL"]`, REAL round-1 positive/negative sets.

**Outputs:** `output/61_docking_top100_diverse_selection/top100_diverse_per_pocket.csv` (`compound_id, smiles, score, pocket_name`), plus `top100_diverse_score_boxplots.png`, `top100_diverse_physchem_profile.png`, `top100_diverse_tsne.png`. Prints a warning per pocket if fewer than `TARGET_N` compounds could be selected.

## Hit aggregation & full 12-pocket cross-matrix docking

### `62_aggregate_hits.py`
Merges the 5 independently-generated hit lists from scripts 52–61 into one deduplicated table — not itself a final hit list (script 61's REAL-10B branch only covers a partial screen of the full non-catalytic search space). For each compound, tracks which source(s) it came from (semicolon-joined, e.g. `"cat_promiscuous;cat_selective"`) and warns if the same compound ID shows conflicting SMILES across sources (keeps the first-seen SMILES). Computes physicochemical/PAINS properties for the aggregated set (`compute_properties`) and produces profiling and t-SNE (ECFP4) plots split by source (`SOURCE_ORDER = ["cat_promiscuous", "cat_selective", "noncat_promiscuous", "noncat_top100_10m", "noncat_top100_10b"]`, one fixed color each) — unlike scripts 52–54/56/61, there's no single shared background library across all 5 sources, so sources are plotted directly against each other instead.

**Inputs:** `output/52_CAT_promiscuous/alaS_aspS_lysS_pheS_REAL_CAT_multi_target_top1000.csv`, `output/53_CAT_selective/alaS_aspS_lysS_pheS_REAL_CAT.csv`, `output/54_NONCAT_promiscuous/alaS_aspS_lysS_pheST_REAL_NONCAT_multi_target_top1000.csv`, `output/56_NONCAT_top100_selection/top100_per_noncat_pocket.csv`, `output/61_docking_top100_diverse_selection/top100_diverse_per_pocket.csv`.

**Outputs**, `output/62_aggregate_hits/`: `aggregated_hits.csv` (`compound_id, smiles, source, MW, cLogP, TPSA, HBD, HBA, RotBonds, AromaticRings, QED, is_pains` — 2,923 rows, the final aggregated set every downstream script builds on), `aggregated_hits_physchem_profile.png`, `aggregated_hits_tsne.png`.

### `63_aggregated_conformations.py`
Generates a 3D conformer (RDKit ETKDGv3 + UFF, same mechanism as scripts 44/58) for every compound in script 62's aggregated hit list, using `N_WORKERS = 16` parallel workers and seed 42.

**Inputs:** `output/62_aggregate_hits/aggregated_hits.csv` (`compound_id, smiles`).

**Outputs:** `output/63_aggregated_conformations/conformations/{compound_id}.sdf` (written atomically via a `.part` temp file + `os.replace`; resumable, skips already-present/non-empty SDFs, reports the embedding-failure count).

### `64_aggregated_ligandprep.py`
Identical wrapper pattern to script 59, applied to script 63's conformers: builds a ligand index, runs `unidocktools ligandprep` (`--batch_size 100`, `unidock_tools` conda env), then rewrites the index to only successfully-prepared compounds.

**Inputs:** `output/63_aggregated_conformations/conformations/*.sdf`.

**Outputs:** `output/64_aggregated_ligandprep/conformations_prepared/`, `output/64_aggregated_ligandprep/input_ligands.txt`.

### `65_aggregated_docking.py`
Docks script 64's prepared ligands as one shared set against **all 12** curated pockets from `output/selected_pockets.csv` (4 CAT + 8 NON-CAT, including the dimer pocket) — a full cross-matrix, unlike script 60's per-pocket subset docking. Each pocket is docked `N_REPLICATES = 5` times with a different seed each run (42–46) to measure Uni-Dock's own run-to-run scoring variance; replicate scores are then averaged (mean + std) per compound. Receptor `.pdbqt` files are reused from wherever they already exist (script 60's, the round-2 REAL docking's, or script 49's outputs); for the dimer pocket specifically, if no prepared receptor exists yet, the stripped structure is protonated with PDB2PQR (AMBER force field, pH 7.0, `adda4tb` conda env) before `unidocktools proteinprep`, mirroring script 49's approach. Box size (22.5 Å), search mode (`fast`) and scoring (`vina`) match the rest of the pipeline. Resumable per replicate (checks `results_{r}.csv` exists with the expected row count) and per receptor. Must run inside the `unidock_tools` conda env on a GPU machine.

**Inputs:** `output/selected_pockets.csv`; `output/pocket_detection_data.csv`; `output/48_detect_pocket_multimers/pocket_detection_data.csv` (dimer centroid); `output/64_aggregated_ligandprep/conformations_prepared/`; receptor sources listed above; `output/48_detect_pocket_multimers/stripped_structures/` (dimer fallback); `output/aligned_relaxed_structures/` (monomer fallback).

**Outputs:** `output/65_aggregated_docking/docking_results/{pocket_name}/results/results_{r}.csv` (per replicate) and `results.csv` (`compound, score, score_std` — mean/std across the 5 replicates; verified 2,923 rows for every one of the 12 pockets, i.e. no compound is lost through conformer generation or ligand prep), plus `docking/docking_{r}.tar.gz`, `logs/logs_{r}.tar.gz`.

### `66_merge_docking_scores.py`
Merges script 65's 12 per-pocket `results.csv` files into one wide table keyed by compound, adding both the mean docking score and a rank for each pocket. The rank background is the 99,105-compound Enamine REAL round-2 library (`output/unidock_REAL_docking_2/docking_results/`) for monomeric pockets, or script 50's multimer docking output for the dimer pocket — computed via `np.searchsorted` of each compound's score into that pocket's sorted background-score array (lower rank = better/more negative score; a compound with no score for that pocket is ranked `+inf` then masked back to `pd.NA`). Also prints (console-only, not saved) summary diagnostics: the fraction of compounds with `QED > 0.5`, the fraction with `MW` in (300, 500), and the fraction of all compound×pocket score "endpoints" clearing various score/std cutoffs.

**Inputs:** `output/62_aggregate_hits/aggregated_hits.csv` (`compound_id, smiles, source, QED, MW`); `output/selected_pockets.csv`; `output/65_aggregated_docking/docking_results/{pocket}/results.csv`; `output/unidock_REAL_docking_2/docking_results/{pocket}/report.csv`; `output/50_unidock_docking_multimers/{pocket}/report.csv` (dimer background).

**Outputs:** `output/66_merge_docking_scores/merged_docking_scores.csv` — `compound_id, smiles, source`, then `{pocket_name}` (mean score) and `{pocket_name}_rank` per pocket, for all 12 pockets (27 columns total, verified). An internal `assert` guards that each per-pocket merge doesn't change row count (catching duplicate `compound_id`s).

## ADMET characterization, plotting & final filtering

### `67_ersilia_characterization.py`
Runs a caller-specified list of Ersilia Model Hub models (`--models`, comma-separated) over script 62's 2,923 aggregated hits, shelling out through the `ersilia` CLI per model in sequence: `catalog --local` check → `fetch` (skipped if already present) → `serve` → `run` → `close` → `delete` (always, even if the model pre-existed). Resumable at the per-model level (skips a model entirely if `{model_id}.csv` already exists). Models used so far (per the script's own usage help, not hardcoded elsewhere): `eos12x7` (Spatial Score: `sps_score`, `nsps_score`), `eos5jv3` (MycoPermeNet: `mycomembrane_permeation`), `eos2xeq` (structural alerts: `has_pains`, `has_brenk`, `is_sim_known_ab`, plus 4 antibiotic-class motif flags — nitrofuran/fluoroquinolone/carbapenem/beta-lactam), `eos42ez` (HepG2/HSkMC/IMR90 cytotoxicity), `eos3ujl`/`eos8d8a`/`eos1lb5` (three independent Mtb permeability models). Must run inside the `ersilia` conda env, on `herbert`.

**Inputs:** `output/62_aggregate_hits/aggregated_hits.csv` (`smiles` column only).

**Outputs:** `output/67_ersilia_characterization/smiles_input.csv` (re-saved in Ersilia's own input format) and one `{model_id}.csv` per requested model (`key, input`, plus that model's own output columns) — no merging across models happens here (scripts 68/70 do that separately).

### `68_plot_results.py`
Despite the "plot results" name, this script does three largely independent things, gated by two flags (`--only-upset` skips the first two; `--annotate` only affects the third):

1. **Ersilia output distributions** — histograms for `eos12x7` (`sps_score`, `nsps_score`), `eos5jv3` (`mycomembrane_permeation`), `eos42ez` (3 cytotoxicity endpoints); a permeability comparison merging `eos3ujl`/`eos8d8a`/`eos1lb5`'s `perm_proba_lepori_mtb` column (distribution + pairwise Spearman scatterplots); console-only reports (% with `nsps_score` in [10,40], % with all 3 cytotoxicity scores < 0.3, per-flag counts for `eos2xeq`'s 7 structural-alert columns — confirming `is_sim_known_ab` is 0 for all 2,923 compounds).
2. **Docking score boxplots** across the 12 curated pockets: Hit Locator, REAL round-1/round-2 pre-screened, "all selected" (script 65) and "top-`TOP_N=10` selected."
3. **Protein-level UpSet plots** for score thresholds `DOCKING_SCORE_THRESHOLDS = [-12, -11, -10, -9]` and rank thresholds `RANK_THRESHOLDS = [100, 500, 1000]` (rank vs. script 66's REAL round-2 background, ~99,105 compounds) — each computed for all pockets together and separately for CAT-only/NON-CAT-only subsets (pheS+pheT merged into "pheST").

**Inputs:** `output/67_ersilia_characterization/*.csv`, `output/selected_pockets.csv`, `output/65_aggregated_docking/docking_results/`, `output/66_merge_docking_scores/merged_docking_scores.csv`.

**Outputs**, `output/68_plot_results/`: `eos12x7_distributions.png`, `eos5jv3_distribution.png`, `eos42ez_distributions.png`, `permeability_model_comparison.png`, `aggregated_docking_score_boxplots.png`, `upset_score_{9,10,11,12}[_CAT|_NONCAT].png`, `upset_rank_{100,500,1000}[_CAT|_NONCAT].png`.

### `69_pymol_visualizations.py`
For each of the 5 curated genes (`--genes`, default all: pheS, pheT, aspS, lysS, alaS), builds one PyMOL session merging every curated pocket's structure + pocket-centroid sphere (`POCKET_SPHERE_SCALE = 10`) into one object per pocket, plus that pocket's top-`--top-n` (default 5) best-scoring compounds from script 65's replicate-averaged `results.csv` — poses are always taken from replicate 1's `docking_1.tar.gz` specifically, since the mean score across 5 replicates doesn't correspond to any single replicate's actual geometry, and replicate 1's pose is used as a representative one. Branches per pocket on whether the pocket name contains `"_model_"`, same monomeric-vs-multimeric resolution logic as scripts 51/54. `--no-surface` swaps to cartoon for faster rendering; no `zoom` call, since a gene's own pockets come from mutually unaligned structures.

**Inputs:** `output/selected_pockets.csv`, `output/65_aggregated_docking/docking_results/`, `output/aligned_relaxed_structures/`, `output/detected_pockets/`, `output/48_detect_pocket_multimers/{stripped_structures,detected_pockets}/`.

**Outputs:** `output/69_pymol_visualizations/session_<uniprot_ac>_<gene>.pse`.

### `70_filtering.py`
Joins scripts 62 (physchem), 66 (docking scores/ranks) and 67 (Ersilia outputs) into one wide table over the 2,923 aggregated hits, then applies a 7-rule sequential filter, printing the running pass count/percentage after each rule:

1. `QED > 0.5`
2. `300 < MW < 500`
3. `nsps_score` in `[10, 40]`
4. `has_pains == 0` (`eos2xeq`)
5. all 3 `eos42ez` cytotoxicity endpoints `< 0.3`
6. `mycomembrane_permeation` at or below its own 80th percentile (`EOS5JV3_TOP_PCT = 0.80`, computed from this run's own data, not a fixed value)
7. (≥2 CAT pockets with rank ≤10,000 **and** score ≤−10) **or** (≥1 NON-CAT pocket with rank ≤10,000 **and** score ≤−8) — rank against the 99,105-compound REAL round-2 background (script 66)

**Inputs:** `output/62_aggregate_hits/aggregated_hits.csv`, `output/66_merge_docking_scores/merged_docking_scores.csv`, `output/67_ersilia_characterization/{eos12x7,eos5jv3,eos2xeq,eos42ez}.csv`, `output/selected_pockets.csv`.

**Outputs**, `output/70_filtering/`: `merged_all_results.csv` (all 2,923 rows, the full joined table) and `filtered_hits.csv` (1,095 rows surviving all 7 rules; none flagged `is_sim_known_ab`, verified).

### `70b_compound_gifs.py`
Splits script 70's 1,095 filtered hits into `N_GROUPS = 6` sequential, roughly-equal, non-reordered chunks (`np.array_split`), then renders each chunk as an animated GIF (one molecule per frame) via the `chemgifs` CLI (`conda run -n chemgifs`, `github.com/ersilia-os/chemical-library-gifs`) — requires a not-yet-merged PR (`--style rdkit`) and a dedicated `chemgifs` conda env with an editable install of that branch.

**Inputs:** `output/70_filtering/filtered_hits.csv`.

**Outputs:** chunk CSVs at `processed/70b_compound_gifs/group_{i}of6.csv` (the only script in this batch writing under `processed/` rather than `output/`), and GIFs at `output/70b_compound_gifs/filtered_hits_group{i}of6.gif`, one fixed brand color per group (`GROUP_COLORS = ["purple","mint","blue","yellow","pink","orange"]`).

## Boltz-2 co-folding validation

### `71_boltz2_prepare_inputs.py`
Builds the two core Boltz-2 inputs for all 12 curated pockets (11 single-chain + the `7K98_pocket_6` pheS–pheT dimer pocket). For each pocket, parses the relevant PDB with Biopython and builds a `{pdb_resnum: 1-indexed_seq_position}` map by enumerating residues in file order (needed because numbering conventions differ across structure sources — SwissModel starts at residue 2/3, AF2/AF3/Chai1 at 1, the 7K98 crystal at −1/−2), then translates each pocket's raw PDB residue numbers into 1-indexed sequence positions (Boltz-2's own pocket-constraint format). For the dimer pocket specifically, both chains (`DIMER_CHAIN_A = "A"` pheS, `DIMER_CHAIN_B = "B"` pheT) are additionally trimmed around their pocket-contact window (`DIMER_TRIM_MARGIN_A = 30`, `DIMER_TRIM_MARGIN_B = 100` residues), because the untrimmed 1,177-residue complex OOMs Boltz-2 on a 24GB GPU — these exact margins were validated empirically (`ligand_iptm=0.984`, `complex_plddt=0.954`, RMSD 2.24 Å vs. the crystal structure on 587 matched Cα atoms, nearest excluded residue 10.57 Å from the ligand).

**Inputs:** `output/selected_pockets.csv`, `output/pocket_detection_data.csv`, `output/aligned_relaxed_structures/`, `output/48_detect_pocket_multimers/{stripped_structures,detected_pockets/7K98/7K98.pdb_predictions.csv}`, `output/70_filtering/filtered_hits.csv`.

**Outputs:** `output/71_boltz2_prepare_inputs/pocket_sequences.csv` (12 rows — the dimer row has extra `_b`-suffixed columns for chain B, plus `chain_a_trim_start/end`, `chain_b_trim_start/end`) and `output/71_boltz2_prepare_inputs/compounds.csv` (1,095 rows: `compound_id, smiles`).

### `72_boltz2_yaml_generation.py`
For every (pocket, compound) pair — **12 pockets × 1,095 compounds = 13,140 total** — writes one Boltz-2 YAML (`version: 1`, protein sequence + ligand SMILES + a `pocket` constraint listing contact residues + an `affinity` property request). The dimer pocket gets a 2-protein-chain YAML with the ligand as chain `C` instead of `B`. Each YAML's `msa:` field is baked in as an **absolute path on nebula** (`NEBULA_REPO_ROOT = "/home/admin/mtb-targeted-protein-degradation"`) rather than relative, since Boltz-2 resolves `msa:` relative to the process's CWD, not the YAML's own location — the referenced MSA cache file doesn't exist yet at generation time, script 73 creates it later. Skips YAMLs that already exist on disk (resumable).

**Inputs:** `output/71_boltz2_prepare_inputs/{pocket_sequences,compounds}.csv`.

**Outputs:** `output/72_boltz2_yaml_generation/input_yamls/<pocket_name>/<compound_id>.yaml` (verified: all 12 pocket subdirectories have exactly 1,095 YAMLs each).

### `73_boltz2_docking.py`
For each target pocket (`--pockets`, default: all 12 in `pocket_sequences.csv`): (1) bootstraps that pocket's MSA if not already cached, by writing a temporary msa-less YAML and calling `boltz predict --use_msa_server` on the shortest-SMILES candidate compounds (retrying up to `MAX_BOOTSTRAP_ATTEMPTS = 3` candidates per round, across `MAX_BOOTSTRAP_ROUNDS = 8` rounds, with a `MSA_REQUEST_BACKOFF_S = 120`-second sleep before every request to respect ColabFold's rate limit); the dimer pocket bootstraps both chains from one joint 2-chain call. (2) runs `boltz predict` in directory mode over script 72's YAML directory (or, with `--max-compounds`, a symlinked subset). (3) parses each compound's `affinity_<id>.json` and `confidence_<id>_model_0.json` into one row, re-aggregating the whole run's results into a CSV after every pocket (not just at the end), so an interrupted multi-day run always has an up-to-date summary. A pocket-level failure is caught and logged, not fatal to the rest of the batch. Must run inside the `boltz` conda env on a GPU machine (`nebula`).

**Inputs:** `output/72_boltz2_yaml_generation/input_yamls/`, `output/71_boltz2_prepare_inputs/{pocket_sequences,compounds}.csv`.

**Outputs:** `output/73_boltz2_docking/msa_cache/<pocket>.csv` (+ `_chainB.csv` for the dimer), `output/73_boltz2_docking/{out_subdir}/<pocket>/boltz_results_<pocket>/predictions/<compound_id>/{affinity,confidence}_*.json`, and `output/73_boltz2_docking/{out_subdir}_affinity_results.csv` (`pocket_name, compound_id, affinity_pred_value, affinity_probability_binary, affinity_pred_value1/2, affinity_probability_binary1/2, confidence_score, ptm, iptm, ligand_iptm, complex_plddt`).

### `74_boltz2_monitor.py`
Read-only progress-reporting utility: prints a per-pocket table (MSA cached yes/no, `n_structures/1095`, `n_affinities/1095`, status DONE/in-progress/not-started) by inspecting files already on disk under `output/73_boltz2_docking/{out_subdir}/`. Writes nothing; meant to be re-run repeatedly during the multi-day run on `nebula`. `N_COMPOUNDS = 1095` is a magic number, not derived from `compounds.csv` — must be kept in sync manually with scripts 71/73/75 if the compound set ever changes.

**Inputs:** `output/71_boltz2_prepare_inputs/pocket_sequences.csv`, `output/73_boltz2_docking/{msa_cache, <out_subdir>/<pocket>/boltz_results_<pocket>/predictions/}`.

**Outputs:** none (stdout only).

### `75_boltz2_collect_affinities.py`
Reads script 73's `{out_subdir}_affinity_results.csv`, joins in `gene_name`/`site_type` (from `pocket_sequences.csv`) and `smiles` (from `compounds.csv`), sorts, and writes one long-format CSV — no aggregation across pheS/pheT into "pheST," left to downstream consumers. Sanity-checks that `selected_pockets.csv`'s pocket list matches `pocket_sequences.csv`'s, reporting (not dropping) any unmatched pockets/compounds; prints a per-pocket completion-percentage table. Per pocket with data: an `affinity_pred_value` histogram (x-axis relabeled to IC50 µM via `10**log10_value` ticks) and an `affinity_pred_value` vs. `affinity_probability_binary` scatter annotated with Pearson r / Spearman ρ, plus a summary CSV of those correlations across all pockets.

**Inputs:** `output/73_boltz2_docking/{out_subdir}_affinity_results.csv`, `output/selected_pockets.csv`, `output/71_boltz2_prepare_inputs/{pocket_sequences,compounds}.csv`.

**Outputs:** `output/75_boltz2_collect_affinities/affinity_results.csv` (`gene_name, site_type, pocket_name, compound_id, smiles, affinity_pred_value, affinity_probability_binary, affinity_pred_value1/2, affinity_probability_binary1/2, confidence_score, ptm, iptm, ligand_iptm, complex_plddt`), `output/75_boltz2_collect_affinities/affinity_probability_correlations.csv` (`pocket_name, gene_name, site_type, n, pearson_r, spearman_rho`), `output/75_boltz2_collect_affinities/plots/<pocket_name>.png`.

## Pending review (no README text yet)

The following script exists on disk but hasn't been individually reviewed against documentation yet: `76_boltz2_docking_supervisor.sh`.
