# Targeted protein degradation in Mycobacterium tuberculosis
Discovery of potential degraders (BacPROTACS) for essential tRNA synthetases in _Mycobacterium tuberculosis_. This repository is work in progress.

## Table of Contents
- [Background](#background-)
- [Domain-specific requirements](#domain-specific-requirements-)
- [Progress reporting](#progress-reporting-)
  - [Sequence and structure annotation](#sequence-and-structure-annotation-of-trna-synthetases-)
    - [Pocket detection](#pocket-detection-)
    - [Pymol visualization](#pymol-visualization-)
    - [Data organization](#data-organization-)
    - [Additional analyses](#additional-analyses-)
    - [Pocket characterization](#pocket-characterization-)
    - [Protein prioritization](#protein-prioritization-)
- [TL;DR](#tldr-)
- [About the Ersilia Open Source Initiative](#about-the-ersilia-open-source-initiative-)

## Background 📚

This project is part of the [GC-ADDA4TB project](https://www.lifearc.org/project/grand-challenges-programme/), led by [Prof. Erick Strauss](https://scholar.google.com/citations?user=zK9kCVUAAAAJ&hl=en).
The project builds upon the [BacPROTACs technology](https://pubmed.ncbi.nlm.nih.gov/35662409/) in which small-molecule bifunctional degraders are designed to harness the proteolytic machinery for targeted protein degradation (TPD) of essential proteins in Mtb. This involves linking (a) a ligand that binds to a protein of interest (POI) to (b) a molecule that recruits the mycobacterial protease machinery, such as the ClpC:ClpP complex (ClpCP).

In a [CRISPR genetic screening](https://pubmed.ncbi.nlm.nih.gov/34297925/), several _Mtb_ tRNA synthetases were identified as highly vulnerable in _Mtb_. Here, targeting these tRNA synthetases with TPD is proposed as a novel therapeutic strategy.

In this project, the main objective is to prioritise a set of purchasable or easily synthesizable compounds with multi-target properties against the tRNA synthetases family. This will be achieved in three consecutive steps:
1. Structural annotation and binding pocket comparison across the tRNA synthetase family.
2. Large-scale virtual screening for compounds with strong predicted binding scores across multiple tRNA synthetase binding sites.
3. Final shortlisting using additional criteria, such as the ligand's amenability to being extended with a linker to dCymM without disrupting the binding pose.

## Domain-specific requirements 🛠️

We recommend creating a Conda environment to run this code. 🐍
```bash
conda create -n adda4tb python=3.10
conda activate adda4tb
pip install -r requirements.txt
```

In addition, the following tools are required:

* [Open source PyMol](https://github.com/schrodinger/pymol-open-source) for protein visualization.
* [PDB2PQR](https://pdb2pqr.readthedocs.io/en/latest/) for protonation at pH 7.0.
* [PyRosetta](https://www.pyrosetta.org/) for protein structure relaxation. PyRosetta can be installed with the [PyRosetta Installer](https://www.pyrosetta.org/downloads).
* [P2Rank](https://github.com/rdk/p2rank) for pocket detection. 
* [PocketVec](https://github.com/sbnb-irb/pocketvec) for pocket characterization.
* [Uni-Dock](https://github.com/dptech-corp/Uni-Dock) for protein-small molecule flexible docking.
* [LazyQSAR](https://github.com/ersilia-os/lazy-qsar) for ML bioactivity modelling. 


To run P2Rank, Java is required. Additionally, Openbabel needs to be installed for file format conversion:

```bash
conda activate adda4tb
conda install -c conda-forge openjdk
conda install -c conda-forge openbabel
```

## Progress reporting 📈

This repository is work in progress, as summarized in the following progress report meetings:

* [Check-in meeting 1](https://docs.google.com/presentation/d/1a7K4EkecYM63CPa7QRw1SOEGv2AEA3FYjpLJPYPpIkY/edit?usp=drive_link) (2025/02/04). Selection and structural annotation of tRNA synthetases.
* [Check-in meeting 2](https://docs.google.com/presentation/d/13RxQsi4-3t9LYxYGGvtfwdeYd3KJTor_07sc7zx3pwM/edit?usp=sharing) (2025/03/19). Binding site detection and visualization.
* [Check-in meeting 3](https://docs.google.com/presentation/d/1N6U4t6mNQmerMz0nD9AVmwABik2z512OxZ6imQbzguY/edit?usp=drive_link) (2025/04/17). Binding site characterization and comparison. Exploration of oligomerization.
* [Check-in meeting 4](https://docs.google.com/presentation/d/1ernke0imNeMticCVAVFo9MbsAVVPie80M7O59HP3DHY/edit?usp=sharing) (2025/05/12). Large-scale protein comparisons: sequence, structure and druggable pockets. 
* [Check-in meeting 5](https://docs.google.com/presentation/d/1o3fwydNJ1JIlEGcfe1XXN29SGv75FAlaJ-uKpIR_9GQ/edit?usp=sharing) (2025/06/18). Protein prioritization. An in-person [follow-up meeting](https://docs.google.com/presentation/d/1LGjTsAx_hhvZWtOQJLbE6Wa53H3UQWUpu0kLTxtnVlg/edit?usp=sharing) was held in Stellenbosch. 
* [Check-in meeting 6](https://docs.google.com/presentation/d/1newAljRfMfg3Jir2Jd1Uzkupof5a3rzpqgZf5P7-Fw4/edit?slide=id.g304d3f7c91c_0_552#slide=id.g304d3f7c91c_0_552) (2025/09/19) Docking Enamine DL HLL-100 (100k compounds), training surrogate models and screening Enamine REAL 9.56M. 
* [Check-in meeting 7](https://docs.google.com/presentation/d/1yv2c_GnoGhCxWjBDK28UrZYc5-cIy0aKm86_I_GSWXg/edit?slide=id.g39dda5591f2_0_552#slide=id.g39dda5591f2_0_552) (2025/11/21) Docking prioritized compounds from Enamine REAL 9.56M. 
* [Check-in meeting 8](https://docs.google.com/presentation/d/1feS_DK45iCSHywRhv-zUl2ittRGeEcKpsDj_YTkfGqw/edit?slide=id.g39dda5591f2_0_552#slide=id.g39dda5591f2_0_552) (2025/12/18) Surrogate modelling of Enamine REAL 9.56M. Moving towards billion scale inferences. 
* [Check-in meeting 9](https://docs.google.com/presentation/d/1EXYev6f3sdS1xrC2bDNniU_oZptOdp2Iqw4PyK-6L8s/edit?slide=id.g3b496d7a12e_0_552&pli=1#slide=id.g3b496d7a12e_0_552) (2026/03/04) Happening soon. Docking prioritized compounds from Enamine REAL 10B. Selecting compounds for experimental validation. 

Below, we explain the progress made chronologically.

### Sequence and structure annotation of tRNA synthetases 🧩

Based on the results of the CRISPR genetic screen ([Bosch et al, 2021; Figure 5](assets/bosch_2021_figure_5.jpg)), we have selected [21 essential tRNA synthetases](data/mtb_trna_synthetases_bosch_2021_fig5_annotated.csv). The UniProt AC, name, protein sequence and EC number have been obtained from UniProt ([Mtb H37RV proteome](data/mtb_h37rv_proteome.tsv)).

#### Protein structures

We have used the following servers or databases to obtain structural data for the selected tRNA synthetases. To ease the query of some resources, we have generated FASTA files for each tRNA synthetase using the `scripts/00_generate_fasta_files.py` script.

* [PDBe](https://www.ebi.ac.uk/pdbe/): Experimental structures (when available) can be found in `/data/structures/pdbe_database`. Note that these structures are often presented in multimeric form, and do not necessarily have full sequence coverage.
* [AlphaFold2 database](https://alphafold.ebi.ac.uk/): Predicted structures with AF2 were downloaded from the AF2-EBI database and stored in `/data/structures/alphafold2_database`. All structures in AF2 had 100% sequence coverage and are monomeric. Only one model was downloaded.
* [AlphaFold3 server](https://alphafoldserver.com/): We predicted structures with the AF3 server and downloaded them into `/data/structures/alphafold3_webserver`. Five models are available per query.
* [Chai-1 server](https://lab.chaidiscovery.com/dashboard): Likewise, we predicted structures with the Chai-1 server, ticking the MSA option. Results are stored in `/data/structures/chai1_server`. Five models per query were returned.
* [AlphaFill](https://alphafill.eu/): The AlphaFill resource was used to obtain AF2 structures along with ligands. We used the `/scripts/01_download_alphafill.py` script to programmatically download the structures and store them in `/data/structures/alphafill_database/`.
* [Swiss-Model](https://swissmodel.expasy.org/): The Swiss-Model server was used to obtain homology models for each sequence. They can be found in `/data/structures/swissmodel`. Note that full coverage is not guaranteed, and that we required a minimum of 80%. A variable number of models per query were returned.

The multiple structure files were organized in the `processed/structures` directory and stored both in `.cif` and `.pdb` formats. This was done with the `/scripts/02_organize_structures.py` script. This scripts ensures that only one chain is saved for each file, and that sequences are not chunked. Note that we omitted the PDBe files in this automated processing. Moreover, to simplify visualization, we aligned all structures using the `/scripts/03_align_structures.py` script. Based on RMSD, we removed structures that seemed to be far apart from the rest.

Then, we prepared these structures for docking with protein protonation with PDB2PQR and relaxation with PyRosetta using the `scripts/04_relax_structures.py` script. This procedure is computationally intensive. 

Afterwards, structures are aligned again with the `scripts/05_align_relaxed_structures.py` script, using their unrelaxed counterparts as reference for all the alignments. At this stage, no structures were removed, even those with high RMSD against the unrelaxed structures. 

#### Sequence data

We downloaded protein family and domain annotations from [InterPro](https://www.ebi.ac.uk/interpro/). Files can be found  in `data/sequences/interpro`. Sequence annotation data was processed using the `scripts/06_sequence_annotation.py` script.

#### Ligand data

In a first instance, we fetched data from [ChEMBL](https://www.ebi.ac.uk/chembl/) using the UniProt AC identifiers. This was done with the `scripts/07_fetch_from_chembl.py` script. We only found data for 3 of the 21 tRNA synthetases.

#### Aggregated data 📊

An aggregated file containing one row per processed structure is available in `/processed/trna_synthetases_data.csv`. This file contains the following information:

| **Field**              | **Description**                                                                                         |
|-------------------------|---------------------------------------------------------------------------------------------------------|
| `file_name`            | Name of the processed PDB structure file                                                               |
| `uniprot_ac`           | Uniprot AC identifier                                                                                  |
| `n_residues`           | Number of residues                                                                                     |
| `start_resid`          | First residue number (first residue is 1) of the sequence available in the PDB file, with respect to the Uniprot full sequence |
| `end_resid`            | Last residue number (first residue is 1) of the sequence available in the PDB file, with respect to the Uniprot full sequence |
| `coverage`             | Percentage sequence coverage                                                                           |
| `sequence_structure`   | Sequence found in the PDB file                                                                         |
| `full_sequence`        | Sequence found in UniProt                                                                              |

#### Pocket detection 🔍🎯

After the complete structural and sequential characterization of tRNA synthetases, we detected pockets in AF2, AF3, Chai-1 and SwissModel predicted models with `scripts/08_detect_pockets.py`. Following the [authors recommendations](https://github.com/rdk/p2rank/issues/76): for each structure, we considered detected pockets as those with a probability (K) > 0.2, but at least Top-3 (N) per structure. After that, we filtered those pockets having at least one residue with pLDDT < 65 (AF2, AF3, Chai-1) or QSQE < 0.65 (SwissModel), discarding about 25-30% of the pockets. Cut-offs are arbitrary; usual recommendations are 70 & 0.7, we’ve been slightly less restrictive.

A summary file containing one row per detected pocket and structure is available in `/processed/pocket_detection_data.csv`. This file contains the following information:

| **Field**                          | **Description**                                                    |
|-------------------------------------|--------------------------------------------------------------------|
| `Uniprot_AC`                       | Uniprot AC identifier                                             |
| `File name`                         | PDB file in which pockets have been detected                     |
| `Prediction type`                   | The method used for protein structure prediction                 |
| `Full path`                         | The full directory path where the PDB file is stored             |
| `Pocket number`                     | The identified pocket number within the structure                 |
| `Pocket score`                      | The score assigned to the detected pocket                        |
| `Pocket probability`                | The probability value indicating confidence in pocket detection  |
| `Pocket centroid coordinate (x y z)` | The (x, y, z) coordinates of the pocket’s centroid               |
| `Pocket residues (chain_resn)`      | List of residues forming the pocket, with chain and residue number |
| `B-factors`                         | Confidence measures: pLDDT (AF2, AF3, Chai) or QSQE (SM)           |

#### Pymol visualization 👀

We prepared PyMOL session files (`.pse`) to facilitate the visualization of detected pockets and their corresponding residues. These were generated using the `scripts/09_prepare_pymol_visualizations.py` script as step 09 in the pipeline.  

Each PyMOL session (one per protein) includes the following elements:  

| **Element**                                      | **Description**                                                                                  | **Displayed?** |
|-------------------------------------------------|----------------------------------------------------------------------------------------------|--------------|
| **Reference structure (AF2)**                    | Wheat color, surface + cartoon representation                                                | ✅ Yes       |
| **Pockets detected in reference structure (AF2)** | Sky-blue spheres with arbitrary size (pocket detection provides a single 3D coordinate)      | ✅ Yes       |
| **Residues defining each pocket in AF2**         | Orange color, surface + cartoon representation                                              | ✅ Yes       |
| **Aligned structures (all but AF2)**             | Gray color, cartoon representation                                                          | ❌ No        |
| **Pockets detected in aligned structures**       | Gray-colored points                                                                         | ❌ No        |
| **InterPro annotations**                         | Includes conserved sites, domains, families, homologous superfamilies etc (red color, surface representation) | ❌ No        |


#### Data organization 🗂️

We then organized sequence and protein information using scripts `scripts/10_organize_sequence_info.py` and `scripts/11_organize_pocket_info.py`, respectively. Sequence summary information can be found at `processed/sequences/interpro_summary.tsv`, and include the following features:

| **Field**                  | **Description**                                                                            |
|----------------------------|--------------------------------------------------------------------------------------------|
| `Accession`                | InterPro Accession identifier                                                              |
| `Name`                     | Name of the InterPro entry                                                                 |
| `Type`                     | Type of annotation (e.g., domain, conserved site, homologous superfamily)                 |
| `Number of proteins`       | Number of proteins associated with the InterPro entry                                     |
| `Average Coverage`         | Average sequence coverage across the associated proteins                                  |

Additionally, the file contains final binary columns (1/0) indicating the presence or absence of each annotation across different proteins. An additional file including manually curated InterPro annotations (i.e. Catalytic domain (ATP binding site), Anticodon Binding Domain, etc.) is found in `processed/sequences/interpro_summary_curated.tsv`.

Pocket summary information together with InterPro annotations can be found at `processed/pocket_detection_data_interpro.csv`, and includes the following features.

| **Field**                          | **Description**                                                    |
|-------------------------------------|--------------------------------------------------------------------|
| `Uniprot_AC`                       | Uniprot AC identifier                                             |
| `File name`                         | PDB file in which pockets have been detected                     |
| `Prediction type`                   | The method used for protein structure prediction (e.g., AlphaFold2, AlphaFold3) |
| `Full path`                         | The full directory path where the PDB file is stored             |
| `Pocket number`                     | The identified pocket number within the structure                 |
| `Pocket score`                      | The score assigned to the detected pocket                        |
| `Pocket probability`                | The probability value indicating confidence in pocket detection  |
| `Pocket centroid coordinate (x y z)` | The (x, y, z) coordinates of the pocket’s centroid               |
| `Pocket residues (chain_resn)`      | List of residues forming the pocket, with chain and residue number |
| `Confidence score`                  | Confidence measures: pLDDT (AF2, AF3, Chai) or QSQE (SM)          |
| `Interpro ID`                       | Identifier of the matched InterPro domain                        |
| `Interpro name`                     | Name of the matched InterPro domain                              |
| `Interpro Matches`                  | Residue ranges corresponding to the InterPro domain             |
| `Residues in pocket`                | Number of residues forming the pocket                           |
| `Residues in Interpro`              | Number of residues forming the InterPro domain                  |
| `Residues overlap`                  | Number of residues shared between the pocket and the InterPro domain |
| `Coverage pocket`                   | Fraction of pocket residues that overlap with an InterPro domain |
| `Coverage domain`                   | Fraction of InterPro domain residues overlapping with the pocket |

For the sake of simplicity, those pocket-InterPro pairs having no overlapping residues have been omitted in the file. 

#### Additional analyses 🧪

In `scripts/12_align_PDB_structures.py`, experimental structures (e.g., from PDBe) are aligned to predicted models to evaluate spatial coherence between structure sources. The script also checks for overlaps between detected pockets and known ligand binding sites, with results summarized in `processed/pdbe_annotation_report.csv`. Meanwhile, `scripts/14_calculate_SeqId.py` computes pairwise sequence identities between the 21 tRNA synthetases using global alignment (Needleman–Wunsch algorithm). At the structural level, `scripts/15_calculate_StSim.py` evaluates the structural similarity between tRNA synthetases, using Pymol. Main outputs are a matrix of sequence identities (`processed/sequences/NW_SeqAlign/SeqId_matrix.tsv` with their corresponding aligned proportions in `processed/sequences/NW_SeqAlign/Prop_matrix.tsv`) and set of CSV files including RMSDs between all structure pairs among all 21 tRNA synthetases (`processed/structural_comparisons`). 

#### Pocket characterization 📐

Pocket characterization was performed using PocketVec descriptors ([Comajuncosa-Creus et al., Nat Commun 2024](https://www.nature.com/articles/s41467-024-52146-3)). Docking calculations were performed in an HPC cluster allowing CPU parallelization and PocketVec descriptors were calculated from raw docking scores in `scripts/16_calculate_PocketVec.py`. Several PocketVec descriptors (12/276) were filtered out due to the excessive presence of outlier compounds (>80), following the authors' recommendations. In line with these, we used a PocketVec distance threshold of 0.17 to classify any pocket pair of interest as similar. PocketVec descriptors are found in `/processed/pocketvec_RUN/fps_rank.pkl`. 

#### Protein prioritization ⚠️

tRNA synthetases have been prioritized taking multiple factors into account. In brief, we exhaustively enumerated protein pairs (210) and triplets (1,330) and classified them according to PocketVec distance, Global Similarity (both structural and sequential) and the number of pockets mapped to the catalytic domain (see `notebooks/17a_Protein_prioritization.ipynb`). Global similarity was established at a 35% sequence identity and 10Å RMSD cut-offs for sequential and structural similarity, respectively. We extended comparisons at the pocket level (using PocketVec descriptors), accounting for 32,561 pairs and 2,499,258 triplets, and collected lenient (PocketVec distance < 0.17) and stringent (PocketVec distance < 0.14) sets of PocketVec descriptor pairs (1,481 and 76) and triplets (3,880 and 807) in `notebooks/17b_Protein_prioritization_pairs.ipynb` and `notebooks/17b_Protein_prioritization_triplets.ipynb`. Finally, we performed an intra-set normalization to derive a 'Prioritization score' per protein, considering how many times the protein appeared in the collected sets, therefore indicating similarity to other tRNA synthetases and potential polypharmacology (see `notebooks/17c_Protein_prioritization_individual.ipynb`). Final results are summarized in `processed/protein_prioritization/final_results.tsv` and include the following information:

| **Field**                          | **Description**                                                    |
|-------------------------------------|--------------------------------------------------------------------|
| `Uniprot_AC`                       | Uniprot AC identifier                                             |
| `Gene name`                         | Standard gene name associated with the Uniprot AC identifier     |
| `Vulnerability score`                   | Vulnerability score derived from [Bosch et al, 2021](https://pubmed.ncbi.nlm.nih.gov/34297925/) |
| `Score`                         | Prioritizazion Score             |
| `Tier`                     | Protein-associated tier                 |
| `sim_Tier1-5`                     | Number of proteins in the tier N having a PocketVec distance < 0.14      |

## High-throughput virtual screening ⚡

After the structural characterization of 21 essential tRNA synthetases, we search for potential small-molecule binders in an active learning fashion. First, we dock a chemically diverse set of 100k compounds against each pocket structure and we then train a surrogate model with [LazyQSAR](https://github.com/ersilia-os/lazy-qsar). 

Computational resources:

- `herbert`
- `norrsken-gpu-wsl` (NVIDIA GeForce RTX 4090)
- `SBNB-IRB cluster` (HPC cluster with /aloy/home)

By default, scripts were executed at `herbert` unless specified otherwise. 

### Characterization of ENAMINE Diversity Library HLL-100 100k compounds ⚗️

We selected the ENAMINE Diversity Library HLL-100 (Sublibrary of HLL-460) as our starting set of small molecules. First, we calculated Morgan Count Fingerprints and [CheMeleon](https://github.com/JacksonBurns/chemeleon) embeddings for 100,157 compounds (see `scripts/18_enamine_mfps.py` and `scripts/24_enamine_chemeleon.py`, respectively. ⚠️ CheMeleon embeddings were calculated using Ersilia’s NVIDIA GeForce RTX 4090. Fingerprints and embeddings are stored in `/processed/enamine_characterization`). We then generated 3D conformations for these compounds using the ETKDGv3 protocol together with UFF energy minimization, totalling 100,154 conformations (see `scripts/19_enamine_conformations.py`, conformations also located in `/processed/enamine_characterization/conformations.tar.gz`). 90% of these small molecules were in the [270,400] Da range.

### Docking ENAMINE Diversity Library HLL-100 (100k compounds)

#### Ligands

100,154 compounds with 3D conformations were prepared using the `ligandprep` functionality of `unidocktools` (see `scripts/20_unidock_ligandprep.py`, executed in `norrsken-gpu-wsl` with a `unidock` conda environment). Results are found in `processed/unidock_docking/conformations_prepared`.

#### Protein structures

178 protein structures (totalling 276 pocket structures, see `processed/pocket_detection_data.csv`) were prepared using the `proteinprep` functionality of `unidocktools` (see `scripts/21_unidock_proteinprep.py`, executed in `norrsken-gpu-wsl` with a `unidock` conda environment having `ambertools` and `reduce` installed as well). 8 of them required manual intervention to fix PDB2PQR/Pyrosetta protonation issues. Results are stored in `processed/unidock_docking/structures_prepared.tar.gz`.

#### Uni-Dock docking (I)

Protein-small molecule docking was performed with [Uni-Dock](https://github.com/dptech-corp/Uni-Dock) with search mode _fast_ and _vina_ scoring function using Ersilia’s NVIDIA GeForce RTX 4090 (276 pocket structures x 100k compounds ~ 1 week of computation time at `norrsken-gpu-wsl`, see `scripts/22_unidock_docking.py`). Docking results are stored in `processed/unidock_docking/docking_results`, inside the `docking.tar.gz` file within each pocket directory (276). A systematic search of errors and warnings is perfomed once the docking is completed, scanning the logs (`logs.tar.gz` files) with `scripts/23_unidock_checks.py`.

#### Training surrogate model (I)

For each pocket structure (x276) a ML model was trained using [LazyQSAR](https://github.com/ersilia-os/lazy-qsar) and [CheMeleon](https://github.com/JacksonBurns/chemeleon) descriptors at 3 binarization percentiles (0.1%, 0.5% and 1%, see `scripts/25_enamine_binarization.py` and `scripts/26_train_models`). Binarized reports mapping compound IDs, docking scores and 3 bin flags are available at `processed/unidock_docking/binarized_reports`. Models were trained at the `SBNB-IRB cluster`, allowing for intensive parallelization (3 bins x 276 pockets = 828 jobs), and saved in `processed/unidock_docking/models` (joblib files). After analyzing the results, only models trained on 0.1% binarized data were used in subsequent steps. 

### Characterization, screening and prioritization of Enamine REAL 9.56M compounds

We selected Enamine REAL 9.56M as the library to screen using surrogate models previously trained on Enamine Diversity Library (100k compounds). First, we calculated [CheMeleon](https://github.com/JacksonBurns/chemeleon) embeddings for 9.56M compounds using Ersilia’s NVIDIA GeForce RTX 4090 at `norrsken-gpu-wsl` (see `scripts/27_enamine_REAL_chemeleon.py`. Embeddings were stored at `/processed/enamine_REAL_characterization/embeddings`, 96 files having ~100k compounds each, ~25GB). Embeddings were then transferred to the `/aloy/home` partition at the `SBNB-IRB cluster` through `dante` (see `scripts/28_transfer_embeddings.sh`). Finally, leveraging the parallelization capacity at the `SBNB-IRB cluster`, surrogate models were used to screen Enamine REAL 9.56M compounds (see `scripts/29_inference_enamine_REAL`, results stored in `processed/unidock_docking/inferences`). A systematic search of errors and warnings is perfomed once the inference is completed, scanning the logs (`processed/unidock_docking/inferences_output` directory) with `scripts/32_enamine_REAL_checks_inferences.py`. Additionally, 
only the probabilities associated to those compounds that did not fail in subsequent steps (scripts `30_enamine_REAL_conformations.py` and `31_enamine_REAL_check_conformations.py`) were kept and stored in `processed/unidock_docking/inference_probs`. This is the reason behind moving from script `29_inference_enamine_REAL` to script `32_enamine_REAL_checks_inferences.py` in this section. Finally, a curated set of 100k compounds was prioritized (highest-scoring by that pocket’s inference probabilities) per each pocket structure (x276), while diversifying by synthons, unweighting those compounds having high inference probabilities systematically (see `scripts/33_enamine_REAL_selection.py`). In addition, we defined a global set of “inactive” molecules (~13k) whose synthons never appear in any pocket’s active set. The list of compounds associated to each pocket structure (100k+13k) is stored in `processed/unidock_REAL_docking/input_ligands`. 

### Characterizing Enamine REAL 9.56M compounds

We first generated conformations for all compounds in Enamine REAL 9.56M using the ETKDGv3 protocol together with UFF energy minimization (see `scripts/30_enamine_REAL_conformations.py`, conformations stored in pkl format for each chunk under `/processed/enamine_REAL_characterization/conformations`). After that, we checked for which compounds the conformation generation protocol was not successful (see `scripts/31_enamine_REAL_checks_conformations.py`) and stored their IDs in `processed/enamine_REAL_characterization/failed_REAL.csv`. Additionally, we sampled 25k compounds and stored them in `processed/enamine_REAL_characterization/random_sample_REAL.csv`, to be used in subsequent analyses. 

We then calculated Morgan Count Fingerprints for all compounds in Enamine REAL 9.56M (see `scripts/34_enamine_REAL_mfps`). Fingerprints were stored in `processed/enamine_REAL_characterization/enamine_REAL_ECFP6.h5` for all molecules. It is crucial that those fingerprints are generated with RDKit version 2025.9.1, mathing the version used in a [related reporsitory](https://github.com/ersilia-os/ready-to-screen-enamine-real) to characterize 10 bilion compounds (see steps below). 
 
### Docking Enamine REAL 9.56M (100k prioritized compounds per pocket)

#### Ligands

Successfully processed Enamine REAL 9.56M compounds with 3D conformations were prepared using the `ligandprep` functionality of `unidocktools` (see `scripts/35_unidock_REAL_ligandprep.py`, executed in `norrsken-gpu-wsl` with a `unidock` conda environment). Results are found in `processed/unidock_REAL_docking/conformations_prepared`.

#### Protein structures

178 protein structures (totalling 276 pocket structures, see `processed/pocket_detection_data.csv`) were already prepared in previous steps. 

#### Uni-Dock docking (II)

Protein-small molecule docking was performed with [Uni-Dock](https://github.com/dptech-corp/Uni-Dock) with search mode _fast_ and _vina_ scoring function using Ersilia’s NVIDIA GeForce RTX 4090 (276 pocket structures x 113k compounds per pocket ~ 1 week of computation time at `norrsken-gpu-wsl`, see `scripts/36_unidock_REAL_docking.py`). Docking results are stored in `processed/unidock_REAL_docking/docking_results`, inside the `docking.tar.gz` file within each pocket directory (276). 

#### Training surrogate model (II)

For each pocket structure (x276) a Naive Bayes Classifier was trained using Morgan Count Fingerprints (RDKit version 2025.9.1) at a 1% binarization threshold on docking scores (1130 actives per pocket), see `scripts/37_surrogate_model.py`. Models were stored at `processed/unidock_REAL_docking/training_results/models`. 

### Characterization, screening and prioritization of Enamine REAL 10B compounds

The characterization of Enamine REAL 10B compounds was performed in a [related repository](https://github.com/ersilia-os/ready-to-screen-enamine-real/tree/main).

The screening of Enamine REAL 10B compounds was performed in a [related repository](https://github.com/ersilia-os/gcadda4tb-enamine-real-screening). 


## TL;DR ⏱️

We’re developing BacPROTAC-based degraders targeting 21 essential tRNA synthetases in *Mycobacterium tuberculosis*. For each of these tRNA synthetases:

1. Sequence annotation via InterPro.
2. Structural characterization: 
    - Data sources: AF2, AF3, Chai-1, SwissModel, etc.
    - Relaxation: PyRosetta
    - Pocket detection: P2Rank
    - Pocket characterization and clustering: PocketVec
3. Protein prioritization based on potential polypharmacology
4. Docking with Uni-dock (I): 276 pocket structures against 100k compounds from Enamine Chemical Diversity
5. Surrogate modelling (I) to prioritize compounds from Enamine REAL 9.56M
6. Docking with Uni-dock (II): 276 pocket structures against 113k compounds (per pocket) from Enamine REAL 9.56M
7. Surrogate modelling (II) to prioritize compounds from Enamine REAL 10B



## Related repositories

1. [ready-to-screen-enamine-real](https://github.com/ersilia-os/ready-to-screen-enamine-real/tree/main): The Enamine REAL 10 billion space characterized using ECFP6s. 
2. [gcadda4tb-enamine-real-screening](https://github.com/ersilia-os/gcadda4tb-enamine-real-screening):


## About the Ersilia Open Source Initiative 🌍❤️

This repository is developed by the [Ersilia Open Source Initiative](https://ersilia.io). Ersilia develops AI/ML tools to support drug discovery research in the Global South. To learn more about us, please visit our [GitBook Documentation](https://ersilia.gitbook.io) 🌐 and our [GitHub profile](https://github.com/ersilia-os/).
