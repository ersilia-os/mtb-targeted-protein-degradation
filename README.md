# Targeted protein degradation in Mycobacterium tuberculosis
Discovery of potential degraders (BacPROTACS) for essential tRNA synthetases in _Mycobacterium tuberculosis_.

## Background

This project is part of the [GC-ADDA4TB project](https://www.lifearc.org/project/grand-challenges-programme/), led by [Prof. Erick Strauss](https://scholar.google.com/citations?user=zK9kCVUAAAAJ&hl=en).
The project builds upon the [BacPROTACs technology](https://pubmed.ncbi.nlm.nih.gov/35662409/) in which small-molecule bifunctional degraders are designed to harness the proteolytic machinery for targeted protein degradataion (TPD) of essential proteins in Mtb. This involves linking (a) a ligand that binds to a protein of interest (POI) to (b) a molecule that recruits the mycobacterial protease machinery, such as the ClpC:ClpP complex (ClpCP).

In a [CRISPR genetic screening](https://pubmed.ncbi.nlm.nih.gov/34297925/), several _Mtb_ tRNA synthetases were identified as highly vulnerable in _Mtb_. Here, targeting these tRNA synthetases with TPD is proposed as a novel therapeutic strategy.

In this project, the main objective is to prioritise a set of purchasable or easily synthesizable compounds with multi-target properties against the tRNA synthetases family. This will be achieved in three consecutive steps:
1. Structural annotation and binding pocket comparison across the tRNA synthetase family.
2. Large-scale virtual screening for compounds with strong predicted binding scores across multiple tRNA synthetase binding sites.
3. Final shortlisting using additional criteria, such as the ligand's amenability to being extended with a linker to dCymM without disrupting the binding pose.

## Domain-specific requirements

We recommend to create a Conda environment to run this code.
```bash
conda create -n adda4tb python=3.10
conda activate adda4tb
pip install -r requirements.txt
```

In addition, the following tools are required:

* [Open source PyMol](https://github.com/schrodinger/pymol-open-source) for protein visualization.
* [PDB2PQR](https://pdb2pqr.readthedocs.io/en/latest/) for protonation at pH 7.0.
* [PyRosetta](https://www.pyrosetta.org/) for protein structure relaxation. PyRosetta can be installed with the [PyRosetta Installer](https://www.pyrosetta.org/downloads).
* [P2RANK](https://github.com/rdk/p2rank) for pocket detection. 


To run P2RANK, Java is required:

```bash
conda activate adda4tb
conda install -c conda-forge openjdk
```

## Progress reporting

This repository is work in progress, as summarized in the following progress report meetings:

* [Check-in meeting 1](https://docs.google.com/presentation/d/1a7K4EkecYM63CPa7QRw1SOEGv2AEA3FYjpLJPYPpIkY/edit?usp=drive_link) (2025/02/04). Selection and structural annotation of tRNA synthetases.
* Check-in meeting 2 (scheduled mid March 2025)

Below, we explain the progress made chronologically.

### Sequence and structure annotation of tRNA synthetases

Based on the results of the CRISPR genetic screen ([Bosch et al, 2021; Figure 5](assets/bosch_2021_figure_5.jpg)), we have selected [21 essential tRNA synthetases](data/mtb_trna_synthetases_bosch_2021_fig5_annotated.csv). The UniProt AC, name, protein sequence and EC number have been obtained from UniProt ([Mtb H37RV proteome](data/mtb_h37rv_proteome.tsv)).

#### Protein structures

We have used the following servers or databases to obtain structural data for the selected tRNA synthetases. To ease the query of some resources, we have generated FASTA files for each tRNA synthetase using the `scripts/00_generate_fasta_files.py` script.

* [PDBe](https://www.ebi.ac.uk/pdbe/): Experimental structures (when available) can be found in the [this subfolder](/data/structures/pdbe_database). Note that these structures are often presented in multimeric form, and do not necessarily have full sequence coverage.
* [AlphaFold2 database](https://alphafold.ebi.ac.uk/): Predicted structures with AF2 were downloaded from the AF2-EBI database and stored in [this subfolder](/data/structures/alphafold2_database). All structures in AF2 had 100% sequence coverage and are monomeric. Only one model was downloaded.
* [AlphaFold3 server](https://alphafoldserver.com/): We predicted structures with the AF3 server and downloaded them into [this subfolder](/data/structures/alphafold3_webserver). Five models are available per query.
* [Chai-1 server](https://lab.chaidiscovery.com/dashboard): Likewise, we predicted structures with the Chai-1 server, ticking the MSA option. Results are stored in [this subfolder](/data/structures/chai1_server). Five models per query were returned.
* [AlphaFill](https://alphafill.eu/): The AlphaFill resource was used to obtain AF2 structures along with ligands. We used the `/scripts/01_download_alphafill.py` script to programmatically download the structures and store them in [this subfolder](/data/structures/alphafill_database/).
* [Swiss-Model](https://swissmodel.expasy.org/): The Swiss-Model server was used to obtain homology models for each sequence. They can be found in [this folder](/data/structures/swissmodel). Note that full coverage is not guaranteed, and that we required a minimum of 80%. A variable number of models per query were returned.

The multiple struture files were organized in the [processed data subfolder](processed/structures) and stored both in `.cif` and `.pdb` formats. This was done with the `/scripts/02_organize_structures.py` script. This scripts ensures that only one chain is saved for each file, and that sequences are not chunked. Note that we ommitted the PDBe files in this automated processing. Moreover, to simplify visualization, we aligned all structures using the `/scripts/03_align_structures.py` script. Based on RMSD, we removed structures that seemed to be far apart from the rest.

Then, we prepared these structures for docking with protein protonation with PDB2PQR and relaxation with PyRosetta using the `scripts/04_relax_structures.py` script. This procedure is computationally intensive. 

#### Sequence data

We downloaded protein family and domain annotations from [InterPro](https://www.ebi.ac.uk/interpro/). Files can be found [here](data/sequences/interpro). Sequence annotation data was processed using the `scripts/05_sequence_annotation.py` script.

#### Ligands data

In a first instance, we fetched data from [ChEMBL](https://www.ebi.ac.uk/chembl/) using the UniProt AC identifiers. This was done with the `scripts/06_fetch_from_chembl.py` script. We only found data for 3 of the 21 tRNA synthetases.

#### Aggregated data

An aggregated file containing one row per processed structure is available [here](/processed/trna_synthetases_data.csv). This file contains the following information:

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


## About the Ersilia Open Source Initiative

This repository is developed by the [Ersilia Open Source Initiative](https://ersilia.io). Ersilia develops AI/ML tools to support drug discovery research in the Global South. To learn more about us, please visit our [GitBook Documentation](https://ersilia.gitbook.io) and our [GitHub profile](https://github.com/ersilia-os/).
