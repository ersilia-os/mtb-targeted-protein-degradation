# Targeted protein degradation in Mycobacterium tuberculosis
Discovery of potential degraders (BacPROTACS) for essential tRNA synthetases in _Mycobacterium tuberculosis_. This repository is work in progress.

## Table of Contents
- [Background](#background-)
- [Domain-specific requirements](#domain-specific-requirements-)
- [Progress reporting](#progress-reporting-)
- [TL;DR](#tldr-)
- [Related repositories](#related-repositories)
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
* [BitBirch](https://github.com/mqcomplab/bitbirch) (`bblean` package) for fingerprint clustering.
* [Ersilia Model Hub](https://github.com/ersilia-os/ersilia) CLI for running pre-trained ADMET/property prediction models (e.g. NSPS, cytotoxicity).

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
* [Check-in meeting 3](https://docs.google.com/presentation/d/11D_Woyr3-8nlitXF87XaACPnfCvDAlQoX372Xf6A24Y/edit?usp=sharing) (2025/04/17). Binding site characterization and comparison. Exploration of oligomerization.
* [Check-in meeting 4](https://docs.google.com/presentation/d/1ernke0imNeMticCVAVFo9MbsAVVPie80M7O59HP3DHY/edit?usp=sharing) (2025/05/12). Large-scale protein comparisons: sequence, structure and druggable pockets. 
* [Check-in meeting 5](https://docs.google.com/presentation/d/1o3fwydNJ1JIlEGcfe1XXN29SGv75FAlaJ-uKpIR_9GQ/edit?usp=sharing) (2025/06/18). Protein prioritization. An in-person [follow-up meeting](https://docs.google.com/presentation/d/1LGjTsAx_hhvZWtOQJLbE6Wa53H3UQWUpu0kLTxtnVlg/edit?usp=sharing) was held in Stellenbosch. 
* [Check-in meeting 6](https://docs.google.com/presentation/d/1newAljRfMfg3Jir2Jd1Uzkupof5a3rzpqgZf5P7-Fw4/edit?slide=id.g304d3f7c91c_0_552#slide=id.g304d3f7c91c_0_552) (2025/09/19) Docking Enamine DL HLL-100 (100k compounds), training surrogate models and screening Enamine REAL 9.56M. 
* [Check-in meeting 7](https://docs.google.com/presentation/d/1yv2c_GnoGhCxWjBDK28UrZYc5-cIy0aKm86_I_GSWXg/edit?slide=id.g39dda5591f2_0_552#slide=id.g39dda5591f2_0_552) (2025/11/21) Docking prioritized compounds from Enamine REAL 9.56M. 
* [Check-in meeting 8](https://docs.google.com/presentation/d/1feS_DK45iCSHywRhv-zUl2ittRGeEcKpsDj_YTkfGqw/edit?slide=id.g39dda5591f2_0_552#slide=id.g39dda5591f2_0_552) (2025/12/18) Surrogate modelling of Enamine REAL 9.56M. Moving towards billion scale inferences. 
* [Check-in meeting 9](https://docs.google.com/presentation/d/1EXYev6f3sdS1xrC2bDNniU_oZptOdp2Iqw4PyK-6L8s/edit?slide=id.g3b496d7a12e_0_552&pli=1#slide=id.g3b496d7a12e_0_552) (2026/03/04) Docking Enamine REAL 10B (99,105 compounds) against the 276 pockets; selecting multi-target hit compounds to prioritize proteins and compounds for experimental validation.
* [Check-in meeting 10](https://docs.google.com/presentation/d/1MJa19RVH7gbu1dTwRH8VSxudkzq9U4kN36GLETrX7JM/edit?slide=id.g3f1bf6e8f71_0_53#slide=id.g3f1bf6e8f71_0_53) (2026/07/16) Prioritizing 4 tRNA synthetases (pheST, aspS, lysS, alaS) for experimental analyses, and final compound prioritization for experimental validation.

For the full script-by-script pipeline documentation — data sources, thresholds, and rationale for each step — see [`scripts/README.md`](scripts/README.md).

## TL;DR ⏱️

We're developing BacPROTAC-based degraders targeting 21 essential tRNA synthetases in *Mycobacterium tuberculosis*. For each of these tRNA synthetases:

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
8. Screening Enamine REAL 10B compounds and reducing hits to a candidate pool (promiscuous + selective selection, ~1M compounds)
9. Drug-likeness filtering and BitBirch clustering to a final ~99k validation set
10. Docking with Uni-dock (III): 276 pocket structures against the ~99k compound set
11. Multi-target and selectivity-driven hit selection against curated reference pockets
12. Cytotoxicity prediction (Ersilia `eos42ez`) and final QED/PAINS/cytotoxicity filtering to a shortlist for experimental validation

## Related repositories

1. [ready-to-screen-enamine-real](https://github.com/ersilia-os/ready-to-screen-enamine-real/tree/main): The Enamine REAL 10 billion space characterized using ECFP6s. 
2. [gcadda4tb-enamine-real-screening](https://github.com/ersilia-os/gcadda4tb-enamine-real-screening): Surrogate-model screening of the Enamine REAL 10 billion compound space against the 276 pocket structures.

## About the Ersilia Open Source Initiative 🌍❤️

This repository is developed by the [Ersilia Open Source Initiative](https://ersilia.io). Ersilia develops AI/ML tools to support drug discovery research in the Global South. To learn more about us, please visit our [GitBook Documentation](https://ersilia.gitbook.io) 🌐 and our [GitHub profile](https://github.com/ersilia-os/).
