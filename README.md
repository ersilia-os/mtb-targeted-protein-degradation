# Targeted protein degradation in Mycobacterium tuberculosis
Discovery of potential degraders (BacPROTACS) for essential Mtb genes

## Background

This project is part of the GC-ADDA4TB project, led by Prof. Erick Strauss. The project builds upon the BacPROTAC technology in which small-molecule bifunctional degraders are designed to harness the proteolytic machinery for TPD of essential proteins in Mtb. This involves linking (a) a ligand that binds to a protein of interest (POI) to (b) a molecule that recruits the mycobacterial protease machinery, such as the ClpC:ClpP complex (ClpCP).

In this project, the main objective is to prioritise a set of purchasable or easily synthesizable compounds with multi-target properties against the tRNA family. This will be achieved in three consecutive steps:
1. Structural annotation and binding pocket comparison across the tRNA synthetase family.
2. Large-scale virtual screening for compounds with strong predicted binding scores across multiple tRNA synthetase binding sites.
3. Final shortlisting using additional criteria, such as the ligand’s amenability to being extended with a linker to dCymM without disrupting the binding pose.

## Methodology

### Annotation of tRNA synthetases

We will start by annotating Mtb tRNA synthetases with several bioinformatics resources. For each tRNA synthetase provided by the GC-ADDA4TB team, besides the essentiality and vulnerability scores, we will collect the protein sequence (UniProt), the domain annotations (InterPro, PFAM, etc.), the structures available (PDB, AlphaFold), and known ligands from ChEMBL, BindingDB, and other databases. All of this information will be compiled in a tRNA synthetase knowledgebase that Ersilia will make available to the GC-ADDA4TB team.

### Binding site identification and comparison

We will then search for putative binding sites across tRNA synthetases and compare them exhaustively to identify pockets with similar characteristics, i.e., those with the potential to bind the same ligand. For this purpose, we will use PocketVec 5, a novel end-to-end approach that has demonstrated robust performance at the proteome-wide level, utilizing both PDB and AlphaFold2 structures. By the end of this analysis, we will generate a network of binding sites linked based on their similarity.

Given that multiple structures or structural models are likely available for a single tRNA synthetase, and that more than one binding site may be detected per structure, the network of binding sites will need to be aggregated to facilitate navigation and clustering. At a high level, the ultimate goal will be to identify subsets of tRNA synthetases that are likely to be engaged by the same ligand. This will be achieved using standard network analysis techniques such as the Markov Clustering Algorithm (MCL).

## Virtual screening for pan-engaging ligands

Broadly speaking, given a binding site of interest, there are two approaches to identify ligands: virtual screening (docking) 6 and de novo design (generative models) 7. In this SOW1, we will focus on the virtual screening of (ultra) large-scale (make-on-demand) libraries for two reasons. First, the hits identified through this approach will already be purchasable (or easily synthesizable), thereby reducing the time required for experimentation. Second, we are addressing a “multi-objective” problem (i.e., binding to multiple binding sites), which is not easily tackled using generative models.

Significant computational workflow design will be required to ensure the efficiency of the virtual screening process and to appropriately rank compounds, as multiple docking scores will be available for each compound. In brief, for each cluster of binding sites of interest, we will develop a composite scoring method to maximize the likelihood that a given ligand will bind to multiple tRNA synthetases.

Regarding the screening library, the goal will be to achieve a billion-scale throughput, in principle using the Enamine REAL collection. This can be further discussed as the project progresses. Acceleration procedures will be necessary at such scales. We will consider one of two approaches: (a) active learning on non-enumerated building blocks 8 or (b) surrogate docking modeling with machine learning 9. An additional filtering based on docking pose will be applied to ensure that the hits can be attached to a linker that does not collide with the binding site.

We aim at providing a list of ~100 virtual hits (i.e. molecules with multi-target potential against tRNA synthetases). Of these, 10-30 will be further shortlisted for experimental testing. The shortlisting approach will depend on the desired diversity of hits, manual inspection, and, if necessary, additional evaluations such as molecular dynamics simulations.


## About the Ersilia Open Source Initiative

This repository is developed by the [Ersilia Open Source Initiative](https://ersilia.io). Ersilia develops AI/ML tools to support drug discovery research in the Global South.
