# 77_pocket_annotation

Corrected redo of domain annotation + pocket-domain overlap (scripts 06/09/10/11/12), plus ligand-evidence scoring they never had. **Coexists, doesn't replace** — old scripts/outputs untouched; script 47 and notebooks reading `pocket_detection_data_interpro.tsv` keep working. Numbered 77 for lack of a free slot, but only depends on scripts 00-08.

Run `./run_all.sh`, or steps individually (`01`-`05`, `09`, `10` = plain python3; `06`-`08`, `11` = `adda4tb` env, needs PyMOL).

## Bugs fixed vs. the old pipeline
- Script 10's cross-protein coverage filter drops correct catalytic annotations on compact proteins (pheS got 0 domain-labeled pockets).
- Keyword curation had substring false-positives (e.g. "Glutamine synthetase" matched a bare `synthetase` regex).
- gatA/gatB use a different catalytic mechanism (transamidation, not aminoacyl-adenylate) and got zero coverage under the old GO-term map.
- Script 11 returns empty for pheS.
- Script 12's whole-chain alignment breaks (10-14 Å RMSD) on proteins with flexible domains — fixed via **per-pocket local alignment**.
- AlphaFill's ligand-provenance JSON (source PDB, `local_rmsd`, identity) is parsed for the first time; the old pipeline only used raw coordinate proximity.

## Steps
1. `fetch_interpro` — live InterPro API pull per protein
2. `build_annotation_table` — dedup + per-protein coverage (informational)
3. `categorize` — GO-term → category, curated overrides, `entry_support_count`
4-5. query + download experimental PDB structures
6. `identify_chains` — chain + numbering-offset auto-detection
7-8. pocket-local-aligned ligand evidence (direct-PDB, then AlphaFill)
9. `cluster_pockets` — pocket-centroid deduplication/clustering (reproducibility, not scored)
10. `assemble_final_table` → `output/77_pocket_annotation/pocket_detection_interpro_updated.csv`
11. `build_pymol_sessions` → `output/77_pymol_sessions/<gene>_pocket_annotation.pse`

## Pocket-centroid clustering (`09_cluster_pockets.py`)
Ported from `scripts/plots/figure_1_calculations.py`'s approach: per protein, sort pocket-instances by `Pocket score` descending, greedily accept a pocket as a new distinct site only if its centroid is farther than a **fixed 6.14 Å** from every already-accepted centroid (`POCKET_DEDUP_DISTANCE_THRESHOLD`, empirically derived in `notebooks/08_coherence_detected_pockets.ipynb`). Every pocket-instance (accepted or not) is then assigned to its nearest accepted site to produce a full `spatial_cluster_id` per instance; site IDs are renumbered by size (largest = 1). Resulting per-protein distinct-pocket counts match `figure_1_calculations.py`'s `gene_to_unique_pocket_count` exactly.

## AlphaFill ligand evidence (`08_extract_alphafill_evidence.py`)
Ligand coordinates are re-fetched immediately after each pocket's own `cmd.align` call, since `cmd.align` mutates the AlphaFill object's real atom positions on every iteration — checking against a coordinate fetched before that pocket's alignment would compare it to the wrong pocket. Matches `07_align_and_extract_ligands.py`'s pattern. `11_build_pymol_sessions.py` computes ligand proximity independently via its own PyMOL alignment and doesn't read these evidence CSVs.

## `has_direct_ligand_evidence`
Pocket-level: `yes` only if *this specific* (protein, structure file, pocket number) has a direct-PDB ligand within the proximity cutoff.

## `catalytic_confidence` (0-4)
**0** whenever the pocket lacks the curated Catalytic Domain label — no ligand evidence can raise it above 0 without that label.
Otherwise: base **1** for the label alone, **+1** (→2) if any weak ligand is nearby.
Ceilings (override the base, still require the label): **4** = strong ligand in a real experimental structure of this protein; **3** = strong ligand only via AlphaFill.
"Strong" = substrate / ATP-ADP-AMP(-APC) cofactor / per-protein reaction-intermediate mimic / curated validated inhibitor (`ligand_classification.py`); everything else weak.

**A low score ≠ invalid pocket** — it means no catalytic-specific evidence was found, not that the pocket isn't real.
