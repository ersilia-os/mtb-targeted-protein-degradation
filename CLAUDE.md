# MTB Targeted Protein Degradation

## Project purpose

Computational discovery of **BacPROTACs** — bacterial PROteolysis TArgeting Chimeras — for targeted protein degradation of essential tRNA synthetases in *Mycobacterium tuberculosis*. Part of the **GC-ADDA4TB** (Grand Challenges) programme led by Prof. Erick Strauss. The 21 target proteins were selected from a CRISPR fitness screen (Bosch 2021). The goal is to identify multi-target virtual screening hits for experimental validation.

---

## Repository layout

```
scripts/          # 47 numbered Python scripts — canonical pipeline (run in order)
notebooks/        # Jupyter notebooks for exploration, QC, and analysis
data/             # Raw inputs: structures, sequences, ligand libraries
output/           # Pipeline outputs: pockets, docking results, models
processed/        # Intermediate compute artefacts: alignments, embeddings, conformations
tools/            # External utility tools
```

---

## Pipeline overview

| Phase | Scripts | What it does |
|-------|---------|--------------|
| Structure preparation | 00–05 | FASTA generation, structure download (6 sources), alignment, relaxation |
| Sequence & pocket analysis | 06–16 | InterPro, ChEMBL, P2Rank pocket detection, PocketVec, SeqId, StSim |
| Protein prioritization | 17 + notebooks 17a–17c | Pair/triplet scoring, 5-tier vulnerability ranking |
| 100k compound screening | 18–26 | Enamine DL 100k MFPs, 3D conformations, Uni-Dock, surrogate models |
| 9.56M inference + docking | 27–37 | CheMeleon embeddings, RF inference, synthon-aware selection, docking, NB models |
| 10B screening & hit selection | 38–42 | External repo inference, A/B set reduction, drug-likeness filtering |
| Final validation docking | 44–47 | Conformations for selected hits, Uni-Dock, summary table |

Script 26 lives at `scripts/26_train_models/main.py`; script 29 at `scripts/29_inference_enamine_REAL/main.py`.

---

## Key data files

| File | Description |
|------|-------------|
| `data/mtb_trna_synthetases_bosch_2021_fig5_annotated.csv` | 21 target proteins — UniProt ACs, gene names, essentiality scores |
| `output/pocket_detection_data.csv` | 276 detected pockets with P2Rank scores and residues |
| `output/pocket_detection_data_interpro.tsv` | Pockets mapped to InterPro functional domains |
| `output/protein_prioritization/final_results.tsv` | 5-tier prioritization ranking (Tier 1 = P9WFU3 pheS) |
| `output/unidock_docking/binarized_reports/` | 100k docking binarized at 0.1%, 0.5%, 1% thresholds |
| `output/unidock_REAL_docking/inference_probs/` | Surrogate model probabilities for 9.56M compounds |
| `data/Enamine/` | Compound libraries: DL 100k, REAL 9.56M, REAL 10B SMILES |
| `data/CheMeleon/` | Pre-trained MPNN embedding weights |
| `data/PocketVec_descriptors/` | Raw PocketVec descriptor files |

---

## External tools

| Tool | Purpose | Notes |
|------|---------|-------|
| **P2Rank** | Pocket detection | Java; run via shell |
| **Uni-Dock** | GPU-accelerated docking | Fast mode + Vina scoring; requires NVIDIA GPU |
| **PyRosetta** | Protein relaxation | Requires academic licence |
| **PDB2PQR** | Protonation state assignment | pH 7.0 |
| **OpenBabel** | Format conversion | PDB → MOL2/SDF for PocketVec |
| **LazyQSAR** | Surrogate QSAR models | Ensemble ML; script 26 |
| **CheMeleon** | MPNN molecular embeddings | 2048-dim; scripts 24, 27 |
| **unidocktools** | Ligand/protein PDBQT preparation | Scripts 20, 21, 35, 45 |

---

## Computational resources

| Host | Role |
|------|------|
| **herbert** | Primary compute server |
| **norrsken-gpu-wsl** | NVIDIA RTX 4090 — docking and CheMeleon embeddings |
| **SBNB-IRB cluster** | HPC with `/aloy/home/` — 10B inference |
| **Google Drive** | Storage for 10B chunk results (downloaded by script 41) |

Script 28 (`28_transfer_embeddings.sh`) uses `rsync` to move embeddings from GPU host to HPC.

---

## External repositories

- **`ersilia-os/gcadda4tb-enamine-real-screening`** — 10B virtual screening pipeline (results referenced by scripts 38–40)
- **`ersilia-os/ready-to-screen-enamine-real`** — 10B compound characterization

---

## Conventions & gotchas

- **Script numbering is execution order.** Run scripts sequentially unless a script is explicitly noted as optional or parallel.
- **RDKit ≥ 2025.9.1 required** for script 34 — the ECFP6 bit ordering changed between versions.
- **Synthon-aware selection** (script 33): max 5 occurrences per synthon ensures chemical diversity in the 9.56M → 113k selection step.
- **Drug-likeness filter** (script 42): MW [250–450], LogP [−1, 5], QED > 0.4, NSPS [10–40].
- **Docking binarization** at 0.1%, 0.5%, 1% of score distribution; three separate model sets trained per threshold.
- **Two-round docking**: 100k Enamine DL for model training → 113k selected from 9.56M for validation.
- **Pocket quality gate**: pLDDT (AlphaFold) or QSQE (homology models) scores used to filter low-confidence structures before docking.
- `output/` holds final artefacts meant to be kept; `processed/` holds intermediate files that can be regenerated.
- 276 pockets across 178 structures from 21 proteins (multiple structure sources per protein).
- Notebooks 17a–17c implement the prioritization logic that script 17 stubs out.

---

## Working with the user

- **Ask, don't assume.** For any non-trivial decision — which approach to take, which dataset to use, what to name something, whether to add a dependency, how to handle an ambiguous case — use the `AskUserQuestion` tool BEFORE editing. Two short questions up front beat a wrong-direction edit.
- **Plans are mandatory.** Anything beyond a one-line fix or a pure read-only investigation MUST go through plan mode. Propose the plan and stop until the user confirms. Never skip the plan step.
- **Surface uncertainty.** When you have multiple reasonable options or are unsure about intent, name them and ask. Don't pick silently.

---

## Hard requirements

- All Python plotting must use the [stylia](https://github.com/ersilia-os/stylia) library. Invoke the `/stylia-plotting` skill for guidance. If the skill is not installed, ask the user to install it.
- Do **not** create new folders at the root level outside the ones listed in the repository layout above.
- Outputs in `output/` should follow the same numbering as the script that produced them.

---

## Ersilia scientific tools

- **Ersilia Model Hub.** For any task involving prediction of small-molecule properties (bioactivity, ADMET, toxicity, target affinity, embeddings, generative chemistry), check the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia) before writing custom code. Browse with `ersilia catalog`; fetch with `ersilia fetch <eos_id>`; serve and run with `ersilia serve <eos_id>` then `ersilia api -i input.csv -o output.csv`. Record the model ID in the script header.
- Do not reimplement a predictor when an Ersilia model already covers it. If no suitable model exists, surface the gap to the user.
- Other Ersilia repositories ([github.com/ersilia-os](https://github.com/ersilia-os)) may contain relevant utilities. Check before writing similar tooling from scratch.

---

## Version control conventions

- **Git** tracks code only: `scripts/`, `notebooks/`, `src/`, `tools/`, `docs/`, `assets/`
- **eosvc** (Ersilia Version Control) tracks data: `data/` and `output/` are linked to an S3 bucket and excluded from git
- `access.json` records whether data/output are public or private
- Empty folders are preserved with `.gitkeep` files; remove `.gitkeep` once a folder has real content
- Do not commit data, outputs, temporary files, secrets, credentials, or API keys

---

## Scientific rigor

- **Citations must be real.** Never invent paper titles, authors, DOIs, journal names, or publication years. Only cite sources the user provided, sources already referenced in the repo, or sources verified through a real web fetch.
- **Claims need sources.** When asserting a scientific fact — a methodology choice, a threshold convention, a biological or chemical mechanism — name the source. Distinguish "the data shows X" (observation) from "X works because Y" (claim that needs a citation).
- **Record dataset provenance.** When pulling from public sources (ChEMBL, PubChem, UniProt, Enamine, etc.), record the version or snapshot date in `scripts/README.md` or as a comment in the downloading script.
- **Set random seeds.** Any script using stochastic methods (train/test split, sampling, model training) must set and record a seed. Use a project-wide `RANDOM_SEED` constant in `src/default.py`.

---

## Python naming conventions

- Project-wide constants reused across scripts must be defined in `src/default.py` and named in `ALL_CAPS`.
- Scripts that import from `src/` must include this path setup at the top, before any `src` imports:

```python
import os
import sys
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
```

- Declare input and output folder paths as module-level variables at the top of the script and ensure they exist:

```python
data_dir = os.path.join(root, "..", "data", "processed")
output_dir = os.path.join(root, "..", "output")
os.makedirs(data_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)
```

---

## README guidelines

- **Root README**: high-level and brief (~50 lines). Cover what the project is, how to get started, main commands, and key outputs. Do not replicate the folder structure.
- **`scripts/README.md`**: one entry per script with a brief description (1–3 sentences). If a script encodes a key decision (a threshold, cutoff, model choice), state the value and its rationale explicitly.
- **`docs/`**: methodology notes, literature summaries, decision logs, and AI-generated reports all belong here, not at the repo root. Use `YYYY-MM-DD_topic.md` or `NN_topic.md` naming.

---

## Human sign-off required

These actions must never be taken autonomously — always explain and ask first:

- **Thresholds and cutoffs:** Never choose, apply, or hardcode a threshold or filtering criterion. Propose options with reasoning and let the user decide.
- **Dropping data:** Never remove data points, even obvious outliers or NaN values. Flag them and ask.
- **Interpreting scientific results:** Do not infer conclusions from analysis outputs. Present what the data shows and ask for the user's interpretation.
- **Deleting files:** Never delete files without explicit confirmation — including old scripts or superseded outputs.
- **Raw data is read-only:** Never modify or overwrite files in `data/raw/`. All transformations produce new files in `data/processed/`.

---

## Available Ersilia skills

Ersilia maintains skills at [ersilia-os/ersilia-skills](https://github.com/ersilia-os/ersilia-skills). Check that repository for the current list. Relevant installed skills include `/stylia-plotting`, `/molecule-auditing`, `/model-incorporation-code`, and `/deep-research`.
