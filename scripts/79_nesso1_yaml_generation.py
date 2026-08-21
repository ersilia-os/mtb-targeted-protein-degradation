#!/usr/bin/env python3
"""
Generates Nesso-1 YAML input files from script 78's structure_sequences.csv and script 71's
compounds.csv -- one YAML per (structure, compound) pair, 11 structures (10 unique single-chain
+ the 7K98 dimer) x 1,095 compounds = 12,045 total.

Nesso-1's YAML schema (docs/prediction.md in recursionpharma/nesso) has no `msa:` key (no MSA
required at all -- ESM-2 embeddings are computed/cached internally per sequence) and no
`constraints:` key (pocket-conditioning isn't supported yet, see script 78's docstring), so each
YAML is just sequences + an affinity property.

Directory layout matters: `nesso predict <dir>` is NOT recursive (only top-level .yaml/.yml
files are picked up, confirmed by reading nesso/main.py's check_inputs()), so each structure's
1,095 compound YAMLs must live in their own flat subdirectory. The YAML filename stem becomes
the output record_id (confirmed in the repo's AGENTS.md and main.py's _process_single_yaml), so
files are named <compound_id>.yaml.

Running Nesso-1 itself is a separate, later script (80).

Usage:
    python 79_nesso1_yaml_generation.py
"""
import os

import pandas as pd
import yaml

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

STRUCTURE_SEQUENCES_CSV = os.path.join(ROOT, "output", "78_nesso1_prepare_inputs", "structure_sequences.csv")
COMPOUNDS_CSV = os.path.join(ROOT, "output", "71_boltz2_prepare_inputs", "compounds.csv")

OUTPUT_DIR = os.path.join(ROOT, "output", "79_nesso1_yaml_generation")
INPUT_YAMLS_DIR = os.path.join(OUTPUT_DIR, "input_yamls")
os.makedirs(INPUT_YAMLS_DIR, exist_ok=True)


def build_yaml(sequence, smiles):
    return {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": sequence}},
            {"ligand": {"id": "B", "smiles": smiles}},
        ],
        "properties": [
            {"affinity": {"binder": "B"}}
        ],
    }


def build_dimer_yaml(seq_a, seq_b, smiles):
    """Same schema as build_yaml(), but for the 7K98 pheS+pheT dimer: two protein chains
    (A=pheS, B=pheT), so the ligand takes id "C" instead of "B"."""
    return {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": seq_a}},
            {"protein": {"id": "B", "sequence": seq_b}},
            {"ligand": {"id": "C", "smiles": smiles}},
        ],
        "properties": [
            {"affinity": {"binder": "C"}}
        ],
    }


def main():
    structures = pd.read_csv(STRUCTURE_SEQUENCES_CSV)
    compounds = pd.read_csv(COMPOUNDS_CSV)
    print(f"Structures: {len(structures)}, compounds: {len(compounds)}, "
          f"target YAMLs: {len(structures) * len(compounds):,}")

    n_written, n_skipped = 0, 0
    for _, srow in structures.iterrows():
        structure_id = srow["structure_id"]
        structure_dir = os.path.join(INPUT_YAMLS_DIR, structure_id)
        os.makedirs(structure_dir, exist_ok=True)

        for _, crow in compounds.iterrows():
            out_path = os.path.join(structure_dir, f"{crow['compound_id']}.yaml")
            if os.path.isfile(out_path):
                n_skipped += 1
                continue

            if srow["is_dimer"]:
                yaml_dict = build_dimer_yaml(srow["sequence"], srow["sequence_b"], crow["smiles"])
            else:
                yaml_dict = build_yaml(srow["sequence"], crow["smiles"])
            with open(out_path, "w") as f:
                yaml.safe_dump(yaml_dict, f, sort_keys=False)
            n_written += 1

        print(f"  {structure_id}: done")

    print(f"\nYAMLs written this run: {n_written:,}")
    print(f"YAMLs already present (skipped): {n_skipped:,}")
    print(f"Total on disk: {n_written + n_skipped:,} (expected {len(structures) * len(compounds):,})")


if __name__ == "__main__":
    main()
