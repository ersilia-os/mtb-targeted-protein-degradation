#!/usr/bin/env python3
"""
Generates Nesso-1 YAML input files from script 78's protein_sequences.csv and script 71's
compounds.csv -- one YAML per (gene, compound) pair, 5 genes x 1,095 compounds = 5,475 total
(revised from the earlier 11-structure design, see script 78's docstring for why).

Nesso-1's YAML schema (docs/prediction.md in recursionpharma/nesso) has no `msa:` key (no MSA
required at all) and no `constraints:` key (no pocket-conditioning), so each YAML is just a
single protein chain + a ligand + an affinity property.

Directory layout matters: `nesso predict <dir>` is NOT recursive (only top-level .yaml/.yml
files are picked up, confirmed by reading nesso/main.py's check_inputs()), so each gene's 1,095
compound YAMLs live in their own flat subdirectory. The YAML filename stem becomes the output
record_id (confirmed in the repo's AGENTS.md and main.py's _process_single_yaml), so files are
named <compound_id>.yaml.

Running Nesso-1 itself is a separate, later script (80).

Usage:
    python 79_nesso1_yaml_generation.py
"""
import os

import pandas as pd
import yaml

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

PROTEIN_SEQUENCES_CSV = os.path.join(ROOT, "output", "78_nesso1_prepare_inputs", "protein_sequences.csv")
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


def main():
    proteins = pd.read_csv(PROTEIN_SEQUENCES_CSV)
    compounds = pd.read_csv(COMPOUNDS_CSV)
    print(f"Genes: {len(proteins)}, compounds: {len(compounds)}, "
          f"target YAMLs: {len(proteins) * len(compounds):,}")

    n_written, n_skipped = 0, 0
    for _, prow in proteins.iterrows():
        gene_name = prow["gene_name"]
        gene_dir = os.path.join(INPUT_YAMLS_DIR, gene_name)
        os.makedirs(gene_dir, exist_ok=True)

        for _, crow in compounds.iterrows():
            out_path = os.path.join(gene_dir, f"{crow['compound_id']}.yaml")
            if os.path.isfile(out_path):
                n_skipped += 1
                continue

            yaml_dict = build_yaml(prow["sequence"], crow["smiles"])
            with open(out_path, "w") as f:
                yaml.safe_dump(yaml_dict, f, sort_keys=False)
            n_written += 1

        print(f"  {gene_name}: done")

    print(f"\nYAMLs written this run: {n_written:,}")
    print(f"YAMLs already present (skipped): {n_skipped:,}")
    print(f"Total on disk: {n_written + n_skipped:,} (expected {len(proteins) * len(compounds):,})")


if __name__ == "__main__":
    main()
