#!/usr/bin/env python3
"""
Generates the actual Boltz-2 YAML input files from script 71's two CSVs -- one YAML per
(pocket, compound) pair, 12,045 total (11 pockets x 1,095 compounds).

Each YAML's `msa:` field is baked in as an ABSOLUTE path on nebula (NEBULA_REPO_ROOT below), not a
relative one: Boltz-2 resolves a YAML's `msa` value via plain `Path(msa_id).exists()`, relative to
the process's current working directory at invocation time, not relative to the YAML file's own
location -- a relative path would silently depend on exactly where `boltz predict` happens to be
invoked from. The referenced MSA cache file does not exist yet at generation time (or at rsync
time) -- it is created later by script 73's MSA-bootstrap step, once per pocket, on nebula.

Running Boltz-2 itself (MSA bootstrapping, `boltz predict`, result aggregation) is a separate,
later script (73).

Usage:
    python 72_boltz2_yaml_generation.py
"""
import os

import pandas as pd
import yaml

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

POCKET_SEQUENCES_CSV = os.path.join(ROOT, "output", "71_boltz2_prepare_inputs", "pocket_sequences.csv")
COMPOUNDS_CSV = os.path.join(ROOT, "output", "71_boltz2_prepare_inputs", "compounds.csv")

OUTPUT_DIR = os.path.join(ROOT, "output", "72_boltz2_yaml_generation")
INPUT_YAMLS_DIR = os.path.join(OUTPUT_DIR, "input_yamls")
os.makedirs(INPUT_YAMLS_DIR, exist_ok=True)

# nebula's actual repo clone path (confirmed via `ssh nebula`) -- NOT the divergent non-git copy
# also present there under ~/Documents/acomajuncosa/. The msa_cache file this points to is created
# by script 73 on nebula, not by this script.
NEBULA_REPO_ROOT = "/home/admin/mtb-targeted-protein-degradation"


def build_yaml(sequence, msa_path, smiles, contacts):
    return {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": sequence, "msa": msa_path}},
            {"ligand": {"id": "B", "smiles": smiles}},
        ],
        "constraints": [
            {"pocket": {"binder": "B", "contacts": [["A", pos] for pos in contacts]}}
        ],
        "properties": [
            {"affinity": {"binder": "B"}}
        ],
    }


def main():
    pockets = pd.read_csv(POCKET_SEQUENCES_CSV)
    compounds = pd.read_csv(COMPOUNDS_CSV)
    print(f"Pockets: {len(pockets)}, compounds: {len(compounds)}, "
          f"target YAMLs: {len(pockets) * len(compounds):,}")

    n_written, n_skipped = 0, 0
    for _, prow in pockets.iterrows():
        pocket_name = prow["pocket_name"]
        sequence = prow["sequence"]
        contacts = [int(p) for p in prow["pocket_contacts"].split()]
        msa_path = f"{NEBULA_REPO_ROOT}/output/73_boltz2_docking/msa_cache/{pocket_name}.csv"

        pocket_dir = os.path.join(INPUT_YAMLS_DIR, pocket_name)
        os.makedirs(pocket_dir, exist_ok=True)

        for _, crow in compounds.iterrows():
            out_path = os.path.join(pocket_dir, f"{crow['compound_id']}.yaml")
            if os.path.isfile(out_path):
                n_skipped += 1
                continue

            yaml_dict = build_yaml(sequence, msa_path, crow["smiles"], contacts)
            with open(out_path, "w") as f:
                yaml.safe_dump(yaml_dict, f, sort_keys=False)
            n_written += 1

        print(f"  {pocket_name}: done")

    print(f"\nYAMLs written this run: {n_written:,}")
    print(f"YAMLs already present (skipped): {n_skipped:,}")
    print(f"Total on disk: {n_written + n_skipped:,} (expected {len(pockets) * len(compounds):,})")


if __name__ == "__main__":
    main()
