#!/usr/bin/env python3
"""
Prepares Uni-Dock receptors (.pdbqt) for script 90's AlphaFold DB monomers -- 38 human aaRS genes
by default, or all 21 Mtb genes via --organism mtb. Mirrors script 49's approach for the unrelaxed
Mtb multimer structures: these are raw AlphaFold DB downloads with no explicit hydrogens (no
relaxation was done for this sub-pipeline, unlike script 04 for the main Mtb structures), so
PDB2PQR protonation (pH 7.0, AMBER force field) runs first, then `unidocktools proteinprep`'s
RDKit-based charge/H-bond typing has hydrogens to work with.

Requires two conda environments, both reached via SimpleConda regardless of which environment this
script itself is launched from:
  - "adda4tb"       for pdb2pqr (same as script 49)
  - "unidock_tools" for `unidocktools proteinprep`

Usage:
    python 94_human_receptor_prep.py [--organism human|mtb] [--genes AARS1,KARS1]
"""
import argparse
import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from utils.conda import SimpleConda  # noqa: E402

import pandas as pd  # noqa: E402

CONDA_ENV_PROTONATION = "adda4tb"
CONDA_ENV_UNIDOCK = "unidock_tools"


def calculate_protonation_states(input_file, output_file):
    """Assigns protonation states at pH 7.0 with PDB2PQR. Verbatim from script 49."""
    cwd = os.getcwd()
    file_name = os.path.basename(input_file).split(".")[0]
    input_file = os.path.abspath(input_file)
    output_file = os.path.abspath(output_file)
    sc = SimpleConda()
    sc.run_commandlines(CONDA_ENV_PROTONATION, [
        "pdb2pqr --ff=AMBER --with-ph=7.0 --pdb-output {0} {1} {2}.pqr".format(output_file, input_file, file_name)
    ])
    os.remove(os.path.join(cwd, file_name + ".pqr"))
    os.remove(os.path.join(cwd, file_name + ".log"))


def run_proteinprep(input_file, output_file):
    """Runs unidocktools proteinprep through the "unidock_tools" conda environment. Verbatim
    approach from script 49 (cd into output_file's directory first so proteinprep's own scratch
    subfolders land next to the prepared structure)."""
    sc = SimpleConda()
    sc.run_commandlines(CONDA_ENV_UNIDOCK, [
        "cd {0}".format(os.path.dirname(os.path.abspath(output_file))),
        "unidocktools proteinprep -r {0} -o {1}".format(os.path.abspath(input_file), os.path.abspath(output_file)),
    ])


def fix_valence_problems(input_file, output_file):
    """Removes PDB2PQR-mislabeled hydrogens that break RDKit's valence sanitization (rare edge
    case found on 4/38 large human aaRS structures -- see src/utils/fix_pdb_valence.py). Always
    run after protonation, before proteinprep: a safe no-op copy when a structure has no valence
    problems, so this doesn't need to be conditional on which genes are affected."""
    sc = SimpleConda()
    sc.run_commandlines(CONDA_ENV_UNIDOCK, [
        "python {0} --input {1} --output {2}".format(
            os.path.join(ROOT, "src", "utils", "fix_pdb_valence.py"),
            os.path.abspath(input_file), os.path.abspath(output_file)),
    ])


def prepare_receptor(uniprot_ac, src_pdb, output_root):
    """Protonates (idempotent), fixes any valence problems (idempotent, safe no-op if none), then
    proteinpreps (idempotent) src_pdb into a receptor .pdbqt. Returns True on success, False on
    any failure (skip-and-continue at the call site)."""
    gene_dir = os.path.join(output_root, uniprot_ac)
    os.makedirs(gene_dir, exist_ok=True)
    protonated_pdb = os.path.join(gene_dir, f"{uniprot_ac}.pdb")
    fixed_pdb = os.path.join(gene_dir, f"{uniprot_ac}_fixed.pdb")
    pdbqt_file = os.path.join(gene_dir, f"{uniprot_ac}.pdbqt")

    if os.path.isfile(pdbqt_file):
        print(f"  Prepared receptor for {uniprot_ac} already exists. Skipping.")
        return True

    try:
        if not os.path.isfile(protonated_pdb):
            print(f"  Protonating {uniprot_ac}...")
            calculate_protonation_states(src_pdb, protonated_pdb)
        else:
            print(f"  Protonated structure for {uniprot_ac} already exists. Skipping protonation.")

        print(f"  Checking/fixing valence problems for {uniprot_ac}...")
        fix_valence_problems(protonated_pdb, fixed_pdb)

        print(f"  Preparing {uniprot_ac} for docking...")
        run_proteinprep(fixed_pdb, pdbqt_file)
    except Exception as e:
        print(f"  WARNING: proteinprep failed for {uniprot_ac} ({e}). Skipping.")
        return False

    if not os.path.isfile(pdbqt_file):
        print(f"  WARNING: proteinprep produced no output for {uniprot_ac}. Skipping.")
        return False
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--organism", choices=["human", "mtb"], default="human",
                         help="Gene set to prepare receptors for (default: human)")
    parser.add_argument("--genes", default=None,
                         help="Comma-separated gene names (default: all in structures_data.csv)")
    args = parser.parse_args()

    structures_data_csv = os.path.join(ROOT, "output", f"90_{args.organism}_download_alphafold", "structures_data.csv")
    output_root = os.path.join(ROOT, "output", f"94_{args.organism}_receptor_prep")
    os.makedirs(output_root, exist_ok=True)

    structures = pd.read_csv(structures_data_csv)
    structures = structures[structures["status"] == "ok"]
    if args.genes:
        wanted = set(args.genes.split(","))
        structures = structures[structures["gene_name"].isin(wanted)]

    n_ok, n_fail = 0, 0
    for _, row in structures.iterrows():
        gene_name, uniprot_ac = row["gene_name"], row["uniprot_ac"]
        src_pdb = os.path.join(ROOT, row["file_path"])
        print(f"--- {gene_name} ({uniprot_ac}) ---")
        if prepare_receptor(uniprot_ac, src_pdb, output_root):
            n_ok += 1
        else:
            n_fail += 1

    print(f"\n{n_ok}/{n_ok + n_fail} receptors prepared successfully.")


if __name__ == "__main__":
    main()
