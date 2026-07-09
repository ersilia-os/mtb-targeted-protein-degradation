#!/usr/bin/env python3
"""
Prepares experimental multimeric PDB structures (from scripts/48_detect_pocket_multimers.py's
stripped_structures/) for Uni-Dock docking.

Mirrors scripts/21_unidock_proteinprep.py's approach for the main single-chain pipeline
(unidocktools proteinprep -> .pdbqt), with one addition: the stripped multimer structures have no
explicit hydrogens (script 48 intentionally skips PDB2PQR protonation/PyRosetta relax for pocket
detection, since P2Rank doesn't need them). Script 21's normal inputs, by contrast, are already
protonated+relaxed (via scripts/04_relax_structures.py) before proteinprep ever sees them. So this
script runs the same PDB2PQR protonation step (pH 7.0, AMBER force field) first, so proteinprep's
RDKit-based charge/H-bond typing has explicit hydrogens to work with, as intended.

Note on script 21's "ambertools and reduce" caution comment: depending on the installed
unidock_tools build, proteinprep may call `pdb4amber` (part of AmberTools) as a subprocess -
confirmed on nebula, resolved via `conda install -c conda-forge ambertools`. `reduce` was not
needed in practice there.

Docking itself runs on a separate machine; this script only produces the .pdbqt receptor files.

Requires two conda environments, both reached via SimpleConda regardless of which environment
this script itself is launched from:
  - "adda4tb"       for pdb2pqr (same as scripts/04_relax_structures.py)
  - "unidock_tools" for `unidocktools proteinprep`
    (pip install git+https://github.com/dptech-corp/Uni-Dock.git#subdirectory=unidock_tools)

Usage:
    python 49_unidock_proteinprep_multimers.py --pdb-codes 7K98,6XYZ
"""
import argparse
import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from utils.conda import SimpleConda

STRIPPED_DIR = os.path.join(ROOT, "output", "48_detect_pocket_multimers", "stripped_structures")
OUTPUT_ROOT = os.path.join(ROOT, "output", "49_unidock_proteinprep_multimers")

CONDA_ENV_PROTONATION = "adda4tb"     # unchanged from scripts/04_relax_structures.py
CONDA_ENV_UNIDOCK = "unidock_tools"   # unchanged from scripts/21_unidock_proteinprep.py


def calculate_protonation_states(input_file, output_file):
    """
    Assigns protonation states at pH 7.0 with PDB2PQR, run through the "adda4tb" conda
    environment. Verbatim from scripts/04_relax_structures.py.
    """
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
    """
    Runs unidocktools proteinprep, through the "unidock_tools" conda environment. Same tool/
    command as scripts/21_unidock_proteinprep.py, invoked via SimpleConda instead of a bare
    os.system call, since (unlike script 21) this script isn't assumed to already be launched
    from inside that environment.

    unidocktools' own "working_dir" (not exposed on the proteinprep CLI, defaults to ".") is
    where it creates its own scratch subfolders (e.g. receptor_reader/, receptor_grid/) - cd into
    output_file's directory first so those land next to the prepared structure instead of
    wherever this script happened to be launched from.
    """
    sc = SimpleConda()
    sc.run_commandlines(CONDA_ENV_UNIDOCK, [
        "cd {0}".format(os.path.dirname(os.path.abspath(output_file))),
        "unidocktools proteinprep -r {0} -o {1}".format(os.path.abspath(input_file), os.path.abspath(output_file)),
    ])


def prepare_structure(pdb_code):
    """
    Protonates the stripped structure at pH 7.0 (idempotent) then runs unidocktools proteinprep
    (idempotent) to produce a receptor .pdbqt ready for Uni-Dock. Returns True on success, False
    on any failure (skip-and-continue at the call site).
    """
    src = os.path.join(STRIPPED_DIR, f"{pdb_code}.pdb")
    if not os.path.isfile(src):
        print(f"Warning: no stripped structure found for {pdb_code} at {src}. Skipping.")
        return False

    pdb_dir = os.path.join(OUTPUT_ROOT, pdb_code)
    os.makedirs(pdb_dir, exist_ok=True)
    protonated_pdb = os.path.join(pdb_dir, f"{pdb_code}.pdb")
    pdbqt_file = os.path.join(pdb_dir, f"{pdb_code}.pdbqt")

    if os.path.isfile(pdbqt_file):
        print(f"Prepared structure for {pdb_code} already exists. Skipping.")
        return True

    try:
        if not os.path.isfile(protonated_pdb):
            print(f"Protonating {pdb_code}...")
            calculate_protonation_states(src, protonated_pdb)
        else:
            print(f"Protonated structure for {pdb_code} already exists. Skipping protonation.")

        print(f"Preparing structure {pdb_code} for docking...")
        run_proteinprep(protonated_pdb, pdbqt_file)
    except Exception as e:
        print(f"Warning: proteinprep failed for {pdb_code} ({e}). Skipping.")
        return False

    if not os.path.isfile(pdbqt_file):
        print(f"Warning: proteinprep produced no output for {pdb_code}. Skipping.")
        return False
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Protonate (PDB2PQR, pH 7.0) and prepare experimental multimeric structures "
                     "(from script 48's stripped_structures/) for Uni-Dock docking via "
                     "unidocktools proteinprep."
    )
    parser.add_argument("--pdb-codes", required=True, help="Comma-separated PDB codes, e.g. 7K98,6XYZ")
    args = parser.parse_args()
    pdb_codes = [c.strip().upper() for c in args.pdb_codes.split(",") if c.strip()]

    os.makedirs(OUTPUT_ROOT, exist_ok=True)

    for pdb_code in pdb_codes:
        print("----------------------------------------")
        prepare_structure(pdb_code)

    print("Done.")


if __name__ == "__main__":
    main()
