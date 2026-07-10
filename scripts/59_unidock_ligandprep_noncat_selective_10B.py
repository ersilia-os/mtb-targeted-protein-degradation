#!/usr/bin/env python3
### CAUTION: THIS SCRIPT NEEDS TO BE RUN
### WITHIN THE unidock_tools CONDA ENVIRONMENT (on nebula)
### using the in-site unidocktools
"""
Prepare script 58's deduplicated conformations
(processed/58_generate_conformations_noncat_selective_10B/conformations/) for Uni-Dock docking, via
unidocktools ligandprep. Mirrors scripts/45_unidock_REAL_2_ligandprep.py /
scripts/XX_unidock_ligandprep_noncat_top10k.py.

Ligand preparation is a property of the molecule alone, not the pocket, so this runs once over ALL
conformations regardless of which of the 7 NON-CAT pockets they'll later be docked against in
scripts/60_unidock_docking_noncat_selective_10B.py - that script builds its own, separate, per-pocket
ligand index files from this script's prepared output; this script's input_ligands.txt is global.

Usage:
    conda activate unidock_tools
    python 59_unidock_ligandprep_noncat_selective_10B.py
"""
import os
import subprocess

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

CONFORMATIONS_DIR = os.path.join(ROOT, "processed", "58_generate_conformations_noncat_selective_10B", "conformations")
OUTPUT_DIR = os.path.join(ROOT, "processed", "59_unidock_ligandprep_noncat_selective_10B")
CONFORMATIONS_PREPARED_DIR = os.path.join(OUTPUT_DIR, "conformations_prepared")
os.makedirs(CONFORMATIONS_PREPARED_DIR, exist_ok=True)

INPUT_LIGANDS_PATH = os.path.join(OUTPUT_DIR, "input_ligands.txt")


def write_input_ligands(source_dir):
    sdfs = sorted(f for f in os.listdir(source_dir) if f.endswith(".sdf"))
    with open(INPUT_LIGANDS_PATH, "w") as outfile:
        for sdf in sdfs:
            outfile.write(os.path.join(source_dir, sdf))
            outfile.write("\n")
    return len(sdfs)


def main():
    print("Creating file with all ligand paths...")
    n_conformations = write_input_ligands(CONFORMATIONS_DIR)
    print(f"  {n_conformations:,} conformations")

    print("Preparing compounds for docking using unidocktools...")
    command = [
        "unidocktools", "ligandprep",
        "--ligand_index", INPUT_LIGANDS_PATH,
        "--savedir", CONFORMATIONS_PREPARED_DIR,
        "--batch_size", "100",
        "--use_file_name",
    ]
    subprocess.run(command, check=True)

    print("Creating file with all ligand paths again, discarding those that failed...")
    n_prepared = write_input_ligands(CONFORMATIONS_PREPARED_DIR)
    print(f"  {n_prepared:,} / {n_conformations:,} prepared successfully")


if __name__ == "__main__":
    main()
