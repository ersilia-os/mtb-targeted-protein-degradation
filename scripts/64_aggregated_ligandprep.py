#!/usr/bin/env python3
### CAUTION: THIS SCRIPT NEEDS TO BE RUN
### WITHIN THE unidock_tools CONDA ENVIRONMENT (on nebula)
### using the in-site unidocktools
"""
Prepare script 63's conformations for docking via unidocktools ligandprep (mirrors script 59/45).

Usage:
    conda activate unidock_tools
    python 64_aggregated_ligandprep.py
"""
import os
import subprocess

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

CONFORMATIONS_DIR = os.path.join(ROOT, "output", "63_aggregated_conformations", "conformations")
OUTPUT_DIR = os.path.join(ROOT, "output", "64_aggregated_ligandprep")
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
