#!/usr/bin/env python3
"""
### CAUTION: THIS SCRIPT NEEDS TO BE RUN
### WITHIN THE unidock_tools CONDA ENVIRONMENT
### ON A GPU MACHINE (the unidock binary requires an NVIDIA GPU)

Runs Uni-Dock docking for a single experimental multimer receptor/pocket (from
scripts/48_detect_pocket_multimers.py and scripts/49_unidock_proteinprep_multimers.py) against
the exact same ~99k Enamine REAL compound set used in scripts/44_generate_conformations.py
through scripts/46_unidock_REAL_2_docking.py (the BitBirch-clustered "Docking (III)" round) -
same box size (22.5 A), seed (42), search mode ("fast") and scoring function ("vina"), so scores
are directly comparable to output/unidock_REAL_docking_2/docking_results/.

Requires the already-prepared ligand set - either already extracted at
output/unidock_REAL_docking_2/conformations_prepared/, or as conformations_prepared.tar in the
same directory (produced by script 45's `unidocktools ligandprep`), which this script will
extract automatically the first time it's needed. This script does not regenerate conformers or
run ligandprep itself, it only builds a fresh ligand index from whatever's in that directory.

Usage:
    python 50_unidock_docking_multimers.py --pdb-code 7K98 --pocket-number 6
"""
import argparse
import csv
import os
import shutil
import subprocess
import tarfile

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
RECEPTOR_DIR = os.path.join(ROOT, "output", "49_unidock_proteinprep_multimers")
POCKET_DATA_CSV = os.path.join(ROOT, "output", "48_detect_pocket_multimers", "pocket_detection_data.csv")
LIGANDS_DIR = os.path.join(ROOT, "output", "unidock_REAL_docking_2", "conformations_prepared")
LIGANDS_TAR = os.path.join(ROOT, "output", "unidock_REAL_docking_2", "conformations_prepared.tar")
OUTPUT_ROOT = os.path.join(ROOT, "output", "50_unidock_docking_multimers")

REPORT_MIN_ROWS = 99105  # same completeness threshold as scripts/46_unidock_REAL_2_docking.py


def ensure_ligands_extracted():
    """
    Extracts conformations_prepared.tar into conformations_prepared/ if the directory doesn't
    already exist. Idempotent - a no-op on subsequent runs.
    """
    if os.path.isdir(LIGANDS_DIR):
        return
    if not os.path.isfile(LIGANDS_TAR):
        raise SystemExit(
            f"Ligand directory not found at {LIGANDS_DIR}, and no {LIGANDS_TAR} to extract it "
            f"from. Copy the prepared ~99k ligand set there first (see script docstring)."
        )
    print(f"Extracting {LIGANDS_TAR} -> {os.path.dirname(LIGANDS_DIR)} ...")
    with tarfile.open(LIGANDS_TAR) as tar:
        tar.extractall(os.path.dirname(LIGANDS_DIR), filter="data")


def build_ligand_index(ligands_dir, index_path):
    """
    Lists every prepared ligand file in ligands_dir into a plain-text index file (absolute paths,
    one per line). Rebuilt fresh each run rather than reusing any old input_ligands.txt, since an
    old index's paths were relative to a different machine's root and won't resolve here.
    """
    files = sorted(
        os.path.join(ligands_dir, f) for f in os.listdir(ligands_dir) if not f.startswith(".")
    )
    with open(index_path, "w") as f:
        f.write("\n".join(files) + "\n")
    return len(files)


def run_unidock(
    receptor,
    ligand_index,
    center_x,
    center_y,
    center_z,
    output_dir,
    log_file,
    size_x=22.5,
    size_y=22.5,
    size_z=22.5,
    seed=42,
    max_gpu_memory=22000,
    num_modes=1,
    verbosity=1,
    cpu=32,
    search_mode="fast",
    scoring="vina",
    batch=500,
):
    """
    Runs the Uni-Dock docking command. Same invocation and defaults as
    scripts/46_unidock_REAL_2_docking.py's run_unidock (box size 22.5 A, seed 42, max_gpu_memory
    22000, num_modes 1, cpu 32, search_mode "fast", scoring "vina", batch 500).
    """
    cmd = [
        "unidock",
        "--receptor", receptor,
        "--ligand_index", ligand_index,
        "--center_x", str(center_x),
        "--center_y", str(center_y),
        "--center_z", str(center_z),
        "--size_x", str(size_x),
        "--size_y", str(size_y),
        "--size_z", str(size_z),
        "--seed", str(seed),
        "--max_gpu_memory", str(max_gpu_memory),
        "--num_modes", str(num_modes),
        "--verbosity", str(verbosity),
        "--cpu", str(cpu),
        "--search_mode", search_mode,
        "--scoring", scoring,
        "--dir", output_dir,
        "--batch", str(batch),
    ]
    print(" ".join(cmd))
    with open(log_file, "w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)


def extract_score_from_sdf(filepath):
    """Extracts the ENERGY value from the <Uni-Dock RESULT> block in an SDF file. Verbatim from
    scripts/46_unidock_REAL_2_docking.py."""
    try:
        with open(filepath, "r") as file:
            for line in file:
                if line.strip() == "> <Uni-Dock RESULT>":
                    result_line = next(file).strip()
                    if "ENERGY=" in result_line:
                        return os.path.basename(filepath), result_line.split("ENERGY=")[1].split()[0]
    except Exception:
        pass
    return os.path.basename(filepath), None


def generate_report(directory, output_csv):
    """Processes all .sdf files in directory, extracts scores, writes compound,score CSV.
    Verbatim from scripts/46_unidock_REAL_2_docking.py."""
    sdf_files = [os.path.join(directory, f) for f in sorted(os.listdir(directory)) if f.endswith(".sdf")]
    results = [extract_score_from_sdf(file) for file in sdf_files]
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["compound", "score"])
        for compound, score in results:
            writer.writerow([compound.replace("_out.sdf", ""), score if score is not None else ""])


def resolve_receptor(pdb_code):
    receptor = os.path.join(RECEPTOR_DIR, pdb_code, f"{pdb_code}.pdbqt")
    if not os.path.isfile(receptor):
        raise SystemExit(
            f"No receptor found for {pdb_code} at {receptor}. "
            f"Run scripts/49_unidock_proteinprep_multimers.py --pdb-codes {pdb_code} first."
        )
    return receptor


def resolve_pocket_centroid(pdb_code, pocket_number):
    df = pd.read_csv(POCKET_DATA_CSV)
    row = df[(df["PDB code"] == pdb_code) & (df["Pocket number"] == pocket_number)]
    if len(row) == 0:
        raise SystemExit(
            f"No pocket {pocket_number} found for {pdb_code} in {POCKET_DATA_CSV}."
        )
    center_x, center_y, center_z = row.iloc[0]["Pocket centroid coordinate (x y z)"].split()
    return center_x, center_y, center_z


def main():
    parser = argparse.ArgumentParser(
        description="Run Uni-Dock docking for one multimer receptor/pocket against the same "
                     "~99k compound set used in scripts 44-46."
    )
    parser.add_argument("--pdb-code", required=True, help="PDB code, e.g. 7K98")
    parser.add_argument("--pocket-number", required=True, type=int, help="Pocket number, e.g. 6")
    args = parser.parse_args()
    pdb_code = args.pdb_code.strip().upper()
    pocket_number = args.pocket_number

    receptor = resolve_receptor(pdb_code)
    center_x, center_y, center_z = resolve_pocket_centroid(pdb_code, pocket_number)
    print(f"Receptor: {receptor}")
    print(f"Pocket {pocket_number} centroid: {center_x} {center_y} {center_z}")

    ensure_ligands_extracted()

    label = f"{pdb_code}_pocket_{pocket_number}"
    outpath = os.path.join(OUTPUT_ROOT, label)
    os.makedirs(outpath, exist_ok=True)

    report_csv = os.path.join(outpath, "report.csv")
    if os.path.isfile(report_csv) and len(pd.read_csv(report_csv)) >= REPORT_MIN_ROWS:
        print(f"{report_csv} already complete. Skipping.")
        return

    ligand_index = os.path.join(outpath, "input_ligands.txt")
    n = build_ligand_index(LIGANDS_DIR, ligand_index)
    print(f"{n} ligands indexed.")

    docking_dir = os.path.join(outpath, "docking")
    os.makedirs(docking_dir, exist_ok=True)
    log_file = os.path.join(outpath, "logs.log")

    print("Running docking!")
    run_unidock(
        receptor=receptor, ligand_index=ligand_index,
        center_x=center_x, center_y=center_y, center_z=center_z,
        output_dir=docking_dir, log_file=log_file,
    )

    print("Generating report!")
    generate_report(docking_dir, report_csv)

    print("Compressing results!")
    with tarfile.open(os.path.join(outpath, "docking.tar.gz"), "w:gz", compresslevel=9) as tar:
        tar.add(docking_dir, arcname="docking")
    with tarfile.open(os.path.join(outpath, "logs.tar.gz"), "w:gz") as tar:
        tar.add(log_file, arcname="logs.log")

    os.remove(log_file)
    shutil.rmtree(docking_dir)

    print("Done.")


if __name__ == "__main__":
    main()
