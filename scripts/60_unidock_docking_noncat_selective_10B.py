#!/usr/bin/env python3
### CAUTION: THIS SCRIPT NEEDS TO BE RUN
### WITHIN THE unidock_tools CONDA ENVIRONMENT (on nebula)
### IN A GPU MACHINE FOR MASSIVE PARALLELIZATION
"""
Dock script 59's prepared ligands against each of the 7 NON-CAT reference pockets independently -
per-compound subset docking (each pocket only docks the compounds it was itself selective for, per
scripts/57_NONCAT_selective_10B.py + scripts/58_generate_conformations_noncat_selective_10B.py's
merged_selective_hits.csv), NOT a full cross-matrix of every compound against every pocket.

run_unidock/extract_score_from_sdf/generate_report are reused verbatim from
scripts/46_unidock_REAL_2_docking.py (seed now sourced from src/default.RANDOM_SEED instead of a
hardcoded 42 - same value).

Receptor structures: 5 of the 7 pockets already have a prepared receptor .pdbqt sitting in
output/unidock_REAL_docking_2/docking_results/{pocket}/{structure}.pdbqt from the prior REAL-2
docking round - reused directly. The other 2 (as of this session: pheT's
alphafold2_P9WFU1_model_0_pocket_1 and alaS's swissmodel_P9WFW7_model_0_pocket_2) get a fresh
unidocktools proteinprep from the already-relaxed/protonated raw PDB in
output/aligned_relaxed_structures/{Uniprot AC}/{File name} (same approach as
scripts/21_unidock_proteinprep.py - no separate PDB2PQR step needed for these monomeric structures).

Docking box center is each pocket's centroid from output/pocket_detection_data.csv (same lookup used
in scripts/39_reduce_n_hits_I.py / 57_NONCAT_selective_10B.py). search_mode="fast" matches script
46's round-2 REAL docking exactly, keeping this round comparable to the existing REAL-1/REAL-2/DL
score distributions scripts/56_NONCAT_top100_selection.py already plots for these same 7 pockets.

Usage:
    conda activate unidock_tools
    python 60_unidock_docking_noncat_selective_10B.py
"""
import csv
import os
import shutil
import subprocess
import sys
import tarfile

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from default import RANDOM_SEED

SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")
POCKET_DETECTION_DATA_CSV = os.path.join(ROOT, "output", "pocket_detection_data.csv")
MERGED_SELECTIVE_HITS_CSV = os.path.join(ROOT, "output", "58_generate_conformations_noncat_selective_10B", "merged_selective_hits.csv")
CONFORMATIONS_PREPARED_DIR = os.path.join(ROOT, "output", "59_unidock_ligandprep_noncat_selective_10B", "conformations_prepared")
REAL2_DOCKING_RESULTS_DIR = os.path.join(ROOT, "output", "unidock_REAL_docking_2", "docking_results")
ALIGNED_RELAXED_STRUCTURES_DIR = os.path.join(ROOT, "output", "aligned_relaxed_structures")

OUTPUT_DIR = os.path.join(ROOT, "output", "60_unidock_docking_noncat_selective_10B")
DOCKING_RESULTS_DIR = os.path.join(OUTPUT_DIR, "docking_results")
os.makedirs(DOCKING_RESULTS_DIR, exist_ok=True)

DIMER_POCKET = "7K98_pocket_6"


def load_noncat_targets():
    """[(gene_name, pocket_name), ...] for the 7 NON-CAT pockets, excluding the dimer pocket."""
    df = pd.read_csv(SELECTED_POCKETS_CSV)
    noncat = df[df["site_type"] == "NON-CAT"]
    noncat = noncat[noncat["pocket_name"] != DIMER_POCKET]
    return list(zip(noncat["gene_name"], noncat["pocket_name"]))


def load_pocket_info():
    """{pocket_name: (uniprot_ac, file_name, (cx, cy, cz))}."""
    df = pd.read_csv(POCKET_DETECTION_DATA_CSV)
    pocket_names = df["File name"].str.replace(".pdb", "", regex=False) + "_pocket_" + df["Pocket number"].astype(str)
    info = {}
    for pocket_name, ac, file_name, centroid in zip(pocket_names, df["Uniprot AC"], df["File name"],
                                                      df["Pocket centroid coordinate (x y z)"]):
        cx, cy, cz = (float(v) for v in centroid.split())
        info[pocket_name] = (ac, file_name, (cx, cy, cz))
    return info


def build_pocket_ligand_index(pocket_name, merged_df, outpath):
    """Filter merged_selective_hits.csv to this pocket's own compounds, point at their prepared
    .sdf files. Returns (ligand_index_path, n_subset, n_found)."""
    subset = merged_df[merged_df["pocket_name"] == pocket_name]
    compound_ids = subset["compound_id"].tolist()
    ligand_index_path = os.path.join(outpath, f"ligand_index_{pocket_name}.txt")
    n_found = 0
    with open(ligand_index_path, "w") as f:
        for cid in compound_ids:
            sdf_path = os.path.join(CONFORMATIONS_PREPARED_DIR, f"{cid}.sdf")
            if os.path.isfile(sdf_path):
                f.write(sdf_path + "\n")
                n_found += 1
    return ligand_index_path, len(compound_ids), n_found


def resolve_receptor(pocket_name, ac, file_name, outpath):
    """Reuse the existing prepared receptor from the prior REAL-2 docking round if present,
    otherwise prepare a fresh one via unidocktools proteinprep. Idempotent: skips if already
    prepared into outpath by a previous run of this script."""
    dst_pdbqt = os.path.join(outpath, file_name.replace(".pdb", ".pdbqt"))
    if os.path.isfile(dst_pdbqt):
        print(f"  Receptor already prepared at {dst_pdbqt}, skipping.")
        return dst_pdbqt

    existing = os.path.join(REAL2_DOCKING_RESULTS_DIR, pocket_name, file_name.replace(".pdb", ".pdbqt"))
    if os.path.isfile(existing):
        print(f"  Reusing existing prepared receptor from {existing}")
        shutil.copyfile(existing, dst_pdbqt)
        return dst_pdbqt

    print(f"  No existing prepared receptor for {pocket_name} - preparing fresh via unidocktools proteinprep...")
    src_pdb = os.path.join(ALIGNED_RELAXED_STRUCTURES_DIR, ac, file_name)
    command = ["unidocktools", "proteinprep", "-r", src_pdb, "-o", dst_pdbqt]
    subprocess.run(command, check=True)
    return dst_pdbqt


def run_unidock(
    receptor: str,
    ligand_index: str,
    center_x: float,
    center_y: float,
    center_z: float,
    size_x: float = 22.5,
    size_y: float = 22.5,
    size_z: float = 22.5,
    seed: int = RANDOM_SEED,
    max_gpu_memory: int = 22000,
    num_modes: int = 1,
    verbosity: int = 1,
    cpu: int = 32,
    search_mode: str = "detail",
    scoring: str = "vina",
    output_dir: str = "output",
    batch: int = 500,
    log_file: str = "log.txt"
):
    """
    Runs the Uni-Dock docking command with the provided parameters. For further information, please
    check https://github.com/dptech-corp/Uni-Dock
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
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)


def extract_score_from_sdf(filepath: str):
    """Extracts the ENERGY value from the <Uni-Dock RESULT> block in an SDF file."""
    try:
        with open(filepath, 'r') as file:
            for line in file:
                if line.strip() == "> <Uni-Dock RESULT>":
                    result_line = next(file).strip()
                    if "ENERGY=" in result_line:
                        return os.path.basename(filepath), result_line.split("ENERGY=")[1].split()[0]
    except Exception:
        pass
    return os.path.basename(filepath), None


def generate_report(directory: str, output_csv: str = "scores.csv"):
    """Processes all .sdf files in a directory, extracts score values, and writes the results as a
    CSV: filename,score."""
    sdf_files = [os.path.join(directory, f) for f in sorted(os.listdir(directory)) if f.endswith(".sdf")]
    results = [extract_score_from_sdf(file) for file in sdf_files]

    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["compound", "score"])
        for compound, score in results:
            writer.writerow([compound.replace("_out.sdf", ""), score if score is not None else ""])


def main():
    targets = load_noncat_targets()
    print(f"NON-CAT targets (excl. dimer pocket {DIMER_POCKET}): {len(targets)}")

    pocket_info = load_pocket_info()
    merged = pd.read_csv(MERGED_SELECTIVE_HITS_CSV)

    for gene, pocket_name in targets:
        print(f"\n\n--- {gene}: {pocket_name} ---\n\n")
        ac, file_name, (cx, cy, cz) = pocket_info[pocket_name]

        outpath = os.path.join(DOCKING_RESULTS_DIR, pocket_name)
        os.makedirs(outpath, exist_ok=True)
        report_path = os.path.join(outpath, "report.csv")

        if os.path.isfile(report_path):
            print(f"  report.csv already exists at {report_path}, skipping.")
            continue

        ligand_index_path, n_subset, n_found = build_pocket_ligand_index(pocket_name, merged, outpath)
        print(f"  Ligand subset: {n_found:,} / {n_subset:,} compounds found prepared")

        receptor = resolve_receptor(pocket_name, ac, file_name, outpath)

        docking_dir = os.path.join(outpath, "docking")
        log_file = os.path.join(outpath, "logs.log")
        os.makedirs(docking_dir, exist_ok=True)

        run_unidock(receptor=receptor, ligand_index=ligand_index_path,
                    center_x=cx, center_y=cy, center_z=cz,
                    search_mode="fast", seed=RANDOM_SEED, cpu=32,
                    output_dir=docking_dir, log_file=log_file)

        print("  Generating report...")
        generate_report(docking_dir, report_path)

        print("  Compressing results...")
        with tarfile.open(os.path.join(outpath, "docking.tar.gz"), "w:gz", compresslevel=9) as tar:
            tar.add(docking_dir, arcname="docking")
        with tarfile.open(os.path.join(outpath, "logs.tar.gz"), "w:gz") as tar:
            tar.add(log_file, arcname="logs.log")

        os.remove(log_file)
        shutil.rmtree(docking_dir)


if __name__ == "__main__":
    main()
