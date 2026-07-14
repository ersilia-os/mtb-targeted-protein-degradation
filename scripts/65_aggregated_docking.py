#!/usr/bin/env python3
### CAUTION: THIS SCRIPT NEEDS TO BE RUN
### WITHIN THE unidock_tools CONDA ENVIRONMENT
### IN A GPU MACHINE FOR MASSIVE PARALLELIZATION
"""
Dock script 64's 2,923 prepared aggregated compounds against each of the 12 curated pockets
(output/selected_pockets.csv, CAT + NON-CAT, dimer included) - one shared ligand set docked
against every pocket, a full cross-matrix (unlike script 60's per-pocket subset). Each pocket is
docked N_REPLICATES times (different seed each time) to measure Uni-Dock's own run-to-run variance.

Usage:
    conda activate unidock_tools
    python 65_aggregated_docking.py
"""
import csv
import glob
import os
import shutil
import subprocess
import sys
import tarfile

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from default import RANDOM_SEED
from utils.conda import SimpleConda

SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")
POCKET_DETECTION_DATA_CSV = os.path.join(ROOT, "output", "pocket_detection_data.csv")
MULTIMER_POCKET_DATA_CSV = os.path.join(ROOT, "output", "48_detect_pocket_multimers", "pocket_detection_data.csv")
LIGANDS_DIR = os.path.join(ROOT, "output", "64_aggregated_ligandprep", "conformations_prepared")

# Existing receptor sources to search (glob across all subfolders - a shared structure's receptor
# may live under a different pocket's own subfolder, e.g. lysS-CAT's receptor is under lysS-NONCAT's).
RECEPTOR_SEARCH_DIRS = [
    os.path.join(ROOT, "output", "60_unidock_docking_noncat_selective_10B", "docking_results"),
    os.path.join(ROOT, "output", "unidock_REAL_docking_2", "docking_results"),
]
MULTIMER_RECEPTOR_DIR = os.path.join(ROOT, "output", "49_unidock_proteinprep_multimers")
ALIGNED_RELAXED_STRUCTURES_DIR = os.path.join(ROOT, "output", "aligned_relaxed_structures")
MULTIMER_STRIPPED_DIR = os.path.join(ROOT, "output", "48_detect_pocket_multimers", "stripped_structures")

OUTPUT_DIR = os.path.join(ROOT, "output", "65_aggregated_docking")
DOCKING_RESULTS_DIR = os.path.join(OUTPUT_DIR, "docking_results")
os.makedirs(DOCKING_RESULTS_DIR, exist_ok=True)

DIMER_POCKET = "7K98_pocket_6"
CONDA_ENV_PROTONATION = "adda4tb"
N_REPLICATES = 5


def load_pockets():
    """[(gene_name, site_type, pocket_name), ...] for all 12 curated pockets."""
    df = pd.read_csv(SELECTED_POCKETS_CSV)
    return list(zip(df["gene_name"], df["site_type"], df["pocket_name"]))


def load_monomeric_pocket_info():
    """{pocket_name: (uniprot_ac, file_name, (cx, cy, cz))}, same lookup as scripts 46/60."""
    df = pd.read_csv(POCKET_DETECTION_DATA_CSV)
    pocket_names = df["File name"].str.replace(".pdb", "", regex=False) + "_pocket_" + df["Pocket number"].astype(str)
    info = {}
    for pocket_name, ac, file_name, centroid in zip(pocket_names, df["Uniprot AC"], df["File name"],
                                                      df["Pocket centroid coordinate (x y z)"]):
        cx, cy, cz = (float(v) for v in centroid.split())
        info[pocket_name] = (ac, file_name, (cx, cy, cz))
    return info


def load_dimer_centroid():
    """Same lookup as script 50."""
    df = pd.read_csv(MULTIMER_POCKET_DATA_CSV)
    row = df[(df["PDB code"] == "7K98") & (df["Pocket number"] == 6)].iloc[0]
    cx, cy, cz = (float(v) for v in row["Pocket centroid coordinate (x y z)"].split())
    return cx, cy, cz


def find_existing_receptor(structure_file):
    for base in RECEPTOR_SEARCH_DIRS:
        matches = glob.glob(os.path.join(base, "*", f"{structure_file}.pdbqt"))
        if matches:
            return matches[0]
    direct = os.path.join(MULTIMER_RECEPTOR_DIR, structure_file, f"{structure_file}.pdbqt")
    if os.path.isfile(direct):
        return direct
    return None


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


def resolve_receptor(pocket_name, structure_file, ac, outpath):
    """Reuse an existing prepared receptor if found anywhere; otherwise prepare fresh (monomeric:
    unidocktools proteinprep directly; dimer: protonate then proteinprep, per script 49) - a
    fallback that should not be needed since all 12 pockets' receptors are already in place."""
    dst_pdbqt = os.path.join(outpath, f"{structure_file}.pdbqt")
    if os.path.isfile(dst_pdbqt):
        print(f"  Receptor already prepared at {dst_pdbqt}, skipping.")
        return dst_pdbqt

    existing = find_existing_receptor(structure_file)
    if existing:
        print(f"  Reusing existing prepared receptor from {existing}")
        shutil.copyfile(existing, dst_pdbqt)
        return dst_pdbqt

    print(f"  No existing prepared receptor for {structure_file} - preparing fresh...")
    if pocket_name == DIMER_POCKET:
        stripped_pdb = os.path.join(MULTIMER_STRIPPED_DIR, f"{structure_file}.pdb")
        protonated_pdb = os.path.join(outpath, f"{structure_file}_protonated.pdb")
        calculate_protonation_states(stripped_pdb, protonated_pdb)
        src_pdb = protonated_pdb
    else:
        src_pdb = os.path.join(ALIGNED_RELAXED_STRUCTURES_DIR, ac, f"{structure_file}.pdb")

    subprocess.run(["unidocktools", "proteinprep", "-r", src_pdb, "-o", dst_pdbqt], check=True)
    return dst_pdbqt


def build_ligand_index(ligands_dir, index_path):
    files = sorted(os.path.join(ligands_dir, f) for f in os.listdir(ligands_dir) if f.endswith(".sdf"))
    with open(index_path, "w") as f:
        f.write("\n".join(files) + "\n")
    return len(files)


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
    search_mode: str = "fast",
    scoring: str = "vina",
    output_dir: str = "output",
    batch: int = 500,
    log_file: str = "log.txt",
):
    """Runs the Uni-Dock docking command. Verbatim from scripts 46/60."""
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


def extract_score_from_sdf(filepath: str):
    """Extracts the ENERGY value from the <Uni-Dock RESULT> block in an SDF file. Verbatim from
    scripts 46/60."""
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


def generate_report(directory: str, output_csv: str):
    """Verbatim from scripts 46/60."""
    sdf_files = [os.path.join(directory, f) for f in sorted(os.listdir(directory)) if f.endswith(".sdf")]
    results = [extract_score_from_sdf(file) for file in sdf_files]
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["compound", "score"])
        for compound, score in results:
            writer.writerow([compound.replace("_out.sdf", ""), score if score is not None else ""])


def average_replicates(results_dir, output_csv):
    """Average score across the N_REPLICATES results_r.csv files for a pocket: compound, score
    (mean), score_std - joined on compound (not row position)."""
    combined = None
    for r in range(1, N_REPLICATES + 1):
        df = pd.read_csv(os.path.join(results_dir, f"results_{r}.csv"))
        df[f"score_{r}"] = pd.to_numeric(df["score"], errors="coerce")
        df = df[["compound", f"score_{r}"]].set_index("compound")
        combined = df if combined is None else combined.join(df, how="outer")

    score_cols = [f"score_{r}" for r in range(1, N_REPLICATES + 1)]
    result = pd.DataFrame({
        "score": combined[score_cols].mean(axis=1).round(3),
        "score_std": combined[score_cols].std(axis=1).round(3),
    })
    result.index.name = "compound"
    result.reset_index().to_csv(output_csv, index=False)


def main():
    pockets = load_pockets()
    print(f"Curated pockets (CAT + NON-CAT, incl. dimer): {len(pockets)}")

    monomeric_info = load_monomeric_pocket_info()
    dimer_centroid = load_dimer_centroid()

    ligand_index_path = os.path.join(OUTPUT_DIR, "input_ligands.txt")
    n_ligands = build_ligand_index(LIGANDS_DIR, ligand_index_path)
    print(f"Ligand set: {n_ligands:,} prepared compounds (shared across all pockets)")

    for gene, site_type, pocket_name in pockets:
        print(f"\n\n--- {gene} ({site_type}): {pocket_name} ---\n\n")

        outpath = os.path.join(DOCKING_RESULTS_DIR, pocket_name)
        results_dir = os.path.join(outpath, "results")
        docking_archive_dir = os.path.join(outpath, "docking")
        logs_archive_dir = os.path.join(outpath, "logs")
        for d in (outpath, results_dir, docking_archive_dir, logs_archive_dir):
            os.makedirs(d, exist_ok=True)

        if pocket_name == DIMER_POCKET:
            structure_file, ac = "7K98", None
            cx, cy, cz = dimer_centroid
        else:
            ac, file_name, (cx, cy, cz) = monomeric_info[pocket_name]
            structure_file = file_name.replace(".pdb", "")

        receptor = resolve_receptor(pocket_name, structure_file, ac, outpath)

        for r in range(1, N_REPLICATES + 1):
            results_path = os.path.join(results_dir, f"results_{r}.csv")
            if os.path.isfile(results_path) and len(pd.read_csv(results_path)) >= n_ligands:
                print(f"  Replicate {r}: already complete ({results_path}), skipping.")
                continue

            print(f"  Replicate {r}/{N_REPLICATES} (seed={RANDOM_SEED + r - 1})...")
            tmp_docking_dir = os.path.join(outpath, "_docking_tmp")
            tmp_log_file = os.path.join(outpath, "_logs_tmp.log")
            os.makedirs(tmp_docking_dir, exist_ok=True)

            run_unidock(receptor=receptor, ligand_index=ligand_index_path,
                        center_x=cx, center_y=cy, center_z=cz, seed=RANDOM_SEED + r - 1,
                        output_dir=tmp_docking_dir, log_file=tmp_log_file)

            print("    Generating report...")
            generate_report(tmp_docking_dir, results_path)

            print("    Compressing results...")
            with tarfile.open(os.path.join(docking_archive_dir, f"docking_{r}.tar.gz"), "w:gz", compresslevel=9) as tar:
                tar.add(tmp_docking_dir, arcname="docking")
            with tarfile.open(os.path.join(logs_archive_dir, f"logs_{r}.tar.gz"), "w:gz") as tar:
                tar.add(tmp_log_file, arcname="logs.log")

            os.remove(tmp_log_file)
            shutil.rmtree(tmp_docking_dir)

        print("  Averaging across replicates...")
        average_replicates(results_dir, os.path.join(outpath, "results.csv"))


if __name__ == "__main__":
    main()
