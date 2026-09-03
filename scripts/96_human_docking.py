#!/usr/bin/env python3
### CAUTION: THIS SCRIPT NEEDS TO BE RUN
### WITHIN THE unidock_tools CONDA ENVIRONMENT
### IN A GPU JOB ON THE IRB CLUSTER (sbnb_gpu_3090 / sbnb_gpu_h200 -- see script 96_human_run_array.sh)
"""
Docks the 1,095 filtered Mtb hits against each of the 389 human pockets (all 38 genes) -- a full
per-pocket run, single shot (no replicates, confirmed with the user: 389 x 1,095 ~= 425,955
endpoints is already a lot of GPU time on cluster nodes we don't have unlimited access to).
Reuses `run_unidock`/`extract_score_from_sdf`/`generate_report` verbatim from scripts 46/60/65
(box 22.5A, seed 42, search_mode=fast, scoring=vina, cpu=32).

Ligands: script 64's already-prepared SDFs (output/64_aggregated_ligandprep/conformations_prepared/),
filtered down to the 1,095 filtered_hits.csv compound IDs -- confirmed directly (not assumed) that
all 1,095 are present there, so no separate ligand-prep script was needed.

Receptors: script 94's per-gene .pdbqt files. Box center per pocket comes from script 91's
pocket_detection_data.csv.

Unlike scripts 80/86 (Nesso-1), there's no --no-aggregate/--aggregate-only split here: each pocket
writes its own independent report.csv (no shared file for concurrent SLURM array tasks to race
on), so plain per-pocket resumability is enough. Script 97 does the cross-pocket merge separately,
once, after the array finishes.

The `unidock` binary is resolved relative to the running interpreter (sys.executable), not looked
up on PATH -- same convention as script 80's NESSO_BIN, for the same reason: this script is always
invoked as `envs/unidock_tools/bin/python scripts/96_human_docking.py` (the SLURM array and the
smoke test alike), never via `conda activate`, so a bare "unidock" is not found on PATH -- its
sibling binary in the same env's bin/ directory is (confirmed the hard way: the first smoke test
run hit FileNotFoundError: 'unidock' before this fix).

Usage:
    # Smoke test: one gene, few compounds, isolated output.
    python 96_human_docking.py --genes AARS1 --max-compounds 5 --out-subdir smoke_test

    # One SLURM array task's worth of work (called from 96_human_run_array.sh).
    python 96_human_docking.py --genes <gene_name> --out-subdir docking_results
"""
import argparse
import csv
import os
import sys
import tarfile

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from default import RANDOM_SEED  # noqa: E402

import subprocess  # noqa: E402

# Resolved relative to sys.executable, not looked up on PATH -- see module docstring for why.
UNIDOCK_BIN = os.path.join(os.path.dirname(sys.executable), "unidock")

POCKET_DETECTION_DATA_CSV = os.path.join(ROOT, "output", "91_human_detect_pockets", "pocket_detection_data.csv")
FILTERED_HITS_CSV = os.path.join(ROOT, "output", "70_filtering", "filtered_hits.csv")
LIGANDS_DIR = os.path.join(ROOT, "output", "64_aggregated_ligandprep", "conformations_prepared")
RECEPTOR_DIR = os.path.join(ROOT, "output", "94_human_receptor_prep")

OUTPUT_DIR = os.path.join(ROOT, "output", "96_human_docking")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_pockets(genes=None):
    """{uniprot_ac: {"gene_name": ..., "pockets": [(pocket_number, file_name, (cx, cy, cz)), ...]}}"""
    df = pd.read_csv(POCKET_DETECTION_DATA_CSV)
    if genes:
        df = df[df["Gene name"].isin(genes)]

    by_gene = {}
    for _, row in df.iterrows():
        ac = row["Uniprot AC"]
        cx, cy, cz = (float(v) for v in row["Pocket centroid coordinate (x y z)"].split())
        entry = by_gene.setdefault(ac, {"gene_name": row["Gene name"], "pockets": []})
        entry["pockets"].append((int(row["Pocket number"]), row["File name"], (cx, cy, cz)))
    return by_gene


def build_ligand_index(index_path, max_compounds=None):
    """Writes an index of script 64's prepared SDF paths, restricted to filtered_hits.csv's 1,095
    compound IDs. Resumable: skipped if index_path already exists (shared across all pockets)."""
    if os.path.isfile(index_path):
        with open(index_path) as f:
            return len(f.readlines())

    compound_ids = pd.read_csv(FILTERED_HITS_CSV)["compound_id"].tolist()
    if max_compounds is not None:
        compound_ids = sorted(compound_ids)[:max_compounds]

    paths = []
    missing = []
    for cid in compound_ids:
        sdf_path = os.path.join(LIGANDS_DIR, f"{cid}.sdf")
        if os.path.isfile(sdf_path):
            paths.append(sdf_path)
        else:
            missing.append(cid)
    if missing:
        print(f"  WARNING: {len(missing)} compound(s) missing from {LIGANDS_DIR}: {missing[:5]}...")

    with open(index_path, "w") as f:
        f.write("\n".join(paths) + "\n")
    return len(paths)


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
    """Runs the Uni-Dock docking command. Adapted from scripts 46/60/65 (those assumed "unidock"
    was already on PATH via `conda activate`; here it's resolved relative to sys.executable, see
    module docstring)."""
    cmd = [
        UNIDOCK_BIN,
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
    scripts 46/60/65."""
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
    """Verbatim from scripts 46/60/65."""
    sdf_files = [os.path.join(directory, f) for f in sorted(os.listdir(directory)) if f.endswith(".sdf")]
    results = [extract_score_from_sdf(file) for file in sdf_files]
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["compound", "score"])
        for compound, score in results:
            writer.writerow([compound.replace("_out.sdf", ""), score if score is not None else ""])


def dock_pocket(uniprot_ac, receptor, pocket_number, file_name, centroid, ligand_index, n_ligands, out_dir):
    outpath = os.path.join(out_dir, uniprot_ac, str(pocket_number))
    os.makedirs(outpath, exist_ok=True)
    report_path = os.path.join(outpath, "report.csv")

    if os.path.isfile(report_path) and len(pd.read_csv(report_path)) >= n_ligands:
        print(f"    Pocket {pocket_number}: already complete, skipping.")
        return

    tmp_docking_dir = os.path.join(outpath, "_docking_tmp")
    tmp_log_file = os.path.join(outpath, "_log_tmp.log")
    os.makedirs(tmp_docking_dir, exist_ok=True)

    cx, cy, cz = centroid
    run_unidock(receptor=receptor, ligand_index=ligand_index, center_x=cx, center_y=cy, center_z=cz,
                output_dir=tmp_docking_dir, log_file=tmp_log_file)

    generate_report(tmp_docking_dir, report_path)

    with tarfile.open(os.path.join(outpath, "docking.tar.gz"), "w:gz", compresslevel=9) as tar:
        tar.add(tmp_docking_dir, arcname="docking")
    with tarfile.open(os.path.join(outpath, "logs.tar.gz"), "w:gz") as tar:
        tar.add(tmp_log_file, arcname="logs.log")

    os.remove(tmp_log_file)
    for f in os.listdir(tmp_docking_dir):
        os.remove(os.path.join(tmp_docking_dir, f))
    os.rmdir(tmp_docking_dir)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genes", default=None, help="Comma-separated gene names (default: all 38)")
    parser.add_argument("--max-compounds", type=int, default=None,
                         help="Limit to the first N compounds, sorted by compound_id (default: all 1,095)")
    parser.add_argument("--out-subdir", default="docking_results",
                         help="Output subdirectory under output/96_human_docking/ (default: docking_results)")
    args = parser.parse_args()

    genes = set(args.genes.split(",")) if args.genes else None
    out_dir = os.path.join(OUTPUT_DIR, args.out_subdir)
    os.makedirs(out_dir, exist_ok=True)

    pockets_by_gene = load_pockets(genes)
    print(f"Genes: {len(pockets_by_gene)}, total pockets: {sum(len(v['pockets']) for v in pockets_by_gene.values())}")

    ligand_index_path = os.path.join(OUTPUT_DIR, "input_ligands.txt" if args.max_compounds is None
                                      else f"input_ligands_max{args.max_compounds}.txt")
    n_ligands = build_ligand_index(ligand_index_path, args.max_compounds)
    print(f"Ligand set: {n_ligands:,} prepared compounds (shared across all pockets)")

    for ac, info in pockets_by_gene.items():
        gene_name = info["gene_name"]
        receptor = os.path.join(RECEPTOR_DIR, ac, f"{ac}.pdbqt")
        if not os.path.isfile(receptor):
            print(f"--- {gene_name} ({ac}): WARNING no receptor at {receptor}, skipping gene ---")
            continue

        print(f"\n--- {gene_name} ({ac}): {len(info['pockets'])} pocket(s) ---")
        for pocket_number, file_name, centroid in sorted(info["pockets"]):
            print(f"  Pocket {pocket_number}...")
            dock_pocket(ac, receptor, pocket_number, file_name, centroid, ligand_index_path, n_ligands, out_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
