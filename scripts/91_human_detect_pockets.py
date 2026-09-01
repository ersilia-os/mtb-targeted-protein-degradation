#!/usr/bin/env python3
"""
Runs P2Rank on script 90's 38 human AlphaFold DB structures and extracts pocket data. Direct port
of script 08's logic, simplified for a single structure source per protein (always AlphaFold2-type
confidence, so no prediction_type branching, and no separate aligned/relaxed vs. original
directory split -- one raw AFDB structure per gene is used throughout).

Keeps every pocket P2Rank reports (no top-K/probability/pLDDT filtering, for the time being --
catalytic-pocket selection is a later step). Confidence is reported as the minimum pLDDT among
residues within RADIUS_A of each pocket's own centroid -- a simple, reproducible geometric
neighborhood, independent of P2Rank's own (SAS-point-based) residue_ids list.

Usage:
    python 91_human_detect_pockets.py
"""
import math
import os
import subprocess

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

STRUCTURES_DATA_CSV = os.path.join(ROOT, "output", "90_human_download_alphafold", "structures_data.csv")
OUTPUT_DIR = os.path.join(ROOT, "output", "91_human_detect_pockets")
DETECTED_POCKETS_DIR = os.path.join(OUTPUT_DIR, "detected_pockets")
os.makedirs(DETECTED_POCKETS_DIR, exist_ok=True)

RADIUS_A = 8.0  # neighborhood radius (Angstrom) around each pocket centroid for the pLDDT check


def detect_pockets(pdb_file, output_dir):
    """Runs P2Rank on pdb_file, saving results to output_dir."""
    path_to_p2rank = "/aloy/home/acomajuncosa/programs/p2rank_2.5/prank"  # Change as needed
    command = [path_to_p2rank, "predict", "-f", pdb_file, "-c", "alphafold", "-o", output_dir,
               "-visualizations", "0"]
    try:
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"Error detecting pockets in {pdb_file}: {e}")


def extract_pocket_centers(csv_file):
    """Pocket rank -> (x, y, z) centroid, from a P2Rank _predictions.csv."""
    df = pd.read_csv(csv_file)
    df.columns = df.columns.str.strip()
    return {int(row["rank"]): (row["center_x"], row["center_y"], row["center_z"])
            for _, row in df.iterrows()}


def extract_pocket_scores(csv_file):
    df = pd.read_csv(csv_file)
    df.columns = df.columns.str.strip()
    return {int(row["rank"]): row["score"] for _, row in df.iterrows()}


def extract_pocket_probabilities(csv_file):
    df = pd.read_csv(csv_file)
    df.columns = df.columns.str.strip()
    return {int(row["rank"]): row["probability"] for _, row in df.iterrows()}


def extract_pocket_residues(csv_file):
    df = pd.read_csv(csv_file)
    df.columns = df.columns.str.strip()
    return {int(row["rank"]): row["residue_ids"].split() for _, row in df.iterrows()}


def parse_residue_geometry(pdb_file):
    """Residue number -> (x, y, z, pLDDT), taken from each residue's CA atom. AlphaFold stores the
    same pLDDT (B-factor) across all of a residue's atoms, so CA's value represents the residue."""
    res_geom = {}
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                resnum = int(line[22:26])
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                bfactor = float(line[60:66])
                res_geom[resnum] = (x, y, z, bfactor)
    return res_geom


def residues_within_radius(res_geom, center, radius=RADIUS_A):
    """[(residue_number, pLDDT), ...] for every residue whose CA falls within radius of center."""
    cx, cy, cz = center
    hits = []
    for resnum, (x, y, z, bfactor) in res_geom.items():
        dist = math.sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2)
        if dist <= radius:
            hits.append((resnum, bfactor))
    return hits


def write_pocket_pdbs(pocket_dict, output_dir):
    """One centroid-marker PDB per pocket."""
    os.makedirs(output_dir, exist_ok=True)
    for pocket_number, (x, y, z) in pocket_dict.items():
        pdb_filename = os.path.join(output_dir, f"pocket_{pocket_number}.pdb")
        with open(pdb_filename, "w") as pdb_file:
            pdb_file.write(f"HEADER    Pocket {pocket_number} Centroid\n")
            pdb_file.write(
                f"HETATM{pocket_number:5d} C   LIG A{pocket_number:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          C\n"
            )
            pdb_file.write("END\n")


def main():
    structures = pd.read_csv(STRUCTURES_DATA_CSV)
    structures = structures[structures["status"] == "ok"]

    report = []
    for _, row in structures.iterrows():
        gene_name, uniprot_ac, file_path = row["gene_name"], row["uniprot_ac"], row["file_path"]  # file_path is repo-relative
        absolute_file_path = os.path.join(ROOT, file_path)
        file_name = os.path.basename(file_path)
        file_stem = file_name.replace(".pdb", "")

        print(f"---------------  Detecting pockets in: {uniprot_ac}, {file_name}  -------------")
        gene_out_dir = os.path.join(DETECTED_POCKETS_DIR, uniprot_ac, file_stem)
        os.makedirs(gene_out_dir, exist_ok=True)

        detect_pockets(absolute_file_path, gene_out_dir)

        predictions_csv = os.path.join(gene_out_dir, f"{file_name}_predictions.csv")
        pocket_centroids = extract_pocket_centers(predictions_csv)
        pocket_scores = extract_pocket_scores(predictions_csv)
        pocket_probabilities = extract_pocket_probabilities(predictions_csv)
        pocket_residues = extract_pocket_residues(predictions_csv)

        res_geom = parse_residue_geometry(absolute_file_path)

        write_pocket_pdbs(pocket_centroids, os.path.join(gene_out_dir, "pockets"))

        for pocket in sorted(pocket_centroids):
            neighborhood = residues_within_radius(res_geom, pocket_centroids[pocket])
            n_within = len(neighborhood)
            min_plddt = min(b for _, b in neighborhood) if neighborhood else None
            mean_plddt = sum(b for _, b in neighborhood) / n_within if neighborhood else None

            report.append([
                gene_name, uniprot_ac, file_name, "alphafold2", file_path, pocket,
                pocket_scores[pocket], pocket_probabilities[pocket],
                " ".join(str(c) for c in pocket_centroids[pocket]),
                n_within, min_plddt, mean_plddt,
                " ".join(pocket_residues[pocket]),
            ])
        print(f"  {len(pocket_centroids)} pocket(s) detected")

    report = pd.DataFrame(report, columns=[
        "Gene name", "Uniprot AC", "File name", "Prediction type", "Full path", "Pocket number",
        "Pocket score", "Pocket probability", "Pocket centroid coordinate (x y z)",
        f"N residues within {RADIUS_A}A", f"Min pLDDT within {RADIUS_A}A",
        f"Mean pLDDT within {RADIUS_A}A", "P2Rank residues (chain_resn)",
    ])
    out_path = os.path.join(OUTPUT_DIR, "pocket_detection_data.csv")
    report.to_csv(out_path, index=False)

    n_genes_with_pockets = report["Uniprot AC"].nunique()
    print(f"\n{len(report)} pockets detected across {n_genes_with_pockets}/{len(structures)} genes.")
    print(f"Saved to {out_path}")


if __name__ == "__main__":
    main()
