#!/usr/bin/env python3
"""
Detect pockets in experimental, potentially multimeric PDB structures.

For each given PDB code: downloads the RCSB-annotated biological assembly (falling back to the
as-deposited asymmetric unit if no assembly is annotated), strips ligands/waters/ions while
keeping all protein chains, and detects pockets with P2Rank (as in
scripts/08_detect_pockets.py), using P2Rank's default config instead of "-c alphafold" and
without script 08's per-residue B-factor/pLDDT confidence gate, since crystallographic B-factors
are not a prediction-confidence score and that gate's semantics don't transfer to experimental
structures. Unlike the main UniProt-AC-keyed pipeline (scripts 00-17), this script is keyed by
PDB code and is not restricted to a single chain, since it is specifically meant to characterize
pockets in real multi-chain biological assemblies (e.g. at subunit interfaces).

Unlike script 04, this script does NOT run PDB2PQR protonation or PyRosetta's FastRelax on the
structure: experimental structures are already validated against experimental data by
crystallographic/cryo-EM refinement, so an unconstrained relax risks moving pocket residues away
from their observed conformation rather than improving them (the opposite of what relax is for
with AlphaFold/homology models, which lack that experimental support). Pocket detection here runs
directly on the ligand-stripped structure.

Usage:
    python 48_detect_pocket_multimers.py --pdb-codes 6XYZ,7ABC
"""
import argparse
import os
import subprocess
import sys
from datetime import date

import pandas as pd
import pymol
import requests
from pymol import cmd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))

OUTPUT_ROOT = os.path.join(ROOT, "output", "48_detect_pocket_multimers")
ASSEMBLIES_DIR = os.path.join(OUTPUT_ROOT, "biological_assemblies")
STRIPPED_DIR = os.path.join(OUTPUT_ROOT, "stripped_structures")
POCKETS_DIR = os.path.join(OUTPUT_ROOT, "detected_pockets")
REPORT_CSV = os.path.join(OUTPUT_ROOT, "pocket_detection_data.csv")
TMP_DIR = os.path.join(ROOT, "tmp")

# Define path to P2RANK binary package - downloaded from https://github.com/rdk/p2rank/releases
PATH_TO_P2RANK = "/aloy/home/acomajuncosa/programs/p2rank_2.5/prank"  # Change as needed

# Pocket selection cut-offs, inherited unchanged from script 08_detect_pockets.py -
# see https://github.com/rdk/p2rank/issues/76
POCKET_PROBABILITY_THRESHOLD = 0.2  # P
POCKET_TOP_K = 3  # K

PDB_DOWNLOAD_URL = "https://files.rcsb.org/download/{}.cif"
PDB_ASSEMBLY_DOWNLOAD_URL = "https://files.rcsb.org/download/{}-assembly{}.cif"
RCSB_ENTRY_API_URL = "https://data.rcsb.org/rest/v1/core/entry/{}"
RCSB_ASSEMBLY_API_URL = "https://data.rcsb.org/rest/v1/core/assembly/{}/{}"


def get_rcsb_entry_info(pdb_id):
    """
    Queries the RCSB entry API for a PDB ID's primary annotated assembly ID(s),
    experimental method and resolution. Returns None on any error.
    """
    try:
        resp = requests.get(RCSB_ENTRY_API_URL.format(pdb_id), timeout=15)
        resp.raise_for_status()
    except requests.RequestException as e:
        print(f"Warning: could not query RCSB entry API for {pdb_id} ({e}).")
        return None

    data = resp.json()
    assembly_ids = data.get("rcsb_entry_container_identifiers", {}).get("assembly_ids", [])
    exptl = data.get("exptl", [])
    method = exptl[0]["method"] if exptl else "N/A"
    resolution_list = data.get("rcsb_entry_info", {}).get("resolution_combined", [])
    resolution = resolution_list[0] if resolution_list else "N/A"

    return {"assembly_ids": assembly_ids, "method": method, "resolution": resolution}


def get_biological_assembly_info(pdb_id, assembly_id):
    """
    Queries the RCSB assembly API for the oligomeric state of a given PDB ID/assembly ID pair.
    Returns a human-readable summary string (e.g. "heterotetrameric (4 chains),
    author_and_software_defined_assembly"), or "N/A" on any error/missing data. The "details"
    field distinguishes author-confirmed assemblies from software-only (PISA) predictions, so it
    is reported verbatim rather than collapsed into a yes/no verdict.
    """
    try:
        resp = requests.get(RCSB_ASSEMBLY_API_URL.format(pdb_id, assembly_id), timeout=15)
        resp.raise_for_status()
    except requests.RequestException as e:
        print(f"Warning: could not query RCSB assembly API for {pdb_id}/{assembly_id} ({e}).")
        return "N/A"

    asm = resp.json().get("pdbx_struct_assembly", {})
    oligomeric = asm.get("oligomeric_details", "unknown oligomeric state")
    count = asm.get("oligomeric_count")
    count_str = f" ({count} chains)" if count else ""
    details = asm.get("details", "unknown method")
    return f"{oligomeric}{count_str}, {details}"


def download_biological_assembly(pdb_code, output_cif, assembly_id):
    """
    Downloads the RCSB biological assembly CIF for a given PDB code and assembly ID.
    Returns True on success, False on any download error.
    """
    try:
        resp = requests.get(PDB_ASSEMBLY_DOWNLOAD_URL.format(pdb_code, assembly_id), timeout=30)
        resp.raise_for_status()
    except requests.RequestException as e:
        print(f"Warning: could not download biological assembly {assembly_id} for {pdb_code} ({e}).")
        return False
    with open(output_cif, "wb") as f:
        f.write(resp.content)
    return True


def download_asymmetric_unit(pdb_code, output_cif):
    """
    Downloads the as-deposited asymmetric unit CIF for a given PDB code.
    Returns True on success, False on any download error (e.g. invalid PDB code).
    """
    try:
        resp = requests.get(PDB_DOWNLOAD_URL.format(pdb_code), timeout=30)
        resp.raise_for_status()
    except requests.RequestException as e:
        print(f"Warning: could not download {pdb_code} ({e}).")
        return False
    with open(output_cif, "wb") as f:
        f.write(resp.content)
    return True


def download_structure(pdb_code, output_cif, entry_info):
    """
    Downloads the RCSB biological assembly for a PDB code if one is annotated in entry_info,
    otherwise falls back to the as-deposited asymmetric unit. Returns (success, assembly_id_used),
    where assembly_id_used is None if the asymmetric-unit fallback was used.
    """
    assembly_ids = entry_info["assembly_ids"] if entry_info else []
    if assembly_ids:
        assembly_id = assembly_ids[0]
        if download_biological_assembly(pdb_code, output_cif, assembly_id):
            return True, assembly_id
        print(f"Warning: falling back to asymmetric unit for {pdb_code}.")
    else:
        print(f"Warning: no annotated biological assembly for {pdb_code}, using asymmetric unit.")
    return download_asymmetric_unit(pdb_code, output_cif), None


def convert_cif_to_pdb(input_cif, output_pdb):
    """
    Converts a CIF file to a PDB file using PyMOL. Verbatim from
    scripts/02_organize_structures.py.
    """
    cmd.reinitialize()
    cmd.load(input_cif)
    cmd.save(output_pdb, "all")
    cmd.delete("all")


def get_chain_ids(input_file):
    """
    Extract all chain IDs from a PDB file. Verbatim from scripts/02_organize_structures.py.
    """
    cmd.load(input_file, "structure")
    chain_ids = cmd.get_chains("structure")
    cmd.delete("all")
    return chain_ids


def strip_heteroatoms_all_chains(pdb_file, output_file):
    """
    Removes ligands, waters and other non-protein heteroatoms from a PDB file, keeping ALL
    protein chains intact. Multi-chain adaptation of scripts/02_organize_structures.py's
    remove_non_sequence_elements, which is scoped to a single chosen chain - here the
    "chain {chain_id} and" restriction is dropped so both the selection and the final save
    act across every chain, since this script is specifically about multimers.
    """
    cmd.load(pdb_file, "structure")
    cmd.select("non_sequence", "not polymer.protein")
    cmd.remove("non_sequence")
    cmd.save(output_file, "all")
    cmd.delete("all")


def detect_pockets(pdb_file, output_dir):
    """
    Runs P2Rank to detect pockets in a given PDB file, using P2Rank's default config (unlike
    scripts/08_detect_pockets.py, which uses "-c alphafold" - not appropriate here since these
    are experimental structures, not AlphaFold-predicted models).
    """
    command = [PATH_TO_P2RANK, "predict", "-f", pdb_file, "-o", output_dir, "-visualizations", "0"]
    try:
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"Error detecting pockets in {pdb_file}: {e}")


def extract_pocket_centers(csv_file):
    """Reads a P2Rank output CSV file and maps pocket number (rank) to its center coordinates."""
    df = pd.read_csv(csv_file)
    df.columns = df.columns.str.strip()
    return {int(row["rank"]): (row["center_x"], row["center_y"], row["center_z"]) for _, row in df.iterrows()}


def extract_pocket_scores(csv_file):
    """Reads a P2Rank output CSV file and maps pocket number (rank) to its score."""
    df = pd.read_csv(csv_file)
    df.columns = df.columns.str.strip()
    return {int(row["rank"]): row["score"] for _, row in df.iterrows()}


def extract_pocket_probabilities(csv_file):
    """Reads a P2Rank output CSV file and maps pocket number (rank) to its probability."""
    df = pd.read_csv(csv_file)
    df.columns = df.columns.str.strip()
    return {int(row["rank"]): row["probability"] for _, row in df.iterrows()}


def extract_pocket_residues(csv_file):
    """Reads a P2Rank output CSV file and maps pocket number (rank) to its residue IDs."""
    df = pd.read_csv(csv_file)
    df.columns = df.columns.str.strip()
    return {int(row["rank"]): row["residue_ids"].split() for _, row in df.iterrows()}


def write_pocket_pdbs(pockets_to_consider, pocket_dict, output_dir):
    """
    Creates a separate PDB file per detected pocket, with a single HETATM pseudo-atom at the
    pocket centroid. Verbatim from scripts/08_detect_pockets.py.
    """
    os.makedirs(output_dir, exist_ok=True)
    for pocket_number, (x, y, z) in pocket_dict.items():
        if pocket_number in pockets_to_consider:
            pdb_filename = os.path.join(output_dir, f"pocket_{pocket_number}.pdb")
            with open(pdb_filename, "w") as pdb_file:
                pdb_file.write(f"HEADER    Pocket {pocket_number} Centroid\n")
                pdb_file.write(
                    f"HETATM{pocket_number:5d} C   LIG A{pocket_number:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          C\n"
                )
                pdb_file.write("END\n")


def main():
    parser = argparse.ArgumentParser(
        description="Download experimental PDB structures, strip ligands/waters while keeping "
                     "all protein chains, and detect pockets (as in script 08, using P2Rank's "
                     "default config and without the pLDDT/QSQE confidence gate, neither of "
                     "which apply to experimental structures - see module docstring)."
    )
    parser.add_argument("--pdb-codes", required=True, help="Comma-separated PDB codes, e.g. 6XYZ,7ABC")
    args = parser.parse_args()

    pdb_codes = [c.strip().upper() for c in args.pdb_codes.split(",") if c.strip()]

    for d in (ASSEMBLIES_DIR, STRIPPED_DIR, POCKETS_DIR, TMP_DIR):
        os.makedirs(d, exist_ok=True)

    pymol.finish_launching(["pymol", "-cq"])

    report_rows = []
    access_date = date.today().isoformat()

    for pdb_code in pdb_codes:
        print("----------------------------------------")
        print(f"Processing {pdb_code}...")
        try:
            entry_info = get_rcsb_entry_info(pdb_code)

            # 1. Download biological assembly, falling back to the asymmetric unit (idempotent)
            assembly_cif = os.path.join(ASSEMBLIES_DIR, f"{pdb_code}.cif")
            if not os.path.exists(assembly_cif):
                success, assembly_id = download_structure(pdb_code, assembly_cif, entry_info)
                if not success:
                    print(f"Warning: could not obtain a structure for {pdb_code}. Skipping.")
                    continue
            else:
                print(f"Structure for {pdb_code} already downloaded. Skipping download.")
                assembly_id = entry_info["assembly_ids"][0] if entry_info and entry_info["assembly_ids"] else None

            assembly_info = (
                get_biological_assembly_info(pdb_code, assembly_id)
                if assembly_id is not None
                else "N/A (asymmetric unit used, no annotated assembly)"
            )

            # 2. Convert to PDB and strip heteroatoms, keeping all protein chains (idempotent)
            stripped_pdb = os.path.join(STRIPPED_DIR, f"{pdb_code}.pdb")
            if not os.path.exists(stripped_pdb):
                tmp_full = os.path.join(TMP_DIR, f"{pdb_code}_full.pdb")
                convert_cif_to_pdb(assembly_cif, tmp_full)
                strip_heteroatoms_all_chains(tmp_full, stripped_pdb)
                os.remove(tmp_full)
            else:
                print(f"Stripped structure for {pdb_code} already exists. Skipping stripping.")

            # 3. Detect pockets (idempotent)
            pocket_dir = os.path.join(POCKETS_DIR, pdb_code)
            os.makedirs(pocket_dir, exist_ok=True)
            predictions_csv = os.path.join(pocket_dir, f"{pdb_code}.pdb_predictions.csv")
            if not os.path.exists(predictions_csv):
                detect_pockets(stripped_pdb, pocket_dir)
            else:
                print(f"P2Rank predictions for {pdb_code} already exist. Skipping P2Rank.")

            if not os.path.exists(predictions_csv):
                print(f"Warning: P2Rank produced no predictions for {pdb_code}. Skipping.")
                continue

            pocket_centroids = extract_pocket_centers(predictions_csv)
            pocket_scores = extract_pocket_scores(predictions_csv)
            pocket_probabilities = extract_pocket_probabilities(predictions_csv)
            pocket_residues = extract_pocket_residues(predictions_csv)

            # Select pockets to consider - same P/K cut-offs as script 08_detect_pockets.py,
            # unchanged - see https://github.com/rdk/p2rank/issues/76.
            # NOTE: script 08_detect_pockets.py additionally drops pockets with any residue's
            # B-factor (read there as pLDDT/QSQE confidence) below a threshold. That gate is
            # SKIPPED here: these are experimental crystal structures, and the PDB B-factor
            # column encodes atomic displacement, not a prediction-confidence score, so 08's
            # threshold semantics do not transfer.
            pockets_to_consider = set(
                i for i in sorted(pocket_centroids)
                if i < POCKET_TOP_K or pocket_probabilities[i] >= POCKET_PROBABILITY_THRESHOLD
            )

            write_pocket_pdbs(pockets_to_consider, pocket_centroids, os.path.join(pocket_dir, "pockets"))

            chains = get_chain_ids(stripped_pdb)

            for ptc in sorted(pockets_to_consider):
                report_rows.append([
                    pdb_code,
                    f"{pdb_code}.pdb",
                    os.path.join("output", "48_detect_pocket_multimers", "stripped_structures", f"{pdb_code}.pdb"),
                    " ".join(chains),
                    ptc,
                    pocket_scores[ptc],
                    pocket_probabilities[ptc],
                    " ".join(str(c) for c in pocket_centroids[ptc]),
                    " ".join(pocket_residues[ptc]),
                    entry_info["method"] if entry_info else "N/A",
                    entry_info["resolution"] if entry_info else "N/A",
                    assembly_info,
                    access_date,
                ])

            print(f"Done with {pdb_code}.")
        except Exception as e:
            print(f"Warning: {pdb_code} failed ({e}). Skipping and continuing.")
            continue

    report_df = pd.DataFrame(report_rows, columns=[
        "PDB code", "File name", "Full path", "Chains", "Pocket number", "Pocket score",
        "Pocket probability", "Pocket centroid coordinate (x y z)", "Pocket residues (chain_resn)",
        "Experimental method", "Resolution (A)", "Biological assembly info", "RCSB access date",
    ])
    report_df.to_csv(REPORT_CSV, index=False)
    print(f"Saved pocket detection report to {REPORT_CSV}")


if __name__ == "__main__":
    main()
