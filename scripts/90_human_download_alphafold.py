#!/usr/bin/env python3
"""
Downloads structures per aminoacyl-tRNA synthetase gene, for the structure-based counter-screen of
the 1,095 filtered Mtb hits: one AlphaFold DB predicted monomer (queried directly, model version
resolved per UniProt AC rather than hardcoded), plus one AlphaFill entry (.cif + .json) per gene --
same source and download mechanism as script 01's Mtb AlphaFill pull, kept as a raw download only
(no PDB conversion, no ligand stripping, no report) for now; a possible later step for
ligand-evidence-based catalytic-pocket annotation, mirroring script 77.

--organism selects the gene set: "human" (default, 38 genes, data/human_trna_synthetases_uniprot.csv)
-- the original counter-screen -- or "mtb" (all 21 CRISPR-screen genes,
data/mtb_trna_synthetases_bosch_2021_fig5_annotated.csv), added to run the exact same
AF2-monomer-only, no-curation recipe against the Mtb targets themselves for a directly comparable
on-target/off-target scan. The two organisms write to entirely separate paths (see
src/utils/counterscreen_targets.py), so running one never touches the other's data.

Resumable: skips a download if its target file already exists on disk.

Usage:
    python 90_human_download_alphafold.py [--organism human|mtb]
"""
import argparse
import os
import sys

import requests
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from utils.counterscreen_targets import load_targets  # noqa: E402

AFDB_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{}"


def fetch_afdb_entry(uniprot_ac):
    """Latest AlphaFold DB prediction entry (dict) for uniprot_ac, or None if unavailable."""
    try:
        resp = requests.get(AFDB_API_URL.format(uniprot_ac), timeout=30)
        resp.raise_for_status()
        entries = resp.json()
    except requests.exceptions.RequestException as e:
        print(f"  WARNING: AFDB API request failed for {uniprot_ac}: {e}")
        return None
    if not entries:
        print(f"  WARNING: no AlphaFold DB entry found for {uniprot_ac}")
        return None
    return entries[0]


def download_pdb(pdb_url, out_path):
    """Downloads pdb_url to out_path; no-op if out_path already exists. Returns success bool."""
    if os.path.isfile(out_path):
        return True
    try:
        resp = requests.get(pdb_url, stream=True, timeout=60)
        resp.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=8192):
                f.write(chunk)
        return True
    except requests.exceptions.RequestException as e:
        print(f"  WARNING: failed to download {pdb_url}: {e}")
        return False


def download_alphafill(uniprot_ac, output_dir):
    """Downloads AlphaFill's .cif (structure + transplanted ligands) and .json (hit/transplant
    metadata) for uniprot_ac into output_dir; no-op per file if it already exists on disk."""
    os.makedirs(output_dir, exist_ok=True)
    base_url = f"https://alphafill.eu/v1/aff/{uniprot_ac}"
    cif_path = os.path.join(output_dir, f"{uniprot_ac}.cif")
    json_path = os.path.join(output_dir, f"{uniprot_ac}.json")

    if not os.path.isfile(cif_path):
        try:
            resp = requests.get(base_url, stream=True, timeout=60)
            resp.raise_for_status()
            with open(cif_path, "wb") as f:
                for chunk in resp.iter_content(chunk_size=8192):
                    f.write(chunk)
        except requests.exceptions.RequestException as e:
            print(f"  WARNING: failed to download AlphaFill cif for {uniprot_ac}: {e}")

    if not os.path.isfile(json_path):
        try:
            resp = requests.get(f"{base_url}/json", timeout=30)
            resp.raise_for_status()
            with open(json_path, "wb") as f:
                f.write(resp.content)
        except requests.exceptions.RequestException as e:
            print(f"  WARNING: failed to download AlphaFill json for {uniprot_ac}: {e}")


def parse_structure_stats(pdb_path):
    """(n_residues, mean_plddt, min_plddt) from a PDB's ATOM lines -- AlphaFold stores per-residue
    pLDDT in the B-factor column, same convention script 08 relies on."""
    res_to_bfactor = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                resnum = int(line[22:26])
                bfactor = float(line[60:66])
                res_to_bfactor[resnum] = bfactor
    if not res_to_bfactor:
        return 0, None, None
    values = list(res_to_bfactor.values())
    return len(res_to_bfactor), sum(values) / len(values), min(values)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--organism", choices=["human", "mtb"], default="human",
                         help="Gene set to download structures for (default: human)")
    args = parser.parse_args()

    structures_dir = os.path.join(ROOT, "data", "structures", f"{args.organism}_alphafold2_database")
    alphafill_dir = os.path.join(ROOT, "data", "structures", f"{args.organism}_alphafill_database")
    output_dir = os.path.join(ROOT, "output", f"90_{args.organism}_download_alphafold")
    os.makedirs(structures_dir, exist_ok=True)
    os.makedirs(alphafill_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    proteins = load_targets(args.organism)
    rows = []

    for _, row in proteins.iterrows():
        gene_name, uniprot_ac = row["gene_name"], row["uniprot_ac"]
        n_residues_uniprot = row["sequence_length"]
        print(f"--- {gene_name} ({uniprot_ac}) ---")

        gene_dir = os.path.join(structures_dir, uniprot_ac)
        os.makedirs(gene_dir, exist_ok=True)

        download_alphafill(uniprot_ac, os.path.join(alphafill_dir, uniprot_ac))

        entry = fetch_afdb_entry(uniprot_ac)
        if entry is None:
            rows.append([gene_name, uniprot_ac, None, n_residues_uniprot, None, None, None, None,
                         "no_afdb_entry"])
            continue

        pdb_url = entry["pdbUrl"]
        file_name = os.path.basename(pdb_url)
        out_path = os.path.join(gene_dir, file_name)

        if not download_pdb(pdb_url, out_path):
            rows.append([gene_name, uniprot_ac, None, n_residues_uniprot, None, None, None, None,
                         "download_failed"])
            continue

        n_residues_pdb, mean_plddt, min_plddt = parse_structure_stats(out_path)
        coverage_pct = round(100 * n_residues_pdb / n_residues_uniprot, 1)
        relative_path = os.path.relpath(out_path, ROOT)

        rows.append([gene_name, uniprot_ac, relative_path, n_residues_uniprot, n_residues_pdb,
                     coverage_pct, mean_plddt, min_plddt, "ok"])
        print(f"  {file_name}: {n_residues_pdb}/{n_residues_uniprot} residues "
              f"({coverage_pct}%), mean pLDDT={mean_plddt:.1f}")

    report = pd.DataFrame(rows, columns=[
        "gene_name", "uniprot_ac", "file_path", "n_residues_uniprot", "n_residues_pdb",
        "coverage_pct", "mean_plddt", "min_plddt", "status",
    ])
    out_path = os.path.join(output_dir, "structures_data.csv")
    report.to_csv(out_path, index=False)

    n_ok = (report["status"] == "ok").sum()
    print(f"\n{n_ok}/{len(report)} structures downloaded successfully.")
    if n_ok < len(report):
        print("Non-ok genes:")
        print(report.loc[report["status"] != "ok", ["gene_name", "uniprot_ac", "status"]]
              .to_string(index=False))
    print(f"Saved report to {out_path}")


if __name__ == "__main__":
    main()
