#!/usr/bin/env python3
"""
Prepares the two core Boltz-2 inputs -- per-pocket protein sequence (with translated pocket-contact
positions) and the compound SMILES list -- for the 11 curated single-chain pockets in
selected_pockets.csv. The 12th, the 7K98_pocket_6 dimer pocket, is deferred: it needs a 2-chain
pheS+pheT complex, unlike the other 11.

Pocket residues in pocket_detection_data.csv are raw PDB residue numbers, not sequence positions,
and numbering differs by structure source (SwissModel chains start at PDB residue 2 or 3;
AlphaFold2/3 and Chai1 chains start at 1). This script builds a per-structure
{pdb_resnum: seq_position} map by enumerating residues in file order, then translates each pocket's
residues through that map -- Boltz-2's pocket-constraint field expects 1-indexed positions in the
sequence string, not raw PDB residue numbers.

YAML generation and running Boltz-2 itself are separate, later scripts.

Usage:
    python 71_boltz2_prepare_inputs.py
"""
import os

import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")
POCKET_DETECTION_DATA_CSV = os.path.join(ROOT, "output", "pocket_detection_data.csv")
ALIGNED_RELAXED_STRUCTURES_DIR = os.path.join(ROOT, "output", "aligned_relaxed_structures")
FILTERED_HITS_CSV = os.path.join(ROOT, "output", "70_filtering", "filtered_hits.csv")

OUTPUT_DIR = os.path.join(ROOT, "output", "71_boltz2_prepare_inputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

DIMER_POCKET = "7K98_pocket_6"


def load_pockets():
    """[(gene_name, site_type, pocket_name), ...] for the 11 curated single-chain pockets
    (excludes the 7K98 dimer pocket)."""
    df = pd.read_csv(SELECTED_POCKETS_CSV)
    df = df[df["pocket_name"] != DIMER_POCKET]
    return list(zip(df["gene_name"], df["site_type"], df["pocket_name"]))


def load_pocket_residues():
    """{pocket_name: (uniprot_ac, file_name, chain, [pdb_resnum, ...])}."""
    df = pd.read_csv(POCKET_DETECTION_DATA_CSV)
    pocket_names = df["File name"].str.replace(".pdb", "", regex=False) + "_pocket_" + df["Pocket number"].astype(str)
    info = {}
    for pocket_name, ac, file_name, residues in zip(
        pocket_names, df["Uniprot AC"], df["File name"], df["Pocket residues (chain_resn)"]
    ):
        tokens = residues.split()
        chain = tokens[0].split("_")[0]
        resnums = [int(tok.split("_")[1]) for tok in tokens]
        info[pocket_name] = (ac, file_name, chain, resnums)
    return info


def extract_sequence_and_map(pdb_path, chain_id):
    """Parses pdb_path's chain_id, returns (1-letter sequence, {pdb_resnum: 1-indexed seq_position}).
    The map is built by enumerating standard residues in file order, so it correctly handles
    per-structure numbering offsets (e.g. SwissModel chains starting at residue 2 or 3) without
    assuming a fixed offset."""
    structure = PDBParser(QUIET=True).get_structure("x", pdb_path)
    residues = [r for r in structure[0][chain_id].get_residues() if r.id[0] == " "]

    sequence_chars = []
    resnum_to_position = {}
    for position, residue in enumerate(residues, start=1):
        one_letter = seq1(residue.get_resname())
        if one_letter == "X":
            print(f"  Warning: unrecognized residue {residue.get_resname()} at "
                  f"{chain_id}_{residue.id[1]} in {pdb_path}, encoded as 'X'.")
        sequence_chars.append(one_letter)
        resnum_to_position[residue.id[1]] = position

    return "".join(sequence_chars), resnum_to_position


def main():
    pockets = load_pockets()
    print(f"Curated single-chain pockets (excl. dimer {DIMER_POCKET}): {len(pockets)}")

    pocket_residues = load_pocket_residues()

    # Keyed by structure file, not pocket, so shared structures (lysS's CAT and NON-CAT pockets
    # both use the same AlphaFold3 model) are only parsed once.
    structure_cache = {}
    rows = []
    for gene_name, site_type, pocket_name in pockets:
        ac, file_name, chain, resnums = pocket_residues[pocket_name]

        if file_name not in structure_cache:
            pdb_path = os.path.join(ALIGNED_RELAXED_STRUCTURES_DIR, ac, file_name)
            structure_cache[file_name] = extract_sequence_and_map(pdb_path, chain)
        sequence, resnum_to_position = structure_cache[file_name]

        missing = [r for r in resnums if r not in resnum_to_position]
        assert not missing, f"{pocket_name}: pocket residue(s) {missing} not found in {file_name} chain {chain}"
        positions = sorted(resnum_to_position[r] for r in resnums)

        rows.append({
            "gene_name": gene_name,
            "site_type": site_type,
            "pocket_name": pocket_name,
            "uniprot_ac": ac,
            "structure_file": file_name,
            "chain": chain,
            "sequence": sequence,
            "sequence_length": len(sequence),
            "pocket_contacts": " ".join(str(p) for p in positions),
            "n_pocket_residues": len(positions),
        })
        print(f"  {gene_name} ({site_type}) {pocket_name}: seq_len={len(sequence)}, "
              f"n_pocket_res={len(positions)}, structure={file_name}")

    pocket_df = pd.DataFrame(rows)
    pocket_out_path = os.path.join(OUTPUT_DIR, "pocket_sequences.csv")
    pocket_df.to_csv(pocket_out_path, index=False)
    print(f"\nSaved {len(pocket_df)} pockets ({len(structure_cache)} unique structures) to {pocket_out_path}")

    compounds = pd.read_csv(FILTERED_HITS_CSV, usecols=["compound_id", "smiles"])
    compounds_out_path = os.path.join(OUTPUT_DIR, "compounds.csv")
    compounds.to_csv(compounds_out_path, index=False)
    print(f"Saved {len(compounds)} compounds to {compounds_out_path}")


if __name__ == "__main__":
    main()
