#!/usr/bin/env python3
"""
Prepares Nesso-1 inputs from Boltz-2's already-vetted script 71 outputs: dedupes the 11 curated
single-chain pockets down to their 10 unique structures (two pockets --
alphafold3_P9WFU9_model_1_pocket_1/_2 for lysS -- share the exact same AlphaFold3 model, so
they'd get byte-identical Nesso-1 results; Nesso-1 has no pocket-conditioning yet to tell them
apart, see docs/prediction.md's "Coming Soon: Pocket Conditioning" in the recursionpharma/nesso
repo), and re-derives the 7K98_pocket_6 pheS+pheT dimer's UNTRIMMED full-length chains
(pocket_sequences.csv only stores the pocket-window-trimmed version Boltz-2 needed to avoid an
OOM -- Nesso-1 is a different, reportedly coarser/faster architecture with an unknown memory
ceiling, so script 80 tries the full 1,177-residue complex first and only falls back to the
pre-trimmed sequences on OOM).

Nesso-1's YAML schema has no pocket-constraint field at all (unlike Boltz-2's
`constraints: pocket: contacts:`), so pocket_contacts are not carried over -- only sequences.

Usage:
    python 78_nesso1_prepare_inputs.py
"""
import os

import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

POCKET_SEQUENCES_CSV = os.path.join(ROOT, "output", "71_boltz2_prepare_inputs", "pocket_sequences.csv")
MULTIMER_STRIPPED_STRUCTURES_DIR = os.path.join(ROOT, "output", "48_detect_pocket_multimers", "stripped_structures")

OUTPUT_DIR = os.path.join(ROOT, "output", "78_nesso1_prepare_inputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

DIMER_POCKET = "7K98_pocket_6"
DIMER_PDB_CODE = "7K98"
DIMER_CHAIN_A = "A"  # pheS
DIMER_CHAIN_B = "B"  # pheT


def extract_sequence(pdb_path, chain_id):
    """1-letter sequence for pdb_path's chain_id, in file order. No residue-number map is built
    here (unlike script 71's extract_sequence_and_map) since Nesso-1 has no pocket-constraint
    field to translate positions for."""
    structure = PDBParser(QUIET=True).get_structure("x", pdb_path)
    residues = [r for r in structure[0][chain_id].get_residues() if r.id[0] == " "]
    return "".join(seq1(r.get_resname()) for r in residues)


def build_single_chain_rows(pocket_sequences):
    """One row per unique structure_file among the 11 single-chain pockets (identified by a
    null sequence_b -- the dimer row is handled separately). Two pockets sharing a
    structure_file are guaranteed identical sequence, since script 71 parses each file once and
    caches it."""
    single = pocket_sequences[pocket_sequences["sequence_b"].isna()].copy()

    rows = []
    for structure_file, group in single.groupby("structure_file", sort=False):
        first = group.iloc[0]
        pocket_names = sorted(group["pocket_name"])
        if len(pocket_names) > 1:
            print(f"  {structure_file}: {len(pocket_names)} pockets share this structure "
                  f"{pocket_names} -- deduplicated to 1 Nesso-1 run.")
        rows.append({
            "structure_id": structure_file.replace(".pdb", ""),
            "pocket_names": " ".join(pocket_names),
            "gene_name": first["gene_name"],
            "structure_file": structure_file,
            "sequence": first["sequence"],
            "sequence_length": first["sequence_length"],
            "is_dimer": False,
            "sequence_b": None,
            "sequence_length_b": None,
            "dimer_trimmed_sequence": None,
            "dimer_trimmed_sequence_b": None,
        })
    return rows


def build_dimer_row(pocket_sequences):
    """Re-derives the FULL, untrimmed pheS (chain A) + pheT (chain B) sequences straight from the
    7K98 structure -- pocket_sequences.csv only has the pocket-window-trimmed version. Also
    carries the trimmed sequences through unchanged, as script 80's OOM-fallback path."""
    trimmed = pocket_sequences.set_index("pocket_name").loc[DIMER_POCKET]

    pdb_path = os.path.join(MULTIMER_STRIPPED_STRUCTURES_DIR, f"{DIMER_PDB_CODE}.pdb")
    seq_a = extract_sequence(pdb_path, DIMER_CHAIN_A)
    seq_b = extract_sequence(pdb_path, DIMER_CHAIN_B)
    print(f"  {DIMER_POCKET}: untrimmed chain A (pheS) len={len(seq_a)}, "
          f"chain B (pheT) len={len(seq_b)}, total={len(seq_a) + len(seq_b)} "
          f"(trimmed Boltz-2 version was {trimmed['sequence_length']} + "
          f"{trimmed['sequence_length_b']}).")

    return {
        "structure_id": f"{DIMER_PDB_CODE}_pheS_pheT_dimer",
        "pocket_names": DIMER_POCKET,
        "gene_name": trimmed["gene_name"],
        "structure_file": f"{DIMER_PDB_CODE}.pdb",
        "sequence": seq_a,
        "sequence_length": len(seq_a),
        "is_dimer": True,
        "sequence_b": seq_b,
        "sequence_length_b": len(seq_b),
        "dimer_trimmed_sequence": trimmed["sequence"],
        "dimer_trimmed_sequence_b": trimmed["sequence_b"],
    }


def main():
    pocket_sequences = pd.read_csv(POCKET_SEQUENCES_CSV)

    rows = build_single_chain_rows(pocket_sequences)
    rows.append(build_dimer_row(pocket_sequences))

    df = pd.DataFrame(rows)
    out_path = os.path.join(OUTPUT_DIR, "structure_sequences.csv")
    df.to_csv(out_path, index=False)

    n_pockets_covered = sum(len(p.split()) for p in df["pocket_names"])
    print(f"\nSaved {len(df)} unique structures ({n_pockets_covered} pockets covered, "
          f"including the dimer) to {out_path}")
    print("Compounds: reused directly from output/71_boltz2_prepare_inputs/compounds.csv, "
          "not copied.")


if __name__ == "__main__":
    main()
