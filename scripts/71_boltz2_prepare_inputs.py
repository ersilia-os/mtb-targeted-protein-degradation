#!/usr/bin/env python3
"""
Prepares the two core Boltz-2 inputs -- per-pocket protein sequence (with translated pocket-contact
positions) and the compound SMILES list -- for all 12 curated pockets in selected_pockets.csv: the
11 single-chain pockets, plus the 7K98_pocket_6 dimer pocket (a pheS+pheT 2-chain complex, since
that pocket sits at the pheS-pheT dimerization interface). The dimer row gets extra "_b"-suffixed
columns (chain_b, uniprot_ac_b, structure_file_b, sequence_b, sequence_length_b,
pocket_contacts_b, n_pocket_residues_b, partner_gene_name) for its second (pheT) chain, plus
chain_a_trim_start/end and chain_b_trim_start/end recording the pocket-window trim applied to
both chains (see DIMER_TRIM_MARGIN_A/B below -- the untrimmed complex OOMs Boltz-2); these are
NaN for the other 11, single-chain, rows.

Pocket residues in pocket_detection_data.csv (single-chain pockets) and in the 7K98 multimer
predictions CSV (the dimer pocket) are raw PDB residue numbers, not sequence positions, and
numbering differs by structure source (SwissModel chains start at PDB residue 2 or 3; AlphaFold2/3
and Chai1 chains start at 1; the 7K98 crystal structure starts at -1/-2). This script builds a
per-structure {pdb_resnum: seq_position} map by enumerating residues in file order, then translates
each pocket's residues through that map -- Boltz-2's pocket-constraint field expects 1-indexed
positions in the sequence string, not raw PDB residue numbers.

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

# The 7K98_pocket_6 dimer pocket lives in the separate multimer-detection pipeline (scripts
# 48-50), not in pocket_detection_data.csv / aligned_relaxed_structures -- 7K98 is a pheS2-pheT2
# heterotetramer (chains A/D = pheS/P9WFU3, chains B/E = pheT/P9WFU1, confirmed by direct sequence
# comparison against the other pheS/pheT pockets); pocket 6 sits at the A/B (one pheS + one pheT)
# interface.
MULTIMER_STRIPPED_STRUCTURES_DIR = os.path.join(ROOT, "output", "48_detect_pocket_multimers", "stripped_structures")
MULTIMER_PREDICTIONS_CSV = os.path.join(
    ROOT, "output", "48_detect_pocket_multimers", "detected_pockets", "7K98", "7K98.pdb_predictions.csv"
)

OUTPUT_DIR = os.path.join(ROOT, "output", "71_boltz2_prepare_inputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

DIMER_POCKET = "7K98_pocket_6"
DIMER_PDB_CODE = "7K98"
DIMER_POCKET_ROW_NAME = "pocket6"
DIMER_CHAIN_A = "A"  # pheS
DIMER_CHAIN_B = "B"  # pheT

# The full complex (343 + 834 = 1,177 residues) reliably OOMs Boltz-2 on a 24GB GPU -- the
# pairwise (residue x residue) representation scales with total complex size, and no exposed
# Boltz-2 CLI flag (MSA depth, recycling steps) reduces that term. Trimming each chain to a
# margin around its own pocket contacts (chain A contacts span 111-312 of 343, already most of
# the chain, hence the smaller margin; chain B contacts span 507-631 of 834, with much more room
# to spare) cuts the total to 587 residues. Validated empirically: a manual boltz predict run on
# exactly these margins completed cleanly (no OOM, ligand_iptm=0.984, complex_plddt=0.954), and
# superposing that prediction back onto the original 7K98 crystal structure (587 matched Ca atoms,
# RMSD 2.24 A) showed the nearest EXCLUDED residue sits 10.57 A from the ligand -- well outside
# protein-ligand contact range (~4-5 A) -- so nothing pocket-relevant is being cut away.
DIMER_TRIM_MARGIN_A = 30
DIMER_TRIM_MARGIN_B = 100


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


def load_dimer_pocket_row():
    """(gene_name, site_type) for the dimer pocket, from selected_pockets.csv (gene_name=pheS,
    site_type=NON-CAT -- the pocket is assigned to pheS's site-type slot even though it's a pheS+
    pheT interface pocket)."""
    df = pd.read_csv(SELECTED_POCKETS_CSV)
    row = df[df["pocket_name"] == DIMER_POCKET].iloc[0]
    return row["gene_name"], row["site_type"]


def load_dimer_contacts():
    """{chain_id: [pdb_resnum, ...]} for the dimer pocket, parsed from the 7K98 multimer
    predictions CSV's residue_ids field (format "A_109 A_113 ... B_504 ..."). Only chains A and B
    are expected for this pocket (the pheS-pheT interface within the pheS2-pheT2 tetramer)."""
    df = pd.read_csv(MULTIMER_PREDICTIONS_CSV)
    df.columns = [c.strip() for c in df.columns]
    row = df[df["name"].str.strip() == DIMER_POCKET_ROW_NAME].iloc[0]
    contacts = {}
    for token in row["residue_ids"].split():
        chain, resnum = token.split("_")
        contacts.setdefault(chain, []).append(int(resnum))
    unexpected = set(contacts) - {DIMER_CHAIN_A, DIMER_CHAIN_B}
    assert not unexpected, f"{DIMER_POCKET}: unexpected chain(s) {unexpected} in pocket contacts"
    return contacts


def trim_to_pocket_window(sequence, contacts, margin):
    """Slices sequence down to [min(contacts) - margin, max(contacts) + margin] (clipped to
    [1, len(sequence)]) and remaps contacts into that window's own 1-indexed coordinates. Returns
    (trimmed_sequence, remapped_contacts, window_start, window_end), where window_start/end are
    in the ORIGINAL sequence's 1-indexed positions (for provenance)."""
    start = max(1, min(contacts) - margin)
    end = min(len(sequence), max(contacts) + margin)
    trimmed_sequence = sequence[start - 1:end]
    remapped_contacts = [p - start + 1 for p in contacts]
    return trimmed_sequence, remapped_contacts, start, end


def build_dimer_row(pocket_residues_by_gene):
    """One row for the 7K98_pocket_6 dimer pocket: chain A (pheS) fields as the normal columns,
    chain B (pheT) fields as "_b"-suffixed columns. Reuses extract_sequence_and_map() unchanged --
    it's already chain/structure-agnostic."""
    gene_name, site_type = load_dimer_pocket_row()
    contacts = load_dimer_contacts()

    pheS_ac = pocket_residues_by_gene["pheS"]
    pheT_ac = pocket_residues_by_gene["pheT"]

    pdb_path = os.path.join(MULTIMER_STRIPPED_STRUCTURES_DIR, f"{DIMER_PDB_CODE}.pdb")
    seq_a, map_a = extract_sequence_and_map(pdb_path, DIMER_CHAIN_A)
    seq_b, map_b = extract_sequence_and_map(pdb_path, DIMER_CHAIN_B)

    missing_a = [r for r in contacts[DIMER_CHAIN_A] if r not in map_a]
    missing_b = [r for r in contacts[DIMER_CHAIN_B] if r not in map_b]
    assert not missing_a, f"{DIMER_POCKET}: chain A residue(s) {missing_a} not found in {pdb_path}"
    assert not missing_b, f"{DIMER_POCKET}: chain B residue(s) {missing_b} not found in {pdb_path}"

    positions_a = sorted(map_a[r] for r in contacts[DIMER_CHAIN_A])
    positions_b = sorted(map_b[r] for r in contacts[DIMER_CHAIN_B])

    seq_a, positions_a, trim_start_a, trim_end_a = trim_to_pocket_window(
        seq_a, positions_a, DIMER_TRIM_MARGIN_A)
    seq_b, positions_b, trim_start_b, trim_end_b = trim_to_pocket_window(
        seq_b, positions_b, DIMER_TRIM_MARGIN_B)

    row = {
        "gene_name": gene_name,
        "site_type": site_type,
        "pocket_name": DIMER_POCKET,
        "uniprot_ac": pheS_ac,
        "structure_file": f"{DIMER_PDB_CODE}.pdb",
        "chain": DIMER_CHAIN_A,
        "sequence": seq_a,
        "sequence_length": len(seq_a),
        "pocket_contacts": " ".join(str(p) for p in positions_a),
        "n_pocket_residues": len(positions_a),
        "chain_a_trim_start": trim_start_a,
        "chain_a_trim_end": trim_end_a,
        "partner_gene_name": "pheT",
        "chain_b": DIMER_CHAIN_B,
        "uniprot_ac_b": pheT_ac,
        "structure_file_b": f"{DIMER_PDB_CODE}.pdb",
        "sequence_b": seq_b,
        "sequence_length_b": len(seq_b),
        "pocket_contacts_b": " ".join(str(p) for p in positions_b),
        "n_pocket_residues_b": len(positions_b),
        "chain_b_trim_start": trim_start_b,
        "chain_b_trim_end": trim_end_b,
    }
    print(f"  {gene_name} ({site_type}) {DIMER_POCKET} [dimer, trimmed]: "
          f"chain A seq_len={len(seq_a)} (window {trim_start_a}-{trim_end_a}) n_pocket_res={len(positions_a)}, "
          f"chain B seq_len={len(seq_b)} (window {trim_start_b}-{trim_end_b}) n_pocket_res={len(positions_b)}, "
          f"structure={DIMER_PDB_CODE}.pdb")
    return row


def main():
    pockets = load_pockets()
    print(f"Curated single-chain pockets (excl. dimer {DIMER_POCKET}): {len(pockets)}")

    pocket_residues = load_pocket_residues()

    # Keyed by structure file, not pocket, so shared structures (lysS's CAT and NON-CAT pockets
    # both use the same AlphaFold3 model) are only parsed once.
    structure_cache = {}
    rows = []
    uniprot_ac_by_gene = {}
    for gene_name, site_type, pocket_name in pockets:
        ac, file_name, chain, resnums = pocket_residues[pocket_name]
        uniprot_ac_by_gene.setdefault(gene_name, ac)

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

    rows.append(build_dimer_row(uniprot_ac_by_gene))

    pocket_df = pd.DataFrame(rows)
    pocket_out_path = os.path.join(OUTPUT_DIR, "pocket_sequences.csv")
    pocket_df.to_csv(pocket_out_path, index=False)
    print(f"\nSaved {len(pocket_df)} pockets ({len(structure_cache)} unique structures + 1 dimer) "
          f"to {pocket_out_path}")

    compounds = pd.read_csv(FILTERED_HITS_CSV, usecols=["compound_id", "smiles"])
    compounds_out_path = os.path.join(OUTPUT_DIR, "compounds.csv")
    compounds.to_csv(compounds_out_path, index=False)
    print(f"Saved {len(compounds)} compounds to {compounds_out_path}")


if __name__ == "__main__":
    main()
