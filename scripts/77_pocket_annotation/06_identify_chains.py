import os
import json
import pandas as pd
from pymol import cmd

root = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.join(root, "..", "..")
pdb_dir = os.path.join(repo_root, "processed", "pocket_annotation", "pdb_structures")
output_dir = os.path.join(repo_root, "output", "77_pocket_annotation")
manifest_dir = os.path.join(repo_root, "processed", "pocket_annotation")
os.makedirs(manifest_dir, exist_ok=True)

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def uniprot_seq_for(uniprot_ac):
    bosch = pd.read_csv(os.path.join(repo_root, "data", "mtb_trna_synthetases_bosch_2021_fig5_annotated.csv"))
    return bosch[bosch["uniprot_ac"] == uniprot_ac].iloc[0]["sequence"]


def chain_seq(obj_name, chain):
    model = cmd.get_model("{} and chain {} and name CA and polymer.protein".format(obj_name, chain))
    seq = {}
    for atom in model.atom:
        aa = AA3TO1.get(atom.resn)
        if aa is not None:
            seq[int(atom.resi)] = aa
    return seq


def best_offset_identity(pdb_seq, uniprot_seq, offset_range=range(-30, 31)):
    """Try PDB_resi + offset == UniProt_1indexed_position over a range of offsets
    (handles N-terminal tags/cloning artifacts shifting numbering), return the
    best-matching offset and its identity."""
    best = (0, 0, 0, 0)  # offset, matched, checked, identity
    for offset in offset_range:
        matched = 0
        checked = 0
        for resi, aa in pdb_seq.items():
            pos = resi + offset
            if pos < 1 or pos > len(uniprot_seq):
                continue
            checked += 1
            if aa == uniprot_seq[pos - 1]:
                matched += 1
        if checked > 0 and matched / checked > best[3]:
            best = (offset, matched, checked, matched / checked)
    return best


if __name__ == "__main__":
    xrefs = json.load(open(os.path.join(output_dir, "pdb_xrefs.json")))

    manifest = []
    for ac, info in xrefs.items():
        if not info["pdb_ids"]:
            continue
        uniprot_seq = uniprot_seq_for(ac)
        for pdb_id in info["pdb_ids"]:
            pdb_path = os.path.join(pdb_dir, "{}.pdb".format(pdb_id))
            cif_path = os.path.join(pdb_dir, "{}.cif".format(pdb_id))
            struct_path = pdb_path if os.path.exists(pdb_path) else cif_path
            fmt = "pdb" if os.path.exists(pdb_path) else "cif"
            if not os.path.exists(struct_path):
                print("MISSING structure file for {}".format(pdb_id))
                continue

            cmd.reinitialize()
            cmd.load(struct_path, "s")
            chains = cmd.get_chains("s")

            chain_seqs = {ch: chain_seq("s", ch) for ch in chains}

            # find the best (chain, offset) pair
            best = None  # (chain, offset, identity, checked)
            for ch, pseq in chain_seqs.items():
                if not pseq:
                    continue
                offset, matched, checked, identity = best_offset_identity(pseq, uniprot_seq)
                if checked == 0:
                    continue
                if best is None or identity > best[2]:
                    best = (ch, offset, identity, checked)

            if best is None:
                print("{} / {}: NO PROTEIN CHAIN FOUND".format(ac, pdb_id))
                continue

            ch, offset, identity, checked = best
            # any other chains that match at the SAME offset (e.g. NCS copies, like pheS's A/D)
            all_good_chains = []
            for ch2, pseq2 in chain_seqs.items():
                if not pseq2:
                    continue
                matched2 = sum(1 for resi, aa in pseq2.items()
                                if 1 <= resi + offset <= len(uniprot_seq) and uniprot_seq[resi + offset - 1] == aa)
                checked2 = sum(1 for resi in pseq2 if 1 <= resi + offset <= len(uniprot_seq))
                if checked2 > 0 and matched2 / checked2 > 0.95:
                    all_good_chains.append(ch2)

            status = "OK" if identity > 0.95 else ("PARTIAL" if identity > 0.5 else "MISMATCH")
            flag = "" if offset == 0 else "  [NUMBERING OFFSET={}]".format(offset)
            print("{:6s} {:6s} fmt={:3s} best_chain={:2s} offset={:3d} identity={:.3f} ({}) all_matching_chains={} -> {}{}".format(
                ac, pdb_id, fmt, ch, offset, identity, checked, all_good_chains, status, flag))

            manifest.append({
                "uniprot_ac": ac, "pdb_id": pdb_id, "format": fmt,
                "chains": all_good_chains if all_good_chains else [ch],
                "offset": offset, "identity": round(identity, 4), "status": status,
            })

    out_path = os.path.join(manifest_dir, "chain_manifest.json")
    with open(out_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print("\nSaved manifest ({} entries) -> {}".format(len(manifest), out_path))
