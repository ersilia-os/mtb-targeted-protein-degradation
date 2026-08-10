import os
import sys
import json
import numpy as np
import pandas as pd
from pymol import cmd

root = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.join(root, "..", "..")
pdb_dir = os.path.join(repo_root, "processed", "pocket_annotation", "pdb_structures")
manifest_dir = os.path.join(repo_root, "processed", "pocket_annotation")
output_dir = os.path.join(repo_root, "output", "77_pocket_annotation")
os.makedirs(output_dir, exist_ok=True)
sys.path.append(root)
from ligand_classification import classify_ligand  # noqa: E402

PROXIMITY_CUTOFF = 10.0  # Angstrom, ligand-to-pocket-centroid
MIN_POCKET_ATOMS_FOR_FIT = 5


def renumber_to_uniprot(obj_name, chain, offset):
    if offset != 0:
        cmd.alter("{} and chain {}".format(obj_name, chain), "resi=str(int(resi)+{})".format(offset))
        cmd.sort()


def parse_pocket_residues(residues_str):
    residues = set()
    for tok in residues_str.split():
        chain, resn = tok.split("_")
        residues.add(int(resn))
    return residues


def pocket_local_align(mobile_obj, chain, ref_obj, pocket_residues):
    """Align mobile chain onto ref using only the pocket's own CA residues
    (rigid local fit) -- avoids needing a single global alignment to hold for
    the whole (possibly multi-domain, independently-flexible) protein. See
    scripts/77_pocket_annotation/README.md for why this replaced a
    whole-chain cmd.align (which broke badly, 10-14 A RMSD, whenever a
    protein has an independently-flexible domain)."""
    resi_str = "+".join(str(r) for r in sorted(pocket_residues))
    mobile_sel = "{} and chain {} and resi {} and name CA and polymer.protein".format(mobile_obj, chain, resi_str)
    ref_sel = "{} and resi {} and name CA and polymer.protein".format(ref_obj, resi_str)
    n_mobile = cmd.select("tmp_check", mobile_sel)
    cmd.delete("tmp_check")
    if n_mobile < MIN_POCKET_ATOMS_FOR_FIT:
        return None
    r = cmd.align(mobile_sel, ref_sel, cycles=0)  # cycles=0: no outlier rejection, use all pocket residues as-is
    return r[0], r[1]


if __name__ == "__main__":
    manifest = json.load(open(os.path.join(manifest_dir, "chain_manifest.json")))
    pockets_df = pd.read_csv(os.path.join(repo_root, "output", "pocket_detection_data.csv"))

    by_protein = {}
    for entry in manifest:
        by_protein.setdefault(entry["uniprot_ac"], []).append(entry)

    all_evidence_rows = []
    fit_log = []

    for ac, entries in by_protein.items():
        ref_path = os.path.join(repo_root, "output", "aligned_relaxed_structures", ac, "alphafold2_{}_model_0.pdb".format(ac))
        if not os.path.exists(ref_path):
            print("WARNING: no AF2 reference for {}, skipping".format(ac))
            continue

        protein_pockets = pockets_df[pockets_df["Uniprot AC"] == ac]

        for entry in entries:
            pdb_id = entry["pdb_id"]
            struct_path = os.path.join(pdb_dir, "{}.{}".format(pdb_id, entry["format"]))

            cmd.reinitialize()
            cmd.load(ref_path, "ref")
            cmd.load(struct_path, "mob")

            offset = entry["offset"]
            for ch in entry["chains"]:
                renumber_to_uniprot("mob", ch, offset)
            main_chain = entry["chains"][0]

            n_fit_ok = 0
            n_fit_skip = 0
            for _, prow in protein_pockets.iterrows():
                pocket_residues = parse_pocket_residues(prow["Pocket residues (chain_resn)"])
                fit = pocket_local_align("mob", main_chain, "ref", pocket_residues)
                if fit is None:
                    n_fit_skip += 1
                    continue
                n_fit_ok += 1
                fit_rmsd, fit_n = fit
                fit_log.append({
                    "uniprot_ac": ac, "pdb_id": pdb_id, "pocket_file_name": prow["File name"],
                    "pocket_number": int(prow["Pocket number"]), "rmsd": round(fit_rmsd, 2), "n_atoms": fit_n,
                })

                pocket_centroid = np.array([float(x) for x in prow["Pocket centroid coordinate (x y z)"].split()])
                lig_sel = "mob and hetatm and not resn HOH and chain {}".format("+".join(entry["chains"]))
                lig_model = cmd.get_model(lig_sel)
                ligs = {}
                for atom in lig_model.atom:
                    key = (atom.resn, atom.chain, atom.resi)
                    ligs.setdefault(key, []).append(atom.coord)

                for (resn, ch, resi), coords in ligs.items():
                    centroid = np.mean(np.array(coords), axis=0)
                    d = np.linalg.norm(centroid - pocket_centroid)
                    if d <= PROXIMITY_CUTOFF:
                        all_evidence_rows.append({
                            "uniprot_ac": ac, "pdb_id": pdb_id, "ligand_resn": resn,
                            "ligand_chain": ch, "ligand_resi": resi,
                            "pocket_file_name": prow["File name"], "pocket_number": int(prow["Pocket number"]),
                            "distance": round(d, 2), "fit_rmsd": round(fit_rmsd, 2),
                            "strength": classify_ligand(resn, ac),
                        })

            print("{} / {}: {} pockets fit OK, {} skipped (too few resolved pocket residues)".format(
                ac, pdb_id, n_fit_ok, n_fit_skip))

    evidence_df = pd.DataFrame(all_evidence_rows)
    evidence_path = os.path.join(output_dir, "direct_pdb_ligand_evidence.csv")
    evidence_df.to_csv(evidence_path, index=False)
    print("\nSaved {} ligand-pocket evidence rows -> {}".format(len(evidence_df), evidence_path))

    fit_log_df = pd.DataFrame(fit_log)
    fit_log_path = os.path.join(output_dir, "pocket_fit_log.csv")
    fit_log_df.to_csv(fit_log_path, index=False)
    print("Saved {} pocket-local fit records -> {}".format(len(fit_log_df), fit_log_path))
    print("\nFit RMSD distribution: mean={:.2f} median={:.2f} max={:.2f} (n>2A: {})".format(
        fit_log_df["rmsd"].mean(), fit_log_df["rmsd"].median(), fit_log_df["rmsd"].max(),
        (fit_log_df["rmsd"] > 2.0).sum()))
