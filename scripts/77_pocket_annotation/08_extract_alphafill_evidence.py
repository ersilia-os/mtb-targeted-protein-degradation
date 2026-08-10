import os
import sys
import json
import numpy as np
import pandas as pd
from pymol import cmd

root = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.join(root, "..", "..")
output_dir = os.path.join(repo_root, "output", "77_pocket_annotation")
os.makedirs(output_dir, exist_ok=True)
alphafill_dir = os.path.join(repo_root, "data", "structures", "alphafill_database")
sys.path.append(root)
from ligand_classification import classify_ligand  # noqa: E402

PROXIMITY_CUTOFF = 10.0
MIN_POCKET_ATOMS_FOR_FIT = 5


def parse_pocket_residues(residues_str):
    residues = set()
    for tok in residues_str.split():
        chain, resn = tok.split("_")
        residues.add(int(resn))
    return residues


def pocket_local_align(mobile_obj, chain, ref_obj, pocket_residues):
    resi_str = "+".join(str(r) for r in sorted(pocket_residues))
    mobile_sel = "{} and chain {} and resi {} and name CA and polymer.protein".format(mobile_obj, chain, resi_str)
    ref_sel = "{} and resi {} and name CA and polymer.protein".format(ref_obj, resi_str)
    n_mobile = cmd.select("tmp_check", mobile_sel)
    cmd.delete("tmp_check")
    if n_mobile < MIN_POCKET_ATOMS_FOR_FIT:
        return None
    r = cmd.align(mobile_sel, ref_sel, cycles=0)
    return r[0], r[1]


if __name__ == "__main__":
    bosch = pd.read_csv(os.path.join(repo_root, "data", "mtb_trna_synthetases_bosch_2021_fig5_annotated.csv"))
    pockets_df = pd.read_csv(os.path.join(repo_root, "output", "pocket_detection_data.csv"))

    all_evidence_rows = []
    fit_log = []

    for _, brow in bosch.iterrows():
        ac = brow["uniprot_ac"]
        cif_path = os.path.join(alphafill_dir, ac, "{}.cif".format(ac))
        json_path = os.path.join(alphafill_dir, ac, "{}.json".format(ac))
        ref_path = os.path.join(repo_root, "output", "aligned_relaxed_structures", ac, "alphafold2_{}_model_0.pdb".format(ac))
        if not (os.path.exists(cif_path) and os.path.exists(json_path) and os.path.exists(ref_path)):
            print("WARNING: missing AlphaFill/reference data for {}, skipping".format(ac))
            continue

        meta = json.load(open(json_path))
        # build transplant lookup: (compound_id) -> list of (pdb_id, local_rmsd, identity)
        # AlphaFill assigns each transplanted ligand its own chain in the CIF;
        # we match by residue/compound identity + proximity, then look up
        # provenance from the JSON's hits[].transplants[] by compound_id.
        transplant_meta = {}
        for hit in meta.get("hits", []):
            pdb_id = hit.get("pdb_id")
            identity = hit.get("alignment", {}).get("identity")
            for t in (hit.get("transplants") or []):
                compound_id = t.get("compound_id") or t.get("analogue_id")
                transplant_meta.setdefault(compound_id, []).append({
                    "pdb_id": pdb_id, "local_rmsd": t.get("local_rmsd"), "identity": identity,
                })

        protein_pockets = pockets_df[pockets_df["Uniprot AC"] == ac]

        cmd.reinitialize()
        cmd.load(ref_path, "ref")
        cmd.load(cif_path, "af")

        n_fit_ok, n_fit_skip = 0, 0
        for _, prow in protein_pockets.iterrows():
            pocket_residues = parse_pocket_residues(prow["Pocket residues (chain_resn)"])
            fit = pocket_local_align("af", "A", "ref", pocket_residues)
            if fit is None:
                n_fit_skip += 1
                continue
            n_fit_ok += 1
            fit_rmsd, fit_n = fit
            fit_log.append({
                "uniprot_ac": ac, "pocket_file_name": prow["File name"],
                "pocket_number": int(prow["Pocket number"]), "rmsd": round(fit_rmsd, 2), "n_atoms": fit_n,
            })

            # re-fetch ligand coordinates AFTER this pocket's own alignment -- cmd.align moves
            # "af"'s real atoms in place, so ligand positions must be read fresh every iteration
            # (caching them once before the loop, as this used to do, compares each pocket's
            # centroid against a stale, essentially-arbitrary ligand position from whichever
            # alignment happened to be active when the cache was built)
            lig_model = cmd.get_model("af and hetatm and not resn HOH")
            ligs = {}
            for atom in lig_model.atom:
                key = (atom.resn, atom.chain, atom.resi)
                ligs.setdefault(key, []).append(atom.coord)

            pocket_centroid = np.array([float(x) for x in prow["Pocket centroid coordinate (x y z)"].split()])
            for (resn, ch, resi), coords in ligs.items():
                centroid = np.mean(np.array(coords), axis=0)
                d = np.linalg.norm(centroid - pocket_centroid)
                if d <= PROXIMITY_CUTOFF:
                    provenance = transplant_meta.get(resn, [{}])[0]
                    all_evidence_rows.append({
                        "uniprot_ac": ac, "ligand_resn": resn, "ligand_chain": ch, "ligand_resi": resi,
                        "pocket_file_name": prow["File name"], "pocket_number": int(prow["Pocket number"]),
                        "distance": round(d, 2), "fit_rmsd": round(fit_rmsd, 2),
                        "source_pdb": provenance.get("pdb_id"), "local_rmsd": provenance.get("local_rmsd"),
                        "homolog_identity": provenance.get("identity"),
                        "strength": classify_ligand(resn, ac),
                    })

        print("{}: {} pockets fit OK, {} skipped".format(ac, n_fit_ok, n_fit_skip))

    evidence_df = pd.DataFrame(all_evidence_rows)
    evidence_path = os.path.join(output_dir, "alphafill_ligand_evidence.csv")
    evidence_df.to_csv(evidence_path, index=False)
    print("\nSaved {} AlphaFill ligand-pocket evidence rows -> {}".format(len(evidence_df), evidence_path))
