import os
import sys
import re
import json
import numpy as np
import pandas as pd
from pymol import cmd

root = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.join(root, "..", "..")
pdb_dir = os.path.join(repo_root, "processed", "pocket_annotation", "pdb_structures")
manifest_dir = os.path.join(repo_root, "processed", "pocket_annotation")
output_dir = os.path.join(repo_root, "output", "77_pocket_annotation")
sessions_dir = os.path.join(repo_root, "output", "77_pymol_sessions")
os.makedirs(sessions_dir, exist_ok=True)
detected_pockets_dir = os.path.join(repo_root, "output", "detected_pockets")
alphafill_dir = os.path.join(repo_root, "data", "structures", "alphafill_database")

sys.path.append(root)
from ligand_classification import classify_ligand  # noqa: E402

COLOR_BACKGROUND = "grey"
CATEGORY_COLOR = {
    "Catalytic Domain (ATP/ligase)": "orange",
    "Anticodon Binding Domain": "blue",
    "Editing Domain": "forest",
    "tRNA Binding Domain": "magenta",
    "RNA Binding Domain": "yellow",
    "Non-catalytic (structural homolog)": "purple",
}
CLUSTER_COLORS = ["red", "cyan", "black", "brown", "teal", "deepsalmon", "lightblue", "olive", "violet", "sand", "pink", "grey70"]
COLOR_EXPERIMENTAL_PROTEIN = "wheat"
COLOR_EXPERIMENTAL_TRNA = "palecyan"
COLOR_LIGAND_WEAK = "green"
COLOR_LIGAND_STRONG = "red"

MIN_POCKET_ATOMS_FOR_FIT = 5


def gene_map():
    bosch = pd.read_csv(os.path.join(repo_root, "data", "mtb_trna_synthetases_bosch_2021_fig5_annotated.csv"))
    return dict(zip(bosch["uniprot_ac"], bosch["gene_name_in_bosch_2021"]))


def safe_object_name(*parts):
    return re.sub(r"[^0-9a-zA-Z]+", "_", "_".join(str(p) for p in parts)).strip("_")


def parse_ranges(ranges_str):
    residues = set()
    for part in ranges_str.split(","):
        s, e = part.split("-")
        residues.update(range(int(s), int(e) + 1))
    return residues


def parse_pocket_residues(residues_str):
    residues = set()
    for tok in residues_str.split():
        _, resn = tok.split("_")
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


def build_category_objects(gene, ac, ref_path):
    cat_path = os.path.join(output_dir, "{}_annotation_table_categorized.csv".format(ac))
    if not os.path.exists(cat_path):
        print("  no categorized table for {}, skipping category objects".format(ac))
        return
    df = pd.read_csv(cat_path, keep_default_na=False)
    category_residues = {}
    for _, row in df.iterrows():
        if row["category"] not in CATEGORY_COLOR:
            continue
        category_residues.setdefault(row["category"], set()).update(parse_ranges(row["ranges"]))

    for category, residues in category_residues.items():
        obj_name = safe_object_name(gene, category)
        cmd.load(ref_path, obj_name)
        cmd.hide("everything", obj_name)
        cmd.show("cartoon", obj_name)
        cmd.show("surface", obj_name)
        cmd.color(COLOR_BACKGROUND, obj_name)
        resi_str = "+".join(str(r) for r in sorted(residues))
        cmd.select("cat_sel", "{} and chain A and resi {}".format(obj_name, resi_str))
        cmd.color(CATEGORY_COLOR[category], "cat_sel")
        cmd.delete("cat_sel")
        print("  category: {} -> {} residues".format(category, len(residues)))


def build_pocket_objects(gene, ac, pockets_df, clusters_df):
    n = 0
    for _, prow in pockets_df.iterrows():
        structure_name = prow["File name"].replace(".pdb", "")
        centroid_path = os.path.join(detected_pockets_dir, ac, structure_name, "pockets", "pocket_{}.pdb".format(prow["Pocket number"]))
        if not os.path.exists(centroid_path):
            continue
        crow = clusters_df[(clusters_df["File name"] == prow["File name"]) & (clusters_df["Pocket number"] == prow["Pocket number"])]
        cluster_id = int(crow["spatial_cluster_id"].iloc[0]) if len(crow) else 0
        color = CLUSTER_COLORS[(cluster_id - 1) % len(CLUSTER_COLORS)]

        obj_name = safe_object_name(gene, "Pocket", structure_name, prow["Pocket number"], "cluster", cluster_id)
        cmd.load(centroid_path, obj_name)
        cmd.show("spheres", obj_name)
        cmd.set("sphere_scale", 1.5, obj_name)
        cmd.color(color, obj_name)
        n += 1
    print("  pockets: {} objects".format(n))


def build_experimental_objects(gene, ac, ref_path, manifest_by_ac, fit_log):
    entries = manifest_by_ac.get(ac, [])
    n = 0
    for entry in entries:
        pdb_id = entry["pdb_id"]
        struct_path = os.path.join(pdb_dir, "{}.{}".format(pdb_id, entry["format"]))
        if not os.path.exists(struct_path):
            continue

        obj_name = safe_object_name(gene, "PDB", pdb_id)
        cmd.load(struct_path, obj_name)
        offset = entry["offset"]
        for ch in entry["chains"]:
            if offset != 0:
                cmd.alter("{} and chain {}".format(obj_name, ch), "resi=str(int(resi)+{})".format(offset))
        if offset != 0:
            cmd.sort()
        main_chain = entry["chains"][0]

        # anchor alignment on whichever pocket this structure fits best (lowest RMSD)
        sub_log = fit_log[(fit_log["uniprot_ac"] == ac) & (fit_log["pdb_id"] == pdb_id)]
        if sub_log.empty:
            cmd.delete(obj_name)
            continue
        best = sub_log.loc[sub_log["rmsd"].idxmin()]
        pockets_row = pd.read_csv(os.path.join(repo_root, "output", "pocket_detection_data.csv"))
        prow = pockets_row[(pockets_row["Uniprot AC"] == ac) & (pockets_row["File name"] == best["pocket_file_name"]) & (pockets_row["Pocket number"] == best["pocket_number"])]
        if prow.empty:
            cmd.delete(obj_name)
            continue
        anchor_residues = parse_pocket_residues(prow.iloc[0]["Pocket residues (chain_resn)"])
        fit = pocket_local_align(obj_name, main_chain, "ref", anchor_residues)
        if fit is None:
            cmd.delete(obj_name)
            continue

        cmd.hide("everything", obj_name)
        cmd.show("cartoon", "{} and polymer.protein".format(obj_name))
        cmd.color(COLOR_EXPERIMENTAL_PROTEIN, "{} and polymer.protein".format(obj_name))
        cmd.show("cartoon", "{} and polymer.nucleic".format(obj_name))
        cmd.color(COLOR_EXPERIMENTAL_TRNA, "{} and polymer.nucleic".format(obj_name))

        lig_sel = "{} and hetatm and not resn HOH and chain {}".format(obj_name, "+".join(entry["chains"]))
        lig_model = cmd.get_model(lig_sel)
        seen = set()
        for atom in lig_model.atom:
            seen.add((atom.resn, atom.chain, atom.resi))
        for resn, ch, resi in seen:
            strength = classify_ligand(resn, ac)
            sel = "{} and chain {} and resi {} and resn {}".format(obj_name, ch, resi, resn)
            cmd.show("sticks", sel)
            color = COLOR_LIGAND_STRONG if strength == "strong" else COLOR_LIGAND_WEAK
            cmd.color(color, "{} and elem C".format(sel))
        n += 1
    print("  experimental PDB structures: {} objects".format(n))


def build_alphafill_object(gene, ac, ref_path, pockets_df):
    cif_path = os.path.join(alphafill_dir, ac, "{}.cif".format(ac))
    if not os.path.exists(cif_path):
        return
    obj_name = safe_object_name(gene, "AlphaFill")
    cmd.load(cif_path, obj_name)

    # anchor on whichever pocket fits best: try each pocket, keep the lowest RMSD
    best_fit = None
    best_pocket_row = None
    for _, prow in pockets_df.iterrows():
        pocket_residues = parse_pocket_residues(prow["Pocket residues (chain_resn)"])
        fit = pocket_local_align(obj_name, "A", "ref", pocket_residues)
        if fit is None:
            continue
        if best_fit is None or fit[0] < best_fit[0]:
            best_fit = fit
            best_pocket_row = prow
    if best_fit is None:
        cmd.delete(obj_name)
        return
    pocket_local_align(obj_name, "A", "ref", parse_pocket_residues(best_pocket_row["Pocket residues (chain_resn)"]))

    cmd.hide("everything", obj_name)
    cmd.show("cartoon", "{} and chain A and polymer.protein".format(obj_name))
    cmd.color(COLOR_EXPERIMENTAL_PROTEIN, "{} and chain A and polymer.protein".format(obj_name))

    lig_model = cmd.get_model("{} and hetatm and not resn HOH".format(obj_name))
    seen = set()
    for atom in lig_model.atom:
        seen.add((atom.resn, atom.chain, atom.resi))
    for resn, ch, resi in seen:
        strength = classify_ligand(resn, ac)
        sel = "{} and chain {} and resi {} and resn {}".format(obj_name, ch, resi, resn)
        cmd.show("sticks", sel)
        color = COLOR_LIGAND_STRONG if strength == "strong" else COLOR_LIGAND_WEAK
        cmd.color(color, "{} and elem C".format(sel))
    print("  alphafill: 1 object, anchored on pocket {} (fit RMSD {:.2f})".format(best_pocket_row["Pocket number"], best_fit[0]))


if __name__ == "__main__":
    genes = gene_map()
    pockets_all = pd.read_csv(os.path.join(repo_root, "output", "pocket_detection_data.csv"))
    clusters_all = pd.read_csv(os.path.join(output_dir, "pocket_clusters.csv"), keep_default_na=False)
    fit_log = pd.read_csv(os.path.join(output_dir, "pocket_fit_log.csv"), keep_default_na=False)

    manifest = json.load(open(os.path.join(manifest_dir, "chain_manifest.json")))
    manifest_by_ac = {}
    for entry in manifest:
        manifest_by_ac.setdefault(entry["uniprot_ac"], []).append(entry)

    target_acs = sys.argv[1:] if len(sys.argv) > 1 else list(genes.keys())

    for ac in target_acs:
        gene = genes.get(ac, ac)
        ref_path = os.path.join(repo_root, "output", "aligned_relaxed_structures", ac, "alphafold2_{}_model_0.pdb".format(ac))
        if not os.path.exists(ref_path):
            print("WARNING: no AF2 reference for {} ({}), skipping".format(gene, ac))
            continue

        print("=== {} ({}) ===".format(gene, ac))
        cmd.reinitialize()
        cmd.load(ref_path, "ref")
        cmd.hide("everything", "ref")

        build_category_objects(gene, ac, ref_path)

        pockets_df = pockets_all[pockets_all["Uniprot AC"] == ac]
        clusters_df = clusters_all[clusters_all["Uniprot AC"] == ac]
        build_pocket_objects(gene, ac, pockets_df, clusters_df)

        build_experimental_objects(gene, ac, ref_path, manifest_by_ac, fit_log)
        build_alphafill_object(gene, ac, ref_path, pockets_df)

        cmd.delete("ref")
        cmd.reset()
        out_path = os.path.join(sessions_dir, "{}_pocket_annotation.pse".format(gene))
        cmd.save(out_path)
        print("  saved -> {}\n".format(out_path))
