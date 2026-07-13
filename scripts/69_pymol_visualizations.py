#!/usr/bin/env python3
"""
One PyMOL session per gene, for the 5 tRNA synthetases with a curated pocket in
output/selected_pockets.csv (pheS, pheT, aspS, lysS, alaS). Each session merges every curated
pocket's structure + centroid spheres into one object (blue=CAT, teal=NON-CAT), plus its top-N
best-scoring compounds from output/65_aggregated_docking as hidden ligand objects, all shown
identically (orange carbons). Adapted from scripts/51_selected_pockets_visualization.py and
scripts/47b_reference_pocket_visualization.py.

Usage:
    python 69_pymol_visualizations.py [--genes pheS,aspS] [--top-n 5] [--no-surface]
"""
import argparse
import os
import sys
import tarfile
import tempfile

import pandas as pd
import pymol
from pymol import cmd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from docking_utils import load_gene_map

SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")
AGGREGATED_DOCKING_DIR = os.path.join(ROOT, "output", "65_aggregated_docking", "docking_results")
OUTPUT_DIR = os.path.join(ROOT, "output", "69_pymol_visualizations")

ALIGNED_DIR = os.path.join(ROOT, "output", "aligned_relaxed_structures")
DETECTED_POCKETS_DIR = os.path.join(ROOT, "output", "detected_pockets")
MULTIMER_STRUCTURES_DIR = os.path.join(ROOT, "output", "48_detect_pocket_multimers", "stripped_structures")
MULTIMER_POCKETS_DIR = os.path.join(ROOT, "output", "48_detect_pocket_multimers", "detected_pockets")

CURATED_GENES = ["pheS", "pheT", "aspS", "lysS", "alaS"]

COLOR_STRUCTURE = [0.7804, 0.8275, 0.8667]
COLOR_LIGAND_DOCKED = [0xF5 / 255, 0xA6 / 255, 0x3A / 255]  # orange
POCKET_COLOR_CAT = "blue"
POCKET_COLOR_NONCAT = "teal"
POCKET_SPHERE_SCALE = 10


def color_ligand(selection, carbon_color, color_name):
    cmd.util.cbag(selection)
    cmd.set_color(color_name, carbon_color)
    cmd.color(color_name, f"{selection} and elem C")


def resolve_pocket(pocket_name, uniprot_ac):
    """(structure_path, pocket_pdb_path) for pocket_name, or None if either file is missing.
    Branches on "_model_" in pocket_name: monomeric (AF2/AF3/SwissModel/Chai1) vs. multimeric
    (experimental PDB, script 48) - same branching as script 51's resolve_pocket."""
    if "_model_" in pocket_name:
        structure_name, pocket_number = pocket_name.rsplit("_pocket_", 1)
        structure_path = os.path.join(ALIGNED_DIR, uniprot_ac, f"{structure_name}.pdb")
        pocket_pdb_path = os.path.join(DETECTED_POCKETS_DIR, uniprot_ac, structure_name, "pockets", f"pocket_{pocket_number}.pdb")
    else:
        pdb_code, pocket_number = pocket_name.rsplit("_pocket_", 1)
        structure_path = os.path.join(MULTIMER_STRUCTURES_DIR, f"{pdb_code}.pdb")
        pocket_pdb_path = os.path.join(MULTIMER_POCKETS_DIR, pdb_code, "pockets", f"pocket_{pocket_number}.pdb")

    if not os.path.isfile(structure_path) or not os.path.isfile(pocket_pdb_path):
        return None
    return structure_path, pocket_pdb_path


def load_and_merge_structure_pocket(structure_path, pocket_pdb_path, final_name, pocket_color, no_surface=False):
    cmd.load(structure_path, "structure")
    cmd.set_color("structure_color", COLOR_STRUCTURE)
    cmd.color("structure_color", "structure")
    cmd.show("cartoon" if no_surface else "surface", "structure")
    cmd.hide("surface" if no_surface else "cartoon", "structure")
    cmd.hide("lines", "structure")
    cmd.hide("sticks", "structure")
    cmd.set("transparency", 0.3, "structure")

    cmd.load(pocket_pdb_path, "pocket")
    cmd.color(pocket_color, "pocket")
    cmd.show("spheres", "pocket")
    cmd.set("sphere_transparency", 0.4, "pocket")
    cmd.set("sphere_scale", POCKET_SPHERE_SCALE, "pocket")

    cmd.create("_tmp_merged", "structure or pocket")
    cmd.delete("structure")
    cmd.delete("pocket")
    cmd.set_name("_tmp_merged", final_name)


def get_top_n_scores(pocket_name, top_n):
    """Top-N (compound_id, score) pairs from output/65_aggregated_docking's report.csv, best first."""
    report_path = os.path.join(AGGREGATED_DOCKING_DIR, pocket_name, "report.csv")
    if not os.path.isfile(report_path):
        return []
    report = pd.read_csv(report_path).dropna(subset=["score"]).sort_values("score", ascending=True)
    return list(report[["compound", "score"]].head(top_n).itertuples(index=False, name=None))


def load_ligand_objects(top_scores, obj_prefix, docking_dir):
    """Loads each top_scores compound as "<obj_prefix>_top<rank>", hidden by default.
    Returns the list of ranks successfully extracted from docking.tar.gz."""
    tar_path = os.path.join(docking_dir, "docking.tar.gz")
    if not os.path.isfile(tar_path):
        print(f"    Warning: {tar_path} not found, no ligand poses loaded for this pocket.")
        return []

    wanted_names = {f"docking/{cid}_out.sdf": rank for rank, (cid, _) in enumerate(top_scores, start=1)}
    sdf_bytes_by_rank = {}
    with tarfile.open(tar_path, "r:gz") as tf:
        for member in tf:
            if member.name in wanted_names:
                sdf_bytes_by_rank[wanted_names[member.name]] = tf.extractfile(member).read()
                if len(sdf_bytes_by_rank) == len(wanted_names):
                    break

    loaded = []
    for rank, sdf_bytes in sorted(sdf_bytes_by_rank.items()):
        obj_name = f"{obj_prefix}_top{rank}"
        with tempfile.NamedTemporaryFile(suffix=".sdf", mode="wb", delete=False) as tmp:
            tmp.write(sdf_bytes)
            tmp_path = tmp.name
        cmd.load(tmp_path, obj_name)
        os.remove(tmp_path)
        cmd.show("sticks", obj_name)
        cmd.hide("lines", obj_name)
        color_ligand(obj_name, COLOR_LIGAND_DOCKED, "ligC_docked")
        cmd.disable(obj_name)
        loaded.append(obj_name)

    if len(loaded) < len(top_scores):
        print(f"    Warning: only found {len(loaded)}/{len(top_scores)} top ligand poses in {tar_path}.")
    return loaded


def build_gene_session(gene, uniprot_ac, gene_pockets, top_n, no_surface):
    pymol.finish_launching(["pymol", "-cq"])
    cmd.reinitialize()

    order_names = []
    for _, row in gene_pockets.sort_values(["site_type", "pocket_name"]).iterrows():
        pocket_name, site_type, comment = row["pocket_name"], row["site_type"], row["comment"]
        print(f"  Pocket: {pocket_name} ({site_type}) - {comment}")

        resolved = resolve_pocket(pocket_name, uniprot_ac)
        if resolved is None:
            print(f"    Warning: structure or pocket centroid file not found for {pocket_name}, skipping.")
            continue
        structure_path, pocket_pdb_path = resolved

        final_name = f"{gene}_{site_type}_{pocket_name}"
        pocket_color = POCKET_COLOR_CAT if site_type == "CAT" else POCKET_COLOR_NONCAT
        load_and_merge_structure_pocket(structure_path, pocket_pdb_path, final_name, pocket_color, no_surface)
        order_names.append(final_name)

        top_scores = get_top_n_scores(pocket_name, top_n)
        if not top_scores:
            print(f"    Warning: no docking report found for {pocket_name}, no ligands loaded.")
            continue

        docking_dir = os.path.join(AGGREGATED_DOCKING_DIR, pocket_name)
        ligand_objs = load_ligand_objects(top_scores, final_name, docking_dir)
        order_names.extend(ligand_objs)
        print(f"    Loaded {len(ligand_objs)} top ligand(s).")

    if order_names:
        cmd.order(" ".join(order_names), location="top", sort=0)

    # No zoom call: a gene's pockets usually come from different, mutually unaligned structures
    # (e.g. pheS spans a SwissModel model, an AlphaFold3 model, and the 7K98 PDB structure), so a
    # single fixed-box zoom across "all" wouldn't be meaningful - inspect one pocket at a time.
    cmd.bg_color("white")
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_trace_fog", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_trace_gain", 0.02)
    cmd.set("spec_reflect", 0)
    cmd.set("specular", 0)
    cmd.set("antialias", 2)

    out_path = os.path.join(OUTPUT_DIR, f"session_{uniprot_ac}_{gene}.pse")
    cmd.save(out_path)
    print(f"  Saved: {out_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genes", default=",".join(CURATED_GENES),
                         help=f"Comma-separated gene name(s) with a curated pocket (default: {','.join(CURATED_GENES)})")
    parser.add_argument("--top-n", type=int, default=5, help="Top-scoring ligands to load per pocket (default: 5)")
    parser.add_argument("--no-surface", action="store_true", help="Cartoon instead of surface (faster)")
    args = parser.parse_args()

    genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    selected = pd.read_csv(SELECTED_POCKETS_CSV)
    curated_genes = sorted(set(selected["gene_name"]))
    invalid = [g for g in genes if g not in curated_genes]
    if invalid:
        print(f"Error: no curated pocket for gene(s): {', '.join(invalid)}. Genes with a curated pocket:")
        print("  " + ", ".join(curated_genes))
        sys.exit(1)

    gene_map = load_gene_map()
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for gene in genes:
        uniprot_ac = gene_map[gene]
        gene_pockets = selected[selected["gene_name"] == gene]
        print(f"\n--- Gene: {gene} ({uniprot_ac}), {len(gene_pockets)} curated pocket(s) ---")
        build_gene_session(gene, uniprot_ac, gene_pockets, args.top_n, args.no_surface)


if __name__ == "__main__":
    main()
