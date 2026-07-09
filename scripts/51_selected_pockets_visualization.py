#!/usr/bin/env python3
"""
Builds a single PyMOL session for all 12 manually-curated pockets in output/selected_pockets.csv
(columns: gene_name, site_type, pocket_name, comment). For each pocket, adds 1+N objects to the
same session:
  - "<gene_name>_<site_type>_<pocket_name>": the structure (surface representation, matching
    notebooks/46_pocket_visualisation.ipynb's convention) merged with the pocket centroid
    (spheres) into one object - same merge pattern as scripts/47_docking_summary.py's
    build_pocket_overview_session / scripts/47b_reference_pocket_visualization.py's
    _load_structure_and_pocket + _merge_structure_and_pocket, adapted from cartoon to surface.
  - "<gene_name>_<site_type>_<pocket_name>_top<rank>" (rank 1..N): the top-N best-scoring docked
    poses from the second Enamine REAL screening (scripts 44-46), loaded but disabled (hidden) by
    default - adapted from 47b's load_top_ligands, which shows them by default.

11 of the 12 pockets follow the standard monomeric convention (resolved via
output/aligned_relaxed_structures/, output/detected_pockets/, output/pocket_detection_data.csv,
output/unidock_REAL_docking_2/docking_results/). Any pocket_name NOT containing "_model_" is
treated as a multimeric PDB-code pocket instead (e.g. "7K98_pocket_6"), resolved via
output/48_detect_pocket_multimers/ and output/50_unidock_docking_multimers/ (scripts 48-50).

Usage:
    python 51_selected_pockets_visualization.py --top-n 5
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
from docking_utils import LIBRARIES, load_gene_map

SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")
OUTPUT_DIR = os.path.join(ROOT, "output", "51_selected_pockets_visualization")
OUTPUT_PSE = os.path.join(OUTPUT_DIR, "session_selected_pockets.pse")

# Monomeric (AF2/AF3/SwissModel/Chai1) pipeline paths - scripts 04-08, 44-46
ALIGNED_DIR = os.path.join(ROOT, "output", "aligned_relaxed_structures")
DETECTED_POCKETS_DIR = os.path.join(ROOT, "output", "detected_pockets")
POCKET_DETECTION_CSV = os.path.join(ROOT, "output", "pocket_detection_data.csv")
REAL_DOCKING_DIR = LIBRARIES["REAL"]  # output/unidock_REAL_docking_2/docking_results

# Multimeric (experimental PDB) pipeline paths - scripts 48-50
MULTIMER_STRUCTURES_DIR = os.path.join(ROOT, "output", "48_detect_pocket_multimers", "stripped_structures")
MULTIMER_POCKETS_DIR = os.path.join(ROOT, "output", "48_detect_pocket_multimers", "detected_pockets")
MULTIMER_POCKET_CSV = os.path.join(ROOT, "output", "48_detect_pocket_multimers", "pocket_detection_data.csv")
MULTIMER_DOCKING_DIR = os.path.join(ROOT, "output", "50_unidock_docking_multimers")

# Same structure coloring convention as notebooks/46_pocket_visualisation.ipynb and
# notebooks/46_protein_visualization.ipynb (verbatim from scripts/47b_reference_pocket_visualization.py).
COLOR_STRUCTURE = [0.7804, 0.8275, 0.8667]
COLOR_LIGAND_DOCKED = [0xF5 / 255, 0xA6 / 255, 0x3A / 255]  # orange
POCKET_SPHERE_COLOR = "blue"
POCKET_SPHERE_SCALE = 10


def color_ligand(selection, carbon_color, color_name):
    """Color-by-element with a custom carbon color. Verbatim from scripts/47b_reference_pocket_visualization.py."""
    cmd.util.cbag(selection)
    cmd.set_color(color_name, carbon_color)
    cmd.color(color_name, f"{selection} and elem C")


def resolve_pocket(gene_name, pocket_name, gene_map):
    """
    Resolves (structure_path, pocket_centroid_path, probability, docking_dir) for a pocket_name,
    branching on the monomeric ("_model_" in pocket_name) vs multimeric convention. Returns None
    if the structure or pocket centroid file can't be found (a hard requirement for the merged
    object); docking_dir is returned regardless of whether it actually contains report.csv/
    docking.tar.gz - that's checked separately, since missing docking data should not prevent the
    structure+pocket object itself from being added.
    """
    if "_model_" in pocket_name:
        structure_name, pocket_number_str = pocket_name.rsplit("_pocket_", 1)
        pocket_number = int(pocket_number_str)
        uniprot_ac = gene_map[gene_name]

        structure_path = os.path.join(ALIGNED_DIR, uniprot_ac, f"{structure_name}.pdb")
        pocket_pdb_path = os.path.join(DETECTED_POCKETS_DIR, uniprot_ac, structure_name, "pockets", f"pocket_{pocket_number}.pdb")
        if not os.path.isfile(structure_path) or not os.path.isfile(pocket_pdb_path):
            return None

        pocket_data = pd.read_csv(POCKET_DETECTION_CSV)
        row = pocket_data[
            (pocket_data["Uniprot AC"] == uniprot_ac)
            & (pocket_data["File name"] == f"{structure_name}.pdb")
            & (pocket_data["Pocket number"] == pocket_number)
        ]
        probability = float(row.iloc[0]["Pocket probability"]) if len(row) else None

        docking_dir = os.path.join(REAL_DOCKING_DIR, pocket_name)
    else:
        pdb_code, pocket_number_str = pocket_name.rsplit("_pocket_", 1)
        pocket_number = int(pocket_number_str)

        structure_path = os.path.join(MULTIMER_STRUCTURES_DIR, f"{pdb_code}.pdb")
        pocket_pdb_path = os.path.join(MULTIMER_POCKETS_DIR, pdb_code, "pockets", f"pocket_{pocket_number}.pdb")
        if not os.path.isfile(structure_path) or not os.path.isfile(pocket_pdb_path):
            return None

        pocket_data = pd.read_csv(MULTIMER_POCKET_CSV)
        row = pocket_data[
            (pocket_data["PDB code"] == pdb_code)
            & (pocket_data["Pocket number"] == pocket_number)
        ]
        probability = float(row.iloc[0]["Pocket probability"]) if len(row) else None

        docking_dir = os.path.join(MULTIMER_DOCKING_DIR, pocket_name)

    return structure_path, pocket_pdb_path, probability, docking_dir


def load_and_merge_structure_pocket(structure_path, pocket_pdb_path, final_name, no_surface=False):
    """
    Loads structure_path as object "structure" (surface representation by default, matching
    notebooks/46_pocket_visualisation.ipynb's convention, or cartoon if no_surface=True - surface
    calculation is expensive, so --no-surface lets it be skipped) and pocket_pdb_path as object
    "pocket" (spheres), then merges them into a single object named final_name. Adapted from
    scripts/47b_reference_pocket_visualization.py's _load_structure_and_pocket +
    _merge_structure_and_pocket, with cartoon swapped for surface.
    """
    cmd.load(structure_path, "structure")
    cmd.set_color("structure_color", COLOR_STRUCTURE)
    cmd.color("structure_color", "structure")
    if no_surface:
        cmd.show("cartoon", "structure")
        cmd.hide("surface", "structure")
    else:
        cmd.show("surface", "structure")
        cmd.hide("cartoon", "structure")
    cmd.hide("lines", "structure")
    cmd.hide("sticks", "structure")
    cmd.set("transparency", 0.3, "structure")

    cmd.load(pocket_pdb_path, "pocket")
    cmd.color(POCKET_SPHERE_COLOR, "pocket")
    cmd.show("spheres", "pocket")
    cmd.set("sphere_transparency", 0.4, "pocket")
    cmd.set("sphere_scale", POCKET_SPHERE_SCALE, "pocket")

    tmp_merged = "_tmp_merged"
    cmd.create(tmp_merged, "structure or pocket")
    cmd.delete("structure")
    cmd.delete("pocket")
    cmd.set_name(tmp_merged, final_name)


def get_top_n_scores(docking_dir, top_n):
    """Returns the top_n (compound_id, score) pairs from docking_dir/report.csv, sorted by score
    ascending (best first), or [] if report.csv doesn't exist."""
    report_path = os.path.join(docking_dir, "report.csv")
    if not os.path.isfile(report_path):
        return []
    report = pd.read_csv(report_path).dropna(subset=["score"]).sort_values("score", ascending=True)
    return list(report[["compound", "score"]].head(top_n).itertuples(index=False, name=None))


def load_ligand_objects(docking_dir, top_ids, obj_prefix):
    """
    Extracts docking/<compound_id>_out.sdf for each id in top_ids from docking_dir/docking.tar.gz,
    loads each as its own object named "<obj_prefix>_top<rank>" (rank 1 = best), colored with
    orange carbons, and disabled (hidden) by default. Adapted from scripts/47b_reference_pocket_
    visualization.py's load_top_ligands, with two changes: cmd.disable() so ligands start hidden
    rather than shown, and objects are created strictly in rank order (1, 2, 3, ...) regardless of
    the tar's internal member order, so top1..top5 always appear in that order in the object panel.
    Returns the list of loaded object names.
    """
    tar_path = os.path.join(docking_dir, "docking.tar.gz")
    if not os.path.isfile(tar_path):
        print(f"  Warning: {tar_path} not found, no ligand poses loaded for this pocket.")
        return []

    wanted_names = {f"docking/{cid}_out.sdf": rank for rank, (cid, _) in enumerate(top_ids, start=1)}
    sdf_bytes_by_rank = {}
    with tarfile.open(tar_path, "r:gz") as tf:
        for member in tf:
            if member.name not in wanted_names:
                continue
            sdf_bytes_by_rank[wanted_names[member.name]] = tf.extractfile(member).read()
            if len(sdf_bytes_by_rank) == len(wanted_names):
                break

    loaded = []
    for rank in range(1, len(top_ids) + 1):
        if rank not in sdf_bytes_by_rank:
            continue
        obj_name = f"{obj_prefix}_top{rank}"
        with tempfile.NamedTemporaryFile(suffix=".sdf", mode="wb", delete=False) as tmp:
            tmp.write(sdf_bytes_by_rank[rank])
            tmp_path = tmp.name
        cmd.load(tmp_path, obj_name)
        os.remove(tmp_path)
        cmd.show("sticks", obj_name)
        cmd.hide("lines", obj_name)
        color_ligand(obj_name, COLOR_LIGAND_DOCKED, "ligC_docked")
        cmd.disable(obj_name)
        loaded.append(obj_name)

    if len(loaded) < len(top_ids):
        print(f"  Warning: only found {len(loaded)}/{len(top_ids)} top ligand poses in {tar_path}.")
    return loaded


def main():
    parser = argparse.ArgumentParser(
        description="Build one PyMOL session containing all pockets in output/selected_pockets.csv, "
                     "each with its top-N best-scoring docked ligand poses (hidden by default)."
    )
    parser.add_argument("--top-n", type=int, default=5, help="Number of top-scoring ligands per pocket (default: 5)")
    parser.add_argument("--no-surface", action="store_true",
                         help="Use cartoon instead of surface for structures, skipping surface calculation (faster)")
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    selected = pd.read_csv(SELECTED_POCKETS_CSV)
    gene_map = load_gene_map()

    pymol.finish_launching(["pymol", "-cq"])
    cmd.reinitialize()

    for _, row in selected.iterrows():
        gene_name, site_type, pocket_name, comment = row["gene_name"], row["site_type"], row["pocket_name"], row["comment"]
        final_name = f"{gene_name}_{site_type}_{pocket_name}"

        print("---------------------------------------------")
        print(f"Pocket: {pocket_name} (gene: {gene_name}, {site_type})")
        print(f"Comment: {comment}")

        resolved = resolve_pocket(gene_name, pocket_name, gene_map)
        if resolved is None:
            print(f"  Warning: structure or pocket centroid file not found for {pocket_name}. Skipping.")
            print("---------------------------------------------")
            continue
        structure_path, pocket_pdb_path, probability, docking_dir = resolved

        print(f"Pocket probability: {probability}")

        load_and_merge_structure_pocket(structure_path, pocket_pdb_path, final_name, no_surface=args.no_surface)

        top_scores = get_top_n_scores(docking_dir, args.top_n)
        if top_scores:
            scores = [score for _, score in top_scores]
            print(f"Top-{args.top_n} docking scores: min={min(scores)}, avg={sum(scores) / len(scores):.3f}")
            load_ligand_objects(docking_dir, top_scores, final_name)
        else:
            print(f"  Warning: no report.csv found at {docking_dir}, no docking scores/ligands for this pocket.")

        print("---------------------------------------------")

    cmd.bg_color("white")
    cmd.save(OUTPUT_PSE)
    print(f"\nSaved: {OUTPUT_PSE}")


if __name__ == "__main__":
    main()
