#!/usr/bin/env python3
"""
PyMOL visualization of the manually-curated reference pocket for a gene.

Loads ONLY the reference structure (output/aligned_relaxed_structures/<uniprot_ac>/<structure>.pdb)
recorded for the gene in output/reference_pocket.csv, and highlights ONLY that gene's
reference pocket (exact P2Rank residues + centroid spheres). Saves a .pse session for
visual QC of the curation done via scripts/47_docking_summary.py.

Usage:
    python 47b_reference_pocket_visualization.py --trna pheS
    python 47b_reference_pocket_visualization.py --trna pheS,aspS
"""
import argparse
import os
import sys
import tarfile
import tempfile
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*converter.*already registered.*")

import matplotlib.pyplot as plt
import pandas as pd
import pymol
import requests
from matplotlib.colors import to_rgb
from pymol import cmd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from docking_utils import LIBRARIES, load_gene_map, load_reference_pockets

POCKET_DETECTION_CSV = os.path.join(ROOT, "output", "pocket_detection_data.csv")
ALIGNED_DIR = os.path.join(ROOT, "output", "aligned_relaxed_structures")
DETECTED_POCKETS_DIR = os.path.join(ROOT, "output", "detected_pockets")
ALPHAFILL_DIR = os.path.join(ROOT, "data", "structures", "alphafill_database")
OUTPUT_DIR = os.path.join(ROOT, "output", "47b_reference_pocket_visualization")

# Same coloring convention as notebooks/46_pocket_visualisation.ipynb and
# notebooks/46_protein_visualization.ipynb, so the same protein always gets the same
# color across the project's PyMOL visualizations.
COLOR_STRUCTURE = [0.7804, 0.8275, 0.8667]
COLOR_LIGAND_DOCKED = [0xF5 / 255, 0xA6 / 255, 0x3A / 255]  # orange, matches the notebook convention
COLOR_LIGAND_PDB = [1.0, 0.0, 1.0]  # magenta, to visually distinguish experimental from docked ligands
COLOR_LIGAND_ALPHAFILL = [0.0, 1.0, 1.0]  # cyan, for AlphaFill-transplanted ligands
UNIPROT_IDS_ORDERED = [
    "P9WFW5", "P9WFW7", "P9WFW3", "P9WQA1", "P9WN61", "P9WFT5", "P9WFV3", "P9WFS9",
    "P9WFV1", "P9WFV9", "P9WFT9", "P9WFV7", "P9WFT7", "P9WFW1", "P9WFU5", "P9WFU9",
    "P9WFV5", "P9WFT3", "P9WFU3", "P9WFU1", "P9WFT1",
]

TOP_N_LIGANDS = 10
UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/{}.json"
PDB_DOWNLOAD_URL = "https://files.rcsb.org/download/{}.cif"


def get_gene_color(uniprot_ac):
    """Per-protein color from the same tab20/tab20b palette + fixed UNIPROT_IDS_ORDERED
    indexing used in notebooks/46_pocket_visualisation.ipynb and 46_protein_visualization.ipynb."""
    palette = list(plt.get_cmap("tab20").colors) + [plt.get_cmap("tab20b").colors[0]]
    return list(to_rgb(palette[UNIPROT_IDS_ORDERED.index(uniprot_ac)]))


def color_ligand(selection, carbon_color, color_name):
    """Color-by-element with a custom carbon color, so ligands from different sources
    (docked vs. experimental) can be told apart at a glance."""
    cmd.util.cbag(selection)
    cmd.set_color(color_name, carbon_color)
    cmd.color(color_name, f"{selection} and elem C")


def get_pdb_structures(uniprot_ac):
    """Query UniProt for all experimental PDB structures cross-referenced for this accession.
    Returns a list of {id, method, resolution, chains} dicts (empty on error/no hits)."""
    try:
        resp = requests.get(UNIPROT_API_URL.format(uniprot_ac), timeout=15)
        resp.raise_for_status()
    except requests.RequestException as e:
        print(f"  Warning: could not query UniProt for {uniprot_ac} ({e}).")
        return []

    pdb_refs = [x for x in resp.json().get("uniProtKBCrossReferences", []) if x["database"] == "PDB"]
    results = []
    for ref in pdb_refs:
        props = {p["key"]: p["value"] for p in ref.get("properties", [])}
        results.append({
            "id": ref["id"],
            "method": props.get("Method", "N/A"),
            "resolution": props.get("Resolution", "N/A"),
            "chains": props.get("Chains", "N/A"),
        })
    return results


def print_pdb_structures(uniprot_ac, pdb_refs):
    if not pdb_refs:
        print(f"  No PDB structures found in UniProt for {uniprot_ac}.")
        return
    print(f"  PDB structures for {uniprot_ac} ({len(pdb_refs)}):")
    for ref in pdb_refs:
        print(f"    {ref['id']}: {ref['method']}, {ref['resolution']}, chains {ref['chains']}")


def parse_poi_chains(chains_str):
    """Parse UniProt's 'Chains' PDB cross-reference property (e.g. 'A/D=1-341', possibly
    comma-separated for multiple ranges) into the list of chain IDs belonging to the
    protein of interest (as opposed to other chains present in the same deposited entry)."""
    chains = []
    if not chains_str or chains_str == "N/A":
        return chains
    for part in chains_str.split(","):
        chain_part = part.split("=")[0] if "=" in part else part
        chains.extend(c.strip() for c in chain_part.split("/") if c.strip())
    return sorted(set(chains))


def download_and_load_pdb_structures(uniprot_ac, pdb_refs, pocket_resi_str):
    """Download (if not already cached) each cross-referenced PDB entry, then split it into
    one object per protein-of-interest chain (as annotated by UniProt's 'Chains' property,
    via parse_poi_chains). Only the chain's bound ligands are shown (as sticks); its own
    cartoon/backbone and solvent are hidden so only the ligands are visible.
    Chains not belonging to the protein of interest (e.g. a complex partner) are dropped.

    Each chain is aligned onto "structure" using ONLY the reference pocket's own residues
    (pocket_resi_str, a PyMOL "resi" selector string e.g. "148+154+158+..."), not the whole
    chain — this protein is multi-domain, so a whole-chain rigid-body fit is dominated by
    domains irrelevant to the pocket and produces a poor local fit around the pocket itself.

    Returns the list of loaded "<pdb_id>_<chain>" object names."""
    if not pdb_refs:
        return []

    cache_dir = os.path.join(OUTPUT_DIR, "pdb_structures", uniprot_ac)
    os.makedirs(cache_dir, exist_ok=True)

    loaded = []
    for ref in pdb_refs:
        pdb_id = ref["id"]
        poi_chains = parse_poi_chains(ref["chains"])
        if not poi_chains:
            print(f"  Warning: could not parse POI chain(s) from '{ref['chains']}' for {pdb_id}, skipping.")
            continue

        cif_path = os.path.join(cache_dir, f"{pdb_id}.cif")
        if not os.path.isfile(cif_path):
            try:
                resp = requests.get(PDB_DOWNLOAD_URL.format(pdb_id), timeout=30)
                resp.raise_for_status()
            except requests.RequestException as e:
                print(f"  Warning: could not download PDB {pdb_id} ({e}), skipping.")
                continue
            with open(cif_path, "wb") as f:
                f.write(resp.content)

        cmd.load(cif_path, pdb_id)
        available_chains = set(cmd.get_chains(pdb_id))
        for chain in poi_chains:
            if chain not in available_chains:
                print(f"  Warning: POI chain {chain} not found in {pdb_id} "
                      f"(has {sorted(available_chains)}), skipping.")
                continue
            obj_name = f"{pdb_id}_{chain}"
            cmd.create(obj_name, f"{pdb_id} and chain {chain}")
            try:
                rms, n_aligned, *_ = cmd.align(
                    f"{obj_name} and resi {pocket_resi_str}",
                    f"structure and resi {pocket_resi_str}",
                )
                print(f"    Aligned {obj_name} to structure (pocket residues only): "
                      f"RMSD={rms:.3f} A over {n_aligned} atoms")
            except pymol.CmdException as e:
                print(f"  Warning: could not align {obj_name} to structure ({e}).")

            # Only the bound ligands matter here: hide the chain's own cartoon/backbone
            # and solvent, showing just its non-solvent heteroatoms as sticks.
            cmd.hide("everything", obj_name)
            ligand_sel = f"{obj_name} and hetatm and not solvent"
            cmd.show("sticks", ligand_sel)
            color_ligand(ligand_sel, COLOR_LIGAND_PDB, "ligC_pdb")

            loaded.append(obj_name)
        cmd.delete(pdb_id)

    return loaded


def load_alphafill_ligands(uniprot_ac, pocket_resi_str):
    """Load the local AlphaFill structure (data/structures/alphafill_database/<uniprot_ac>/<uniprot_ac>.cif)
    — AlphaFold's predicted structure with ligands transplanted onto it from homologous PDB
    entries, one per chain (the protein itself is also just one of the chains; no fixed
    chain letter is assumed). Align it onto "structure" using the pocket residues only
    (AlphaFill is anchored to AlphaFold2, not necessarily the same frame as our reference
    structure), then show only its transplanted (non-solvent) heteroatoms as sticks, across
    all chains. Returns True if loaded, False otherwise."""
    cif_path = os.path.join(ALPHAFILL_DIR, uniprot_ac, f"{uniprot_ac}.cif")
    if not os.path.isfile(cif_path):
        print(f"  Warning: AlphaFill structure not found at {cif_path}, skipping.")
        return False

    cmd.load(cif_path, "alphafill")
    try:
        rms, n_aligned, *_ = cmd.align(
            f"alphafill and polymer and resi {pocket_resi_str}",
            f"structure and resi {pocket_resi_str}",
        )
        print(f"    Aligned alphafill to structure (pocket residues only): "
              f"RMSD={rms:.3f} A over {n_aligned} atoms")
    except pymol.CmdException as e:
        print(f"  Warning: could not align alphafill to structure ({e}).")

    cmd.hide("everything", "alphafill")
    ligand_sel = "alphafill and hetatm and not solvent"
    cmd.show("sticks", ligand_sel)
    color_ligand(ligand_sel, COLOR_LIGAND_ALPHAFILL, "ligC_alphafill")
    return True


def load_top_ligands(pocket_name):
    """Load the TOP_N_LIGANDS best-scoring (lowest score) REAL-library docked poses
    for this pocket as PyMOL objects named "top_<rank>_<compound_id>" (rank 1 = best),
    colored by element with orange carbons (distinct from the magenta used for
    PDB-derived experimental ligands, so the two sources are easy to tell apart).
    Returns the list of loaded object names."""
    results_dir = os.path.join(LIBRARIES["REAL"], pocket_name)

    report_path = os.path.join(results_dir, "report.csv")
    if not os.path.isfile(report_path):
        print(f"  Warning: {report_path} not found, skipping ligands.")
        return []
    report = pd.read_csv(report_path).sort_values("score", ascending=True)
    top_ids = report["compound"].head(TOP_N_LIGANDS).tolist()

    tar_path = os.path.join(results_dir, "docking.tar.gz")
    if not os.path.isfile(tar_path):
        print(f"  Warning: {tar_path} not found, skipping ligands.")
        return []

    wanted = {f"docking/{cid}_out.sdf": f"top_{rank}_{cid}" for rank, cid in enumerate(top_ids, start=1)}
    loaded = []
    with tarfile.open(tar_path, "r|gz") as tf:
        for member in tf:
            if member.name not in wanted:
                continue
            obj_name = wanted[member.name]
            data = tf.extractfile(member).read()
            with tempfile.NamedTemporaryFile(suffix=".sdf", mode="wb", delete=False) as tmp:
                tmp.write(data)
                tmp_path = tmp.name
            cmd.load(tmp_path, obj_name)
            os.remove(tmp_path)
            cmd.show("sticks", obj_name)
            cmd.hide("lines", obj_name)
            color_ligand(obj_name, COLOR_LIGAND_DOCKED, "ligC_docked")
            loaded.append(obj_name)
            if len(loaded) == len(wanted):
                break

    if len(loaded) < len(top_ids):
        print(f"  Warning: only found {len(loaded)}/{len(top_ids)} top ligand poses in {tar_path}.")
    return loaded


def zoom_to_fixed_box(selection="structure", box_size=100.0, box_name="_zoom_box"):
    """Zoom to a fixed-size cube centered on `selection`, for consistent framing."""
    (xmin, ymin, zmin), (xmax, ymax, zmax) = cmd.get_extent(selection)
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    cz = 0.5 * (zmin + zmax)
    h = box_size / 2.0

    cmd.delete(box_name)
    corners = [
        (cx - h, cy - h, cz - h), (cx - h, cy - h, cz + h),
        (cx - h, cy + h, cz - h), (cx - h, cy + h, cz + h),
        (cx + h, cy - h, cz - h), (cx + h, cy - h, cz + h),
        (cx + h, cy + h, cz - h), (cx + h, cy + h, cz + h),
    ]
    for i, (x, y, z) in enumerate(corners):
        cmd.pseudoatom(f"{box_name}_{i}", pos=[x, y, z])
    cmd.group(box_name, f"{box_name}_*")
    cmd.zoom(box_name, buffer=0.0, complete=1)
    cmd.delete(box_name)
    cmd.delete(f"{box_name}_*")


def build_session(gene, uniprot_ac, pocket_name, pocket_data, pdb_refs):
    structure_name, pocket_number = pocket_name.rsplit("_pocket_", 1)

    structure_path = os.path.join(ALIGNED_DIR, uniprot_ac, f"{structure_name}.pdb")
    if not os.path.isfile(structure_path):
        print(f"  Error: structure file not found: {structure_path}")
        return False

    pocket_pdb_path = os.path.join(
        DETECTED_POCKETS_DIR, uniprot_ac, structure_name, "pockets", f"pocket_{pocket_number}.pdb"
    )
    if not os.path.isfile(pocket_pdb_path):
        print(f"  Error: pocket centroid file not found: {pocket_pdb_path}")
        return False

    row = pocket_data[
        (pocket_data["Uniprot AC"] == uniprot_ac)
        & (pocket_data["File name"] == f"{structure_name}.pdb")
        & (pocket_data["Pocket number"] == int(pocket_number))
    ]
    if row.empty:
        print(f"  Error: no entry in {POCKET_DETECTION_CSV} for "
              f"Uniprot AC={uniprot_ac}, File name={structure_name}.pdb, Pocket number={pocket_number}")
        return False
    pocket_residues = row.iloc[0]["Pocket residues (chain_resn)"].split()

    pymol.finish_launching(["pymol", "-cq"])
    cmd.reinitialize()

    # Load reference structure only
    cmd.load(structure_path, "structure")
    cmd.set_color("structure_color", COLOR_STRUCTURE)
    cmd.color("structure_color", "structure")
    cmd.show("surface", "structure")
    cmd.hide("cartoon", "structure")
    cmd.hide("lines", "structure")
    cmd.hide("sticks", "structure")
    cmd.set("transparency", 0.3, "structure")

    # Per-protein color (same tab20/tab20b convention as the 46_*.ipynb notebooks),
    # used for both the pocket centroid spheres and the pocket residues below.
    gene_color = get_gene_color(uniprot_ac)
    cmd.set_color("pocket_color", gene_color)

    # Pocket centroid spheres
    cmd.load(pocket_pdb_path, "pocket")
    cmd.color("pocket_color", "pocket")
    cmd.show("spheres", "pocket")
    cmd.set("sphere_transparency", 0.4, "pocket")
    cmd.set("sphere_scale", 6, "pocket")

    # Exact P2Rank residues for this pocket
    selection_str = " or ".join(
        f"resi {res.split('_')[1]} and chain {res.split('_')[0]} and structure"
        for res in pocket_residues
    )
    cmd.select("pocket_residues", selection_str)
    cmd.color("pocket_color", "pocket_residues")
    cmd.delete("pocket_residues")

    # Residue numbers of the pocket itself, used below to align experimental structures
    # on the pocket's local domain only (this protein is multi-domain; a whole-chain fit
    # is dominated by domains irrelevant to the pocket).
    pocket_resnums = sorted({int(res.split("_")[1]) for res in pocket_residues})
    pocket_resi_str = "+".join(str(n) for n in pocket_resnums)

    # Top-10 best REAL-library docked ligands
    ligands = load_top_ligands(pocket_name)
    if ligands:
        print(f"  Loaded {len(ligands)} top ligand(s): {', '.join(ligands)}")

    # Every POI chain of every experimental PDB structure cross-referenced in UniProt
    # (each chain its own object, showing only its bound ligands)
    pdb_chains = download_and_load_pdb_structures(uniprot_ac, pdb_refs, pocket_resi_str)
    if pdb_chains:
        print(f"  Loaded {len(pdb_chains)} experimental PDB chain(s): {', '.join(pdb_chains)}")

    # AlphaFill-transplanted ligands (local, not downloaded)
    load_alphafill_ligands(uniprot_ac, pocket_resi_str)

    # Fancy visualization
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
    zoom_to_fixed_box("structure", box_size=100)

    out_path = os.path.join(OUTPUT_DIR, f"session_{uniprot_ac}_{gene}.pse")
    cmd.save(out_path)
    print(f"  Saved: {out_path}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="PyMOL visualization of the manually-curated reference pocket."
    )
    parser.add_argument("--trna", required=True,
                        help="Comma-separated gene name(s) (e.g. pheS or pheS,aspS)")
    args = parser.parse_args()

    genes = [g.strip() for g in args.trna.split(",")]
    gene_map = load_gene_map()
    ref_pockets = load_reference_pockets()

    invalid = [g for g in genes if g not in gene_map]
    if invalid:
        print(f"Error: unknown gene(s): {', '.join(invalid)}. Available genes:")
        print("  " + ", ".join(sorted(gene_map)))
        sys.exit(1)

    missing_ref = [g for g in genes if g not in ref_pockets]
    for g in missing_ref:
        print(f"Skipping '{g}': no reference pocket defined. Run scripts/47_docking_summary.py "
              f"and record one in output/reference_pocket.csv (columns: gene_name, pocket_name).")
    genes = [g for g in genes if g not in missing_ref]

    if not genes:
        sys.exit(0)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    pocket_data = pd.read_csv(POCKET_DETECTION_CSV)

    for gene in genes:
        uniprot_ac = gene_map[gene]
        pocket_name = ref_pockets[gene]
        print(f"Gene: {gene}  (pocket: {pocket_name})")
        pdb_refs = get_pdb_structures(uniprot_ac)
        print_pdb_structures(uniprot_ac, pdb_refs)
        build_session(gene, uniprot_ac, pocket_name, pocket_data, pdb_refs)


if __name__ == "__main__":
    main()
