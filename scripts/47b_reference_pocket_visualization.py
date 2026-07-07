#!/usr/bin/env python3
"""
PyMOL visualization of the manually-curated reference pocket for a gene.

Loads the reference structure (output/aligned_relaxed_structures/<uniprot_ac>/<structure>.pdb)
recorded for the gene in output/reference_pocket_catalytic.csv and/or
output/reference_pocket_noncatalytic.csv — both are checked automatically, and a session is
built for each one that's defined for the gene — and highlights that pocket (exact P2Rank
residues + centroid spheres). On top of it, overlays the top-10 best-scoring REAL-library
docked poses, every experimental PDB structure cross-referenced in UniProt for that protein
(kept as a full assembly, aligned via one POI chain's pocket residues), any AlphaFill-
transplanted ligands, and optionally a user-supplied list of homolog PDB structures from
another protein/organism (loaded as "H_<pdb_id>", aligned via whole-chain sequence alignment
since residue numbering won't match across organisms). Saves one .pse session per pocket,
named after the pocket itself (e.g. session_swissmodel_P9WFU3_model_0_pocket_1.pse) for visual
QC of the curation done via scripts/47_docking_summary.py.

Usage:
    python 47b_reference_pocket_visualization.py --trna pheS
    python 47b_reference_pocket_visualization.py --trna pheS,aspS --homologs output/homolog_pdb_ids.csv
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
from docking_utils import LIBRARIES, load_gene_map

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
COLOR_LIGAND_HOMOLOG = [1.0, 1.0, 0.0]  # yellow, for user-supplied homolog structures
UNIPROT_IDS_ORDERED = [
    "P9WFW5", "P9WFW7", "P9WFW3", "P9WQA1", "P9WN61", "P9WFT5", "P9WFV3", "P9WFS9",
    "P9WFV1", "P9WFV9", "P9WFT9", "P9WFV7", "P9WFT7", "P9WFW1", "P9WFU5", "P9WFU9",
    "P9WFV5", "P9WFT3", "P9WFU3", "P9WFU1", "P9WFT1",
]

TOP_N_LIGANDS = 5
UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/{}.json"
PDB_DOWNLOAD_URL = "https://files.rcsb.org/download/{}.cif"
RCSB_ENTRY_API_URL = "https://data.rcsb.org/rest/v1/core/entry/{}"
RCSB_ASSEMBLY_API_URL = "https://data.rcsb.org/rest/v1/core/assembly/{}/{}"


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


# def add_hbond_dashes(receptor_sel, ligand_sel, dash_name, cutoff=3.5, angle=45):
#     """Detect donor/acceptor H-bond contacts between ligand_sel and receptor_sel within
#     `cutoff` Angstrom / `angle` degrees (same criteria as notebooks/46_docking_exploration_IIa.ipynb
#     and IIb.ipynb), drawn as a black dashed-line object. Returns the number of contacts found."""
#     nearby = f"({receptor_sel} within 6 of {ligand_sel}) and not solvent"
#     pairs = []
#     pairs += cmd.find_pairs(f"{ligand_sel} and donor", f"{nearby} and acceptor", mode=1, cutoff=cutoff, angle=angle)
#     pairs += cmd.find_pairs(f"{ligand_sel} and acceptor", f"{nearby} and donor", mode=1, cutoff=cutoff, angle=angle)
#     pairs = list(dict.fromkeys(pairs))
#     if not pairs:
#         return 0
#     cmd.delete(dash_name)
#     for a1, a2 in pairs:
#         cmd.distance(dash_name, a1, a2)
#     cmd.hide("labels", dash_name)
#     cmd.color("black", dash_name)
#     cmd.set("dash_width", 2, dash_name)
#     return len(pairs)


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


def get_biological_assembly_info(pdb_id):
    """Query RCSB for pdb_id's primary annotated biological assembly (the first of
    rcsb_entry_container_identifiers.assembly_ids) and return a short human-readable summary,
    e.g. "heterotetrameric (4 chains), author_and_software_defined_assembly", or None on any
    error/missing data. The 'details' field distinguishes author-confirmed assemblies (real
    experimental evidence, e.g. SEC-MALS) from software-only (PISA) predictions — reported
    verbatim rather than collapsed into a yes/no verdict, since that distinction is the point."""
    try:
        entry_resp = requests.get(RCSB_ENTRY_API_URL.format(pdb_id), timeout=15)
        entry_resp.raise_for_status()
        assembly_ids = entry_resp.json().get("rcsb_entry_container_identifiers", {}).get("assembly_ids", [])
        if not assembly_ids:
            return None

        asm_resp = requests.get(RCSB_ASSEMBLY_API_URL.format(pdb_id, assembly_ids[0]), timeout=15)
        asm_resp.raise_for_status()
        asm = asm_resp.json().get("pdbx_struct_assembly", {})
    except requests.RequestException:
        return None

    oligomeric = asm.get("oligomeric_details", "unknown oligomeric state")
    count = asm.get("oligomeric_count")
    count_str = f" ({count} chains)" if count else ""
    details = asm.get("details", "unknown method")
    extra = f", {len(assembly_ids)} assemblies reported" if len(assembly_ids) > 1 else ""
    return f"{oligomeric}{count_str}, {details}{extra}"


def print_pdb_structures(uniprot_ac, pdb_refs):
    if not pdb_refs:
        print(f"  No PDB structures found in UniProt for {uniprot_ac}.")
        return
    print(f"  PDB structures for {uniprot_ac} ({len(pdb_refs)}):")
    for ref in pdb_refs:
        print(f"    {ref['id']}: {ref['method']}, {ref['resolution']}, chains {ref['chains']}")


def load_reference_pocket_csv(path):
    """Load a gene_name -> pocket_name CSV (output/reference_pocket_catalytic.csv or
    output/reference_pocket_noncatalytic.csv, selected via --catalytic / --noncatalytic).
    Rows with an empty/NaN pocket_name (not yet curated for that gene) are dropped rather
    than erroring, since that's the expected state for genes not yet assigned a pocket of
    this kind. Exits with a clear error if the file itself is missing."""
    if not os.path.isfile(path):
        print(f"Error: {path} not found.")
        print("Run script 47 and manually populate it (columns: gene_name, pocket_name).")
        sys.exit(1)
    df = pd.read_csv(path)
    df = df[df["pocket_name"].notna() & (df["pocket_name"].astype(str).str.strip() != "")]
    return dict(zip(df["gene_name"], df["pocket_name"]))


def load_homolog_pdb_ids(path, gene_map):
    """Load a CSV mapping gene_name -> ';'-separated homolog PDB IDs (see --homologs help
    for the expected format and where to source candidates from). Returns
    {gene_name: [pdb_id, ...]}; genes not present in the file simply get no homolog overlay.
    Warns (does not fail) on gene names not recognized in gene_map, to catch typos."""
    if not os.path.isfile(path):
        print(f"Error: {path} not found.")
        print("Expected a CSV with columns: gene_name, homolog_pdb_ids (';'-separated PDB IDs)")
        sys.exit(1)

    df = pd.read_csv(path)
    homolog_map = {}
    for _, row in df.iterrows():
        gene = row["gene_name"]
        if gene not in gene_map:
            print(f"  Warning: '{gene}' in {path} is not a recognized gene, skipping.")
            continue
        ids = [pid.strip() for pid in str(row["homolog_pdb_ids"]).split(";") if pid.strip()]
        if ids:
            homolog_map[gene] = ids
    return homolog_map


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
    """Download (if not already cached) each cross-referenced PDB entry and keep the full
    deposited assembly intact — all chains, not just the protein-of-interest chain(s). This
    lets a second POI copy or a complex partner subunit show up in its natural position
    relative to the reference pocket, so you can tell whether the pocket actually sits at a
    protein-protein interface (and whether a bound ligand is really a PPI-site binder rather
    than sitting in an isolated pocket).

    Alignment onto "structure" uses ONLY the reference pocket's own residues (pocket_resi_str,
    a PyMOL "resi" selector string e.g. "148+154+158+...") on ONE representative POI chain
    (as annotated by UniProt's 'Chains' property, via parse_poi_chains) — not the whole chain,
    since this protein is multi-domain and a whole-chain rigid-body fit is dominated by domains
    irrelevant to the pocket. PyMOL's cmd.align moves the entire mobile object, so every other
    chain in the entry rides along rigidly with the aligned POI chain.

    Returns the list of loaded "<pdb_id>" object names."""
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
        align_chain = next((c for c in poi_chains if c in available_chains), None)
        if align_chain is None:
            print(f"  Warning: none of POI chain(s) {poi_chains} found in {pdb_id} "
                  f"(has {sorted(available_chains)}), skipping.")
            cmd.delete(pdb_id)
            continue

        try:
            rms, n_aligned, *_ = cmd.align(
                f"{pdb_id} and chain {align_chain} and resi {pocket_resi_str}",
                f"structure and resi {pocket_resi_str}",
            )
            print(f"    Aligned {pdb_id} to structure (chain {align_chain} pocket residues only): "
                  f"RMSD={rms:.3f} A over {n_aligned} atoms")
        except pymol.CmdException as e:
            print(f"  Warning: could not align {pdb_id} ({e}).")

        assembly_info = get_biological_assembly_info(pdb_id)
        if assembly_info:
            print(f"    Biological assembly: {assembly_info}")

        # Keep the whole assembly as cartoon (all chains) plus every non-solvent heteroatom
        # as sticks, so ligands bound to non-POI chains stay visible too.
        cmd.hide("everything", pdb_id)
        cmd.show("cartoon", pdb_id)
        ligand_sel = f"{pdb_id} and hetatm and not solvent"
        cmd.show("sticks", ligand_sel)
        color_ligand(ligand_sel, COLOR_LIGAND_PDB, "ligC_pdb")

        loaded.append(pdb_id)

    return loaded


def download_and_load_homolog_structures(uniprot_ac, homolog_ids):
    """Download (if not already cached) each user-supplied homolog PDB entry — a structure
    of a DIFFERENT protein/organism than the target, picked by the user (e.g. from AlphaFill's
    own donor-structure list in data/structures/alphafill_database/<uniprot_ac>/<uniprot_ac>.json,
    or an independent homology search) — and load the full deposited assembly (all chains),
    named "H_<pdb_id>" so it's visually distinguishable from same-organism experimental
    structures (loaded as plain "<pdb_id>" by download_and_load_pdb_structures).

    Since these come from a different protein/organism, residue numbering does not correspond
    to the reference structure's numbering, so alignment can't use pocket_resi_str like
    download_and_load_pdb_structures does. Instead, cmd.super is run on the first available
    chain against the whole reference structure. cmd.super (not cmd.align or cmd.cealign) was
    picked after empirically comparing all three across every homolog gathered so far: cmd.align
    (pure sequence-based) is thrown off by divergent domains (e.g. aspS's 3G1Z: RMSD=6.43 A;
    lysS's 7MPP: RMSD=12.63 A), and cmd.cealign (CE, fragment-chaining) is worse than plain align
    for several close orthologs despite being better on those same divergent cases (e.g. pheS's
    3PCO: align=1.34 A but cealign=6.43 A). cmd.super's structure-guided sequence alignment matched
    or beat cmd.align on every close-ortholog case AND fixed every divergent case (3G1Z: 1.60 A;
    7MPP: 3.55 A; aspS's 5YAL: 22.35->3.52 A), while still correctly reporting a bad RMSD (16.9 A)
    for a genuinely unrelated donor (alaS's spurious AlphaFill hit 3TZI, murine COX-2) rather than
    forcing a false-looking-good fit.

    Returns the list of loaded "H_<pdb_id>" object names."""
    if not homolog_ids:
        return []

    cache_dir = os.path.join(OUTPUT_DIR, "homolog_structures", uniprot_ac)
    os.makedirs(cache_dir, exist_ok=True)

    loaded = []
    for pdb_id in homolog_ids:
        obj_name = f"H_{pdb_id}"
        cif_path = os.path.join(cache_dir, f"{pdb_id}.cif")
        if not os.path.isfile(cif_path):
            try:
                resp = requests.get(PDB_DOWNLOAD_URL.format(pdb_id), timeout=30)
                resp.raise_for_status()
            except requests.RequestException as e:
                print(f"  Warning: could not download homolog PDB {pdb_id} ({e}), skipping.")
                continue
            with open(cif_path, "wb") as f:
                f.write(resp.content)

        cmd.load(cif_path, obj_name)
        available_chains = sorted(cmd.get_chains(obj_name))
        if not available_chains:
            print(f"  Warning: no chains found in homolog {pdb_id}, skipping.")
            cmd.delete(obj_name)
            continue
        align_chain = available_chains[0]

        try:
            rms, n_aligned, *_ = cmd.super(
                f"{obj_name} and chain {align_chain} and polymer",
                "structure and polymer",
            )
            print(f"    Aligned homolog {pdb_id} to structure (chain {align_chain}, "
                  f"super): RMSD={rms:.3f} A over {n_aligned} atoms")
        except pymol.CmdException as e:
            print(f"  Warning: could not align homolog {pdb_id} ({e}).")

        assembly_info = get_biological_assembly_info(pdb_id)
        if assembly_info:
            print(f"    Biological assembly: {assembly_info}")

        # Keep the whole assembly as cartoon (all chains) plus every non-solvent heteroatom
        # as sticks, same as native PDB structures, but in a distinct color (yellow).
        cmd.hide("everything", obj_name)
        cmd.show("cartoon", obj_name)
        ligand_sel = f"{obj_name} and hetatm and not solvent"
        cmd.show("sticks", ligand_sel)
        color_ligand(ligand_sel, COLOR_LIGAND_HOMOLOG, "ligC_homolog")

        loaded.append(obj_name)

    return loaded


def load_alphafill_ligands(uniprot_ac, pocket_resi_str):
    """Load the local AlphaFill structure (data/structures/alphafill_database/<uniprot_ac>/<uniprot_ac>.cif)
    — AlphaFold's predicted structure with ligands transplanted onto it from homologous PDB
    entries, one per chain (the protein itself is also just one of the chains; no fixed
    chain letter is assumed). Align it onto "structure" using the pocket residues only
    (AlphaFill is anchored to AlphaFold2, not necessarily the same frame as our reference
    structure), then show it as cartoon plus its transplanted (non-solvent) heteroatoms as
    sticks, across all chains. Returns True if loaded, False otherwise."""
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
    cmd.show("cartoon", "alphafill")
    ligand_sel = "alphafill and hetatm and not solvent"
    cmd.show("sticks", ligand_sel)
    color_ligand(ligand_sel, COLOR_LIGAND_ALPHAFILL, "ligC_alphafill")
    return True


def load_top_ligands(pocket_name, kind):
    """Load the TOP_N_LIGANDS best-scoring (lowest score) REAL-library docked poses
    for this pocket as PyMOL objects named "top_<rank>_<kind>" (kind = "CAT"/"NONCAT",
    rank 1 = best), colored by element with orange carbons (distinct from the magenta
    used for PDB-derived experimental ligands, so the two sources are easy to tell apart).
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

    wanted = {f"docking/{cid}_out.sdf": f"top_{rank}_{kind}" for rank, cid in enumerate(top_ids, start=1)}
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
            # cmd.h_add(obj_name)
            # n_hbonds = add_hbond_dashes("pocket_residues", obj_name, f"hbonds_{obj_name}")
            # print(f"    {obj_name}: {n_hbonds} H-bond(s) with pocket_residues")
            loaded.append(obj_name)
            if len(loaded) == len(wanted):
                break

    if len(loaded) < len(top_ids):
        print(f"  Warning: only found {len(loaded)}/{len(top_ids)} top ligand poses in {tar_path}.")
    # cmd.hide("everything", "hydro")
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


def _load_structure_and_pocket(pocket_name, uniprot_ac, pocket_data):
    """Load <structure_name>.pdb as object "structure" and pocket_<n>.pdb as object
    "pocket" for the given pocket_name. Returns (True, pocket_resi_str) on success,
    (False, None) if the structure file, pocket centroid file, or pocket_data row is
    missing (each case prints its own error)."""
    structure_name, pocket_number = pocket_name.rsplit("_pocket_", 1)

    structure_path = os.path.join(ALIGNED_DIR, uniprot_ac, f"{structure_name}.pdb")
    if not os.path.isfile(structure_path):
        print(f"  Error: structure file not found: {structure_path}")
        return False, None

    pocket_pdb_path = os.path.join(
        DETECTED_POCKETS_DIR, uniprot_ac, structure_name, "pockets", f"pocket_{pocket_number}.pdb"
    )
    if not os.path.isfile(pocket_pdb_path):
        print(f"  Error: pocket centroid file not found: {pocket_pdb_path}")
        return False, None

    row = pocket_data[
        (pocket_data["Uniprot AC"] == uniprot_ac)
        & (pocket_data["File name"] == f"{structure_name}.pdb")
        & (pocket_data["Pocket number"] == int(pocket_number))
    ]
    if row.empty:
        print(f"  Error: no entry in {POCKET_DETECTION_CSV} for "
              f"Uniprot AC={uniprot_ac}, File name={structure_name}.pdb, Pocket number={pocket_number}")
        return False, None
    pocket_residues = row.iloc[0]["Pocket residues (chain_resn)"].split()

    # Load reference structure only
    cmd.load(structure_path, "structure")
    cmd.set_color("structure_color", COLOR_STRUCTURE)
    cmd.color("structure_color", "structure")
    cmd.show("cartoon", "structure")
    cmd.hide("surface", "structure")
    cmd.hide("lines", "structure")
    cmd.hide("sticks", "structure")
    cmd.set("transparency", 0.3, "structure")

    # Per-protein color (same tab20/tab20b convention as the 46_*.ipynb notebooks),
    # used for the pocket centroid spheres below.
    gene_color = get_gene_color(uniprot_ac)
    cmd.set_color("pocket_color", gene_color)

    # Pocket centroid spheres
    cmd.load(pocket_pdb_path, "pocket")
    cmd.color("pocket_color", "pocket")
    cmd.show("spheres", "pocket")
    cmd.set("sphere_transparency", 0.4, "pocket")
    cmd.set("sphere_scale", 6, "pocket")

    # Residue numbers of the pocket itself, used to align experimental structures on the
    # pocket's local domain only (this protein is multi-domain; a whole-chain fit is
    # dominated by domains irrelevant to the pocket).
    pocket_resnums = sorted({int(res.split("_")[1]) for res in pocket_residues})
    pocket_resi_str = "+".join(str(n) for n in pocket_resnums)
    return True, pocket_resi_str


def _merge_structure_and_pocket(final_name):
    """Merge the current "structure" + "pocket" objects into one, renamed to final_name
    (e.g. "swissmodel_P9WFU3_model_0_pocket_1") — same merge pattern used by
    scripts/47_docking_summary.py's build_pocket_overview_session."""
    tmp_merged = "_tmp_merged"
    cmd.create(tmp_merged, "structure or pocket")
    cmd.delete("structure")
    cmd.delete("pocket")
    cmd.set_name(tmp_merged, final_name)


def build_session(gene, uniprot_ac, cat_pocket_name, noncat_pocket_name, pocket_data, pdb_refs, homolog_ids=None):
    """Build ONE session for this gene containing both the catalytic and non-catalytic
    reference pocket (whichever are defined). Each pocket's structure + its own centroid
    sphere are merged into a single object named after that pocket. Native PDB structures,
    --homologs, and AlphaFill ligands are loaded ONCE, aligned via the "primary" pocket
    (catalytic if defined, else non-catalytic) — not duplicated for the secondary pocket."""
    pymol.finish_launching(["pymol", "-cq"])
    cmd.reinitialize()

    if cat_pocket_name:
        primary_kind, primary_name = "CAT", cat_pocket_name
        secondary_kind, secondary_name = "NONCAT", noncat_pocket_name
    else:
        primary_kind, primary_name = "NONCAT", noncat_pocket_name
        secondary_kind, secondary_name = "CAT", cat_pocket_name

    ok, pocket_resi_str = _load_structure_and_pocket(primary_name, uniprot_ac, pocket_data)
    if not ok:
        return False

    ligands = load_top_ligands(primary_name, primary_kind)
    if ligands:
        print(f"  Loaded {len(ligands)} top ligand(s) ({primary_kind}): {', '.join(ligands)}")
    cat_ligands, noncat_ligands = (ligands, []) if primary_kind == "CAT" else ([], ligands)

    # Every experimental PDB structure cross-referenced in UniProt, kept as a full assembly
    # (all chains) aligned via one representative POI chain's pocket residues
    pdb_structures = download_and_load_pdb_structures(uniprot_ac, pdb_refs, pocket_resi_str)
    if pdb_structures:
        print(f"  Loaded {len(pdb_structures)} experimental PDB structure(s): {', '.join(pdb_structures)}")

    # AlphaFill-transplanted ligands (local, not downloaded)
    alphafill_loaded = load_alphafill_ligands(uniprot_ac, pocket_resi_str)

    # User-supplied homolog PDB structures (different protein/organism), kept as full
    # assemblies named "H_<pdb_id>" and aligned via cmd.super (structure-guided, robust to
    # inter-domain conformational differences between the homolog and the reference)
    homolog_structures = download_and_load_homolog_structures(uniprot_ac, homolog_ids)
    if homolog_structures:
        print(f"  Loaded {len(homolog_structures)} homolog structure(s): {', '.join(homolog_structures)}")

    _merge_structure_and_pocket(primary_name)

    if secondary_name:
        ok, _ = _load_structure_and_pocket(secondary_name, uniprot_ac, pocket_data)
        if ok:
            ligands = load_top_ligands(secondary_name, secondary_kind)
            if ligands:
                print(f"  Loaded {len(ligands)} top ligand(s) ({secondary_kind}): {', '.join(ligands)}")
            if secondary_kind == "CAT":
                cat_ligands = ligands
            else:
                noncat_ligands = ligands
            _merge_structure_and_pocket(secondary_name)

    # Object panel order: CAT structure, NONCAT structure, native PDBs, AlphaFill,
    # homologs, top-CAT ligands, top-NONCAT ligands.
    order_names = []
    if cat_pocket_name:
        order_names.append(cat_pocket_name)
    if noncat_pocket_name:
        order_names.append(noncat_pocket_name)
    order_names += pdb_structures
    if alphafill_loaded:
        order_names.append("alphafill")
    order_names += homolog_structures
    order_names += cat_ligands
    order_names += noncat_ligands
    cmd.order(" ".join(order_names), location="top", sort=0)

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
    zoom_to_fixed_box("all", box_size=100)

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
    parser.add_argument("--homologs",
                        help="Path to a CSV mapping gene_name -> ';'-separated homolog PDB ID(s) "
                             "from a different protein/organism to overlay, e.g. a row "
                             "'alaS,3HXU;3HXV;3HXW' for E. coli AlaRS structures when the target "
                             "has no native PDB coverage. Loaded as 'H_<pdb_id>' in the session. "
                             "Genes not listed in the file simply get no homolog overlay. A good "
                             "source of candidates is the donor PDB list AlphaFill already used, "
                             "in data/structures/alphafill_database/<uniprot_ac>/<uniprot_ac>.json "
                             "('hits' field) — pick the ones with low global_rmsd and a "
                             "biologically sensible donor (see 47b's own history for the "
                             "3TZI/COX-2 counterexample of a spurious hit to exclude).")
    args = parser.parse_args()

    catalytic_path = os.path.join(ROOT, "output", "reference_pocket_catalytic.csv")
    noncatalytic_path = os.path.join(ROOT, "output", "reference_pocket_noncatalytic.csv")
    ref_pockets_catalytic = load_reference_pocket_csv(catalytic_path)
    ref_pockets_noncatalytic = load_reference_pocket_csv(noncatalytic_path)

    genes = [g.strip() for g in args.trna.split(",")]
    gene_map = load_gene_map()

    invalid = [g for g in genes if g not in gene_map]
    if invalid:
        print(f"Error: unknown gene(s): {', '.join(invalid)}. Available genes:")
        print("  " + ", ".join(sorted(gene_map)))
        sys.exit(1)

    missing_ref = [g for g in genes if g not in ref_pockets_catalytic and g not in ref_pockets_noncatalytic]
    for g in missing_ref:
        print(f"Skipping '{g}': no catalytic or non-catalytic reference pocket defined. "
              f"Run scripts/47_docking_summary.py and record one in {catalytic_path} "
              f"and/or {noncatalytic_path} (columns: gene_name, pocket_name).")
    genes = [g for g in genes if g not in missing_ref]

    if not genes:
        sys.exit(0)

    homolog_map = load_homolog_pdb_ids(args.homologs, gene_map) if args.homologs else {}

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    pocket_data = pd.read_csv(POCKET_DETECTION_CSV)

    for gene in genes:
        uniprot_ac = gene_map[gene]
        print(f"Gene: {gene}")
        pdb_refs = get_pdb_structures(uniprot_ac)
        print_pdb_structures(uniprot_ac, pdb_refs)

        # Drop any "homolog" IDs that are actually already covered by this gene's own
        # native PDB structures (e.g. AlphaFill's donor list for a well-studied target
        # often includes the target's own solved structures) — loading them a second
        # time via the homolog path only gives a redundant, worse-aligned duplicate.
        homolog_ids = homolog_map.get(gene)
        if homolog_ids:
            native_ids = {ref["id"] for ref in pdb_refs}
            overlap = [h for h in homolog_ids if h in native_ids]
            if overlap:
                print(f"  Excluding {len(overlap)} homolog ID(s) already covered by native "
                      f"PDB structures: {', '.join(overlap)}")
                homolog_ids = [h for h in homolog_ids if h not in native_ids]

        cat_pocket_name = ref_pockets_catalytic.get(gene)
        noncat_pocket_name = ref_pockets_noncatalytic.get(gene)
        print(f"  Pocket (catalytic): {cat_pocket_name or 'not defined'}")
        print(f"  Pocket (non-catalytic): {noncat_pocket_name or 'not defined'}")
        build_session(gene, uniprot_ac, cat_pocket_name, noncat_pocket_name, pocket_data, pdb_refs, homolog_ids)


if __name__ == "__main__":
    main()
