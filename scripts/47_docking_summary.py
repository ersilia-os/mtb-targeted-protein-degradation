#!/usr/bin/env python3
"""
Print a docking score summary table for all pocket structures of a given gene.

Usage:
    python 47_docking_summary.py --trna pheS
    python 47_docking_summary.py --trna pheS,aspS --lib DL
    python 47_docking_summary.py --trna pheS --lib REAL
    python 47_docking_summary.py --trna pheS --pymol
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, Superimposer

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

# Same structure color as scripts/47b_reference_pocket_visualization.py, for visual
# consistency across the project's PyMOL outputs.
COLOR_STRUCTURE = [0.7804, 0.8275, 0.8667]
COLOR_POCKET = [1.0, 0.0, 1.0]  # magenta, for pocket centroid spheres
ALIGNED_DIR = os.path.join(ROOT, "output", "aligned_relaxed_structures")
DETECTED_POCKETS_DIR = os.path.join(ROOT, "output", "detected_pockets")
PYMOL_OUTPUT_DIR = os.path.join(ROOT, "output", "47_docking_summary")

LIBRARIES = {
    "DL":   os.path.join(ROOT, "output", "unidock_docking",        "docking_results"),
    "REAL": os.path.join(ROOT, "output", "unidock_REAL_docking_2", "docking_results"),
}

PERCENTILES = [0.01, 0.1, 1]

SWISSMODEL_DIR = os.path.join(ROOT, "data", "structures", "swissmodel")
PLDDT_TYPES = {"alphafold2", "alphafold3", "chai1"}


def load_gene_map():
    path = os.path.join(ROOT, "data", "mtb_trna_synthetases_bosch_2021_fig5.csv")
    df = pd.read_csv(path)
    return {row["gene_name_in_bosch_2021"]: row["uniprot_ac"] for _, row in df.iterrows()}


def load_vulnerability(gene):
    path = os.path.join(ROOT, "data", "bosch_2021_DataS2.xlsx")
    df = pd.read_excel(path, sheet_name="(1) Mtb H37Rv")
    row = df[df["name"] == gene]
    if row.empty:
        return None, None, None
    row = row.iloc[0]
    return row["Vulnerability Index"], row["VI Lower Bound"], row["VI Upper Bound"]


DOMAIN_TO_CODE = {
    "Catalytic Domain (ATP Binding Site)":          "CAT",
    "Anticodon Binding Domain":                     "ABD",
    "Editing Domain":                               "ED",
    "tRNA Binding Domain":                          "tRNA-BD",
    "Other too broad/unspecified functional entities": "Other",
    "Other":                                        "Other",
    "NA":                                           "NA",
}


def load_pocket_domains():
    path = os.path.join(ROOT, "output", "pocket_detection_data_interpro.tsv")
    df = pd.read_csv(path, sep="\t")
    domains = {}
    for _, row in df.iterrows():
        key = f"{row['File name'].replace('.pdb', '')}_pocket_{row['Pocket number']}"
        code = DOMAIN_TO_CODE.get(row["Interpro curated annotation"], "Other")
        domains.setdefault(key, set()).add(code)
    return {k: ";".join(sorted(v - {"Other", "NA"}) or sorted(v)) for k, v in domains.items()}


def load_pocket_data():
    path = os.path.join(ROOT, "output", "pocket_detection_data.csv")
    df = pd.read_csv(path)
    probs, coords, pred_types, min_plddt = {}, {}, {}, {}
    for _, row in df.iterrows():
        key = f"{row['File name'].replace('.pdb', '')}_pocket_{row['Pocket number']}"
        probs[key] = row["Pocket probability"]
        coords[key] = np.array(row["Pocket centroid coordinate (x y z)"].split(), dtype=float)
        pred_types[key] = row["Prediction type"]
        if row["Prediction type"] in PLDDT_TYPES:
            bfactors = [float(x) for x in str(row["B-factors"]).split()]
            min_plddt[key] = round(min(bfactors), 1) if bfactors else None
        else:
            min_plddt[key] = None
    return probs, coords, pred_types, min_plddt


def load_gmqe(uniprot_ac, model_idx):
    """Return GMQE from REMARK 3 of the raw Swiss-Model PDB, or None if not found."""
    try:
        entries = [e for e in os.listdir(SWISSMODEL_DIR) if uniprot_ac in e]
    except FileNotFoundError:
        return None
    if not entries:
        return None
    pdb_path = os.path.join(SWISSMODEL_DIR, entries[0], "models", f"{model_idx + 1:02d}", "model.pdb")
    if not os.path.isfile(pdb_path):
        return None
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("REMARK   3  GMQE"):
                try:
                    return round(float(line.split()[-1]), 3)
                except (ValueError, IndexError):
                    return None
    return None


SOLVENT = {"HOH", "WAT", "DOD"}
ALPHAFILL_THR = 6.0


def alphafill_ligands(uniprot_ac, pocket_coords):
    """
    Align the AlphaFill CIF for uniprot_ac onto the AF3 aligned structure,
    then return a dict {pocket_name: "Yes"/"No"} for HETATM proximity within ALPHAFILL_THR Å.
    Returns None if either file is missing.
    """
    af3_path = os.path.join(ROOT, "output", "aligned_relaxed_structures", uniprot_ac,
                            f"alphafold3_{uniprot_ac}_model_0.pdb")
    cif_path = os.path.join(ROOT, "data", "structures", "alphafill_database",
                            uniprot_ac, f"{uniprot_ac}.cif")
    if not os.path.isfile(af3_path) or not os.path.isfile(cif_path):
        return None

    pdb = PDBParser(QUIET=True).get_structure("ref", af3_path)
    cif = MMCIFParser(QUIET=True).get_structure("af", cif_path)

    # Build Cα maps for chain A
    pdb_ca = {r.id[1]: r["CA"] for r in pdb[0]["A"].get_residues()
              if r.id[0] == " " and "CA" in r}
    cif_ca = {r.id[1]: r["CA"] for r in cif[0]["A"].get_residues()
              if r.id[0] == " " and "CA" in r}
    common = sorted(set(pdb_ca) & set(cif_ca))
    if len(common) < 10:
        return None

    sup = Superimposer()
    sup.set_atoms([pdb_ca[i] for i in common], [cif_ca[i] for i in common])
    sup.apply(cif[0].get_atoms())

    # Collect transformed HETATM coordinates per residue
    het_residues = [r for r in cif[0].get_residues()
                    if r.id[0] != " " and r.resname not in SOLVENT]

    result = {}
    for pocket_name, centroid in pocket_coords.items():
        hits = set()
        for r in het_residues:
            for atom in r.get_atoms():
                d = np.linalg.norm(atom.get_vector().get_array() - centroid)
                if d <= ALPHAFILL_THR:
                    hits.add(r.resname)
                    break
        result[pocket_name] = "Yes" if hits else "No"
    return result


def collect_rows(uniprot_ac, lib_paths):
    # One row per unique pocket name
    all_pockets = set()
    for results_dir in lib_paths.values():
        if os.path.isdir(results_dir):
            all_pockets.update(p for p in os.listdir(results_dir) if uniprot_ac in p)

    pocket_probs, pocket_coords, pred_types, pocket_min_plddt = load_pocket_data()
    pocket_domains = load_pocket_domains()

    # AlphaFill: run once per gene (uniprot_ac), keyed by pocket name
    gene_coords = {p: pocket_coords[p] for p in pocket_coords if uniprot_ac in p}
    af_ligands = alphafill_ligands(uniprot_ac, gene_coords)  # None if files missing

    # GMQE cache: (uniprot_ac, model_idx) -> value
    gmqe_cache = {}

    rows = []
    for pocket in sorted(all_pockets):
        ptype = pred_types.get(pocket)

        # GMQE: only for swissmodel structures
        if ptype == "swissmodel":
            structure_name = pocket.rsplit("_pocket_", 1)[0]
            try:
                model_idx = int(structure_name.split("_model_")[-1])
            except (ValueError, IndexError):
                model_idx = 0
            cache_key = (uniprot_ac, model_idx)
            if cache_key not in gmqe_cache:
                gmqe_cache[cache_key] = load_gmqe(uniprot_ac, model_idx)
            gmqe = gmqe_cache[cache_key]
        else:
            gmqe = None

        row = {
            "pocket":     pocket,
            "prob":       pocket_probs.get(pocket, None),
            "domain":     pocket_domains.get(pocket, "NA"),
            "alphafill":  af_ligands[pocket] if af_ligands is not None else "N/A",
            "min_plddt":  pocket_min_plddt.get(pocket),
            "gmqe":       gmqe,
        }
        for lib_name, results_dir in lib_paths.items():
            report = os.path.join(results_dir, pocket, "report.csv")
            if os.path.isfile(report):
                scores = pd.read_csv(report)["score"].values
                row[f"{lib_name}_n"] = len(scores)
                for p in PERCENTILES:
                    row[f"{lib_name}_p{p}"] = round(np.percentile(scores, p), 3)
            else:
                row[f"{lib_name}_n"] = None
                for p in PERCENTILES:
                    row[f"{lib_name}_p{p}"] = None
        rows.append(row)
    return rows


def fmt(val, default="N/A"):
    return str(val) if val is not None else default


def print_table(rows, gene, lib_names, vi=None, vi_lo=None, vi_hi=None):
    if not rows:
        print(f"No docking results found for gene '{gene}'.")
        return

    rows = sorted(rows, key=lambda r: r["prob"] if r["prob"] is not None else 0, reverse=True)

    best_row = rows[0]

    # Column widths
    col_pocket   = "Structure"
    col_prob     = "Prob"
    col_domain   = "Domain"
    col_alphafill = "AlphaFill"
    col_plddt    = "Min pLDDT"
    col_gmqe     = "GMQE"
    w_pocket    = max(len(col_pocket),    max(len(r["pocket"])                    for r in rows))
    w_prob      = max(len(col_prob),      max(len(fmt(r["prob"]))                 for r in rows))
    w_domain    = max(len(col_domain),    max(len(fmt(r["domain"]))               for r in rows))
    w_alphafill = max(len(col_alphafill), max(len(fmt(r["alphafill"]))            for r in rows))
    w_plddt     = max(len(col_plddt),     max(len(fmt(r["min_plddt"], "-"))       for r in rows))
    w_gmqe      = max(len(col_gmqe),      max(len(fmt(r["gmqe"], "-"))            for r in rows))

    # Per-library column specs: (header, field, width) — N excluded from table
    lib_spec = {}
    for lib in lib_names:
        lib_spec[lib] = [
            (f"{lib} p0.01",  f"{lib}_p0.01", max(len(f"{lib} p0.01"), max(len(fmt(r[f"{lib}_p0.01"]))  for r in rows))),
            (f"{lib} p0.1",   f"{lib}_p0.1",  max(len(f"{lib} p0.1"),  max(len(fmt(r[f"{lib}_p0.1"]))   for r in rows))),
            (f"{lib} p1",     f"{lib}_p1",    max(len(f"{lib} p1"),    max(len(fmt(r[f"{lib}_p1"]))      for r in rows))),
        ]

    sep_parts = [f"+-{'-'*w_pocket}-+-{'-'*w_prob}-+-{'-'*w_domain}-+-{'-'*w_alphafill}-+-{'-'*w_plddt}-+-{'-'*w_gmqe}-+"]
    hdr_parts = [f"| {col_pocket:<{w_pocket}} | {col_prob:<{w_prob}} | {col_domain:<{w_domain}} | {col_alphafill:<{w_alphafill}} | {col_plddt:<{w_plddt}} | {col_gmqe:<{w_gmqe}} |"]
    for lib in lib_names:
        for lbl, _, w in lib_spec[lib]:
            sep_parts.append(f"-{'-'*w}-+")
            hdr_parts.append(f" {lbl:<{w}} |")
    sep = "".join(sep_parts)
    hdr = "".join(hdr_parts)

    title_line = f"  Gene: {gene}  "
    print(sep)
    print(f"|{title_line:^{len(sep)-2}}|")
    print(sep)
    print(hdr)
    print(sep)

    for r in rows:
        parts = [f"| {r['pocket']:<{w_pocket}} | {fmt(r['prob']):<{w_prob}} | {fmt(r['domain']):<{w_domain}} | {fmt(r['alphafill']):<{w_alphafill}} | {fmt(r['min_plddt'], '-'):<{w_plddt}} | {fmt(r['gmqe'], '-'):<{w_gmqe}} |"]
        for lib in lib_names:
            for _, field, w in lib_spec[lib]:
                parts.append(f" {fmt(r[field]):<{w}} |")
        print("".join(parts))

    print(sep)
    print(f"  {len(rows)} structures")
    for lib in lib_names:
        n = next((r[f"{lib}_n"] for r in rows if r.get(f"{lib}_n") is not None), None)
        if n is not None:
            print(f"  {lib}: {n} compounds")
    # Consistency check: is the highest-prob pocket also the best scorer across all lib/percentile combos?
    consistent = all(
        best_row.get(f"{lib}_p{p}") == min((r[f"{lib}_p{p}"] for r in rows if r.get(f"{lib}_p{p}") is not None), default=None)
        for lib in lib_names
        for p in PERCENTILES
    )
    print(f"  Highest prob = best scores: {'Yes' if consistent else 'No'}")
    if vi is not None:
        print(f"  Vulnerability Index: {round(vi, 3)} [{round(vi_lo, 3)}, {round(vi_hi, 3)}]")


def build_pocket_overview_session(gene, uniprot_ac, rows):
    """Build a PyMOL session with one merged object per POCKET — a copy of that pocket's
    structure cartoon plus its own centroid sphere, combined into a single object named
    after the pocket (e.g. "swissmodel_P9WFU3_model_0_pocket_1"). A structure with multiple
    reported pockets is loaded once and copied per pocket, so the structure is intentionally
    duplicated across pockets rather than shared by a single merged object."""
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*converter.*already registered.*")
    import pymol
    from pymol import cmd

    os.makedirs(PYMOL_OUTPUT_DIR, exist_ok=True)

    pymol.finish_launching(["pymol", "-cq"])
    cmd.reinitialize()
    cmd.feedback("disable", "executive", "actions")  # silence "Executive: object X created." spam from cmd.copy
    cmd.set_color("structure_color", COLOR_STRUCTURE)
    cmd.set_color("pocket_color", COLOR_POCKET)

    cached_structures = set()

    for r in rows:
        pocket_name = r["pocket"]
        structure_name, pocket_number = pocket_name.rsplit("_pocket_", 1)

        structure_path = os.path.join(ALIGNED_DIR, uniprot_ac, f"{structure_name}.pdb")
        if not os.path.isfile(structure_path):
            print(f"  Warning: structure file not found: {structure_path}, skipping.")
            continue
        pocket_pdb_path = os.path.join(
            DETECTED_POCKETS_DIR, uniprot_ac, structure_name, "pockets", f"pocket_{pocket_number}.pdb"
        )
        if not os.path.isfile(pocket_pdb_path):
            print(f"  Warning: pocket file not found: {pocket_pdb_path}, skipping.")
            continue

        cache_name = f"_cache_{structure_name}"
        if cache_name not in cached_structures:
            cmd.load(structure_path, cache_name)
            cached_structures.add(cache_name)
        struct_part = f"_tmp_struct_{pocket_name}"
        cmd.copy(struct_part, cache_name)
        cmd.show("cartoon", struct_part)
        cmd.color("structure_color", struct_part)
        cmd.set("transparency", 0.3, struct_part)

        pocket_part = f"_tmp_pocket_{pocket_name}"
        cmd.load(pocket_pdb_path, pocket_part)
        cmd.color("pocket_color", pocket_part)
        cmd.show("spheres", pocket_part)
        cmd.set("sphere_transparency", 0.4, pocket_part)
        cmd.set("sphere_scale", 6, pocket_part)

        cmd.create(pocket_name, f"{struct_part} or {pocket_part}")
        cmd.delete(struct_part)
        cmd.delete(pocket_part)

    for cache_name in cached_structures:
        cmd.delete(cache_name)

    cmd.bg_color("white")
    cmd.set("internal_gui_width", 400)  # wide enough for full object names, e.g. alphafold3_P9WFU3_model_2
    out_path = os.path.join(PYMOL_OUTPUT_DIR, f"session_{uniprot_ac}_{gene}.pse")
    cmd.save(out_path)
    print(f"  Saved pocket overview session: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Docking score table by gene.")
    parser.add_argument("--trna", required=True,
                        help="Comma-separated gene name(s) (e.g. pheS or pheS,aspS)")
    parser.add_argument(
        "--lib",
        choices=["DL", "REAL", "both"],
        default="both",
        help="Library to show: DL, REAL, or both (default: both)",
    )
    parser.add_argument("--pymol", action="store_true",
                        help="Also build a PyMOL session (output/47_docking_summary/session_<AC>_<gene>.pse) "
                             "with one merged object per structure (cartoon + its own pocket spheres) "
                             "for every reported structure/pocket.")
    args = parser.parse_args()

    genes = [g.strip() for g in args.trna.split(",")]
    gene_map = load_gene_map()
    lib_names = ["DL", "REAL"] if args.lib == "both" else [args.lib]
    lib_paths = {lib: LIBRARIES[lib] for lib in lib_names}

    invalid = [g for g in genes if g not in gene_map]
    if invalid:
        print(f"Error: unknown gene(s): {', '.join(invalid)}. Available genes:")
        print("  " + ", ".join(sorted(gene_map)))
        sys.exit(1)

    for gene in genes:
        uniprot_ac = gene_map[gene]
        rows = collect_rows(uniprot_ac, lib_paths)
        vi, vi_lo, vi_hi = load_vulnerability(gene)
        print_table(rows, gene, lib_names, vi, vi_lo, vi_hi)
        if args.pymol:
            build_pocket_overview_session(gene, uniprot_ac, rows)
        print()

    print("NOTE: manually define one CATALYTIC and (optionally) one NON-CATALYTIC reference pocket")
    print("  per tRNA synthetase and record them in:")
    print("    output/reference_pocket_catalytic.csv     (columns: gene_name, pocket_name)")
    print("    output/reference_pocket_noncatalytic.csv  (columns: gene_name, pocket_name)")
    print("Script 47b automatically visualizes whichever of the two is defined for each gene.")


if __name__ == "__main__":
    main()
