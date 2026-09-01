#!/usr/bin/env python3
### CAUTION: Step 3 (AlphaFill ligand evidence) needs PyMOL -- run in the adda4tb conda env.
"""
Annotates the human pockets (output/91_human_detect_pockets/pocket_detection_data.csv) for the 38
human aaRS genes, mimicking script 77's Mtb approach but condensed into one script and scoped to
exactly three steps -- no catalytic_confidence scoring, no direct-PDB evidence, no clustering.

1) Download InterPro data per gene (skip if already cached).
2) Map pockets to InterPro domain categories (Catalytic Domain / Anticodon Binding / Editing /
   tRNA Binding / RNA Binding / Other), reusing script 77's GO-term map and name-fallback patterns.
   Script 77's CURATED_OVERRIDES/STRONG_NAME_OVERRIDE_PATTERNS are NOT reused here -- those were
   literature-verified fixes for specific Mtb InterPro entries and don't apply to human proteins.
3) Map AlphaFill-transplanted ligands onto pockets, via a per-pocket local PyMOL alignment (cmd.align,
   restricted to that pocket's own residues) between the AlphaFold DB structure used for pocket
   detection and AlphaFill's structure for the same protein. This local alignment is required: a
   whole-chain fit between the two disagrees by ~34-46A for a real human aaRS (AlphaFold-version
   domain motion), confirmed empirically before writing this script, while a pocket-local region
   agrees to <0.3A -- the same phenomenon script 77 solved for Mtb via per-pocket local alignment.
   No strong/weak ligand classification (script 77's ligand_classification.py is curated per-Mtb-
   protein and has no human equivalent yet) -- this stops at raw evidence.

Usage:
    python 92_human_pocket_annotation.py
"""
import json
import os
import re

import numpy as np
import pandas as pd
import requests
from pymol import cmd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

HUMAN_AARS_CSV = os.path.join(ROOT, "data", "human_trna_synthetases_uniprot.csv")
POCKETS_CSV = os.path.join(ROOT, "output", "91_human_detect_pockets", "pocket_detection_data.csv")
ALPHAFILL_DIR = os.path.join(ROOT, "data", "structures", "human_alphafill_database")

INTERPRO_RAW_DIR = os.path.join(ROOT, "processed", "human_pocket_annotation", "interpro_raw")
OUTPUT_DIR = os.path.join(ROOT, "output", "92_human_pocket_annotation")
os.makedirs(INTERPRO_RAW_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

INTERPRO_BASE_URL = "https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{}/?page_size=200"

PROXIMITY_CUTOFF = 10.0
MIN_POCKET_ATOMS_FOR_FIT = 5


# ---------------------------------------------------------------------------
# Step 1: InterPro download (skip if cached)
# ---------------------------------------------------------------------------

def fetch_interpro(uniprot_ac):
    """All InterPro entries for uniprot_ac, cached to disk; skips the API call if already cached."""
    out_path = os.path.join(INTERPRO_RAW_DIR, f"{uniprot_ac}_interpro.json")
    if os.path.isfile(out_path):
        with open(out_path) as f:
            return json.load(f)

    url = INTERPRO_BASE_URL.format(uniprot_ac)
    results = []
    while url:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        results.extend(data["results"])
        url = data.get("next")

    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    return results


# ---------------------------------------------------------------------------
# Step 2: InterPro domain categorization + pocket mapping
# ---------------------------------------------------------------------------

# Same GO-id -> category map as script 77's 03_categorize.py (organism-independent: aminoacylation
# GO terms are conserved). GO:0000166 ("nucleotide binding") deliberately excluded, same rationale.
GO_ID_CATEGORY = {
    "GO:0004812": "Catalytic Domain (ATP/ligase)",   # aminoacyl-tRNA ligase activity
    "GO:0005524": "Catalytic Domain (ATP/ligase)",   # ATP binding
    "GO:0043039": "Catalytic Domain (ATP/ligase)",   # tRNA aminoacylation
    "GO:0006418": "Catalytic Domain (ATP/ligase)",   # tRNA aminoacylation for protein translation
    "GO:0043771": "Anticodon Binding Domain",        # anticodon binding
    "GO:0002161": "Editing Domain",                  # aminoacyl-tRNA editing activity
    "GO:0000049": "tRNA Binding Domain",              # tRNA binding
    "GO:0003723": "RNA Binding Domain",               # RNA binding
    "GO:0050567": "Catalytic Domain (ATP/ligase)",   # glutaminyl-tRNA synthase (glutamine-hydrolyzing) activity
    "GO:0016884": "Catalytic Domain (ATP/ligase)",   # carbon-nitrogen ligase activity
    "GO:0016874": "Catalytic Domain (ATP/ligase)",   # ligase activity
}
CATALYTIC_NAME_PATTERN = re.compile(r"-tRNA (ligase|synthetase) activity", re.IGNORECASE)
CATEGORY_PRIORITY = [
    "Catalytic Domain (ATP/ligase)",
    "Anticodon Binding Domain",
    "Editing Domain",
    "tRNA Binding Domain",
    "RNA Binding Domain",
]
# Broad, last-resort name fallback (only used when an entry has no informative GO term at all) --
# same as script 77, requires "synthetase"/"ligase" to co-occur with "tRNA"/"transfer RNA".
NAME_FALLBACK_PATTERNS = [
    (re.compile(r"[a-z]_core\b", re.IGNORECASE), "Catalytic Domain (ATP/ligase)"),
    (re.compile(r"(tRNA|transfer RNA)[\s-]*(synthetase|ligase)|aminoacyl.{0,3}tRNA", re.IGNORECASE), "Catalytic Domain (ATP/ligase)"),
]
AARS_CLASS_II_PATTERN = re.compile(r"\bclass[- ]?II\b", re.IGNORECASE)
AARS_CLASS_I_PATTERN = re.compile(r"\bclass[- ]?I\b(?!\s*and\s*II)", re.IGNORECASE)

CATEGORY_COLUMNS = {
    "Catalytic Domain (ATP/ligase)": "Catalytic_Domain_support",
    "Anticodon Binding Domain": "Anticodon_Binding_support",
    "Editing Domain": "Editing_Domain_support",
    "tRNA Binding Domain": "tRNA_Binding_support",
    "RNA Binding Domain": "RNA_Binding_support",
}

ANNOTATION_COLUMNS = [
    "accession", "name", "source_database", "type", "specificity",
    "start", "end", "ranges", "coverage", "protein_length",
    "member_db_accession", "member_db_name", "go_terms",
]


def locations_str(entry_protein_locations):
    ranges = []
    for loc in entry_protein_locations:
        for frag in loc["fragments"]:
            ranges.append((frag["start"], frag["end"]))
    return sorted(set(ranges))


def range_str(ranges):
    return ",".join(f"{s}-{e}" for s, e in ranges)


def go_terms_str(go_terms):
    if not go_terms:
        return ""
    return ";".join(f"{g['identifier']}|{g['name']}" for g in go_terms)


def build_annotation_table(results):
    """Per-InterPro-entry table (accession, name, GO terms, residue span) from raw InterPro JSON,
    same fields/logic as script 77's 02_build_annotation_table.py."""
    parsed = {}
    for r in results:
        m = r["metadata"]
        prot = r["proteins"][0]
        ranges = locations_str(prot["entry_protein_locations"])
        if not ranges:
            continue
        parsed[m["accession"]] = {
            "accession": m["accession"], "name": m["name"] or "",
            "source_database": m["source_database"], "type": m["type"],
            "integrated_into": m["integrated"], "go_terms": go_terms_str(m["go_terms"]),
            "start": ranges[0][0], "end": ranges[-1][1], "ranges": range_str(ranges),
            "protein_length": prot["protein_length"],
            "coverage": round(sum(e - s + 1 for s, e in ranges) / prot["protein_length"], 3),
            "member_db_accession": [], "member_db_name": [],
        }

    rows = []
    for acc, entry in parsed.items():
        parent_acc = entry["integrated_into"]
        if parent_acc and parent_acc in parsed:
            parsed[parent_acc]["member_db_accession"].append(acc)
            parsed[parent_acc]["member_db_name"].append(entry["name"])
            continue
        rows.append(entry)

    if not rows:
        return pd.DataFrame(columns=ANNOTATION_COLUMNS)

    for entry in rows:
        entry["member_db_accession"] = ";".join(entry["member_db_accession"])
        entry["member_db_name"] = ";".join(n or "" for n in entry["member_db_name"])
        entry["specificity"] = "specific" if entry["type"] in {"domain", "family", "site"} else "structural/broad"

    return pd.DataFrame(rows)[ANNOTATION_COLUMNS].sort_values("start").reset_index(drop=True)


def categorize_go_terms(go_terms_str_val):
    if not go_terms_str_val:
        return None, []
    hits = []
    for term in go_terms_str_val.split(";"):
        go_id, go_name = term.split("|", 1)
        if go_id in GO_ID_CATEGORY:
            hits.append((go_id, go_name, GO_ID_CATEGORY[go_id]))
        elif CATALYTIC_NAME_PATTERN.search(go_name):
            hits.append((go_id, go_name, "Catalytic Domain (ATP/ligase)"))
    if not hits:
        return None, []
    categories_hit = {h[2] for h in hits}
    for cat in CATEGORY_PRIORITY:
        if cat in categories_hit:
            return cat, hits
    return None, hits


def categorize_by_name_fallback(name):
    for pattern, category in NAME_FALLBACK_PATTERNS:
        if pattern.search(name):
            return category
    return None


def categorize(name, go_terms_val):
    cat, hits = categorize_go_terms(go_terms_val)
    if cat is not None:
        return cat, hits, "GO term"
    cat = categorize_by_name_fallback(name)
    if cat is not None:
        return cat, [], "name fallback"
    return "Other/Unclassified", [], "none"


def detect_aars_class(entry_names):
    is_class_i = any(AARS_CLASS_I_PATTERN.search(n) for n in entry_names)
    is_class_ii = any(AARS_CLASS_II_PATTERN.search(n) for n in entry_names)
    if is_class_i and is_class_ii:
        return "Ambiguous (both class I and II entries found)"
    if is_class_i:
        return "Class I"
    if is_class_ii:
        return "Class II"
    return "N/A (not a classical aaRS)"


def compute_residue_support(categorized_df):
    """(residue, category) -> entry_support_count: how many independent InterPro entries assign
    that category to that residue. Same weakest-link logic as script 77."""
    support_rows = []
    for _, row in categorized_df.iterrows():
        if row["category"] == "Other/Unclassified":
            continue
        for r in range(row["start"], row["end"] + 1):
            support_rows.append((r, row["category"], row["accession"]))
    if not support_rows:
        return pd.DataFrame(columns=["residue", "category", "entry_support_count"])
    support_df = pd.DataFrame(support_rows, columns=["residue", "category", "accession"])
    return (support_df.groupby(["residue", "category"])["accession"]
            .nunique().reset_index().rename(columns={"accession": "entry_support_count"}))


def parse_pocket_residues(residues_str):
    residues = set()
    for tok in residues_str.split():
        _, resn = tok.split("_")
        residues.add(int(resn))
    return residues


def pocket_labels_and_support(pocket_residues, support_by_category):
    labels = []
    support_col_values = {}
    for category, col in CATEGORY_COLUMNS.items():
        cat_support = support_by_category.get(category, {})
        hits = [cat_support[r] for r in pocket_residues if r in cat_support]
        if hits:
            labels.append(category)
            support_col_values[col] = min(hits)
        else:
            support_col_values[col] = ""
    return labels, support_col_values


# ---------------------------------------------------------------------------
# Step 3: AlphaFill ligand evidence via per-pocket local alignment
# ---------------------------------------------------------------------------

def load_alphafill_transplant_meta(json_path):
    """analogue_id -> [{pdb_id, local_rmsd, identity}, ...], from AlphaFill's own JSON."""
    with open(json_path) as f:
        meta = json.load(f)
    transplant_meta = {}
    for hit in meta.get("hits") or []:
        pdb_id = hit.get("pdb_id")
        identity = (hit.get("alignment") or {}).get("identity")
        for t in (hit.get("transplants") or []):
            analogue_id = t.get("analogue_id")
            if not analogue_id:
                continue
            transplant_meta.setdefault(analogue_id, []).append({
                "pdb_id": pdb_id, "local_rmsd": t.get("local_rmsd"), "identity": identity,
            })
    return transplant_meta


def pocket_local_align(mobile_sel_base, ref_sel_base, pocket_residues):
    """Superposition-only (cycles=0) cmd.align restricted to this pocket's own CA atoms. Returns
    (rmsd, n_atoms), or None if fewer than MIN_POCKET_ATOMS_FOR_FIT shared CA atoms exist."""
    resi_str = "+".join(str(r) for r in sorted(pocket_residues))
    mobile_sel = f"{mobile_sel_base} and resi {resi_str} and name CA and polymer.protein"
    ref_sel = f"{ref_sel_base} and resi {resi_str} and name CA and polymer.protein"
    n_mobile = cmd.select("tmp_check", mobile_sel)
    cmd.delete("tmp_check")
    if n_mobile < MIN_POCKET_ATOMS_FOR_FIT:
        return None
    r = cmd.align(mobile_sel, ref_sel, cycles=0)
    return r[0], r[1]


def extract_alphafill_evidence(pockets_df):
    evidence_rows = []
    for uniprot_ac, protein_pockets in pockets_df.groupby("Uniprot AC"):
        cif_path = os.path.join(ALPHAFILL_DIR, uniprot_ac, f"{uniprot_ac}.cif")
        json_path = os.path.join(ALPHAFILL_DIR, uniprot_ac, f"{uniprot_ac}.json")
        ref_path = os.path.join(ROOT, protein_pockets["Full path"].iloc[0])
        if not (os.path.isfile(cif_path) and os.path.isfile(json_path) and os.path.isfile(ref_path)):
            print(f"  WARNING: missing AlphaFill/reference data for {uniprot_ac}, skipping")
            continue

        transplant_meta = load_alphafill_transplant_meta(json_path)

        cmd.reinitialize()
        cmd.load(ref_path, "ref")
        cmd.load(cif_path, "af")

        n_fit_ok, n_fit_skip = 0, 0
        for _, prow in protein_pockets.iterrows():
            pocket_residues = parse_pocket_residues(prow["P2Rank residues (chain_resn)"])
            fit = pocket_local_align("af and chain A", "ref and chain A", pocket_residues)
            if fit is None:
                n_fit_skip += 1
                continue
            n_fit_ok += 1
            fit_rmsd, fit_n = fit

            # re-fetch ligand coordinates AFTER this pocket's own alignment -- cmd.align moves
            # "af"'s real atoms in place (same gotcha as script 77's 08_extract_alphafill_evidence.py)
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
                    provenance = (transplant_meta.get(resn) or [{}])[0]
                    evidence_rows.append({
                        "uniprot_ac": uniprot_ac, "ligand_analogue_id": resn, "ligand_chain": ch,
                        "ligand_resi": resi, "pocket_file_name": prow["File name"],
                        "pocket_number": int(prow["Pocket number"]), "distance": round(d, 2),
                        "fit_rmsd": round(fit_rmsd, 2), "fit_n_atoms": fit_n,
                        "source_pdb": provenance.get("pdb_id"), "local_rmsd": provenance.get("local_rmsd"),
                        "homolog_identity": provenance.get("identity"),
                    })

        print(f"  {uniprot_ac}: {n_fit_ok} pockets fit OK, {n_fit_skip} skipped")

    return pd.DataFrame(evidence_rows)


# ---------------------------------------------------------------------------

def main():
    proteins = pd.read_csv(HUMAN_AARS_CSV)
    pockets_df = pd.read_csv(POCKETS_CSV)

    print("=== Step 1+2: InterPro download & domain-pocket mapping ===")
    domain_rows = []
    for _, prow in proteins.iterrows():
        gene_name, uniprot_ac = prow["gene_name"], prow["uniprot_ac"]
        print(f"--- {gene_name} ({uniprot_ac}) ---")

        results = fetch_interpro(uniprot_ac)
        annotation_df = build_annotation_table(results)
        annotation_df.to_csv(os.path.join(OUTPUT_DIR, f"{uniprot_ac}_annotation_table.csv"), index=False)

        categories, matched_go, sources = [], [], []
        for _, row in annotation_df.iterrows():
            cat, hits, source = categorize(row["name"], row["go_terms"])
            categories.append(cat)
            matched_go.append(";".join(f"{h[0]}|{h[1]}" for h in hits))
            sources.append(source)
        annotation_df["category"] = categories
        annotation_df["category_matched_go_terms"] = matched_go
        annotation_df["category_source"] = sources
        aars_class = detect_aars_class(annotation_df["name"].tolist())
        annotation_df["aars_class"] = aars_class
        annotation_df.to_csv(os.path.join(OUTPUT_DIR, f"{uniprot_ac}_annotation_table_categorized.csv"), index=False)

        support_counts = compute_residue_support(annotation_df)
        support_counts.to_csv(os.path.join(OUTPUT_DIR, f"{uniprot_ac}_residue_support.csv"), index=False)

        support_by_category = {}
        for _, r in support_counts.iterrows():
            support_by_category.setdefault(r["category"], {})[r["residue"]] = r["entry_support_count"]

        gene_pockets = pockets_df[pockets_df["Uniprot AC"] == uniprot_ac]
        for _, prow2 in gene_pockets.iterrows():
            pocket_residues = parse_pocket_residues(prow2["P2Rank residues (chain_resn)"])
            labels, support_vals = pocket_labels_and_support(pocket_residues, support_by_category)
            row = {
                "Uniprot AC": uniprot_ac, "Gene name": gene_name, "File name": prow2["File name"],
                "Pocket number": prow2["Pocket number"], "curated_labels": "|".join(labels),
                "aars_class": aars_class,
            }
            row.update(support_vals)
            domain_rows.append(row)

    domain_df = pd.DataFrame(domain_rows)
    domain_out = os.path.join(OUTPUT_DIR, "pocket_domain_labels.csv")
    domain_df.to_csv(domain_out, index=False)
    n_catalytic = domain_df["curated_labels"].str.contains("Catalytic Domain", na=False).sum()
    print(f"\nSaved {len(domain_df)} pocket domain-label rows -> {domain_out}")
    print(f"{n_catalytic}/{len(domain_df)} pockets carry a Catalytic Domain (ATP/ligase) label")

    print("\n=== Step 3: AlphaFill ligand evidence ===")
    evidence_df = extract_alphafill_evidence(pockets_df)
    evidence_out = os.path.join(OUTPUT_DIR, "alphafill_ligand_evidence.csv")
    evidence_df.to_csv(evidence_out, index=False)
    print(f"\nSaved {len(evidence_df)} AlphaFill ligand-pocket evidence rows -> {evidence_out}")


if __name__ == "__main__":
    main()
