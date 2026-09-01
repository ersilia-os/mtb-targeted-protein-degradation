#!/usr/bin/env python3
"""
Merges the human pocket-detection table (script 91) with its InterPro domain labels (script 92)
and AlphaFill ligand evidence (script 92) into one wide table, one row per pocket. Still no derived
"is catalytic" verdict column from the InterPro label alone -- it's broad (388/389 pockets), so
collapsing it into a single flag would be an unreviewed threshold decision. What IS added here is a
strong/weak classification per AlphaFill ligand (human equivalent of script 77's
ligand_classification.py, which only covers the 18 Mtb proteins), since not every ligand within the
10A cutoff is equally meaningful evidence -- distinguishing them lets the merged table carry real
signal instead of just raw ligand IDs.

Classification tiers per (gene, ligand), each verified against RCSB's Chemical Component Dictionary
or a literature search rather than guessed (see comments below for what was checked and how):

1. STRONG -- self-structure: the ligand's source_pdb is one of this exact gene's own experimental
   PDB entries (fetched live from UniProt, same mechanism as script 77's 04_query_pdb_xrefs.py).
   Chemistry doesn't need separate classification here: if it's IN that protein's own solved
   structure, it's real self-evidence by construction.
2. STRONG -- universal cofactor/byproduct: ATP/ADP/AMP/APC (same as script 77) plus the
   aminoacylation reaction's actual pyrophosphate byproduct and its non-hydrolyzable mimics
   (POP, PPV = pyrophosphate; 2PN = imidodiphosphoric acid) -- relevant regardless of which aaRS.
3. STRONG -- matching reaction-intermediate mimic: an aminoacyl-sulfamoyl-adenosine/aminoacyl-AMP
   class compound (RCSB-confirmed identity) whose amino acid matches THIS gene's own substrate --
   e.g. A5A (alanyl-...) is strong for AARS1/AARS2, but NOT for GARS1 even though it's the same
   compound class, since alanine isn't glycine's business.
4. STRONG -- literature-documented near-cognate/editing substrate: VRT (norvaline-adenosine
   conjugate) for VARS1/VARS2 -- ValRS is well documented to misactivate the non-proteinogenic
   near-cognate norvaline and rely on its editing domain to hydrolyze it (Nureki/Yokoyama-era ValRS
   fidelity literature); GGB (L-canavanine) for RARS1/RARS2 -- canavanine is a well-documented
   arginine mimic that gets misactivated by ArgRS. Confirmed independently for KRS (cladosporin,
   PDB 4YCU = cladosporin bound to HUMAN KARS1 specifically -- also caught by tier 1 here since
   4YCU is in KARS1's own PDB list) and 2CR (borrelidin, PDB 4P3N-4P3P = borrelidin bound to human
   TARS1 -- also caught by tier 1, TARS1's own PDB list includes 4P3N).
Everything else defaults to WEAK (metal ions, cryoprotectants/buffers, dyes, unrelated cofactors
like NAD/FAD/heme/Fe-S clusters, wrong-amino-acid-specificity homology-transplanted mimics) --
same "default weak unless proven otherwise" convention as script 77's classify_ligand.

Usage:
    python 93_human_merge_pocket_annotations.py
"""
import json
import os

import pandas as pd
import requests

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

HUMAN_AARS_CSV = os.path.join(ROOT, "data", "human_trna_synthetases_uniprot.csv")
POCKETS_CSV = os.path.join(ROOT, "output", "91_human_detect_pockets", "pocket_detection_data.csv")
DOMAIN_LABELS_CSV = os.path.join(ROOT, "output", "92_human_pocket_annotation", "pocket_domain_labels.csv")
LIGAND_EVIDENCE_CSV = os.path.join(ROOT, "output", "92_human_pocket_annotation", "alphafill_ligand_evidence.csv")

PDB_XREF_CACHE_DIR = os.path.join(ROOT, "processed", "human_pocket_annotation", "pdb_xrefs")
OUTPUT_DIR = os.path.join(ROOT, "output", "93_human_merge_pocket_annotations")
os.makedirs(PDB_XREF_CACHE_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

MERGE_KEYS = ["Uniprot AC", "File name", "Pocket number"]
LIGAND_EVIDENCE_GROUP_KEYS = ["uniprot_ac", "pocket_file_name", "pocket_number"]

AGGREGATED_LIGAND_COLUMNS = {
    "ligand_analogue_id": "alphafill_ligand_ids",
    "distance": "alphafill_ligand_distances",
    "fit_rmsd": "alphafill_ligand_fit_rmsd",
    "source_pdb": "alphafill_ligand_sources",
    "local_rmsd": "alphafill_ligand_local_rmsd",
    "homolog_identity": "alphafill_ligand_homolog_identity",
    "strength": "alphafill_ligand_strengths",
    "strength_basis": "alphafill_ligand_strength_basis",
}


# ---------------------------------------------------------------------------
# PDB cross-references (skip if cached) -- same mechanism as script 77's 04_query_pdb_xrefs.py
# ---------------------------------------------------------------------------

def fetch_pdb_xrefs(uniprot_ac):
    """Set of this protein's own experimental PDB IDs, cached to disk."""
    cache_path = os.path.join(PDB_XREF_CACHE_DIR, f"{uniprot_ac}.json")
    if os.path.isfile(cache_path):
        with open(cache_path) as f:
            return set(json.load(f))

    resp = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}.json", timeout=30)
    resp.raise_for_status()
    d = resp.json()
    pdb_ids = sorted({x["id"] for x in d.get("uniProtKBCrossReferences", []) if x["database"] == "PDB"})

    with open(cache_path, "w") as f:
        json.dump(pdb_ids, f, indent=2)
    return set(pdb_ids)


# ---------------------------------------------------------------------------
# Ligand strong/weak classification (human equivalent of script 77's ligand_classification.py)
# ---------------------------------------------------------------------------

COFACTOR_CODES = {"ATP", "ADP", "AMP", "APC"}
REACTION_BYPRODUCT_CODES = {"POP", "PPV", "2PN"}  # pyrophosphate and non-hydrolyzable PPi mimics

# gene -> free amino acid 3-letter PDB code(s) this enzyme charges. EPRS1 is the bifunctional
# cytoplasmic Glu-Pro ligase; EARS2/PARS2 are its separate monofunctional mitochondrial
# counterparts. FARSA/FARSB are the cytoplasmic PheRS alpha/beta heterotetramer subunits.
SUBSTRATE_BY_GENE = {
    "AARS1": {"ALA"}, "AARS2": {"ALA"},
    "CARS1": {"CYS"}, "CARS2": {"CYS"},
    "DARS1": {"ASP"}, "DARS2": {"ASP"},
    "EARS2": {"GLU"}, "EPRS1": {"GLU", "PRO"},
    "FARS2": {"PHE"}, "FARSA": {"PHE"}, "FARSB": {"PHE"},
    "GARS1": {"GLY"},
    "HARS1": {"HIS"}, "HARS2": {"HIS"},
    "IARS1": {"ILE"}, "IARS2": {"ILE"},
    "KARS1": {"LYS"},
    "LARS1": {"LEU"}, "LARS2": {"LEU"},
    "MARS1": {"MET"}, "MARS2": {"MET"},
    "NARS1": {"ASN"}, "NARS2": {"ASN"},
    "PARS2": {"PRO"},
    "QARS1": {"GLN"},
    "RARS1": {"ARG"}, "RARS2": {"ARG"},
    "SARS1": {"SER"}, "SARS2": {"SER"},
    "TARS1": {"THR"}, "TARS2": {"THR"}, "TARS3": {"THR"},
    "VARS1": {"VAL"}, "VARS2": {"VAL"},
    "WARS1": {"TRP"}, "WARS2": {"TRP"},
    "YARS1": {"TYR"}, "YARS2": {"TYR"},
}

# Reaction-intermediate-mimic ligand -> amino acid, confirmed via RCSB Chemical Component
# Dictionary lookups this session (aminoacyl-sulfamoyl-adenosine / aminoacyl-AMP class).
INTERMEDIATE_ANALOG_TO_AA = {
    "A5A": "ALA",   # 5'-O-(N-(L-alanyl)-sulfamoyl)adenosine
    "DSZ": "ASP",   # 5'-O-(L-alpha-aspartylsulfamoyl)adenosine
    "G5A": "GLY",   # 5'-O-(glycylsulfamoyl)adenosine
    "KAA": "LYS",   # 5'-O-[(L-lysylamino)sulfonyl]adenosine
    "LSS": "LEU",   # 5'-O-(L-leucylsulfamoyl)adenosine
    "SSA": "SER",   # 5'-O-(N-(L-seryl)-sulfamoyl)adenosine
    "A3S": "SER",   # SERINE-3'-AMINOADENOSINE
    "TSB": "THR",   # 5'-O-(N-(L-threonyl)-sulfamoyl)adenosine
    "TYM": "TRP",   # tryptophanyl-5'-AMP
    "YSA": "TYR",   # 5'-O-[N-(L-tyrosyl)sulfamoyl]adenosine
    "HAM": "HIS",   # HISTIDYL-ADENOSINE MONOPHOSPHATE
}

# Literature-documented near-cognate/editing-site substrates -- confirmed via web search this
# session, not assumed from chemical similarity alone.
NEAR_COGNATE_ANALOG_BY_GENE = {
    "VRT": {"VARS1", "VARS2"},  # 2'-(L-norvalyl)amino-2'-deoxyadenosine: ValRS misactivates/edits
                                  # the non-proteinogenic near-cognate norvaline (well documented in
                                  # aaRS fidelity/editing-domain literature).
    "GGB": {"RARS1", "RARS2"},  # L-canavanine: a plant-derived arginine mimic well documented to be
                                  # misactivated by ArgRS.
}


def classify_ligand(gene_name, ligand_code, source_pdb, own_pdb_ids):
    source_pdb_upper = str(source_pdb).upper() if pd.notna(source_pdb) else None
    if source_pdb_upper and source_pdb_upper in own_pdb_ids:
        return "strong", "self-structure"
    if ligand_code in COFACTOR_CODES:
        return "strong", "universal cofactor"
    if ligand_code in REACTION_BYPRODUCT_CODES:
        return "strong", "reaction byproduct (pyrophosphate)"
    if INTERMEDIATE_ANALOG_TO_AA.get(ligand_code) in SUBSTRATE_BY_GENE.get(gene_name, set()):
        return "strong", "matching reaction-intermediate mimic"
    if gene_name in NEAR_COGNATE_ANALOG_BY_GENE.get(ligand_code, set()):
        return "strong", "literature-documented near-cognate/editing substrate"
    return "weak", "default (buffer/crystallization/unrelated/wrong-specificity)"


# ---------------------------------------------------------------------------
# Merge
# ---------------------------------------------------------------------------

def classify_ligand_evidence(evidence_df, pdb_xrefs_by_ac):
    evidence_df = evidence_df.copy()
    strengths, bases = [], []
    for _, row in evidence_df.iterrows():
        gene_name = row["gene_name"]
        own_pdb_ids = pdb_xrefs_by_ac.get(row["uniprot_ac"], set())
        strength, basis = classify_ligand(gene_name, row["ligand_analogue_id"], row["source_pdb"], own_pdb_ids)
        strengths.append(strength)
        bases.append(basis)
    evidence_df["strength"] = strengths
    evidence_df["strength_basis"] = bases
    return evidence_df


def aggregate_ligand_evidence(evidence_df):
    """One row per pocket, pipe-joining every ligand match within that pocket's evidence, plus a
    strong/weak count summary."""
    out_columns = (MERGE_KEYS + ["n_alphafill_ligands", "n_strong_ligands", "n_weak_ligands"]
                   + list(AGGREGATED_LIGAND_COLUMNS.values()))
    if evidence_df.empty:
        return pd.DataFrame(columns=out_columns)

    rows = []
    for (uniprot_ac, file_name, pocket_number), sub in evidence_df.groupby(LIGAND_EVIDENCE_GROUP_KEYS):
        row = {
            "Uniprot AC": uniprot_ac, "File name": file_name, "Pocket number": pocket_number,
            "n_alphafill_ligands": len(sub),
            "n_strong_ligands": int((sub["strength"] == "strong").sum()),
            "n_weak_ligands": int((sub["strength"] == "weak").sum()),
        }
        for src_col, out_col in AGGREGATED_LIGAND_COLUMNS.items():
            row[out_col] = "|".join(sub[src_col].astype(str))
        rows.append(row)
    return pd.DataFrame(rows, columns=out_columns)


def main():
    proteins = pd.read_csv(HUMAN_AARS_CSV)
    pockets = pd.read_csv(POCKETS_CSV)
    domain_labels = pd.read_csv(DOMAIN_LABELS_CSV)
    ligand_evidence = pd.read_csv(LIGAND_EVIDENCE_CSV)

    print("Fetching PDB cross-references (skip if cached)...")
    pdb_xrefs_by_ac = {}
    for _, row in proteins.iterrows():
        pdb_xrefs_by_ac[row["uniprot_ac"]] = fetch_pdb_xrefs(row["uniprot_ac"])

    gene_by_ac = dict(zip(proteins["uniprot_ac"], proteins["gene_name"]))
    ligand_evidence["gene_name"] = ligand_evidence["uniprot_ac"].map(gene_by_ac)
    ligand_evidence = classify_ligand_evidence(ligand_evidence, pdb_xrefs_by_ac)
    ligand_evidence.to_csv(os.path.join(OUTPUT_DIR, "alphafill_ligand_evidence_classified.csv"), index=False)

    n_before = len(pockets)

    merged = pockets.merge(
        domain_labels.drop(columns=["Gene name"]), on=MERGE_KEYS, how="left", validate="one_to_one",
    )
    assert len(merged) == n_before, "domain-label merge changed row count"

    ligand_agg = aggregate_ligand_evidence(ligand_evidence)
    merged = merged.merge(ligand_agg, on=MERGE_KEYS, how="left", validate="one_to_one")
    assert len(merged) == n_before, "ligand-evidence merge changed row count"
    for col in ["n_alphafill_ligands", "n_strong_ligands", "n_weak_ligands"]:
        merged[col] = merged[col].fillna(0).astype(int)

    out_path = os.path.join(OUTPUT_DIR, "merged_pocket_data.csv")
    merged.to_csv(out_path, index=False)

    n_with_catalytic_label = merged["curated_labels"].str.contains("Catalytic Domain", na=False).sum()
    n_with_ligand_evidence = (merged["n_alphafill_ligands"] > 0).sum()
    n_with_strong_ligand = (merged["n_strong_ligands"] > 0).sum()
    n_with_label_and_strong = (merged["curated_labels"].str.contains("Catalytic Domain", na=False)
                                & (merged["n_strong_ligands"] > 0)).sum()

    print(f"\nMerged {len(merged)} pockets -> {out_path}")
    print(f"  {n_with_catalytic_label}/{len(merged)} have a Catalytic Domain (ATP/ligase) label")
    print(f"  {n_with_ligand_evidence}/{len(merged)} have >=1 AlphaFill ligand within 10A of centroid")
    print(f"  {n_with_strong_ligand}/{len(merged)} have >=1 STRONG AlphaFill ligand")
    print(f"  {n_with_label_and_strong}/{len(merged)} have BOTH a Catalytic Domain label AND >=1 strong ligand")

    n_strong_total = (ligand_evidence["strength"] == "strong").sum()
    print(f"\nLigand evidence rows: {len(ligand_evidence)} total, {n_strong_total} classified strong")
    print("Strong-evidence basis breakdown:")
    print(ligand_evidence.loc[ligand_evidence["strength"] == "strong", "strength_basis"].value_counts())


if __name__ == "__main__":
    main()
