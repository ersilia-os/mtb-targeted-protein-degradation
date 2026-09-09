"""
Extended Data Table 3: one row per compound in the 1,095-compound filtered-hits set (script 70) -
identity, physchem/druglikeness, liabilities, on-target docking/Boltz-2 scores against the 12
hand-curated Mtb pockets, off-target Nesso-1 co-folding affinity against the full human (38-gene)
and Mtb (21-gene) counter-screens, and selection provenance.

Sources, all already computed elsewhere in the project - this script joins/reformats them, it
doesn't recompute anything:
- output/70_filtering/filtered_hits.csv: compound_id, smiles, source, and every physchem/ADMET/
  liability column (MW/cLogP/TPSA/HBD/HBA/RotBonds/AromaticRings/QED, PAINS/Brenk/known-antibiotic-
  scaffold-motif flags, SPS/NSPS, mycomembrane permeation, cytotoxicity).
- "Source library": derived from filtered_hits.csv's own "source" column (which encodes
  SELECTION METHOD, e.g. "CAT-promiscuous", not library) rather than read directly - Hit Locator
  (100k) compounds never appear in this set at all (HL was only used for surrogate-model training,
  never final hit selection - confirmed), so exactly 2 libraries are possible: "REAL 9.56M" for
  source == "NONCAT-top100-10M" (scripts 56/61 explicitly split 10M vs 10B for that one selection
  method), "REAL 9.92B" for every other source value (CAT-promiscuous/CAT-selective/NONCAT-
  promiscuous/NONCAT-top100-10B all hardcode LIB = "REAL", scripts 52/53/54, which is the round-2/
  10B library - confirmed by reading those scripts directly).
- output/98_compound_docking_summary/compound_docking_summary_condensed.csv: the 12 RAW (not
  target-group-aggregated) curated-pocket docking_<gene>_<CAT|NONCAT>[_n] / boltz2_<gene>_
  <CAT|NONCAT>[_n] columns - each docking value is already a mean across script 65's 5 Uni-Dock
  replicates (user-confirmed as the intended granularity, over the coarser 8-column pheST/aspS/
  lysS/alaS x CAT/NONCAT group-level file). Boltz-2 values already have script 98's low-confidence
  gating applied upstream (NaN wherever affinity_probability_binary < BOLTZ2_MIN_PROBABILITY,
  0.4 - a project-wide policy decided from the real probability distribution, reused here
  unchanged rather than picking a new cutoff for this table specifically - user-confirmed).
- output/98_compound_docking_summary/compound_docking_summary.csv: the full (non-condensed)
  per-gene nesso1_human_<GENE> (38 columns) / nesso1_mtb_<GENE> (21 columns) Nesso-1-predicted
  IC50 (nM) columns - user-confirmed full per-gene detail over the compact top1/5/10 summary.
- output/101_merge_selections/audit_input.csv: "Prioritized" (bool) and "Origins" (catalytic/
  non-catalytic/both, blank if not prioritized). "Type" (single/dual/multi) is derived here from
  that file's own "selected_targets" (pipe-joined unique gene list): 1 distinct target -> single,
  2 -> dual, 3+ -> multi; blank for a non-prioritized compound (no targets to count).

Rows are sorted by Compound ID for a stable, citable order.

Usage:
    python ExtendedDataTable3.py
"""
import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "..", "src"))

import pandas as pd

from xlsx_utils import save_table_with_legend

output_dir = os.path.join(root, "..", "..", "..", "output")
table_dir = os.path.join(output_dir, "plots", "FigSupp", "ExtendedDataTable3")
os.makedirs(table_dir, exist_ok=True)

FILTERED_HITS_CSV = os.path.join(output_dir, "70_filtering", "filtered_hits.csv")
CONDENSED_SUMMARY_CSV = os.path.join(output_dir, "98_compound_docking_summary",
                                      "compound_docking_summary_condensed.csv")
FULL_SUMMARY_CSV = os.path.join(output_dir, "98_compound_docking_summary",
                                 "compound_docking_summary.csv")
AUDIT_INPUT_CSV = os.path.join(output_dir, "101_merge_selections", "audit_input.csv")

PHYSCHEM_COLS = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotBonds", "AromaticRings", "QED"]
LIABILITY_COLS = ["is_pains", "has_pains", "has_brenk", "is_sim_known_ab"]
MOTIF_COLS = ["nitrofuran_motif", "fluoroquinolone_motif", "carbepenem_motif", "betalactam_motif"]
ADMET_COLS = ["sps_score", "nsps_score", "mycomembrane_permeation",
              "cytotoxicity_hepg2", "cytotoxicity_hskmc", "cytotoxicity_imr90"]

CURATED_POCKET_COLS = [
    "pheS_CAT", "pheS_NONCAT_1", "pheS_NONCAT_2", "pheT_NONCAT",
    "aspS_CAT", "aspS_NONCAT_1", "aspS_NONCAT_2",
    "lysS_CAT", "lysS_NONCAT",
    "alaS_CAT", "alaS_NONCAT_1", "alaS_NONCAT_2",
]
DOCKING_COLS = [f"docking_{c}" for c in CURATED_POCKET_COLS]
BOLTZ2_COLS = [f"boltz2_{c}" for c in CURATED_POCKET_COLS]


def source_library(source):
    return "REAL 9.56M" if source == "NONCAT-top100-10M" else "REAL 9.92B"


def target_type(selected_targets):
    if pd.isna(selected_targets) or not str(selected_targets).strip():
        return ""
    n = len(str(selected_targets).split("|"))
    if n == 1:
        return "single"
    if n == 2:
        return "dual"
    return "multi"


def pocket_phrase(curated_pocket_col):
    """"pheS_CAT" -> "catalytic pocket of pheS"; "pheS_NONCAT_2" -> "non-catalytic pocket #2 of
    pheS" - for the docking/Boltz-2 legend text, one of the 12 hand-curated pheST/aspS/lysS/alaS
    pockets (pheS+pheT sharing the "pheST" complex but scored under their own pheS_/pheT_ prefix
    here, same as compound_docking_summary_condensed.csv's own column names)."""
    parts = curated_pocket_col.split("_")
    gene, site = parts[0], parts[1]
    site_phrase = "catalytic" if site == "CAT" else "non-catalytic"
    pocket_number = f" #{parts[2]}" if len(parts) == 3 else ""
    return f"{site_phrase} pocket{pocket_number} of {gene}"


def build_legend(nesso1_human_cols, nesso1_mtb_cols):
    explanations = {
        "Compound ID": "Enamine REAL compound identifier (source_key).",
        "SMILES": "Compound structure, SMILES notation.",
        "Source library": "Which Enamine REAL screening round this compound was drawn from - "
            "REAL 9.56M (round-1) or REAL 9.92B (round-2). Hit Locator (100k) compounds never "
            "appear in this set (used only for surrogate-model training, not final hit selection).",
        "MW": "Molecular weight.",
        "cLogP": "Calculated logP (octanol-water partition coefficient).",
        "TPSA": "Topological polar surface area.",
        "HBD": "Number of hydrogen bond donors.",
        "HBA": "Number of hydrogen bond acceptors.",
        "RotBonds": "Number of rotatable bonds.",
        "AromaticRings": "Number of aromatic rings.",
        "QED": "Quantitative Estimate of Drug-likeness (0-1, higher is more drug-like).",
        "is_pains": "Whether the compound matches a PAINS (Pan-Assay Interference) substructure "
            "filter.",
        "has_pains": "Whether the compound has any PAINS alert (may differ from is_pains by the "
            "specific filter set applied).",
        "has_brenk": "Whether the compound matches a Brenk structural-alert filter (reactive/"
            "unstable/toxic substructures).",
        "is_sim_known_ab": "Whether the compound is structurally similar to a known antibiotic.",
        "nitrofuran_motif": "Whether the compound contains a nitrofuran substructure (a known "
            "antibiotic scaffold).",
        "fluoroquinolone_motif": "Whether the compound contains a fluoroquinolone substructure.",
        "carbepenem_motif": "Whether the compound contains a carbapenem substructure.",
        "betalactam_motif": "Whether the compound contains a beta-lactam substructure.",
        "sps_score": "Synthetic Proximity Score (raw).",
        "nsps_score": "Normalized Synthetic Proximity Score.",
        "mycomembrane_permeation": "Predicted mycomembrane permeation (M. tuberculosis outer "
            "membrane penetration).",
        "cytotoxicity_hepg2": "Predicted cytotoxicity against the HepG2 (liver) cell line.",
        "cytotoxicity_hskmc": "Predicted cytotoxicity against the HSkMC (skeletal muscle) cell "
            "line.",
        "cytotoxicity_imr90": "Predicted cytotoxicity against the IMR-90 (lung fibroblast) cell "
            "line.",
        "Prioritized": "Whether this compound was one of the final prioritized hits (scripts "
            "99/100's catalytic/non-catalytic selection, plus 4 hand-picked near-misses).",
        "Origins": "For a prioritized compound, whether it qualified via catalytic selection, "
            "non-catalytic selection, or both. Empty if not prioritized.",
        "Type": "For a prioritized compound, how many distinct Mtb targets (genes) it was "
            "selected for: single (1), dual (2), or multi (3+). Empty if not prioritized.",
    }
    for col in DOCKING_COLS:
        pocket_col = col.replace("docking_", "", 1)
        explanations[col] = (
            f"Uni-Dock docking score (mean of 5 replicates) against the {pocket_phrase(pocket_col)}, "
            "one of the 12 hand-curated Mtb pockets."
        )
    for col in BOLTZ2_COLS:
        pocket_col = col.replace("boltz2_", "", 1)
        explanations[col] = (
            f"Boltz-2-predicted binding affinity (IC50, nM) against the {pocket_phrase(pocket_col)}, "
            "one of the 12 hand-curated Mtb pockets. NaN where Boltz-2's own predicted confidence "
            "(affinity_probability_binary) is below 0.4 (low-confidence gate, project-wide policy)."
        )
    for col in nesso1_human_cols:
        gene = col.replace("nesso1_human_", "", 1)
        explanations[col] = (
            f"Nesso-1-predicted binding affinity (IC50, nM) against the human off-target {gene}, "
            "protein-level (no pocket conditioning) co-folding counter-screen."
        )
    for col in nesso1_mtb_cols:
        gene = col.replace("nesso1_mtb_", "", 1)
        explanations[col] = (
            f"Nesso-1-predicted binding affinity (IC50, nM) against the Mtb off-target {gene} "
            "(AF2-monomer counter-screen), protein-level (no pocket conditioning)."
        )
    return explanations


def main():
    filtered_hits = pd.read_csv(FILTERED_HITS_CSV)
    assert len(filtered_hits) == 1095, f"Expected 1,095 compounds, got {len(filtered_hits)}."

    condensed = pd.read_csv(CONDENSED_SUMMARY_CSV).set_index("compound_id")
    full_summary = pd.read_csv(FULL_SUMMARY_CSV).set_index("compound_id")
    nesso1_human_cols = [c for c in full_summary.columns
                         if c.startswith("nesso1_human_") and "top" not in c]
    nesso1_mtb_cols = [c for c in full_summary.columns
                       if c.startswith("nesso1_mtb_") and "top" not in c]
    assert len(nesso1_human_cols) == 38 and len(nesso1_mtb_cols) == 21, (
        f"Expected 38 human/21 Mtb Nesso-1 columns, got {len(nesso1_human_cols)}/{len(nesso1_mtb_cols)}.")

    audit_input = pd.read_csv(AUDIT_INPUT_CSV).set_index("compound_id")

    ids = filtered_hits["compound_id"]
    columns = {
        "Compound ID": ids,
        "SMILES": filtered_hits["smiles"],
        "Source library": filtered_hits["source"].map(source_library),
    }
    for col in PHYSCHEM_COLS + LIABILITY_COLS + MOTIF_COLS + ADMET_COLS:
        columns[col] = filtered_hits[col].values
    for col in DOCKING_COLS + BOLTZ2_COLS:
        columns[col] = ids.map(condensed[col]).values
    for col in nesso1_human_cols + nesso1_mtb_cols:
        columns[col] = ids.map(full_summary[col]).values
    columns["Prioritized"] = ids.map(audit_input["prioritized"]).fillna(False).values
    columns["Origins"] = ids.map(audit_input["origins"]).fillna("").values
    columns["Type"] = ids.map(audit_input["selected_targets"]).apply(target_type).values

    table = pd.DataFrame(columns)

    table = table.sort_values("Compound ID").reset_index(drop=True)

    output_path = os.path.join(table_dir, "ExtendedDataTable3.xlsx")
    save_table_with_legend(table, build_legend(nesso1_human_cols, nesso1_mtb_cols), output_path)
    print(f"Saved {len(table)} row(s), {len(table.columns)} column(s), to {output_path}")
    print(f"Source library counts:\n{table['Source library'].value_counts()}")
    print(f"Prioritized: {table['Prioritized'].sum()} of {len(table)}")
    print(f"Type counts (prioritized only):\n{table.loc[table['Prioritized'], 'Type'].value_counts()}")


if __name__ == "__main__":
    main()
