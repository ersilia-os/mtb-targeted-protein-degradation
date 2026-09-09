"""
Extended Data Table 2: one row per detected pocket (276 total - one per P2Rank pocket per
structure, NOT deduplicated to canonical/physical sites the way figures 2-4 group pockets by
spatial_cluster_id; every raw detection gets its own row here).

Sources, all already computed elsewhere in the project - this script joins/reformats them, it
doesn't recompute pocket detection, domain annotation, or docking itself:
- output/77_pocket_annotation/pocket_detection_interpro_updated.csv: identity (Uniprot AC, Gene,
  File name, Pocket number), P2Rank score/probability, curated domain label(s) (pipe-delimited
  "curated_labels", parsed into 4 bool columns - same 4 domains as ExtendedDataFigure1.py/
  ExtendedDataTable1.py, ignoring this file's own extra RNA_Binding/NonCatalytic support columns),
  catalytic_confidence (0-4), n_models_in_cluster (how many raw detections, across this protein's
  structures, collapsed into this pocket's own canonical/spatial cluster), and direct-PDB/AlphaFill
  ligand evidence - reformatted here from that file's raw (duplicated, unsorted) pipe-delimited
  strings into deduplicated, alphabetically-sorted ones (user request).
- output/pocket_detection_data.csv: "Pocket residues (chain_resn)" (space-delimited), joined in
  for "Number of residues" - not carried into the interpro-annotated file above.
- output/selected_pockets.csv: the hand-curated pheS/pheT/aspS/lysS/alaS pockets figures 3/4
  showcase - "Selected pocket" (bool) plus that curation's own site type (CAT/NON-CAT) and
  comment (user-confirmed: keep the curator's rationale, not just a flag). 7K98_pocket_6 (the
  one selected_pockets.csv row from the experimental-structure pipeline, not the monomeric
  P2Rank pipeline this table covers) has no matching row here and is simply absent, not an error.
- src/docking_utils.py's LIBRARIES/load_scores/load_real_positive_scores: per-pocket docking score
  summaries (n scored, min, 1st/10th percentile, median) for all 3 screening rounds - HLL (Enamine
  Hit Locator 100K), REAL 9.56M (round-1 Enamine REAL screening, PRIORITIZED top-100k "active" set
  only, matching figure_2_calculations.py's own convention over the full screened set - user-
  confirmed), and REAL 9.92B (round-2 Enamine REAL 10B screening). NaN/n=0 for a pocket with no
  report.csv in a given library, not an error - not every pocket was screened in every round.

Rows are sorted by (Gene name, Structure name, Pocket number) for a stable, citable order.

Usage:
    python ExtendedDataTable2.py
"""
import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "..", "src"))

import numpy as np
import pandas as pd

from docking_utils import LIBRARIES, load_scores, load_real_positive_scores
from xlsx_utils import save_table_with_legend

output_dir = os.path.join(root, "..", "..", "..", "output")
table_dir = os.path.join(output_dir, "plots", "FigSupp", "ExtendedDataTable2")
os.makedirs(table_dir, exist_ok=True)

ANNOTATION_CSV = os.path.join(output_dir, "77_pocket_annotation", "pocket_detection_interpro_updated.csv")
POCKET_DETECTION_CSV = os.path.join(output_dir, "pocket_detection_data.csv")
SELECTED_POCKETS_CSV = os.path.join(output_dir, "selected_pockets.csv")

# Same 4 domains as ExtendedDataFigure1.py/ExtendedDataTable1.py's DOMAIN_LABELS - full curated
# label text (as it appears, pipe-delimited, in curated_labels) -> short column name.
DOMAIN_LABELS = {
    "Catalytic Domain (ATP/ligase)": "Catalytic",
    "tRNA Binding Domain": "tRNA binding",
    "Editing Domain": "Editing",
    "Anticodon Binding Domain": "Anticodon",
}

# Docking round -> (library dir for HLL/REAL 9.92B; None for REAL 9.56M, which uses its own
# prioritized-set loader instead of a plain report.csv path) and the 4 summary stats requested.
DOCKING_ROUNDS = ["HLL", "REAL 9.56M", "REAL 9.92B"]
STATS = ["n", "min", "perc1", "perc10", "median"]


def domain_flags(curated_labels):
    labels = set(str(curated_labels).split("|")) if pd.notna(curated_labels) else set()
    return {short: (full in labels) for full, short in DOMAIN_LABELS.items()}


def clean_ligand_list(raw):
    """Pipe-delimited ligand codes as stored (e.g. "TYM|ATP|MG|5BX|ATP|MG|...", with real
    duplicates from multiple co-crystallized copies) -> deduplicated, alphabetically sorted,
    pipe-joined (user request) - empty string (not NaN) when there's no evidence at all."""
    if pd.isna(raw) or not str(raw).strip():
        return ""
    return "|".join(sorted(set(str(raw).split("|"))))


def docking_stats(scores):
    scores = np.asarray(scores, dtype=float)
    if len(scores) == 0:
        return {"n": 0, "min": np.nan, "perc1": np.nan, "perc10": np.nan, "median": np.nan}
    return {
        "n": len(scores),
        "min": round(float(np.min(scores)), 3),
        "perc1": round(float(np.percentile(scores, 1)), 3),
        "perc10": round(float(np.percentile(scores, 10)), 3),
        "median": round(float(np.median(scores)), 3),
    }


def safe_load_scores(library_dir, pocket_name):
    report_path = os.path.join(library_dir, pocket_name, "report.csv")
    if not os.path.isfile(report_path):
        return pd.Series(dtype=float)
    return load_scores(report_path)


STAT_EXPLANATIONS = {
    "n": "number of compounds scored",
    "min": "best (most negative) docking score",
    "perc1": "1st percentile docking score",
    "perc10": "10th percentile docking score",
    "median": "median docking score",
}
ROUND_EXPLANATIONS = {
    "HLL": "Enamine Hit Locator 100K library",
    "REAL 9.56M": "round-1 Enamine REAL screening, prioritized top-100k \"active\" set only",
    "REAL 9.92B": "round-2 Enamine REAL 10B screening",
}


def build_legend():
    explanations = {
        "Pocket name": "Structure name + pocket number, e.g. \"swissmodel_P9WFU3_model_0_pocket_1\" "
                       "- unique identifier for this raw P2Rank pocket detection.",
        "Gene name": "Gene symbol used throughout this project.",
        "Uniprot AC": "UniProt accession code for the M. tuberculosis protein.",
        "Structure name": "Structure this pocket was detected on (source model + model index), "
                           "e.g. \"swissmodel_P9WFU3_model_0\".",
        "Pocket number": "P2Rank's own pocket index on this structure (1-based, ranked by "
                          "P2Rank score).",
        "Number of residues": "Number of residues P2Rank assigned to this pocket.",
        "P2Rank probability": "P2Rank's own predicted probability that this is a true ligand-"
                               "binding pocket.",
        "P2Rank score": "P2Rank's own pocket score (not a probability - a ranking score).",
        "Number of pockets for the same canonical pocket": "How many raw detections (across this "
            "protein's structures) collapsed into this pocket's own canonical/spatial cluster "
            "(6.14 A greedy centroid dedup) - reproducibility across independently-modeled "
            "structures.",
        "Catalytic": "Whether this pocket carries the curated Catalytic Domain (ATP/ligase) label.",
        "tRNA binding": "Whether this pocket carries the curated tRNA Binding Domain label.",
        "Editing": "Whether this pocket carries the curated Editing Domain label.",
        "Anticodon": "Whether this pocket carries the curated Anticodon Binding Domain label.",
        "Catalytic confidence": "0-4 score: 0 whenever the pocket lacks the curated Catalytic "
            "Domain label; otherwise 1 for the label alone, +1 for any weak ligand evidence, "
            "ceiling 4 for a strong ligand found directly in an experimental structure of this "
            "protein (ceiling 3 if the only strong evidence is AlphaFill-transplanted).",
        "Ligand evidence PDB": "Ligand codes (deduplicated, alphabetically sorted, pipe-separated) "
            "found directly in an experimental PDB structure near this pocket. Empty if none.",
        "Ligand evidence Alphafill": "Ligand codes (deduplicated, alphabetically sorted, pipe-"
            "separated) transplanted onto this structure by AlphaFill from a homologous "
            "experimental structure. Empty if none.",
        "Selected pocket": "Whether this pocket is one of the hand-curated pockets figures 3/4 "
            "showcase (output/selected_pockets.csv).",
        "Site type": "For a selected pocket, the curator's site-type call: CAT (catalytic) or "
            "NON-CAT (non-catalytic). Empty if not a selected pocket.",
        "Comment": "For a selected pocket, the curator's free-text rationale. Empty if not a "
            "selected pocket.",
    }
    for round_name in DOCKING_ROUNDS:
        for stat_name in STATS:
            explanations[f"{round_name} {stat_name}"] = (
                f"{STAT_EXPLANATIONS[stat_name].capitalize()} for this pocket, {ROUND_EXPLANATIONS[round_name]}. "
                "NaN (n=0) if this pocket has no docking data for this round."
            )
    return explanations


def main():
    annotation = pd.read_csv(ANNOTATION_CSV)
    assert len(annotation) == 276, f"Expected 276 pockets in {ANNOTATION_CSV}, got {len(annotation)}."

    pocket_detection = pd.read_csv(POCKET_DETECTION_CSV)
    annotation = annotation.merge(
        pocket_detection[["Uniprot AC", "File name", "Pocket number", "Pocket residues (chain_resn)"]],
        on=["Uniprot AC", "File name", "Pocket number"], how="left", validate="one_to_one")

    selected = pd.read_csv(SELECTED_POCKETS_CSV).set_index("pocket_name")

    rows = []
    for _, row in annotation.iterrows():
        structure_name = row["File name"].replace(".pdb", "")
        pocket_name = f"{structure_name}_pocket_{int(row['Pocket number'])}"

        sel = selected.loc[pocket_name] if pocket_name in selected.index else None

        record = {
            "Pocket name": pocket_name,
            "Gene name": row["Gene"],
            "Uniprot AC": row["Uniprot AC"],
            "Structure name": structure_name,
            "Pocket number": int(row["Pocket number"]),
            "Number of residues": len(str(row["Pocket residues (chain_resn)"]).split()),
            "P2Rank probability": row["Pocket probability"],
            "P2Rank score": row["Pocket score"],
            "Number of pockets for the same canonical pocket": int(row["n_models_in_cluster"]),
            **domain_flags(row["curated_labels"]),
            "Catalytic confidence": int(row["catalytic_confidence"]),
            "Ligand evidence PDB": clean_ligand_list(row["direct_pdb_ligands"]),
            "Ligand evidence Alphafill": clean_ligand_list(row["alphafill_ligands"]),
            "Selected pocket": sel is not None,
            "Site type": sel["site_type"] if sel is not None else "",
            "Comment": sel["comment"] if sel is not None else "",
        }

        hll_scores = safe_load_scores(LIBRARIES["DL"], pocket_name)
        real9_56m_scores = load_real_positive_scores(pocket_name)
        real9_92b_scores = safe_load_scores(LIBRARIES["REAL"], pocket_name)
        for round_name, scores in zip(DOCKING_ROUNDS, [hll_scores, real9_56m_scores, real9_92b_scores]):
            for stat_name, value in docking_stats(scores).items():
                record[f"{round_name} {stat_name}"] = value

        rows.append(record)

    table = pd.DataFrame(rows)
    table = table.sort_values(["Gene name", "Structure name", "Pocket number"]).reset_index(drop=True)

    output_path = os.path.join(table_dir, "ExtendedDataTable2.xlsx")
    save_table_with_legend(table, build_legend(), output_path)
    print(f"Saved {len(table)} row(s) to {output_path}")
    print(f"Selected pockets found: {table['Selected pocket'].sum()} of {len(selected)} in {SELECTED_POCKETS_CSV}")
    for round_name in DOCKING_ROUNDS:
        n_missing = (table[f"{round_name} n"] == 0).sum()
        print(f"  {round_name}: {n_missing} of {len(table)} pocket(s) with no docking data")


if __name__ == "__main__":
    main()
