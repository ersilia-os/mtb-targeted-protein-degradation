#!/usr/bin/env python3
"""
Merges script 99's catalytic and script 100's non-catalytic compound selections, and prepares the
input the molecule-auditing skill (~/github/claude-ersilia-skills/skills/molecule-auditing) needs
to build an interactive explorer over the result.

Writes two files:

1. `merged_selections.csv` -- scripts 99 and 100's output concatenated as-is (both already share
   the exact same schema: compound_id, smiles, hit_type, method, rank, targets, scores, diff_top1/
   5/10, diff_nesso1_top1/5/10, diff_curated), plus one new `origin` column ("catalytic" /
   "non-catalytic") so provenance survives the merge. One row per selection event -- a compound
   selected under multiple (targets, method) combinations still appears multiple times here.

2. `audit_input.csv` -- one row per UNIQUE compound (the actual audit-skill input), joining:
   * physchem/ADMET/liability columns from script 70's filtered_hits.csv (MW, cLogP, QED, PAINS/
     Brenk flags, cytotoxicity, mycomembrane permeation, etc.) -- notably, NONE of these ever made
     it into script 70's own audit explorer (confirmed by inspecting its config.json), even though
     they were sitting right there in the same input file.
   * on-target scores from script 98's compound_docking_summary.csv: the 8 curated
     docking_<gene>_CAT/NONCAT and boltz2_<gene>_CAT/NONCAT columns, plus the diff_* off-target-
     selectivity columns (already present per-row in merged_selections.csv too, but a single value
     per compound here).
   * selection provenance, summarized from merged_selections.csv: `origins` (catalytic/non-
     catalytic/both), `n_selections` (row count for this compound across both files), `best_rank`,
     `hit_types`/`selected_targets` (sorted, semicolon-joined unique values), and `starred_columns`
     -- the set of on-target audit columns (docking_<gene>_<CAT|NONCAT> / boltz2_<gene>_<CAT|NONCAT>)
     that actually qualified this compound for inclusion, across every selection event. Script 99's
     `targets` are bare gene names (always CAT); script 100's are individual pocket labels (e.g.
     "alaS_NONCAT_1", "pheS_NONCAT_2") that map back to the single gene-level NONCAT column the
     audit explorer displays (pheS/pheT pockets both map to "pheST_NONCAT", the same merge used
     throughout scripts 98-100). The explorer stars these cells so it's visible at a glance which
     specific score(s) put a compound on the list, not just that it's on the list.

No column curation, renaming, badges, or filters are applied here -- that interview-driven step is
the molecule-auditing skill's own Step 1-2 job, run interactively afterward.

Usage:
    python 101_merge_selections.py
"""
import os
import re
from collections import Counter

import numpy as np
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

CATALYTIC_CSV = os.path.join(ROOT, "output", "99_catalytic_selection", "catalytic_hits.csv")
NONCATALYTIC_CSV = os.path.join(ROOT, "output", "100_noncatalytic_selection", "noncatalytic_hits.csv")
FILTERED_HITS_CSV = os.path.join(ROOT, "output", "70_filtering", "filtered_hits.csv")
SUMMARY_CSV = os.path.join(ROOT, "output", "98_compound_docking_summary", "compound_docking_summary.csv")
HUMAN_GENE_MIN_CSV = os.path.join(ROOT, "output", "97_human_merge_docking_scores", "gene_min_scores.csv")

OUTPUT_DIR = os.path.join(ROOT, "output", "101_merge_selections")
os.makedirs(OUTPUT_DIR, exist_ok=True)

ADMET_COLS = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotBonds", "AromaticRings", "QED", "is_pains",
              "sps_score", "nsps_score", "mycomembrane_permeation", "has_pains", "has_brenk",
              "is_sim_known_ab", "nitrofuran_motif", "fluoroquinolone_motif", "carbepenem_motif",
              "betalactam_motif", "cytotoxicity_hepg2", "cytotoxicity_hskmc", "cytotoxicity_imr90"]

GENES = ["alaS", "aspS", "lysS", "pheST"]
ONTARGET_COLS = [f"{method}_{gene}_{site}" for method in ("docking", "boltz2")
                  for gene in GENES for site in ("CAT", "NONCAT")]
DIFF_COLS = [f"diff_top{n}" for n in (1, 5, 10)] + [f"diff_nesso1_top{n}" for n in (1, 5, 10)] + ["diff_curated"]
SELECTIVITY_NS = [1, 5, 10]
CURATED_GROUP_COLS = [f"docking_{gene}_{site}" for gene in GENES for site in ("CAT", "NONCAT")]

# Manually added near-misses (human-reviewed, 2026-09-07): each ranked 11th or 12th -- just outside
# script 99/100's top-10 cutoff -- on one on-target column, at a razor-thin margin from the actual
# cutoff value (see the conversation this came from: 3 of the 4 are within ~2% of the 10th-place
# score; the 4th, boltz2_lysS_NONCAT, is the weakest case at ~9% but still the same single-rank-11
# story). `targets` uses each script's own token convention (bare gene for catalytic, pocket label
# for non-catalytic) so target_to_gene_site() resolves them identically to a real selection row.
MANUAL_ADDITIONS = [
    {"compound_id": "s_22____25571760____15235754", "origin": "catalytic", "method": "docking", "targets": "aspS", "rank": 11},
    {"compound_id": "s_282770____26484422____26578612", "origin": "catalytic", "method": "docking", "targets": "pheST", "rank": 11},
    {"compound_id": "s_2708____26075032____24002312", "origin": "catalytic", "method": "boltz2", "targets": "alaS", "rank": 12},
    {"compound_id": "m_279130____28136586____28146358____28610778", "origin": "non-catalytic", "method": "boltz2", "targets": "lysS_NONCAT", "rank": 11},
]


def manual_addition_rows(summary):
    """Builds rows for MANUAL_ADDITIONS in the exact schema scripts 99/100 produce, so they merge
    into `merged` and flow through starred_columns/provenance identically to a real selection --
    smiles/scores/diff_* all read from script 98's own summary table, same source everything else
    uses."""
    summary = summary.set_index("compound_id")
    rows = []
    for add in MANUAL_ADDITIONS:
        row = summary.loc[add["compound_id"]]
        gene_site = target_to_gene_site(add["targets"], add["origin"])
        col = f"{add['method']}_{gene_site}"
        diffs = {f"diff_top{n}": round(row[f"mtb_top{n}"] - row[f"human_top{n}"], 3) for n in SELECTIVITY_NS}
        diffs.update({f"diff_nesso1_top{n}": round(row[f"nesso1_mtb_top{n}"] - row[f"nesso1_human_top{n}"], 3) for n in SELECTIVITY_NS})
        diffs["diff_curated"] = round(row[CURATED_GROUP_COLS].min() - row["human_top1"], 3)
        rows.append({
            "compound_id": add["compound_id"],
            "smiles": row["smiles"],
            "hit_type": "single",
            "method": add["method"],
            "rank": add["rank"],
            "targets": add["targets"],
            "scores": f"{add['targets']}:{row[col]}",
            "origin": add["origin"],
            **diffs,
        })
    return pd.DataFrame(rows)


def merged_selections():
    """Concatenates scripts 99 and 100's output (tagging each with an `origin` column) plus
    MANUAL_ADDITIONS' hand-picked near-miss rows."""
    cat = pd.read_csv(CATALYTIC_CSV)
    cat.insert(0, "origin", "catalytic")
    noncat = pd.read_csv(NONCATALYTIC_CSV)
    noncat.insert(0, "origin", "non-catalytic")
    manual = manual_addition_rows(pd.read_csv(SUMMARY_CSV))
    return pd.concat([cat, noncat, manual], ignore_index=True)


def target_to_gene_site(target, origin):
    """Maps one script 99/100 `targets` token to the "<gene>_<CAT|NONCAT>" suffix used by the
    audit's docking_<gene>_<site> / boltz2_<gene>_<site> columns. Script 99 (origin="catalytic")
    tokens are bare gene names, always CAT. Script 100 (origin="non-catalytic") tokens are
    individual pocket labels ("<gene>_NONCAT[_n]") -- pheS/pheT pockets both collapse to
    "pheST_NONCAT" (same merge as curated_pocket_groups() throughout scripts 98-100)."""
    if origin == "catalytic":
        return f"{target}_CAT"
    if target.startswith("pheS_NONCAT") or target == "pheT_NONCAT":
        return "pheST_NONCAT"
    return re.match(r"^([A-Za-z0-9]+)_(CAT|NONCAT)", target).group(0)


def starred_columns_for_row(row):
    """[(<method>_<gene>_<site>, hit_type), ...] -- one (column, hit_type) pair per target token
    in this selection-event row (a list, not a set: a compound qualified by the same (column,
    hit_type) pair via N different selection events should count N times, so the explorer can show
    N stars / let the Type filter know how many Single/Dual/Multi events landed on that column)."""
    return [(f"{row['method']}_{target_to_gene_site(t, row['origin'])}", row["hit_type"])
            for t in row["targets"].split("|")]


def format_starred_columns(lists):
    """"<col>:<hit_type>:<count>|..." -- how many selection events (across all of this compound's
    rows) named each (column, hit_type) pair, sorted for determinism. hit_type travels alongside
    the column so the explorer's Type (Single/Dual/Multi) filter can be scoped to only the
    currently-selected proteins/method, not applied globally."""
    counts = Counter(pair for lst in lists for pair in lst)
    return "|".join(f"{col}:{ht}:{n}" for (col, ht), n in sorted(counts.items()))


HIT_TYPE_ORDER = ["single", "dual", "multi"]


def format_hit_types(values):
    """"Single, dual, multi" -- fixed order (not alphabetical), comma-separated, only the first
    word capitalized (sentence-style, not each word), including only the hit_types this compound
    actually has. When all three are present (the longest, most cramped combination), "Single" is
    abbreviated to "Sing." to save card space."""
    present = [t for t in HIT_TYPE_ORDER if t in set(values)]
    if not present:
        return ""
    first = "Sing." if present == HIT_TYPE_ORDER else present[0].capitalize()
    return ", ".join([first] + present[1:])


def hit_types_by_origin(merged, origin, index):
    """format_hit_types(), computed separately for just this origin's rows and reindexed to the
    full compound index -- real NaN (not the string "None") for a compound with no selection under
    this origin at all, so Catalytic/Non-catalytic can be shown as two independent rows rather than
    one merged Origin + Type pair. NaN, not "None": pandas' read_csv treats the literal string
    "None" as a missing value by default, so it would silently round-trip back to NaN on the next
    read (e.g. inside build_table.py) -- the "None" label is applied at display time instead (see
    the explorer's config.json `na_display`)."""
    sub = merged[merged["origin"] == origin].groupby("compound_id")["hit_type"].agg(format_hit_types)
    return sub.reindex(index)


def selection_provenance(merged):
    """One row per unique compound_id: origins (catalytic/non-catalytic/both, "|"-joined sorted
    unique values), n_selections (row count across both files), best_rank (min rank achieved by
    any of its selections), hit_types_catalytic/hit_types_noncatalytic (format_hit_types() per
    origin, "None" if this compound has no selection under that origin -- see
    hit_types_by_origin()) and selected_targets (sorted, "|"-joined unique values), plus the diff_*
    columns -- identical across every row for a given compound_id (computed purely from
    compound_id in scripts 99/100), so `.first()` is safe, not an arbitrary pick."""
    merged = merged.copy()
    merged["starred"] = merged.apply(starred_columns_for_row, axis=1)

    grouped = merged.groupby("compound_id")
    compound_index = grouped.size().index
    out = pd.DataFrame({
        "origins": grouped["origin"].agg(lambda s: "|".join(sorted(s.unique()))),
        "n_selections": grouped.size(),
        "best_rank": grouped["rank"].min(),
        "hit_types_catalytic": hit_types_by_origin(merged, "catalytic", compound_index),
        "hit_types_noncatalytic": hit_types_by_origin(merged, "non-catalytic", compound_index),
        "selected_targets": grouped["targets"].agg(lambda s: "|".join(sorted(set("|".join(s).split("|"))))),
        "starred_columns": grouped["starred"].agg(format_starred_columns),
    })
    return out.join(grouped[DIFF_COLS].first())


def top2_human_offtargets(compound_index):
    """{compound_id: human_offtarget_1, human_offtarget_2} -- the compound's two most favorable
    (lowest) docking scores among the 38 human off-target genes, as plain numbers (equal to
    human_top1/human_top2 -- kept here, rather than just reusing those, so the exact same
    per-gene lookup is available if a gene-name label is wanted again later)."""
    wide = pd.read_csv(HUMAN_GENE_MIN_CSV).set_index("compound_id").reindex(compound_index)
    sorted_vals = np.sort(wide.to_numpy(dtype="float64", na_value=np.inf), axis=1)
    out = pd.DataFrame(index=wide.index)
    out["human_offtarget_1"] = sorted_vals[:, 0]
    out["human_offtarget_2"] = sorted_vals[:, 1]
    return out


def main():
    merged = merged_selections()
    merged_path = os.path.join(OUTPUT_DIR, "merged_selections.csv")
    merged.to_csv(merged_path, index=False)
    print(f"Saved {len(merged):,} rows -> {merged_path}")
    print(f"  catalytic: {(merged.origin == 'catalytic').sum():,} rows, "
          f"non-catalytic: {(merged.origin == 'non-catalytic').sum():,} rows")

    provenance = selection_provenance(merged)
    unique_ids = provenance.index
    print(f"\n{len(unique_ids):,} unique compounds "
          f"({(provenance['origins'] == 'catalytic').sum():,} catalytic-only, "
          f"{(provenance['origins'] == 'non-catalytic').sum():,} non-catalytic-only, "
          f"{(provenance['origins'] == 'catalytic|non-catalytic').sum():,} both)")

    admet = pd.read_csv(FILTERED_HITS_CSV)[["compound_id", "smiles"] + ADMET_COLS].set_index("compound_id")
    ontarget = pd.read_csv(SUMMARY_CSV)[["compound_id"] + ONTARGET_COLS].set_index("compound_id")
    human_offtargets = top2_human_offtargets(unique_ids)

    audit_input = provenance.join([admet, ontarget, human_offtargets], how="left")
    audit_input = audit_input.reset_index()

    audit_path = os.path.join(OUTPUT_DIR, "audit_input.csv")
    audit_input.to_csv(audit_path, index=False)
    print(f"\nSaved {len(audit_input):,} rows x {len(audit_input.columns)} columns -> {audit_path}")
    print(f"Missing cells (ADMET/on-target join): {audit_input[ADMET_COLS + ONTARGET_COLS].isna().sum().sum()}")


if __name__ == "__main__":
    main()
