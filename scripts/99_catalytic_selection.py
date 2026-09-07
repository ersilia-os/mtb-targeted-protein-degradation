#!/usr/bin/env python3
"""
Selects, out of the 1,095 filtered hits (script 70), compounds of interest against the 4 curated
catalytic (CAT) Mtb pockets -- alaS, aspS, lysS, pheST (pheS+pheT merged, same convention as
script 98) -- each of which has exactly one CAT pocket (output/selected_pockets.csv), so script
98's docking_<gene>_CAT / boltz2_<gene>_CAT columns are already single-pocket values, not an
aggregate across several.

Single hits: for each of the 4 genes, independently by each method:
* docking -- top `TOP_N` (10) compounds by docking_<gene>_CAT (ascending, most negative = best)
* boltz2  -- top `TOP_N` (10) compounds by boltz2_<gene>_CAT (ascending, lowest predicted IC50 nM = best)
A compound can rank top-N by both methods for the same gene -- kept as two separate rows (no
dedup), since the two methods are independent evidence and the overlap itself is a signal.

Dual hits: for each of the C(4,2) = 6 unordered gene pairs, independently by each method, the two
genes' scores are averaged per compound (mean of docking_<geneA>_CAT and docking_<geneB>_CAT, or
of the two boltz2_ columns) and the top `TOP_N` (10) compounds by that average are kept -- a
joint-potency read, not a top-N-per-gene-then-intersect.

Multi hits: same averaging idea, but over all 4 genes at once (one group, not a combinatorial
choice like dual) -- per method, the top `TOP_N_MULTI = 50` compounds (5x single/dual's top 10,
since there's only one 4-gene group to select from rather than 4 genes or 6 pairs) by the mean of
all 4 docking_<gene>_CAT (or boltz2_<gene>_CAT) columns.

All three sections share one `hit_type` column ("single"/"dual"/"multi") and one `targets` /
`scores` schema: `targets` is the "|"-joined gene name(s) (one for single, two for dual, all four
for multi), and `scores` is "|"-joined "<gene>:<value>" pairs holding each gene's own raw score
(NOT the average dual/multi hits are ranked by -- that average isn't itself stored anywhere, only
used to pick the top N).

Also reports each selected compound's Mtb-vs-human off-target selectivity, for both counter-screen
methods (script 98's own off-target columns -- NOT the 4 curated on-target CAT pockets this script
otherwise deals with):

* `diff_top1/5/10` = `mtb_top<n>` minus `human_top<n>` (Uni-Dock docking, scripts 90-97).
* `diff_nesso1_top1/5/10` = `nesso1_mtb_top<n>` minus `nesso1_human_top<n>` (Nesso-1 co-folding,
  scripts 78-88, predicted IC50 in nM).

Both are mtb-minus-human (not human-minus-mtb) so that a NEGATIVE value means mtb_top<n> is the
more favorable (lower/more negative docking score, or lower predicted IC50) of the two -- i.e. the
compound's off-target profile genuinely favors engaging Mtb targets over human ones at that N. A
positive value means the opposite: the human off-target score is the more favorable one, a
selectivity liability.

One more, on-target-vs-off-target this time (single number, not split by N or pocket):
`diff_curated` = the compound's best (min) docking score across all 12 curated on-target pockets
minus `human_top1` -- min-of-the-8-already-aggregated docking_<gene>_CAT/_NONCAT columns equals
the true min across all 12 raw pockets, since a min-of-mins is still the overall min. Same sign
convention: negative means the curated on-target score is the more favorable one (good -- potency
beats the worst-case human off-target liability); positive means the compound's single strongest
human off-target hit outscores even its best on-target pocket.

Usage:
    python 99_catalytic_selection.py
"""
import itertools
import os

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

SUMMARY_CSV = os.path.join(ROOT, "output", "98_compound_docking_summary", "compound_docking_summary.csv")

OUTPUT_DIR = os.path.join(ROOT, "output", "99_catalytic_selection")
os.makedirs(OUTPUT_DIR, exist_ok=True)

GENES = ["alaS", "aspS", "lysS", "pheST"]
PAIRS = list(itertools.combinations(GENES, 2))
ALL_TARGET_GROUPS = [(g,) for g in GENES] + PAIRS + [tuple(GENES)]
HIT_TYPE_BY_SIZE = {1: "single", 2: "dual", 4: "multi"}
TOP_N = 10
TOP_N_MULTI = 50
METHODS = {"docking": "docking_{gene}_CAT", "boltz2": "boltz2_{gene}_CAT"}
SELECTIVITY_NS = [1, 5, 10]


CURATED_GROUP_COLS = [f"docking_{gene}_{site}" for gene in GENES for site in ("CAT", "NONCAT")]


def off_target_selectivity(df):
    """{compound_id: diff_top<n>, diff_nesso1_top<n>, diff_curated} -- mtb_top<n> minus human_top<n>
    (docking) and nesso1_mtb_top<n> minus nesso1_human_top<n> (Nesso-1), both from script 98's own
    off-target counter-screen columns, for n in SELECTIVITY_NS; plus diff_curated = the compound's
    best (min) score across all 12 curated on-target pockets (CURATED_GROUP_COLS's 8
    already-aggregated columns -- min-of-mins is still the overall min) minus human_top1. Negative
    always means the Mtb/on-target side is the more favorable (lower) value."""
    out = df[["compound_id"]].copy()
    for n in SELECTIVITY_NS:
        out[f"diff_top{n}"] = df[f"mtb_top{n}"] - df[f"human_top{n}"]
        out[f"diff_nesso1_top{n}"] = df[f"nesso1_mtb_top{n}"] - df[f"nesso1_human_top{n}"]
    out["diff_curated"] = df[CURATED_GROUP_COLS].min(axis=1) - df["human_top1"]
    diff_cols = [c for c in out.columns if c.startswith("diff_")]
    out[diff_cols] = out[diff_cols].round(3)
    return out


def hits_for_targets(df, targets, method, col_template, top_n):
    """Top-N rows (by column average across `targets`, ascending -- lower is always better) for
    one (targets, method) combination -- shared by single hits (targets = one gene, average of one
    column is just that column), dual hits (targets = a gene pair) and multi hits (targets = all
    4 genes). skipna=False: a dual/multi group is only ranked if EVERY target in it has a real
    (non-gated) score -- otherwise pandas' default skipna=True would silently average over fewer
    targets than claimed whenever the Boltz-2 low-confidence gating (script 98) NaN's one of them,
    misreporting a single-target result as a dual/multi hit."""
    cols = [col_template.format(gene=gene) for gene in targets]
    avg = df[cols].mean(axis=1, skipna=False)
    top = df.loc[avg.nsmallest(top_n).index]
    rows = []
    for rank, (_, row) in enumerate(top.iterrows(), start=1):
        rows.append({
            "compound_id": row["compound_id"],
            "smiles": row["smiles"],
            "hit_type": HIT_TYPE_BY_SIZE[len(targets)],
            "method": method,
            "rank": rank,
            "targets": "|".join(targets),
            "scores": "|".join(f"{gene}:{row[col]}" for gene, col in zip(targets, cols)),
        })
    return rows


def catalytic_hits(df):
    """One row per (targets, method, rank) -- see module docstring for the single/dual/multi
    schema. A compound can rank top-N by both methods for the same targets, or independently for
    more than one gene/pair -- kept as separate rows throughout, no dedup: each (targets, method)
    selection is independent evidence."""
    rows = []
    for targets in ALL_TARGET_GROUPS:
        top_n = TOP_N_MULTI if len(targets) == len(GENES) else TOP_N
        for method, col_template in METHODS.items():
            rows.extend(hits_for_targets(df, targets, method, col_template, top_n))
    return pd.DataFrame(rows)


def main():
    df = pd.read_csv(SUMMARY_CSV)

    out = catalytic_hits(df)
    out = out.merge(off_target_selectivity(df), on="compound_id", how="left", validate="many_to_one")

    out_path = os.path.join(OUTPUT_DIR, "catalytic_hits.csv")
    out.to_csv(out_path, index=False)
    n_by_type = out["hit_type"].value_counts()
    print(f"Saved {len(out):,} rows ({', '.join(f'{n} {t}' for t, n in n_by_type.items())}) -> {out_path}")

    n_both = (out.groupby(["hit_type", "targets", "compound_id"])["method"].nunique() == len(METHODS)).sum()
    print(f"Compounds ranking top-N by both methods for the same targets: {n_both}")

    print(f"Selected compounds (rows): {len(out):,}")
    print(f"Selected unique compounds: {out['compound_id'].nunique():,}")


if __name__ == "__main__":
    main()
