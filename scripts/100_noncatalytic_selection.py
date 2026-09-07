#!/usr/bin/env python3
"""
Selects, out of the 1,095 filtered hits (script 70), compounds of interest against the 8
individual curated non-catalytic (NON-CAT) Mtb pockets -- unlike script 99's catalytic selection,
these are NOT aggregated per gene: pheS and aspS each have 2 raw NON-CAT pockets (labeled _1/_2,
same convention as script 98's curated_pocket_labels()), alaS has 2, pheT and lysS have 1 each --
8 pockets total, read directly from script 98's own condensed summary
(output/98_compound_docking_summary/compound_docking_summary_condensed.csv), which already carries
these un-aggregated raw scores.

Single hits: for each of the 8 pockets, independently by each method (docking, boltz2), the top `TOP_N`
compounds by that pocket's own raw score (ascending -- more negative docking score, or lower
predicted IC50 nM, is always better).

Dual hits: for each INTER-PROTEIN pocket pair -- i.e. excluding pairs where both pockets belong to
the same protein, where pheS+pheT count as one protein (same "pheST" merge as script 99, since
they're an obligate heterodimer) -- the two pockets' scores are averaged per compound and the top `TOP_N`
compounds by that average are kept. Of the 8 pockets grouped into 4 proteins (pheST: 3 pockets --
2 pheS + 1 pheT, aspS: 2, alaS: 2, lysS: 1), C(8,2) = 28 total pairs minus 5 intra-protein pairs
(pheST: C(3,2)=3, aspS: C(2,2)=1, alaS: C(2,2)=1, lysS: 0) leaves 23 inter-protein pairs.

Also reports each selected compound's Mtb-vs-human off-target selectivity, same three comparisons
and same mtb/on-target-minus-human sign convention as script 99 (negative = Mtb-favoring):
`diff_top1/5/10` (docking counter-screen), `diff_nesso1_top1/5/10` (Nesso-1 counter-screen), and
`diff_curated` (best score across all 12 curated on-target pockets, CAT and NON-CAT alike, minus
human_top1) -- all read from the condensed summary's own already-present columns.

Usage:
    python 100_noncatalytic_selection.py
"""
import itertools
import os

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

CONDENSED_CSV = os.path.join(ROOT, "output", "98_compound_docking_summary", "compound_docking_summary_condensed.csv")

OUTPUT_DIR = os.path.join(ROOT, "output", "100_noncatalytic_selection")
os.makedirs(OUTPUT_DIR, exist_ok=True)

NONCAT_POCKETS = ["pheS_NONCAT_1", "pheS_NONCAT_2", "pheT_NONCAT", "aspS_NONCAT_1",
                  "aspS_NONCAT_2", "lysS_NONCAT", "alaS_NONCAT_1", "alaS_NONCAT_2"]
CAT_POCKETS = ["pheS_CAT", "aspS_CAT", "lysS_CAT", "alaS_CAT"]
POCKET_PROTEIN = {  # pheS+pheT count as one protein (obligate heterodimer, same merge as script 99)
    "pheS_NONCAT_1": "pheST", "pheS_NONCAT_2": "pheST", "pheT_NONCAT": "pheST",
    "aspS_NONCAT_1": "aspS", "aspS_NONCAT_2": "aspS",
    "lysS_NONCAT": "lysS",
    "alaS_NONCAT_1": "alaS", "alaS_NONCAT_2": "alaS",
}
INTER_PROTEIN_PAIRS = [(p1, p2) for p1, p2 in itertools.combinations(NONCAT_POCKETS, 2)
                        if POCKET_PROTEIN[p1] != POCKET_PROTEIN[p2]]
TOP_N = 10
METHODS = ["docking", "boltz2"]
SELECTIVITY_NS = [1, 5, 10]


def hits_for_pockets(df, pockets, method):
    """Top-N rows (by column average across `pockets`, ascending -- lower is always better) for
    one (pockets, method) combination -- shared by single hits (pockets = one pocket, average of
    one column is just that column) and dual hits (pockets = an inter-protein pocket pair).
    skipna=False: a dual pair is only ranked if BOTH pockets have a real (non-gated) score --
    otherwise pandas' default skipna=True would silently average over just one pocket whenever the
    Boltz-2 low-confidence gating (script 98) NaN's the other, misreporting a single-pocket result
    as a dual hit."""
    cols = [f"{method}_{pocket}" for pocket in pockets]
    avg = df[cols].mean(axis=1, skipna=False)
    top = df.loc[avg.nsmallest(TOP_N).index]
    rows = []
    for rank, (_, row) in enumerate(top.iterrows(), start=1):
        rows.append({
            "compound_id": row["compound_id"],
            "smiles": row["smiles"],
            "hit_type": "single" if len(pockets) == 1 else "dual",
            "method": method,
            "rank": rank,
            "targets": "|".join(pockets),
            "scores": "|".join(f"{pocket}:{row[col]}" for pocket, col in zip(pockets, cols)),
        })
    return rows


def noncatalytic_hits(df):
    """One row per (pockets, method, rank) -- see module docstring for the single/dual schema. No
    dedup: each (pockets, method) selection is independent evidence, so a compound can rank top-5
    in more than one pocket/pair/method."""
    rows = []
    pocket_groups = [(p,) for p in NONCAT_POCKETS] + INTER_PROTEIN_PAIRS
    for pockets in pocket_groups:
        for method in METHODS:
            rows.extend(hits_for_pockets(df, pockets, method))
    return pd.DataFrame(rows)


def off_target_selectivity(df):
    """{compound_id: diff_top<n>, diff_nesso1_top<n>, diff_curated} -- same three comparisons as
    script 99's off_target_selectivity(), read from the condensed summary's own columns: mtb_top<n>
    minus human_top<n> (docking), nesso1_mtb_top<n> minus nesso1_human_top<n> (Nesso-1), and
    diff_curated = min across all 12 curated on-target pockets (CAT_POCKETS + NONCAT_POCKETS's raw
    docking_<pocket> columns -- min-of-mins is still the overall min) minus human_top1. Negative
    always means the Mtb/on-target side is the more favorable (lower) value."""
    out = df[["compound_id"]].copy()
    for n in SELECTIVITY_NS:
        out[f"diff_top{n}"] = df[f"mtb_top{n}"] - df[f"human_top{n}"]
        out[f"diff_nesso1_top{n}"] = df[f"nesso1_mtb_top{n}"] - df[f"nesso1_human_top{n}"]
    curated_cols = [f"docking_{p}" for p in CAT_POCKETS + NONCAT_POCKETS]
    out["diff_curated"] = df[curated_cols].min(axis=1) - df["human_top1"]
    diff_cols = [c for c in out.columns if c.startswith("diff_")]
    out[diff_cols] = out[diff_cols].round(3)
    return out


def main():
    df = pd.read_csv(CONDENSED_CSV)

    out = noncatalytic_hits(df)
    out = out.merge(off_target_selectivity(df), on="compound_id", how="left", validate="many_to_one")

    out_path = os.path.join(OUTPUT_DIR, "noncatalytic_hits.csv")
    out.to_csv(out_path, index=False)
    n_by_type = out["hit_type"].value_counts()
    print(f"Saved {len(out):,} rows ({', '.join(f'{n} {t}' for t, n in n_by_type.items())}) -> {out_path}")
    print(f"  single: {len(NONCAT_POCKETS)} pockets x {len(METHODS)} methods x top {TOP_N}")
    print(f"  dual: {len(INTER_PROTEIN_PAIRS)} inter-protein pairs x {len(METHODS)} methods x top {TOP_N}")
    print(f"Selected compounds (rows): {len(out):,}")
    print(f"Selected unique compounds: {out['compound_id'].nunique():,}")


if __name__ == "__main__":
    main()
