#!/usr/bin/env python3
"""
Multi-target hit-overlap analysis for the NON-CAT pockets against the second Enamine REAL
screening (output/unidock_REAL_docking_2, ~99,105 compounds) - the NON-CAT counterpart of
scripts/52_CAT_promiscuous.py.

Structurally more complex than CAT: output/selected_pockets.csv has 8 NON-CAT pocket rows, and
pheS and pheT (subunits of the same PheRS heterodimer) are combined into a single target unit
("pheST") so the analysis has 4 targets, matching CAT's target count:

    pheST: alphafold3_P9WFU3_model_2_pocket_2, 7K98_pocket_6, alphafold2_P9WFU1_model_0_pocket_1
           (2 pockets from pheS incl. the multimer, 1 from pheT)
    aspS:  alphafold3_P9WFW3_model_3_pocket_2, chai1_P9WFW3_model_0_pocket_2
    lysS:  alphafold3_P9WFU9_model_1_pocket_2
    alaS:  swissmodel_P9WFW7_model_0_pocket_2, alphafold2_P9WFW7_model_0_pocket_3

output/reference_pocket_noncatalytic.csv is stale (predates 7K98_pocket_6), so the target->pockets
mapping is read from output/selected_pockets.csv (site_type == "NON-CAT") instead.

A "top-N hit" for a target means the compound is a top-N hit in AT LEAST ONE of that target's
pockets (union across its pockets) - so a target's observed hit-set size S_target ranges from N
(full overlap between its pockets) up to N * n_pockets (no overlap). Expected-by-chance for a
specific combination of targets uses each target's ACTUAL observed S_target, not N:
expected(combo) = product(S_t for t in combo) / M**(len(combo) - 1)
reported per exact combination (not a k-level aggregate), since S_target varies by target.

For each of top-100 and top-1,000 (top-10,000 is not run), the final selected/saved compound set
is: >=2 of the 4 NON-CAT targets hit, AND NOT CAT-promiscuous - using one threshold level up from
the NON-CAT pass on script 52's CAT multi-target lists (NON-CAT top-100 excludes CAT-promiscuous-
at-top-1,000; NON-CAT top-1,000 excludes CAT-promiscuous-at-top-10,000, ~10,900 compounds). A
NON-CAT multi-target compound that merely binds ONE CAT pocket well is not concerning on its own;
only ALSO being CAT-promiscuous (>=2 CAT targets, suggesting a non-specific, dock-well-everywhere
artifact rather than genuine site selectivity) triggers exclusion.

For each of top-100 and top-1,000:
  - per-target hit-set sizes (S_target)
  - an UpSet plot of hit overlap across the 4 targets
  - observed vs. expected-by-chance for every combination of 2, 3, or 4 targets
  - a console count of the CAT-promiscuous exclusion actually applied (see criterion above), plus
    an overlap count vs. script 53's ~500 CAT-selective compounds, purely informational
  - a CSV of the final selected compounds (per the criterion above), with per-POCKET (not
    per-target) score/rank for all 8 pockets, target-level n_targets/targets_hit, SMILES and
    physchem properties

Then, for the >=2-target compounds from the top-1,000 CSV vs. a random 20,000-compound REAL
background: score boxplots (per pocket, vs. DL/REAL-round-1-negative-set references), physchem
profiling, and a t-SNE (ECFP4) plot.

Usage:
    python 54_NONCAT_promiscuous.py
"""
import os
import sys
from collections import Counter
from itertools import combinations

import numpy as np
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from default import RANDOM_SEED
from docking_utils import (
    LIBRARIES, compute_properties, load_real_negative_scores, load_scores, lookup_smiles,
    plot_profiling, plot_score_boxplots, plot_tsne, sample_prescreened_smiles,
)

SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")
OUTPUT_DIR = os.path.join(ROOT, "output", "54_NONCAT_promiscuous")
MULTIMER_DOCKING_DIR = os.path.join(ROOT, "output", "50_unidock_docking_multimers")
CAT_MULTI_TARGET_TOP1000_CSV = os.path.join(
    ROOT, "output", "52_CAT_promiscuous", "alaS_aspS_lysS_pheS_REAL_CAT_multi_target_top1000.csv"
)
CAT_MULTI_TARGET_TOP10000_CSV = os.path.join(
    ROOT, "output", "52_CAT_promiscuous", "alaS_aspS_lysS_pheS_REAL_CAT_multi_target_top10000.csv"
)
CAT_SELECTIVE_CSV = os.path.join(
    ROOT, "output", "53_CAT_selective", "alaS_aspS_lysS_pheS_REAL_CAT.csv"
)

LIB = "REAL"
TOP_NS = [100, 1_000]
# NON-CAT top_n -> the CAT-promiscuous reference threshold used to exclude likely artifacts,
# one level up (NON-CAT top-100 checked against CAT top-1,000; NON-CAT top-1,000 against CAT
# top-10,000). CAT_EXCLUSION_CSV maps each to the corresponding script-52 output path.
CAT_EXCLUSION_CSV = {100: CAT_MULTI_TARGET_TOP1000_CSV, 1_000: CAT_MULTI_TARGET_TOP10000_CSV}
PROFILING_TOP_N = 1_000
BG_SAMPLE_SIZE = 20_000

PROP_COLUMNS = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotBonds", "AromaticRings", "QED", "is_pains"]


def load_noncat_targets():
    """target_name -> [pocket_name, ...] for the NON-CAT pockets, with pheS+pheT merged into one
    "pheST" target unit (1-3 pockets per target)."""
    df = pd.read_csv(SELECTED_POCKETS_CSV)
    df = df[df["site_type"] == "NON-CAT"].copy()
    df["gene_name"] = df["gene_name"].replace({"pheS": "pheST", "pheT": "pheST"})
    return df.groupby("gene_name")["pocket_name"].apply(list).to_dict()


def resolve_pocket_results_dir(pocket_name):
    """Same "_model_" heuristic as scripts/51_selected_pockets_visualization.py: monomeric
    (AF2/AF3/SwissModel/Chai1) pockets live under the shared REAL round-2 results dir, the
    multimeric one (7K98_pocket_6) under its own script-50 output tree."""
    if "_model_" in pocket_name:
        return LIBRARIES[LIB]
    return MULTIMER_DOCKING_DIR


def load_pocket_scores(pocket_names):
    """{pocket_name: pd.Series} - per-pocket equivalent of docking_utils.build_matrix, which
    assumes one shared results_dir for every pocket (not true here: 7K98_pocket_6 lives elsewhere)."""
    scores = {}
    for pocket_name in pocket_names:
        report = os.path.join(resolve_pocket_results_dir(pocket_name), pocket_name, "report.csv")
        if not os.path.isfile(report):
            print(f"  Warning: report not found for pocket '{pocket_name}', skipping.")
            continue
        scores[pocket_name] = load_scores(report)
    return scores


def target_hit_set(target_pockets, pocket_top, target, top_n):
    """Union of each of target's pockets' top-top_n compound sets."""
    return set().union(*(set(pocket_top[p][:top_n]) for p in target_pockets[target]))


def save_target_upset(target_hit_sets, top_n, output_dir, lib, trna_tag):
    """Same as docking_utils.save_upset, but takes pre-computed hit-sets directly instead of full
    sorted lists to slice - needed because a multi-pocket target's "top-N" set is a union across
    pockets, not a simple prefix of one ranked list."""
    import warnings
    import matplotlib.pyplot as plt
    from upsetplot import from_contents, plot as upset_plot
    warnings.filterwarnings("ignore", category=FutureWarning, module="upsetplot")

    data = from_contents(target_hit_sets)
    upset_plot(data)
    plt.suptitle(f"{lib} — top {top_n:,}")
    fname = f"{trna_tag}_{lib}_upset_top{top_n}.png"
    path = os.path.join(output_dir, fname)
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved: {path}")


def build_multi_target_table(targets, target_pockets, hit_sets, pocket_rank, pocket_scores, all_pockets):
    """
    Builds a DataFrame of every compound hitting >=2 of the 4 targets' hit-sets, with per-pocket
    score/rank (NaN if that specific pocket didn't have the compound in its own top-N),
    target-level n_targets/targets_hit, SMILES and physchem properties. Also returns
    {target: S_target} for the observed-vs-expected reporting.
    """
    n_targets_of, targets_hit_of = {}, {}
    all_in_top = set.union(*hit_sets.values())
    for cid in all_in_top:
        hit = [t for t in targets if cid in hit_sets[t]]
        n_targets_of[cid] = len(hit)
        targets_hit_of[cid] = ";".join(hit)

    multi_ids = [cid for cid in all_in_top if n_targets_of[cid] >= 2]
    smiles_map = lookup_smiles(multi_ids, LIB)

    rows = []
    for cid in multi_ids:
        row = {
            "compound": cid, "smiles": smiles_map.get(cid, ""),
            "n_targets": n_targets_of[cid], "targets_hit": targets_hit_of[cid],
        }
        for p in all_pockets:
            in_top = p in pocket_rank and cid in pocket_rank[p]
            row[f"score_{p}"] = pocket_scores[p].get(cid, np.nan) if in_top else np.nan
            row[f"rank_{p}"] = pocket_rank[p].get(cid, np.nan) if in_top else np.nan
        rows.append(row)
    multi_df = pd.DataFrame(rows)

    props = compute_properties(smiles_map)
    multi_df = multi_df.merge(props, left_on="compound", right_index=True, how="left")
    return multi_df


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    target_pockets = load_noncat_targets()
    targets = sorted(target_pockets)
    trna_tag = "_".join(targets)
    all_pockets = [p for pockets in target_pockets.values() for p in pockets]
    print(f"NON-CAT targets: {target_pockets}")

    pocket_scores = load_pocket_scores(all_pockets)
    M = max(len(s) for s in pocket_scores.values())

    pocket_top = {p: s.sort_values(ascending=True).index.tolist() for p, s in pocket_scores.items()}
    pocket_rank = {p: {cid: i + 1 for i, cid in enumerate(ids)} for p, ids in pocket_top.items()}

    def load_cat_ids(path):
        if os.path.isfile(path):
            return set(pd.read_csv(path)["compound"])
        print(f"  Warning: {path} not found, can't cross-check against it.")
        return set()

    cat_exclusion_ids = {top_n: load_cat_ids(path) for top_n, path in CAT_EXCLUSION_CSV.items()}
    cat_selective_ids = load_cat_ids(CAT_SELECTIVE_CSV)

    multi_target_dfs = {}
    for top_n in TOP_NS:
        print(f"\n--- Top-{top_n} ---")
        hit_sets = {t: target_hit_set(target_pockets, pocket_top, t, top_n) for t in targets}
        for t in targets:
            print(f"  {t}: S_target={len(hit_sets[t]):,} (from {len(target_pockets[t])} pocket(s))")

        save_target_upset(hit_sets, top_n, OUTPUT_DIR, LIB, f"{trna_tag}_NONCAT")

        print("  Observed vs. expected-by-chance, per combination:")
        for k in range(2, len(targets) + 1):
            for combo in combinations(targets, k):
                observed = len(set.intersection(*(hit_sets[t] for t in combo)))
                expected = 1.0
                for t in combo:
                    expected *= len(hit_sets[t])
                expected /= M ** (k - 1)
                print(f"    {'+'.join(combo)}: observed={observed:,}  expected={expected:.2f}")

        multi_df = build_multi_target_table(targets, target_pockets, hit_sets, pocket_rank, pocket_scores, all_pockets)
        cat_excl_top_n = {100: 1_000, 1_000: 10_000}[top_n]
        in_cat_exclusion = multi_df["compound"].isin(cat_exclusion_ids[top_n])
        in_cat_selective = multi_df["compound"].isin(cat_selective_ids)
        print(f"  Compounds by number of targets hit, after excluding likely artifacts\n"
              f"  (also CAT-promiscuous at top-{cat_excl_top_n:,}); "
              f"CAT-selective overlap is informational, not an exclusion:")
        for k in range(2, len(targets) + 1):
            at_k = multi_df["n_targets"] == k
            print(f"    {k}/{len(targets)} targets: {at_k.sum():,} total"
                  f" | excl. CAT-promiscuous top-{cat_excl_top_n:,}: {(at_k & ~in_cat_exclusion).sum():,}"
                  f" | overlap w/ 500 CAT-selective: {(at_k & in_cat_selective).sum():,}")

        # Actual selection: >=2 targets (already guaranteed by build_multi_target_table) AND not
        # CAT-promiscuous at the paired threshold (one level up from this NON-CAT top_n).
        multi_df = multi_df[~in_cat_exclusion].copy()

        out_csv = os.path.join(OUTPUT_DIR, f"{trna_tag}_{LIB}_NONCAT_multi_target_top{top_n}.csv")
        multi_df = multi_df.sort_values("n_targets", ascending=False)
        multi_df.to_csv(out_csv, index=False)
        print(f"  Saved: {out_csv} ({len(multi_df):,} selected)")
        multi_target_dfs[top_n] = multi_df

    sel_df = multi_target_dfs[PROFILING_TOP_N]
    sel_ids = set(sel_df["compound"])

    print("\n--- Score boxplots ---")
    scores1 = pd.concat(pocket_scores, axis=1)
    dl_scores = pd.concat(
        {p: load_scores(os.path.join(ROOT, "output", "unidock_docking", "docking_results", p, "report.csv"))
         for p in all_pockets if os.path.isfile(os.path.join(ROOT, "output", "unidock_docking", "docking_results", p, "report.csv"))},
        axis=1,
    )
    real_negative_scores = {p: load_real_negative_scores(p) for p in all_pockets}
    scores_path = os.path.join(OUTPUT_DIR, f"{trna_tag}_{LIB}_NONCAT_scores.png")
    plot_score_boxplots(scores1, sel_ids, all_pockets, scores_path, dl_scores=dl_scores, real_negative_scores=real_negative_scores, xtick_rotation=20)
    print(f"Saved: {scores_path}")

    print("\n--- Physchem profiling ---")
    sel_props = sel_df.set_index("compound")[PROP_COLUMNS]
    sel_smiles = dict(zip(sel_df["compound"], sel_df["smiles"]))

    bg_smiles = sample_prescreened_smiles(LIB, sel_ids, BG_SAMPLE_SIZE, RANDOM_SEED)
    bg_props = compute_properties(bg_smiles)

    profiling_path = os.path.join(OUTPUT_DIR, f"{trna_tag}_{LIB}_NONCAT_profiling.png")
    plot_profiling(sel_props, bg_props, profiling_path)
    print(f"Saved: {profiling_path}")

    print("\n--- t-SNE (ECFP4) ---")
    tsne_path = os.path.join(OUTPUT_DIR, f"{trna_tag}_{LIB}_NONCAT_tsne.png")
    plot_tsne(sel_smiles, bg_smiles, tsne_path, seed=RANDOM_SEED)
    print(f"Saved: {tsne_path}")


if __name__ == "__main__":
    main()
