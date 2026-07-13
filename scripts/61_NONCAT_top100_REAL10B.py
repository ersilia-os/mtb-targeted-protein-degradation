#!/usr/bin/env python3
"""
Per NON-CAT pocket: sort script 60's docking report by score, greedily select the top TARGET_N
compounds with no synthon reused (MAX_SYN=1). Output: top100_diverse_per_pocket.csv (compound_id,
smiles, score, pocket_name).

Also produces a 6-group score boxplot (script 56's groups + this script's own "REAL 2 - selected"),
plus a physchem profile + t-SNE vs a background sampled evenly per pocket from script 57's raw
selective-hit pools (not the capped output/58 pool).

Usage:
    python 61_NONCAT_top100_REAL10B.py
"""
import glob
import os
import sys

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from default import RANDOM_SEED
from docking_utils import (
    LIBRARIES, compute_properties, load_real_negative_scores, load_real_positive_scores, load_scores,
    plot_profiling, plot_score_boxplots_multi, plot_tsne,
)

SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")
SCRIPT57_OUTPUT_DIR = os.path.join(ROOT, "output", "57_NONCAT_selective_10B")
DOCKING_RESULTS_DIR = os.path.join(ROOT, "output", "60_unidock_docking_noncat_selective_10B", "docking_results")
MERGED_SELECTIVE_HITS_CSV = os.path.join(ROOT, "output", "58_generate_conformations_noncat_selective_10B", "merged_selective_hits.csv")
NONCAT_TOP100_REAL10M_CSV = os.path.join(ROOT, "output", "56_NONCAT_top100_selection", "top100_per_noncat_pocket.csv")
OUTPUT_DIR = os.path.join(ROOT, "output", "61_docking_top100_diverse_selection")
os.makedirs(OUTPUT_DIR, exist_ok=True)

DIMER_POCKET = "7K98_pocket_6"
TARGET_N = 100
MAX_SYN = 1
BG_SAMPLE_PER_POCKET = 2_000


def load_noncat_targets():
    """[(gene_name, pocket_name), ...] for the 7 NON-CAT pockets, excluding the dimer pocket."""
    df = pd.read_csv(SELECTED_POCKETS_CSV)
    noncat = df[df["site_type"] == "NON-CAT"]
    noncat = noncat[noncat["pocket_name"] != DIMER_POCKET]
    return list(zip(noncat["gene_name"], noncat["pocket_name"]))


def get_synthons(compound_id):
    return compound_id.split("____")[1:]


def select_diverse_top_n(df, max_syn=MAX_SYN, target_n=TARGET_N):
    """df must be sorted ascending by score. Greedily picks target_n compounds, skipping any that
    reuse a synthon already selected."""
    used_synthons = set()
    selected_rows = []
    for _, row in df.iterrows():
        synthons = get_synthons(row["compound_id"])
        if any(s in used_synthons for s in synthons):
            continue
        selected_rows.append(row)
        used_synthons.update(synthons)
        if len(selected_rows) == target_n:
            break
    return pd.DataFrame(selected_rows)


def build_score_sources(pockets, result):
    """score_sources for docking_utils.plot_score_boxplots_multi: script 56's 5 groups plus this
    script's own new "REAL 2 - selected"."""
    dl_scores = {p: load_scores(os.path.join(LIBRARIES["DL"], p, "report.csv")) for p in pockets}
    real1_neg_scores = {p: load_real_negative_scores(p) for p in pockets}
    real1_pos_scores = {p: load_real_positive_scores(p) for p in pockets}
    real2_all_scores = {p: load_scores(os.path.join(LIBRARIES["REAL"], p, "report.csv")) for p in pockets}

    round1_selected = pd.read_csv(NONCAT_TOP100_REAL10M_CSV)
    real1_selected_scores = {
        p: round1_selected.loc[round1_selected["pocket_name"] == p].set_index("compound")["score"]
        for p in pockets
    }
    real2_selected_scores = {
        p: result.loc[result["pocket_name"] == p].set_index("compound_id")["score"]
        for p in pockets
    }

    return [
        ("Hit Locator", "mint", dl_scores, True),
        ("REAL 1 - negatives", "blue", real1_neg_scores, True),
        ("REAL 1 - positives", "orange", real1_pos_scores, False),
        ("REAL 2 - all", "yellow", real2_all_scores, True),
        ("REAL 1 - selected", "purple", real1_selected_scores, False),
        ("REAL 2 - selected", "pink", real2_selected_scores, False),
    ]


def sample_background(targets, exclude_ids, n_per_pocket, seed):
    """{compound_id: smiles}, n_per_pocket sampled independently from each pocket's own script-57
    selective-hit pool (equal representation, since pool sizes vary 36,878-559,663), excluding
    exclude_ids. Deduplicates within each pocket first (see script 58's gather_pocket_hits)."""
    sampled_dfs = []
    for gene, pocket in targets:
        chunk_files = sorted(glob.glob(os.path.join(SCRIPT57_OUTPUT_DIR, f"{gene}_{pocket}", "*.csv")))
        dfs = [pd.read_csv(f, usecols=["compound_id", "smiles"]) for f in chunk_files]
        if not dfs:
            continue
        pocket_pool = pd.concat(dfs, ignore_index=True).drop_duplicates(subset="compound_id")
        pocket_pool = pocket_pool[~pocket_pool["compound_id"].isin(exclude_ids)]
        sampled_dfs.append(pocket_pool.sample(n=n_per_pocket, random_state=seed))

    sampled = pd.concat(sampled_dfs, ignore_index=True)
    return dict(zip(sampled["compound_id"], sampled["smiles"]))


def main():
    targets = load_noncat_targets()
    print(f"NON-CAT targets (excl. dimer pocket {DIMER_POCKET}): {len(targets)}")

    smiles_lookup = pd.read_csv(MERGED_SELECTIVE_HITS_CSV, usecols=["compound_id", "smiles"])
    smiles_lookup = smiles_lookup.drop_duplicates(subset="compound_id").set_index("compound_id")["smiles"]

    selected_dfs = []
    for gene, pocket in targets:
        report_path = os.path.join(DOCKING_RESULTS_DIR, pocket, "report.csv")
        if not os.path.isfile(report_path):
            print(f"  {gene} ({pocket}): report.csv not found, skipping.")
            continue

        report = pd.read_csv(report_path)
        report = report.rename(columns={"compound": "compound_id"})
        report = report[report["score"].notna() & (report["score"] != "")]
        report["score"] = report["score"].astype(float)
        report = report.sort_values("score", ascending=True)

        selected = select_diverse_top_n(report)
        selected["pocket_name"] = pocket
        print(f"  {gene} ({pocket}): {len(report):,} candidates -> {len(selected)} selected"
              + (f" (WARNING: fewer than {TARGET_N})" if len(selected) < TARGET_N else ""))
        selected_dfs.append(selected)

    result = pd.concat(selected_dfs, ignore_index=True)
    result["smiles"] = result["compound_id"].map(smiles_lookup)
    result = result[["compound_id", "smiles", "score", "pocket_name"]]

    out_path = os.path.join(OUTPUT_DIR, "top100_diverse_per_pocket.csv")
    result.to_csv(out_path, index=False)
    print(f"\nSaved {len(result)} rows to {out_path}")

    print("\nBuilding score boxplots across 6 reference distributions...")
    pockets = [pocket for _, pocket in targets]
    score_sources = build_score_sources(pockets, result)
    plot_path = os.path.join(OUTPUT_DIR, "top100_diverse_score_boxplots.png")
    plot_score_boxplots_multi(pockets, score_sources, plot_path, xtick_rotation=20)
    print(f"Saved: {plot_path}")

    print("\nSampling background from script 57's raw selective-hit pools "
          f"({BG_SAMPLE_PER_POCKET:,}/pocket)...")
    sel_smiles = dict(zip(result["compound_id"], result["smiles"]))
    bg_smiles = sample_background(targets, set(sel_smiles), BG_SAMPLE_PER_POCKET, RANDOM_SEED)
    print(f"  Background: {len(bg_smiles):,} compounds")

    print("\nComputing physicochemical properties...")
    sel_props = compute_properties(sel_smiles)
    bg_props = compute_properties(bg_smiles)
    profiling_path = os.path.join(OUTPUT_DIR, "top100_diverse_physchem_profile.png")
    plot_profiling(sel_props, bg_props, profiling_path)
    print(f"Saved: {profiling_path}")

    print("\nBuilding t-SNE (ECFP4) projection...")
    tsne_path = os.path.join(OUTPUT_DIR, "top100_diverse_tsne.png")
    plot_tsne(sel_smiles, bg_smiles, tsne_path, seed=RANDOM_SEED)
    print(f"Saved: {tsne_path}")


if __name__ == "__main__":
    main()
