#!/usr/bin/env python3
"""
Top-100 non-promiscuous docking hits per NON-CAT pocket, for the 7 NON-CAT pockets EXCLUDING the
dimer-interface pocket (pheS's 7K98_pocket_6 - output/selected_pockets.csv's comment for that row
literally reads "pocket at dimerization interface"). The pocket list (gene_name -> pocket_name) is
read directly from output/selected_pockets.csv (site_type == "NON-CAT"), not hardcoded, matching
script 54's approach to the same NON-CAT target set:

    pheS: alphafold3_P9WFU3_model_2_pocket_2
    pheT: alphafold2_P9WFU1_model_0_pocket_1
    aspS: alphafold3_P9WFW3_model_3_pocket_2, chai1_P9WFW3_model_0_pocket_2
    lysS: alphafold3_P9WFU9_model_1_pocket_2
    alaS: swissmodel_P9WFW7_model_0_pocket_2, alphafold2_P9WFW7_model_0_pocket_3

For each pocket: load its Enamine REAL round-1 docking report (output/unidock_REAL_docking/
docking_results/{pocket}/report.csv - each pocket's own top-100k prioritized compounds + a
background sample, ~112,959 rows, NOT the full 9.56M), sort ascending by score (Vina/Uni-Dock
convention: more negative = better, so ascending = best-to-worst), discard compounds flagged
promiscuous by script 55 (output/55_identify_promiscuous_enamine_REAL/promiscuous_hits.csv), and
keep the top 100 remaining. The per-pocket 100th-place (cutoff) score is printed as each pocket is
processed.

The 7 per-pocket top-100 tables are concatenated WITHOUT deduplication: if the same compound is a
top-100 hit in more than one pocket, it appears once per pocket (~700 rows total, human-confirmed
choice - a final count of such cross-pocket duplicates is printed, not an error).

Also produces a boxplot (per pocket) of raw docking scores across 5 reference distributions:
Hit Locator (Enamine DL 100k), REAL 1 - negatives (~13k, round-1 background), REAL 1 - positives
(~100k, round-1 surrogate-prioritized set), REAL 2 - all (~99,105, the pre-screened synthon-diverse
redocked set), and REAL 1 - selected (this script's own top-100/pocket). Uses the new
docking_utils.plot_score_boxplots_multi rather than plot_score_boxplots (used by scripts 52-54),
since "REAL 1 - selected" here is drawn from round 1 while "REAL 2 - all" is a separately
synthon-deduped round-2 selection - they are not nested, unlike scripts 52-54's selected-subset-of-
pre-screened assumption.

Usage:
    python 56_NONCAT_top100_selection.py
"""
import os
import sys

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from default import RANDOM_SEED
from docking_utils import (
    LIBRARIES, REAL_ROUND1_RESULTS_DIR, compute_properties, load_real_negative_scores,
    load_real_positive_scores, load_scores, lookup_smiles, plot_profiling,
    plot_score_boxplots_multi, plot_tsne, sample_prescreened_smiles,
)

SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")
PROMISCUOUS_HITS_CSV = os.path.join(ROOT, "output", "55_identify_promiscuous_enamine_REAL",
                                     "promiscuous_hits.csv")
OUTPUT_DIR = os.path.join(ROOT, "output", "56_NONCAT_top100_selection")
os.makedirs(OUTPUT_DIR, exist_ok=True)

DIMER_POCKET = "7K98_pocket_6"
TOP_N = 100
BG_SAMPLE_SIZE = 20_000


def main():
    pockets_df = pd.read_csv(SELECTED_POCKETS_CSV)
    noncat = pockets_df[pockets_df["site_type"] == "NON-CAT"]
    noncat = noncat[noncat["pocket_name"] != DIMER_POCKET]
    print(f"NON-CAT pockets (excl. dimer): {len(noncat)}")
    for gene, pocket in zip(noncat["gene_name"], noncat["pocket_name"]):
        print(f"  {gene}: {pocket}")

    promiscuous_ids = set(pd.read_csv(PROMISCUOUS_HITS_CSV, usecols=["id"])["id"])
    print(f"\nLoaded {len(promiscuous_ids)} promiscuous compound IDs to exclude")

    tables = []
    for gene, pocket in zip(noncat["gene_name"], noncat["pocket_name"]):
        report = pd.read_csv(os.path.join(REAL_ROUND1_RESULTS_DIR, pocket, "report.csv"))
        report = report.sort_values("score", ascending=True)
        report = report[~report["compound"].isin(promiscuous_ids)]
        top = report.head(TOP_N).copy()
        top["gene_name"] = gene
        top["pocket_name"] = pocket
        top["rank"] = range(1, len(top) + 1)
        print(f"{pocket}: 100th-place docking score = {top['score'].iloc[TOP_N - 1]:.3f}")
        tables.append(top)

    result = pd.concat(tables, ignore_index=True)

    print(f"\nLooking up SMILES for {result['compound'].nunique()} unique compounds...")
    smiles_map = lookup_smiles(result["compound"].unique().tolist(), "REAL_ROUND1")
    result["smiles"] = result["compound"].map(smiles_map)

    print("\nComputing physicochemical properties for the selected compounds...")
    sel_smiles = dict(zip(result["compound"], result["smiles"]))
    sel_props = compute_properties(sel_smiles)
    result = result.merge(sel_props, left_on="compound", right_index=True, how="left")

    result = result[[
        "gene_name", "pocket_name", "rank", "compound", "score", "smiles",
        "MW", "cLogP", "TPSA", "HBD", "HBA", "RotBonds", "AromaticRings", "QED", "is_pains",
    ]]

    out_path = os.path.join(OUTPUT_DIR, "top100_per_noncat_pocket.csv")
    result.to_csv(out_path, index=False)
    print(f"\nSaved {len(result)} rows to {out_path}")

    dup_counts = result["compound"].value_counts()
    n_dup_compounds = int((dup_counts > 1).sum())
    n_dup_rows = int(dup_counts[dup_counts > 1].sum())
    print(f"\n{n_dup_compounds} compounds appear in more than one pocket's top-{TOP_N} "
          f"({n_dup_rows} rows total) - expected duplication across pockets, not an error.")

    print("\nBuilding score boxplots across 5 reference distributions...")
    pockets = list(noncat["pocket_name"])
    selected_scores = {
        pocket: result.loc[result["pocket_name"] == pocket, "score"]
        for pocket in pockets
    }
    dl_scores = {pocket: load_scores(os.path.join(LIBRARIES["DL"], pocket, "report.csv"))
                 for pocket in pockets}
    real1_neg_scores = {pocket: load_real_negative_scores(pocket) for pocket in pockets}
    real1_pos_scores = {pocket: load_real_positive_scores(pocket) for pocket in pockets}
    real2_all_scores = {pocket: load_scores(os.path.join(LIBRARIES["REAL"], pocket, "report.csv"))
                         for pocket in pockets}

    # common=True: same compound set scored against every pocket (Hit Locator, REAL 2 - all, and
    # REAL 1 - negatives, which is a single fixed ~12,958-compound background sample shared by
    # every pocket - confirmed identical across pockets). common=False: genuinely per-pocket
    # (REAL 1 - positives is each pocket's own top-100k surrogate ranking; REAL 1 - selected is
    # each pocket's own top-100 post-filter picks).
    score_sources = [
        ("Hit Locator", "mint", dl_scores, True),
        ("REAL 1 - negatives", "blue", real1_neg_scores, True),
        ("REAL 1 - positives", "orange", real1_pos_scores, False),
        ("REAL 2 - all", "yellow", real2_all_scores, True),
        ("REAL 1 - selected", "purple", selected_scores, False),
    ]
    plot_path = os.path.join(OUTPUT_DIR, "noncat_score_boxplots.png")
    plot_score_boxplots_multi(pockets, score_sources, plot_path, xtick_rotation=20)
    print(f"Saved boxplot to {plot_path}")

    print("\nBuilding physicochemical profile vs. a REAL round-1 background sample...")
    bg_smiles = sample_prescreened_smiles("REAL_ROUND1", set(sel_smiles), BG_SAMPLE_SIZE, RANDOM_SEED)
    bg_props = compute_properties(bg_smiles)
    profiling_path = os.path.join(OUTPUT_DIR, "noncat_top100_physchem_profile.png")
    plot_profiling(sel_props, bg_props, profiling_path)
    print(f"Saved physchem profile to {profiling_path}")

    print("\nBuilding t-SNE (ECFP4) projection...")
    tsne_path = os.path.join(OUTPUT_DIR, "noncat_top100_tsne.png")
    plot_tsne(sel_smiles, bg_smiles, tsne_path, seed=RANDOM_SEED)
    print(f"Saved t-SNE figure to {tsne_path}")


if __name__ == "__main__":
    main()
