#!/usr/bin/env python3
"""
Rank-matrix selectivity analysis restricted to the 4 catalytic (CAT) pockets - pheS, aspS, lysS,
alaS (pheT has no CAT entry) - against the second Enamine REAL screening
(output/unidock_REAL_docking_2, ~99,105 compounds).

Uses the shared helpers in src/docking_utils.py. Reference pockets are loaded directly from
output/reference_pocket_catalytic.csv (output/reference_pocket.csv no longer exists, having been
split into output/reference_pocket_catalytic.csv / _noncatalytic.csv), via the same
load_cat_pockets() pattern as scripts/52_CAT_promiscuous.py.

Builds two rank matrices:
  Matrix 1 - compounds x the 4 CAT reference pockets
  Matrix 2 - compounds x all non-target pockets (every pocket not belonging to one of the 4
             CAT genes' UniProt ACs, PLUS pheT - the other subunit of pheS's PheRS heterodimer,
             excluded even though it has no CAT entry of its own - i.e. the full 276-pocket
             library minus our own 4 genes' 40 pockets minus pheT's 13, leaving 223)

Applies five complementary selectivity metrics (m1-m5) to select ~500 hits, with cross-metric
scaffold deduplication:
  m1: max potency, low selectivity bar (nontarget_p50 >= 20k)
  m2: max potency, high selectivity bar (nontarget_p50 >= 50k)
  m3: potency + selectivity gap (nontarget_p10 - max_target_rank)
  m4: potency + max selectivity (nontarget_p1)
  m5: diversity rescue - unique Murcko scaffolds among compounds hitting >=2 targets well

For the ~500 selected compounds vs. a random 20,000-compound REAL background:
  - UpSet plots of target membership (top-100/1,000/10,000)
  - score boxplots (vs. DL/REAL-round-1-negative-set references)
  - physchem profiling
  - a t-SNE projection of ECFP4 fingerprints

Then reports overlap with script 52's CAT multi-target top-1,000 binder list.

Usage:
    python 53_CAT_selective.py
"""
import os
import sys
from collections import Counter
from itertools import combinations

import matplotlib
matplotlib.use("Agg")
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*converter.*already registered.*")

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from default import RANDOM_SEED
from docking_utils import (
    LIBRARIES,
    build_matrix,
    compute_properties,
    load_gene_map,
    load_real_negative_scores,
    lookup_smiles,
    plot_profiling,
    plot_score_boxplots,
    plot_tsne,
    sample_prescreened_smiles,
    save_upset,
)

CAT_POCKETS_CSV = os.path.join(ROOT, "output", "reference_pocket_catalytic.csv")
OUTPUT_DIR = os.path.join(ROOT, "output", "53_CAT_selective")

LIB = "REAL"
BG_SAMPLE_SIZE = 20_000
UPSET_THRESHOLDS = [100, 1_000, 10_000]
BOXPLOT_GENE_ORDER = ["pheS", "aspS", "lysS", "alaS"]


def load_cat_pockets():
    """gene_name -> pocket_name for the 4 catalytic pockets."""
    df = pd.read_csv(CAT_POCKETS_CSV)
    return dict(zip(df["gene_name"], df["pocket_name"]))


def murcko_inchikey(smi):
    if not isinstance(smi, str):
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        return None
    return Chem.MolToInchiKey(scaffold)


def rank_transform(df):
    """Per-column ascending rank: rank 0 = lowest score (best binder)."""
    arr = df.to_numpy(dtype=np.float64)
    ranks = np.argsort(np.argsort(arr, axis=0), axis=0)
    return pd.DataFrame(ranks, index=df.index, columns=df.columns)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    gene_map = load_gene_map()
    target_pocket_map = load_cat_pockets()
    genes = sorted(target_pocket_map)
    trna_tag = "_".join(genes)
    target_uniprot_acs = {gene_map[g] for g in genes}
    if "pheS" in genes:
        # pheT is the other subunit of the same PheRS heterodimer as pheS - exclude its
        # pockets from the non-target background too, not just pheS's own (same convention
        # as scripts/54_NONCAT_promiscuous.py's combined "pheST" target unit).
        target_uniprot_acs.add(gene_map["pheT"])
    results_dir = LIBRARIES[LIB]

    print(f"Library : {LIB}")
    print(f"Targets : {', '.join(genes)}")

    # --- Matrix 2: all pockets from non-target proteins ---
    all_pockets = sorted(
        p for p in os.listdir(results_dir)
        if os.path.isdir(os.path.join(results_dir, p))
    )
    non_target_pocket_map = {
        p: p for p in all_pockets
        if not any(ac in p for ac in target_uniprot_acs)
    }

    print(f"Building Matrix 1 ({len(target_pocket_map)} reference pocket(s))...")
    scores1 = build_matrix(target_pocket_map, results_dir, label="Matrix 1")

    print(f"Building Matrix 2 ({len(non_target_pocket_map)} non-target pocket(s))...")
    scores2 = build_matrix(non_target_pocket_map, results_dir, label="Matrix 2")

    errors_found = False
    for name, df in [("Matrix 1", scores1), ("Matrix 2", scores2)]:
        nan_counts = df.isna().sum()
        bad = nan_counts[nan_counts > 0]
        if not bad.empty:
            print(f"Error: NaN scores found in {name}:")
            for col, n in bad.items():
                print(f"  {col}: {n} NaN value(s)")
            errors_found = True
    if errors_found:
        sys.exit(1)

    # Per-gene compound ranking, reused below (after selection) for the UpSet plot.
    gene_top = {g: scores1[g].sort_values().index.tolist() for g in genes if g in scores1.columns}

    ranks1 = rank_transform(scores1)
    ranks2 = rank_transform(scores2).reindex(ranks1.index)

    # --- Selectivity metrics ---
    max_target_rank = ranks1.to_numpy().max(axis=1)
    p1, p5, p10, p50 = np.percentile(ranks2.to_numpy(), [1, 5, 10, 50], axis=1)

    results = pd.DataFrame({
        "max_target_rank": max_target_rank,
        **{f"rank_{g}": ranks1[g] for g in genes},
        "nontarget_p1":  np.round(p1).astype(int),
        "nontarget_p5":  np.round(p5).astype(int),
        "nontarget_p10": np.round(p10).astype(int),
        "nontarget_p50": np.round(p50).astype(int),
    }, index=ranks1.index)
    results.index.name = "compound"

    # m1: max potency, low selectivity
    results["m1"] = (
        results["max_target_rank"]
        .where(results["nontarget_p50"] >= 20_000)
        .astype(pd.Int64Dtype())
    )
    # m2: max potency, high selectivity
    results["m2"] = (
        results["max_target_rank"]
        .where(results["nontarget_p50"] >= 50_000)
        .astype(pd.Int64Dtype())
    )
    # m3: high potency, selectivity gap (nontarget_p10 − max_target_rank)
    _gap = results["nontarget_p10"] - results["max_target_rank"]
    results["m3"] = (
        _gap
        .where((results["max_target_rank"] <= 20_000) & (_gap > 0) & (results["nontarget_p50"] >= 20_000))
        .astype(pd.Int64Dtype())
    )
    # m4: high potency, max selectivity
    results["m4"] = (
        results["nontarget_p1"]
        .where(results["max_target_rank"] <= 20_000)
        .astype(pd.Int64Dtype())
    )

    # top2_rank: 2nd-best target rank per compound (compound must bind well to ≥2 targets)
    _rank_cols = [f"rank_{g}" for g in genes]
    if len(genes) >= 2:
        _arr = results[_rank_cols].to_numpy()
        _top2 = np.partition(_arr, 1, axis=1)[:, 1]
    else:
        _top2 = results[_rank_cols[0]].to_numpy()
    results.insert(1, "top2_rank", _top2.astype(int))

    # m5: diversity rescue
    # Sequential m1–m4 selection (same order/logic as the METRIC_CONFIGS loop below) —
    # used here to know which compounds/scaffolds m5 must treat as "already selected".
    _dedup_m1_m4 = set()
    for _m, _asc in [("m1", True), ("m2", True), ("m3", False), ("m4", False)]:
        _h = results[results[_m].notna()].sort_values(_m, ascending=_asc)
        _dedup_m1_m4.update(_h[~_h.index.isin(_dedup_m1_m4)].head(50).index)
    _m5_cap = 500 - len(_dedup_m1_m4)

    _dedup_smiles = lookup_smiles(list(_dedup_m1_m4), LIB)
    _dedup_scaffold_set = {
        ik for ik in (murcko_inchikey(smi) for smi in _dedup_smiles.values())
        if ik is not None
    }
    print(f"\nm1–m4 selection: {len(_dedup_m1_m4):,} compounds, "
          f"{len(_dedup_scaffold_set):,} unique Murcko scaffolds")

    _c1 = results[results["top2_rank"] <= 20_000]
    _c2 = _c1[_c1["nontarget_p50"] >= 50_000]
    _rescue = _c2[~_c2.index.isin(_dedup_m1_m4)].sort_values("top2_rank", ascending=True).copy()
    print(f"Rescue pool funnel:")
    print(f"  top2_rank ≤ 20k              : {len(_c1):,}")
    print(f"  + nontarget_p50 ≥ 50k        : {len(_c2):,}")
    print(f"  + not in m1-m4 ({len(_dedup_m1_m4)} compounds): {len(_rescue):,}")

    _rescue_smiles = lookup_smiles(_rescue.index.tolist(), LIB)
    _rescue["_smiles"] = _rescue.index.map(_rescue_smiles)
    _rescue["_scaffold_ik"] = _rescue["_smiles"].map(murcko_inchikey)
    _valid_rescue = _rescue[_rescue["_scaffold_ik"].notna()]
    _n_no_smiles = _rescue["_smiles"].isna().sum()
    _n_bad_scaffold = (_rescue["_smiles"].notna() & _rescue["_scaffold_ik"].isna()).sum()
    print(f"  + valid Murcko scaffold          : {len(_valid_rescue):,}"
          + (f"  ({_n_no_smiles} SMILES not found, {_n_bad_scaffold} invalid scaffold)"
             if _n_no_smiles or _n_bad_scaffold else ""))

    _novel = _valid_rescue[~_valid_rescue["_scaffold_ik"].isin(_dedup_scaffold_set)]
    _m5 = _novel.drop_duplicates(subset="_scaffold_ik")
    print(f"  + scaffold not in m1-m4          : {len(_novel):,}"
          + (f"  ({len(_valid_rescue) - len(_novel)} scaffold(s) already covered by m1-m4)"
             if len(_valid_rescue) > len(_novel) else ""))
    print(f"  + one per scaffold (dedup)       : {len(_m5):,}"
          + (f"  ({len(_novel) - len(_m5)} duplicate scaffold(s) dropped)"
             if len(_novel) > len(_m5) else ""))
    print(f"m5: {len(_m5):,} qualifying compounds  "
          f"(will display top {_m5_cap} = 500 − {len(_dedup_m1_m4)} actual m1–m4 selection)")

    results["m5"] = (
        results["top2_rank"]
        .where(results.index.isin(_m5.index))
        .astype(pd.Int64Dtype())
    )


    out_path = os.path.join(OUTPUT_DIR, f"{trna_tag}_{LIB}_CAT.csv")
    results["selected"] = pd.NA

    METRIC_CONFIGS = [
        ("m1", True,  "max potency, low selectivity → sorted by max_target_rank ascending, excluding nontarget_p50 < 20k",  50),
        ("m2", True,  "max potency, high selectivity → sorted by max_target_rank ascending, excluding nontarget_p50 < 50k",  50),
        ("m3", False, "high potency, selectivity gap → sorted by (nontarget_p10 − max_target_rank) descending, excluding max_target_rank > 20k, nontarget_p10 <= max_target_rank, and nontarget_p50 < 20k", 50),
        ("m4", False, "high potency, max selectivity → sorted by nontarget_p1 descending, excluding max_target_rank > 20k", 50),
        ("m5", True,  "diversity rescue → unique new Murcko scaffolds, top2_rank ≤ 20k, nontarget_p50 ≥ 50k, sorted by top2_rank ascending", _m5_cap),
    ]
    N = len(results)
    _already_selected = set()
    for metric, asc, desc, top_n in METRIC_CONFIGS:
        hits = results[results[metric].notna()].sort_values(metric, ascending=asc)
        hits_new = hits[~hits.index.isin(_already_selected)]
        top_hits = hits_new.head(top_n)
        results.loc[top_hits.index, "selected"] = metric
        _already_selected.update(top_hits.index)
        display = results.loc[top_hits.index[:10]].copy()
        print(f"\n— {metric}: {desc}")
        print(f"Passing filter: {len(hits):,} / {N:,}  ({len(hits_new):,} not yet selected in a prior metric)")
        print(display.to_string())

    # --- UpSet plots: target membership (top 100 / 1,000 / 10,000, as in script 52) among selected ---
    # Slice each gene's GLOBAL ranking to top_n first, then intersect with the selected set —
    # NOT the other way around (filtering to selected first would make the top_n slice a no-op,
    # since the selected set is smaller than every threshold here).
    if len(gene_top) >= 2:
        for top_n in UPSET_THRESHOLDS:
            gene_top_selected = {
                g: [cid for cid in ids[:top_n] if cid in _already_selected]
                for g, ids in gene_top.items()
            }
            top_sets = {g: set(ids) for g, ids in gene_top_selected.items()}
            print(f"\nShared hits among the {len(_already_selected):,} selected compounds "
                  f"— top {top_n:,}:")
            for r in range(2, len(gene_top_selected) + 1):
                for combo in combinations(sorted(gene_top_selected), r):
                    shared = len(set.intersection(*(top_sets[g] for g in combo)))
                    print(f"  {' & '.join(combo)}: {shared:,}")

            # Multi-target binder counts (exactly k targets), same format as script 52.
            # No "expected by chance" figure here (unlike script 52): that estimate assumes a
            # random draw from the full ~M-compound population, but these counts are already
            # capped at len(_already_selected) — a pre-filtered "good binder" pool, not a random
            # sample — so a full-population chance expectation isn't a meaningful comparison.
            all_in_top = set.union(*top_sets.values()) if any(top_sets.values()) else set()
            tally = Counter(
                sum(1 for g in gene_top_selected if cid in top_sets[g])
                for cid in all_in_top
            )
            for k in range(len(gene_top_selected), 1, -1):
                print(f"  {k}/{len(gene_top_selected)} targets: {tally.get(k, 0):,} compound(s)")

            save_upset(gene_top_selected, top_n, OUTPUT_DIR, LIB, f"{trna_tag}_CAT")
        print()

    sel_smiles = lookup_smiles(list(_already_selected), LIB)
    results.insert(0, "smiles", results.index.map(sel_smiles))
    print(f"\n  Computing properties for {len(sel_smiles):,} selected...")
    sel_props = compute_properties(sel_smiles)
    results = results.join(sel_props, how="left")

    results = results[results["selected"].notna()]
    results.to_csv(out_path)
    print(f"Saved results to {out_path} ({len(results):,} selected compounds)")

    print("Loading Enamine DL reference scores...")
    dl_scores = build_matrix(target_pocket_map, LIBRARIES["DL"], label="DL reference")
    print("Loading Enamine REAL 10M round-1 negative-set reference scores...")
    real_negative_scores = {g: load_real_negative_scores(target_pocket_map[g]) for g in genes}

    scores_path = os.path.join(OUTPUT_DIR, f"{trna_tag}_{LIB}_CAT_scores.png")
    plot_score_boxplots(scores1, _already_selected, BOXPLOT_GENE_ORDER, scores_path,
                        sel_label="Selected selective", sel_color="orange",
                        dl_scores=dl_scores, real_negative_scores=real_negative_scores)
    print(f"Saved score boxplots to {scores_path}")

    print("\nPhysicochemical profiling...")
    bg_smiles = sample_prescreened_smiles(LIB, _already_selected, BG_SAMPLE_SIZE, RANDOM_SEED)
    print(f"  Pre-screened: {len(bg_smiles):,} compounds...")
    bg_props  = compute_properties(bg_smiles)
    profiling_path = os.path.join(OUTPUT_DIR, f"{trna_tag}_{LIB}_CAT_profiling.png")
    plot_profiling(sel_props, bg_props, profiling_path, sel_color="orange")
    print(f"  Saved profiling figure to {profiling_path}")

    print("\nt-SNE (ECFP4)...")
    tsne_path = os.path.join(OUTPUT_DIR, f"{trna_tag}_{LIB}_CAT_tsne.png")
    plot_tsne(sel_smiles, bg_smiles, tsne_path, seed=RANDOM_SEED, sel_color="orange")
    print(f"  Saved t-SNE figure to {tsne_path}")

    # --- Overlap with script 52 multi-target binders ---
    raw_path = os.path.join(ROOT, "output", "52_CAT_promiscuous",
                            f"{trna_tag}_{LIB}_CAT_multi_target_top1000.csv")
    if os.path.isfile(raw_path):
        raw_ids = set(pd.read_csv(raw_path, index_col=0).index)
        overlap = _already_selected & raw_ids
        print(f"\nOverlap with script 52 multi-target binders (top 1,000):")
        print(f"  Script 52 selected : {len(raw_ids):,} compounds")
        print(f"  Script 53 selected : {len(_already_selected):,} compounds")
        print(f"  Overlap            : {len(overlap):,} compounds "
              f"({100 * len(overlap) / len(_already_selected):.1f}% of script 53 selection)")
    else:
        print(f"\nSkipping overlap with script 52: {raw_path} not found.")
        print("  Run script 52 first.")


if __name__ == "__main__":
    main()
