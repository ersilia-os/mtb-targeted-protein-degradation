#!/usr/bin/env python3
"""
Rank-matrix selectivity analysis for tRNA synthetase docking hits.

Builds two rank matrices:
  Matrix 1 — compounds × reference pockets (target genes only)
  Matrix 2 — compounds × all non-target pockets

Applies five complementary selectivity metrics (m1–m5) to select ~500 hits,
with cross-metric scaffold deduplication.

Usage:
    python 49_docking_hits_selective.py --trna pheS,aspS,lysS,alaS --lib REAL
    python 49_docking_hits_selective.py --trna pheS,aspS,lysS,alaS --lib REAL
"""

import argparse
import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
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
    load_reference_pockets,
    lookup_smiles,
    plot_profiling,
    plot_score_boxplots,
    sample_background_smiles,
)

BG_SAMPLE_SIZE = 25_000


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
    parser = argparse.ArgumentParser(
        description="Rank-matrix selectivity analysis for tRNA synthetase docking hits."
    )
    parser.add_argument("--trna", required=True,
                        help="Comma-separated gene name(s) (e.g. pheS,aspS,lysS,alaS)")
    parser.add_argument("--lib", choices=["DL", "REAL"], required=True,
                        help="Compound library: DL or REAL")
    args = parser.parse_args()

    genes = [g.strip() for g in args.trna.split(",")]
    results_dir = LIBRARIES[args.lib]

    gene_map = load_gene_map()
    ref_pockets = load_reference_pockets()

    invalid = [g for g in genes if g not in gene_map]
    if invalid:
        print(f"Error: unknown gene(s): {', '.join(invalid)}. Available genes:")
        print("  " + ", ".join(sorted(gene_map)))
        sys.exit(1)

    missing_ref = [g for g in genes if g not in ref_pockets]
    if missing_ref:
        print(f"Error: no reference pocket defined for: {', '.join(missing_ref)}")
        sys.exit(1)

    target_uniprot_acs = {gene_map[g] for g in genes}

    output_dir = os.path.join(ROOT, "output", "49_docking_hits_selective")
    os.makedirs(output_dir, exist_ok=True)
    trna_tag = "_".join(genes)

    # --- Matrix 1: target reference pockets ---
    target_pocket_map = {g: ref_pockets[g] for g in genes}

    # --- Matrix 2: all pockets from non-target proteins ---
    all_pockets = sorted(
        p for p in os.listdir(results_dir)
        if os.path.isdir(os.path.join(results_dir, p))
    )
    non_target_pocket_map = {
        p: p for p in all_pockets
        if not any(ac in p for ac in target_uniprot_acs)
    }

    print(f"Library : {args.lib}")
    print(f"Targets : {', '.join(genes)}")
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
    _M_SORT = {"m1": True, "m2": True, "m3": False, "m4": False}
    _union_ids = set()
    for _m, _asc in _M_SORT.items():
        _top50 = results[results[_m].notna()].sort_values(_m, ascending=_asc).head(50)
        _union_ids.update(_top50.index)
    _union_smiles = lookup_smiles(list(_union_ids), args.lib)
    _union_scaffold_set = {
        ik for ik in (murcko_inchikey(smi) for smi in _union_smiles.values())
        if ik is not None
    }
    print(f"\nUnion (top-50 × m1–m4): {len(_union_ids):,} compounds, "
          f"{len(_union_scaffold_set):,} unique Murcko scaffolds")

    _c1 = results[results["top2_rank"] <= 20_000]
    _c2 = _c1[_c1["nontarget_p50"] >= 50_000]
    _rescue = _c2[~_c2.index.isin(_union_ids)].sort_values("top2_rank", ascending=True).copy()
    print(f"Rescue pool funnel:")
    print(f"  top2_rank ≤ 20k              : {len(_c1):,}")
    print(f"  + nontarget_p50 ≥ 50k        : {len(_c2):,}")
    print(f"  + not in union ({len(_union_ids)} compounds): {len(_rescue):,}")

    _rescue_smiles = lookup_smiles(_rescue.index.tolist(), args.lib)
    _rescue["_smiles"] = _rescue.index.map(_rescue_smiles)
    _rescue["_scaffold_ik"] = _rescue["_smiles"].map(murcko_inchikey)
    _valid_rescue = _rescue[_rescue["_scaffold_ik"].notna()]
    _n_no_smiles = _rescue["_smiles"].isna().sum()
    _n_bad_scaffold = (_rescue["_smiles"].notna() & _rescue["_scaffold_ik"].isna()).sum()
    print(f"  + valid Murcko scaffold          : {len(_valid_rescue):,}"
          + (f"  ({_n_no_smiles} SMILES not found, {_n_bad_scaffold} invalid scaffold)"
             if _n_no_smiles or _n_bad_scaffold else ""))

    _dedup_m1_m4 = set()
    for _m, _asc in [("m1", True), ("m2", True), ("m3", False), ("m4", False)]:
        _h = results[results[_m].notna()].sort_values(_m, ascending=_asc)
        _dedup_m1_m4.update(_h[~_h.index.isin(_dedup_m1_m4)].head(50).index)
    _m5_cap = 500 - len(_dedup_m1_m4)
    _novel = _valid_rescue[~_valid_rescue["_scaffold_ik"].isin(_union_scaffold_set)]
    _m5 = _novel.drop_duplicates(subset="_scaffold_ik")
    print(f"  + scaffold not in union          : {len(_novel):,}"
          + (f"  ({len(_valid_rescue) - len(_novel)} scaffold(s) already covered by union)"
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


    out_path = os.path.join(output_dir, f"{trna_tag}_{args.lib}.csv")
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

    sel_smiles = lookup_smiles(list(_already_selected), args.lib)
    results["smiles"] = results.index.map(sel_smiles)
    print(f"\n  Computing properties for {len(sel_smiles):,} selected...")
    sel_props = compute_properties(sel_smiles)
    results = results.join(sel_props, how="left")

    results.to_csv(out_path)
    print(f"Saved results to {out_path}")

    scores_path = os.path.join(output_dir, f"{trna_tag}_{args.lib}_scores.png")
    plot_score_boxplots(scores1, _already_selected, genes, scores_path,
                        sel_label="Selected (m1–m5)", sel_color="orange")
    print(f"Saved score boxplots to {scores_path}")

    print("\nPhysicochemical profiling...")
    bg_smiles = sample_background_smiles(args.lib, _already_selected, BG_SAMPLE_SIZE, RANDOM_SEED)
    print(f"  Background: {len(bg_smiles):,} compounds...")
    bg_props  = compute_properties(bg_smiles)
    profiling_path = os.path.join(output_dir, f"{trna_tag}_{args.lib}_profiling.png")
    plot_profiling(sel_props, bg_props, profiling_path, sel_color="orange")
    print(f"  Saved profiling figure to {profiling_path}")


if __name__ == "__main__":
    main()
