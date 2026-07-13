#!/usr/bin/env python3
"""
Aggregate the 5 compound sets from scripts 52-61 into one CSV: compound_id, smiles, source,
physchem properties. Not a final hit list - script 61's set is still partial.

Also produces a physchem/PAINS profile and a t-SNE plot with the 5 sources overlaid (no single
background library applies to all 5, unlike scripts 52-54/56/61's selected-vs-background versions).

Usage:
    python 62_aggregate_hits.py
"""
import os
import sys

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from default import RANDOM_SEED
from docking_utils import compute_properties, plot_profiling_multi, plot_tsne_multi

CAT_PROMISCUOUS_CSV = os.path.join(ROOT, "output", "52_CAT_promiscuous",
                                    "alaS_aspS_lysS_pheS_REAL_CAT_multi_target_top1000.csv")
CAT_SELECTIVE_CSV = os.path.join(ROOT, "output", "53_CAT_selective", "alaS_aspS_lysS_pheS_REAL_CAT.csv")
NONCAT_PROMISCUOUS_CSV = os.path.join(ROOT, "output", "54_NONCAT_promiscuous",
                                       "alaS_aspS_lysS_pheST_REAL_NONCAT_multi_target_top1000.csv")
NONCAT_TOP100_10M_CSV = os.path.join(ROOT, "output", "56_NONCAT_top100_selection", "top100_per_noncat_pocket.csv")
NONCAT_TOP100_10B_CSV = os.path.join(ROOT, "output", "61_docking_top100_diverse_selection",
                                      "top100_diverse_per_pocket.csv")

OUTPUT_DIR = os.path.join(ROOT, "output", "62_aggregate_hits")
os.makedirs(OUTPUT_DIR, exist_ok=True)

SOURCE_ORDER = ["cat_promiscuous", "cat_selective", "noncat_promiscuous", "noncat_top100_10m", "noncat_top100_10b"]
SOURCE_COLORS = {
    "cat_promiscuous": "blue",
    "cat_selective": "purple",
    "noncat_promiscuous": "orange",
    "noncat_top100_10m": "mint",
    "noncat_top100_10b": "pink",
}


def merge_into(agg, compound_id, smiles, source):
    entry = agg.setdefault(compound_id, {"smiles": smiles, "sources": set()})
    if entry["smiles"] != smiles:
        print(f"Warning: SMILES mismatch for {compound_id}: '{entry['smiles']}' vs '{smiles}' (keeping first)")
    entry["sources"].add(source)


def process_simple(agg, path, compound_col, source_name):
    df = pd.read_csv(path)
    for _, row in df.iterrows():
        merge_into(agg, row[compound_col], row["smiles"], source_name)
    print(f"{source_name}: {len(df)} rows")
    return len(df)


def process_grouped_by_pocket(agg, path, compound_col, source_name):
    """Handles the two per-(compound, pocket)-row sources (56, 61): dedupes by compound_col,
    keeping the first SMILES seen for each compound."""
    df = pd.read_csv(path)
    for compound_id, group in df.groupby(compound_col):
        if group["smiles"].nunique() > 1:
            print(f"Warning: {source_name} has >1 distinct SMILES for {compound_id}, using the first.")
        merge_into(agg, compound_id, group["smiles"].iloc[0], source_name)
    print(f"{source_name}: {len(df)} rows, {df[compound_col].nunique()} unique compounds")
    return len(df)


def main():
    agg = {}
    n_rows = {
        "cat_promiscuous": process_simple(agg, CAT_PROMISCUOUS_CSV, "compound", "cat_promiscuous"),
        "cat_selective": process_simple(agg, CAT_SELECTIVE_CSV, "compound", "cat_selective"),
        "noncat_promiscuous": process_simple(agg, NONCAT_PROMISCUOUS_CSV, "compound", "noncat_promiscuous"),
        "noncat_top100_10m": process_grouped_by_pocket(agg, NONCAT_TOP100_10M_CSV, "compound", "noncat_top100_10m"),
        "noncat_top100_10b": process_grouped_by_pocket(agg, NONCAT_TOP100_10B_CSV, "compound_id", "noncat_top100_10b"),
    }

    rows = []
    n_multi_source = 0
    for compound_id, entry in agg.items():
        sources = sorted(entry["sources"], key=SOURCE_ORDER.index)
        if len(sources) > 1:
            n_multi_source += 1
        rows.append({
            "compound_id": compound_id,
            "smiles": entry["smiles"],
            "source": ";".join(sources),
        })

    result = pd.DataFrame(rows).sort_values("compound_id").reset_index(drop=True)

    print("\nComputing physicochemical properties...")
    props = compute_properties(dict(zip(result["compound_id"], result["smiles"])))
    result = result.merge(props, left_on="compound_id", right_index=True, how="left")

    out_path = os.path.join(OUTPUT_DIR, "aggregated_hits.csv")
    result.to_csv(out_path, index=False)

    print(f"\nRows per source: {n_rows}")
    print(f"Sum of per-source rows (raw, before dedup): {sum(n_rows.values())}")
    print(f"Compounds appearing in more than one source: {n_multi_source}")
    print(f"Aggregated unique compounds: {len(result)}")
    print(f"Saved: {out_path}")

    source_ids = {}
    for source in SOURCE_ORDER:
        in_source = result["source"].apply(lambda s: source in s.split(";"))
        source_ids[source] = result.loc[in_source, "compound_id"]

    print("\nBuilding per-source physchem/PAINS profile...")
    source_props = [(source, SOURCE_COLORS[source], props.loc[props.index.isin(source_ids[source])])
                     for source in SOURCE_ORDER]
    profile_path = os.path.join(OUTPUT_DIR, "aggregated_hits_physchem_profile.png")
    plot_profiling_multi(source_props, profile_path)
    print(f"Saved: {profile_path}")

    print("\nBuilding per-source t-SNE (ECFP4)...")
    smiles_lookup = dict(zip(result["compound_id"], result["smiles"]))
    source_smiles = [(source, SOURCE_COLORS[source],
                       {cid: smiles_lookup[cid] for cid in source_ids[source]})
                      for source in SOURCE_ORDER]
    tsne_path = os.path.join(OUTPUT_DIR, "aggregated_hits_tsne.png")
    plot_tsne_multi(source_smiles, tsne_path, seed=RANDOM_SEED)
    print(f"Saved: {tsne_path}")


if __name__ == "__main__":
    main()
