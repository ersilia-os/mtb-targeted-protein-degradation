#!/usr/bin/env python3
"""
Merges script 96's docking results (pockets x 1,095 compounds) with pocket annotation and script
70's compound info (SMILES, source), into one long table -- one row per (pocket, compound).

--organism human (default, 389 pockets x 1,095 = 425,955 endpoints): merges with script 93's
pocket annotation (InterPro domain labels + classified AlphaFill ligand evidence).

--organism mtb (all 21 Mtb CRISPR-screen genes' AF2-monomer pockets): merges with script 91's raw
pocket_detection_data.csv directly instead -- scripts 92/93 (InterPro/AlphaFill annotation) are
deliberately not mirrored for Mtb, since Mtb already has a much richer, curated version of that via
script 77, so the Mtb output has no InterPro/AlphaFill columns. Expected row count is computed
dynamically (n_pockets_detected * 1,095) rather than hardcoded, since the Mtb pocket count isn't
fixed like the human screen's 389.

Also writes two more per-compound files, both reindexed to the full 1,095-compound list and full
expected gene list (NaN for a gene not detected/docked yet, same partial-run handling as script
98's best_score_per_gene) so their shape is always the final one regardless of how much of the
counter-screen has completed so far:

* `gene_min_scores.csv` -- wide table, one row per compound, one column per gene (alphabetically
  sorted) holding that gene's best (lowest) score across all of its own pockets.
* `top_n_summary.csv` -- top1/top5/top10 = the Nth-most-favorable value among a compound's
  per-gene best scores (from the file above) -- e.g. `top5` is the 5th-most-favored of this
  organism's genes, a promiscuity/robustness read ("even the 5th-most-favored off-target scores
  this well"), not an average. Same convention as script 98's own human_top1/5/10 / mtb_top1/5/10
  columns -- kept here too so this organism's summary is available directly off of script 97,
  without needing script 98's full cross-pipeline join.

Usage:
    python 97_merge_docking_scores.py [--organism human|mtb]
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

TOP_NS = [1, 5, 10]

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

FILTERED_HITS_CSV = os.path.join(ROOT, "output", "70_filtering", "filtered_hits.csv")

N_COMPOUNDS = 1095


def load_docking_scores(docking_results_dir):
    """Long table: Uniprot AC, Pocket number, compound_id, score -- one row per (pocket, compound),
    read directly from script 96's per-pocket report.csv files."""
    report_paths = sorted(glob.glob(os.path.join(docking_results_dir, "*", "*", "report.csv")))
    rows = []
    for path in report_paths:
        pocket_number = int(os.path.basename(os.path.dirname(path)))
        uniprot_ac = os.path.basename(os.path.dirname(os.path.dirname(path)))
        df = pd.read_csv(path).rename(columns={"compound": "compound_id"})
        df.insert(0, "Pocket number", pocket_number)
        df.insert(0, "Uniprot AC", uniprot_ac)
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


def min_score_per_gene(merged, fallback_genes, fallback_compounds):
    """Wide table (compound_id index, one alphabetically-sorted column per gene) -- best (min)
    score per (compound, gene) across all of that gene's own pockets, reindexed to the full
    expected gene list and the full 1,095-compound list (NaN where not docked/detected yet)."""
    per_gene = merged.groupby(["compound_id", "Gene name"])["score"].min().unstack("Gene name")
    return per_gene.reindex(index=fallback_compounds, columns=sorted(fallback_genes))


def top_n_summary(per_gene):
    """{compound_id: top1/top5/top10} -- the n-th smallest (n-th most favorable) value per row of
    min_score_per_gene()'s wide table, for each n in TOP_NS (ascending sort, so index 0 = best/most
    negative score; NaNs -- a gene not yet docked -- sort to the end)."""
    sorted_vals = np.sort(per_gene.to_numpy(), axis=1)
    out = pd.DataFrame(index=per_gene.index)
    for n in TOP_NS:
        out[f"top{n}"] = sorted_vals[:, n - 1]
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--organism", choices=["human", "mtb"], default="human",
                         help="Which counter-screen's docking results to merge (default: human)")
    args = parser.parse_args()

    docking_results_dir = os.path.join(ROOT, "output", f"96_{args.organism}_docking", "docking_results")
    output_dir = os.path.join(ROOT, "output", f"97_{args.organism}_merge_docking_scores")
    os.makedirs(output_dir, exist_ok=True)

    print("Loading docking scores from all pockets...")
    scores = load_docking_scores(docking_results_dir)
    n_pockets = scores.groupby(["Uniprot AC", "Pocket number"]).ngroups
    print(f"  {len(scores):,} (pocket, compound) rows from {n_pockets} pockets")

    if args.organism == "human":
        pocket_annotation_csv = os.path.join(ROOT, "output", "93_human_merge_pocket_annotations", "merged_pocket_data.csv")
    else:
        pocket_annotation_csv = os.path.join(ROOT, "output", "91_mtb_detect_pockets", "pocket_detection_data.csv")
    pocket_annotation = pd.read_csv(pocket_annotation_csv)
    merged = scores.merge(pocket_annotation, on=["Uniprot AC", "Pocket number"],
                           how="left", validate="many_to_one")
    assert len(merged) == len(scores), "pocket-annotation merge changed row count"

    compounds = pd.read_csv(FILTERED_HITS_CSV)[["compound_id", "smiles", "source"]]
    merged = merged.merge(compounds, on="compound_id", how="left", validate="many_to_one")
    assert len(merged) == len(scores), "compound-info merge changed row count"

    out_path = os.path.join(output_dir, "docking_scores.csv")
    merged.to_csv(out_path, index=False)

    n_expected = n_pockets * N_COMPOUNDS
    print(f"\nSaved {len(merged):,} rows x {len(merged.columns)} columns -> {out_path}")
    print(f"Expected {n_expected:,} rows ({n_pockets} pockets x {N_COMPOUNDS:,} compounds): "
          f"{'MATCH' if len(merged) == n_expected else 'MISMATCH'}")
    print(f"\nScore summary:\n{merged['score'].describe()}")
    print(f"\nMissing scores: {merged['score'].isna().sum()}")

    fallback_genes = sorted(pocket_annotation["Gene name"].unique())
    fallback_compounds = compounds["compound_id"]
    per_gene = min_score_per_gene(merged, fallback_genes, fallback_compounds)

    gene_min_path = os.path.join(output_dir, "gene_min_scores.csv")
    per_gene.reset_index(names="compound_id").to_csv(gene_min_path, index=False)
    print(f"\nSaved {len(per_gene):,} rows x {len(per_gene.columns)} genes -> {gene_min_path}")
    print(f"Missing (compound, gene) cells: {per_gene.isna().sum().sum()}")

    topn = top_n_summary(per_gene)
    topn_path = os.path.join(output_dir, "top_n_summary.csv")
    topn.reset_index(names="compound_id").to_csv(topn_path, index=False)
    print(f"\nSaved top-N summary ({', '.join(f'top{n}' for n in TOP_NS)}) -> {topn_path}")


if __name__ == "__main__":
    main()
