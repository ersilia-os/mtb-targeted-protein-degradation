#!/usr/bin/env python3
"""
Builds one summary table over the 1,095 filtered hits (script 70), one row per compound, combining
the two structure-based Uni-Dock counter-screens (scripts 90-97), the 12 hand-curated Mtb pockets
(scripts 62-66), and the two Nesso-1 co-folding counter-screens (scripts 78-88):

1. `human_<GENE>` (38 columns) -- best (lowest/most negative) Uni-Dock score per human gene, taken
   across ALL of that gene's detected pockets. Read directly from script 97's own
   `output/97_human_merge_docking_scores/gene_min_scores.csv` (already computed and reindexed to
   the full 1,095-compound x 38-gene shape there) rather than recomputed from the raw
   `docking_scores.csv` long table -- avoids duplicating that aggregation in two places.
2. `mtb_<GENE>` (21 columns) -- same, for the AF2-monomer-only Mtb counter-screen
   (`--organism mtb`, scripts 90-97), from `output/97_mtb_merge_docking_scores/gene_min_scores.csv`.
3. `docking_<group>_CAT` / `docking_<group>_NONCAT` (8 columns) -- best Uni-Dock score per curated
   Mtb pocket group (pheST/aspS/lysS/alaS x CAT/NON-CAT; pheS+pheT merged into "pheST", same
   convention as script 54), read directly from `output/selected_pockets.csv` rather than
   hardcoded, from `output/66_merge_docking_scores/merged_docking_scores.csv`. Prefixed
   "docking_" (not "mtb_") so it can't be confused with column 2's mtb_<GENE> AF2-only
   counter-screen columns -- same prefix the condensed file's raw per-pocket columns use.
4. `human_top1/top5/top10` and `mtb_top1/top5/top10` -- the Nth-best (Nth most negative) value
   among that organism's per-gene best Uni-Dock scores (columns 1/2 above), read directly from
   script 97's own `top_n_summary.csv` (same reasoning as columns 1/2 -- not recomputed here) --
   e.g. `human_top5` is the 5th most favorable of the 38 human_<GENE> values, a
   robustness/promiscuity read ("even the 5th-most-favored off-target scores this well"), not an
   average.
5. `nesso1_mtb_<GENE>` (21 columns) / `nesso1_human_<GENE>` (38 columns) -- Nesso-1's predicted
   IC50 in **nM (integer)** per gene, converted from its raw `affinity_pred_value` (log10 IC50, uM)
   via `ic50_nm()` -- order-preserving, so lower still means more potent. Nesso-1 is protein-level
   (no pocket conditioning, see script 78's docstring), so each (gene, compound) already has
   exactly one row -- no per-pocket aggregation needed, unlike the Uni-Dock columns above. Plus
   `nesso1_human_top1/5/10` and `nesso1_mtb_top1/5/10` -- same Nth-best-value convention as
   columns 4, computed over these per-gene IC50 (nM) columns instead of the Uni-Dock scores.
6. `boltz2_<group>_CAT` / `boltz2_<group>_NONCAT` (8 columns) -- best (lowest predicted IC50, nM,
   integer, same `ic50_nm()` conversion as column 5) Boltz-2 affinity per curated Mtb pocket group,
   same 12-pocket groups as columns 3 (reuses `curated_pocket_groups()`), from
   `output/75_boltz2_collect_affinities/affinity_results.csv`. Boltz-2 only covers the 12 curated
   Mtb pockets (scripts 71-75) -- no human run, no AF2-only 21-gene Mtb run.

Low-confidence affinity gating (both raw and aggregated Boltz-2/Nesso-1 columns above, columns 5
and 6): before any pivoting/aggregation, `affinity_pred_value` is set to NaN wherever:

* Boltz-2 -- `affinity_probability_binary < BOLTZ2_MIN_PROBABILITY` (0.4). Not a threshold from the
  Boltz-2 paper/docs (no per-prediction confidence cutoff is documented there -- confirmed by
  reading both the bioRxiv preprint and the official repo docs); this is a project-level policy
  choice, decided after inspecting our own 13,140-row distribution (median 0.425, so ~44% of rows
  are affected overall, unevenly across pockets: 22.5% for pheS_CAT vs. 72.2% for pheT_NONCAT).
* Nesso-1 -- `entropy_crop_pl == 0.0`, exactly. This one IS an explicit, author-stated rule from
  `docs/prediction.md` in `recursionpharma/nesso`: "If the value of entropy_crop_pl is 0.0, then it
  means the model wasn't able to confidently place the ligand, the predictions in those cases
  should not be trusted." Affects very few rows (3 of 22,995 Mtb, 0 of 41,610 human).

A NaN'd cell still participates in min()/aggregation as a missing value (skipped, not treated as
0 or dropped from the row entirely) -- a curated pocket GROUP's score is NaN only if every pocket
in that group was itself gated out.

Also saves a condensed second file, `compound_docking_summary_condensed.csv`: the 12 raw
(un-aggregated) curated-pocket scores for docking and for Boltz-2 (24 columns,
`docking_<gene>_<CAT|NONCAT>[_n]` / `boltz2_<gene>_<CAT|NONCAT>[_n]`, via curated_pocket_labels())
plus only the top1/5/10 columns for docking and Nesso-1 (human and mtb, 12 columns) -- everything
else (the 38/21 per-gene docking columns, the 38+21 per-gene Nesso-1 columns, and the aggregated
CAT/NONCAT columns) omitted.

Usage:
    python 98_compound_docking_summary.py
"""
import os

import numpy as np
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

FILTERED_HITS_CSV = os.path.join(ROOT, "output", "70_filtering", "filtered_hits.csv")
SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")
CURATED_DOCKING_CSV = os.path.join(ROOT, "output", "66_merge_docking_scores", "merged_docking_scores.csv")
HUMAN_GENE_MIN_CSV = os.path.join(ROOT, "output", "97_human_merge_docking_scores", "gene_min_scores.csv")
MTB_GENE_MIN_CSV = os.path.join(ROOT, "output", "97_mtb_merge_docking_scores", "gene_min_scores.csv")
HUMAN_TOPN_CSV = os.path.join(ROOT, "output", "97_human_merge_docking_scores", "top_n_summary.csv")
MTB_TOPN_CSV = os.path.join(ROOT, "output", "97_mtb_merge_docking_scores", "top_n_summary.csv")
NESSO1_MTB_CSV = os.path.join(ROOT, "output", "82_nesso1_collect_affinities", "affinity_results.csv")
NESSO1_HUMAN_CSV = os.path.join(ROOT, "output", "88_nesso1_human_collect_affinities", "affinity_results.csv")
BOLTZ2_MTB_CSV = os.path.join(ROOT, "output", "75_boltz2_collect_affinities", "affinity_results.csv")

OUTPUT_DIR = os.path.join(ROOT, "output", "98_compound_docking_summary")
os.makedirs(OUTPUT_DIR, exist_ok=True)

TOP_NS = [1, 5, 10]
BOLTZ2_MIN_PROBABILITY = 0.4  # see module docstring's "Low-confidence affinity gating" section


def ic50_nm(log10_ic50_um):
    """Converts Boltz-2/Nesso-1's raw affinity_pred_value (log10 of predicted IC50, in uM) to
    predicted IC50 in nM, rounded to the nearest integer. Order-preserving (10**x is monotonically
    increasing), so downstream min()/rank aggregation over the converted values is still correct."""
    return (10 ** log10_ic50_um * 1000).round().astype("Int64")


def load_per_gene_csv(csv_path, prefix):
    """Reads one of script 97's own per-organism summary files (gene_min_scores.csv or
    top_n_summary.csv, both indexed by compound_id and already reindexed there to the full
    1,095-compound x full-gene-list shape) and adds the "<prefix>_" column prefix ("human_"/"mtb_")
    script 98's summary table uses -- avoids recomputing the same per-gene min / top-N aggregation
    a second time from the raw docking_scores.csv long table."""
    return pd.read_csv(csv_path).set_index("compound_id").add_prefix(f"{prefix}_")


def curated_pocket_groups():
    """{"docking_<group>_<CAT|NONCAT>": [pocket_name, ...]} read from output/selected_pockets.csv,
    merging pheS+pheT into "pheST" (same convention as script 54 -- pheT has no CAT pocket of its
    own, so pheST's CAT group is pheS's single CAT pocket). Prefixed "docking_" rather than "mtb_"
    so it can't be confused with the mtb_<GENE> AF2-only counter-screen columns (different pipeline,
    different granularity) -- same prefix convention raw_docking_pocket_scores() uses."""
    df = pd.read_csv(SELECTED_POCKETS_CSV)
    df["group_gene"] = df["gene_name"].replace({"pheS": "pheST", "pheT": "pheST"})
    groups = {}
    for (gene, site_type), sub in df.groupby(["group_gene", "site_type"]):
        col = f"docking_{gene}_{site_type.replace('-', '')}"
        groups[col] = sub["pocket_name"].tolist()
    return groups


def curated_pocket_scores():
    """{"docking_<group>_<CAT|NONCAT>": Series indexed by compound_id} -- best (min) mean docking
    score per curated pocket group, from script 66's wide table."""
    wide = pd.read_csv(CURATED_DOCKING_CSV).set_index("compound_id")
    out = pd.DataFrame(index=wide.index)
    for col, pockets in curated_pocket_groups().items():
        out[col] = wide[pockets].min(axis=1)
    return out


def nesso1_scores(affinity_csv, prefix):
    """{"<prefix>_<gene>": Series indexed by compound_id} -- Nesso-1's predicted IC50 (nM, integer,
    see ic50_nm()) per gene. Protein-level (no pocket concept), so each (gene, compound) already
    has exactly one row -- pivoted directly, no aggregation needed. affinity_pred_value is NaN'd
    first wherever entropy_crop_pl == 0.0 (see module docstring's low-confidence gating section)."""
    df = pd.read_csv(affinity_csv, usecols=["gene_name", "compound_id", "affinity_pred_value", "entropy_crop_pl"])
    df.loc[df["entropy_crop_pl"] == 0.0, "affinity_pred_value"] = np.nan
    wide = df.pivot(index="compound_id", columns="gene_name", values="affinity_pred_value")
    return ic50_nm(wide).add_prefix(f"{prefix}_")


def load_boltz2_affinities():
    """Long-format Boltz-2 affinity results (script 75), with affinity_pred_value NaN'd wherever
    affinity_probability_binary < BOLTZ2_MIN_PROBABILITY (see module docstring's low-confidence
    gating section) -- shared by curated_boltz2_scores() and raw_boltz2_pocket_scores() so the
    gating logic lives in one place."""
    df = pd.read_csv(BOLTZ2_MTB_CSV, usecols=["pocket_name", "compound_id", "affinity_pred_value", "affinity_probability_binary"])
    df.loc[df["affinity_probability_binary"] < BOLTZ2_MIN_PROBABILITY, "affinity_pred_value"] = np.nan
    return df


def curated_boltz2_scores():
    """{"boltz2_<group>_<CAT|NONCAT>": Series indexed by compound_id} -- best (lowest predicted
    IC50, nM, integer, see ic50_nm()) Boltz-2 affinity per curated Mtb pocket group, reusing
    curated_pocket_groups() (same 12 pockets as the Uni-Dock curated columns). A group's score is
    NaN only if every one of its pockets was gated out by load_boltz2_affinities()."""
    long_df = load_boltz2_affinities()
    wide = ic50_nm(long_df.pivot(index="compound_id", columns="pocket_name", values="affinity_pred_value"))
    out = pd.DataFrame(index=wide.index)
    for col, pockets in curated_pocket_groups().items():
        out[col.replace("docking_", "boltz2_")] = wide[pockets].min(axis=1)
    return out


def curated_pocket_labels():
    """{pocket_name: "<gene>_<CAT|NONCAT>[_n]"} -- human-readable label per curated pocket, using
    each pocket's own gene_name (pheS/pheT kept separate here, unlike curated_pocket_groups()'s
    "pheST" merge -- that merge only makes sense once scores are aggregated) and site_type. A
    numeric suffix disambiguates genes with more than one pocket of the same site_type (e.g. pheS's
    two NON-CAT pockets -> pheS_NONCAT_1, pheS_NONCAT_2), in `output/selected_pockets.csv` row
    order; genes with a single pocket of that site_type get the plain "<gene>_<CAT|NONCAT>" label."""
    df = pd.read_csv(SELECTED_POCKETS_CSV)
    df["site_type_clean"] = df["site_type"].str.replace("-", "", regex=False)
    df["rank_in_group"] = df.groupby(["gene_name", "site_type_clean"]).cumcount() + 1
    df["group_size"] = df.groupby(["gene_name", "site_type_clean"])["pocket_name"].transform("size")
    labels = {}
    for _, row in df.iterrows():
        label = f"{row['gene_name']}_{row['site_type_clean']}"
        if row["group_size"] > 1:
            label += f"_{row['rank_in_group']}"
        labels[row["pocket_name"]] = label
    return labels


def raw_docking_pocket_scores():
    """{"docking_<gene>_<CAT|NONCAT>[_n]": Series indexed by compound_id} -- un-aggregated
    per-pocket mean Uni-Dock score for each of the 12 curated pockets (script 66's own columns),
    labeled via curated_pocket_labels() instead of the raw pocket/structure filename."""
    wide = pd.read_csv(CURATED_DOCKING_CSV).set_index("compound_id")
    labels = curated_pocket_labels()
    return wide[list(labels.keys())].rename(columns=labels).add_prefix("docking_")


def raw_boltz2_pocket_scores():
    """{"boltz2_<gene>_<CAT|NONCAT>[_n]": Series indexed by compound_id} -- un-aggregated
    per-pocket Boltz-2 predicted IC50 (nM, integer, see ic50_nm()) for each of the 12 curated
    pockets, same labeling as raw_docking_pocket_scores(). NaN wherever load_boltz2_affinities()
    gated that (pocket, compound) out."""
    long_df = load_boltz2_affinities()
    wide = ic50_nm(long_df.pivot(index="compound_id", columns="pocket_name", values="affinity_pred_value"))
    labels = curated_pocket_labels()
    return wide[list(labels.keys())].rename(columns=labels).add_prefix("boltz2_")


def topn_columns(per_gene_df, prefix, top_ns):
    """{"<prefix>_top<n>": Series} -- the n-th best (n-th most negative) value per row across
    per_gene_df's columns, for each n in top_ns. `na_value=np.nan` forces a plain float64 array --
    needed since the low-confidence gating can leave pandas' nullable Int64 dtype with pd.NA cells,
    and np.sort can't compare pd.NA (raises "boolean value of NA is ambiguous")."""
    sorted_vals = np.sort(per_gene_df.to_numpy(dtype="float64", na_value=np.nan), axis=1)  # ascending: index 0 = best (min) score
    out = pd.DataFrame(index=per_gene_df.index)
    for n in top_ns:
        out[f"{prefix}_top{n}"] = sorted_vals[:, n - 1]
    return out


def main():
    hits = pd.read_csv(FILTERED_HITS_CSV)[["compound_id", "smiles"]].set_index("compound_id")

    human_per_gene = load_per_gene_csv(HUMAN_GENE_MIN_CSV, "human")
    mtb_per_gene = load_per_gene_csv(MTB_GENE_MIN_CSV, "mtb")
    curated = curated_pocket_scores()
    human_topn = load_per_gene_csv(HUMAN_TOPN_CSV, "human")
    mtb_topn = load_per_gene_csv(MTB_TOPN_CSV, "mtb")
    nesso1_mtb = nesso1_scores(NESSO1_MTB_CSV, "nesso1_mtb")
    nesso1_human = nesso1_scores(NESSO1_HUMAN_CSV, "nesso1_human")
    nesso1_mtb_topn = topn_columns(nesso1_mtb, "nesso1_mtb", TOP_NS)
    nesso1_human_topn = topn_columns(nesso1_human, "nesso1_human", TOP_NS)
    boltz2_curated = curated_boltz2_scores()

    final = hits.join([human_per_gene, mtb_per_gene, curated, human_topn, mtb_topn,
                        nesso1_mtb, nesso1_human, nesso1_mtb_topn, nesso1_human_topn,
                        boltz2_curated], how="left")
    final = final.reset_index()

    out_path = os.path.join(OUTPUT_DIR, "compound_docking_summary.csv")
    final.to_csv(out_path, index=False)

    print(f"Saved {len(final):,} rows x {len(final.columns)} columns -> {out_path}")
    for name, df in [("human per-gene", human_per_gene), ("mtb per-gene", mtb_per_gene),
                      ("curated 12-pocket", curated), ("nesso1 mtb", nesso1_mtb),
                      ("nesso1 human", nesso1_human), ("boltz2 curated", boltz2_curated)]:
        n_missing = df.reindex(hits.index).isna().sum().sum()
        print(f"  {name}: {df.shape[1]} columns, {n_missing} missing cells")

    n_boltz2_gated = (pd.read_csv(BOLTZ2_MTB_CSV, usecols=["affinity_probability_binary"])["affinity_probability_binary"] < BOLTZ2_MIN_PROBABILITY).sum()
    n_nesso1_mtb_gated = (pd.read_csv(NESSO1_MTB_CSV, usecols=["entropy_crop_pl"])["entropy_crop_pl"] == 0.0).sum()
    n_nesso1_human_gated = (pd.read_csv(NESSO1_HUMAN_CSV, usecols=["entropy_crop_pl"])["entropy_crop_pl"] == 0.0).sum()
    print(f"\nLow-confidence affinity gating (see module docstring):")
    print(f"  Boltz-2: {n_boltz2_gated:,} / 13,140 raw (pocket, compound) rows NaN'd (affinity_probability_binary < {BOLTZ2_MIN_PROBABILITY})")
    print(f"  Nesso-1 mtb: {n_nesso1_mtb_gated:,} rows NaN'd (entropy_crop_pl == 0.0)")
    print(f"  Nesso-1 human: {n_nesso1_human_gated:,} rows NaN'd (entropy_crop_pl == 0.0)")

    docking_raw = raw_docking_pocket_scores()
    boltz2_raw = raw_boltz2_pocket_scores()

    condensed = hits.join([docking_raw, boltz2_raw, human_topn, mtb_topn,
                            nesso1_human_topn, nesso1_mtb_topn], how="left")
    condensed = condensed.reset_index()

    condensed_out_path = os.path.join(OUTPUT_DIR, "compound_docking_summary_condensed.csv")
    condensed.to_csv(condensed_out_path, index=False)
    print(f"\nSaved condensed table: {len(condensed):,} rows x {len(condensed.columns)} columns -> {condensed_out_path}")


if __name__ == "__main__":
    main()
