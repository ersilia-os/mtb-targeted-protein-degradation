"""
Computes the data feeding figure_3_plot.py: multi-target docking hits from the Enamine REAL 10B
round-2 screening (output/unidock_REAL_docking_2, ~99,105 compounds x 276 pockets across all 21
genes).

A compound "hits" a gene's target if its docking score beats a cutoff in ANY of that gene's own
pockets (OR across pockets, i.e. the gene's best/min score vs. the cutoff) - the same pocket->gene
aggregation used by scripts/68_plot_results.py's _score_upsets and
scripts/39_reduce_n_hits_I.py's collapse_per_protein, generalized here from ~4 curated genes to all
21. Target identity doesn't matter for this analysis - only how many distinct genes a compound
hits.

compute_gene_min_scores(): builds the full compound x pocket score matrix (docking_utils.
build_matrix over all 276 pockets) and reduces pockets -> genes via min-per-gene-group. Not cached
to disk - only the final per-cutoff hit CSVs are kept.

compute_multi_target_hits(): for each cutoff in CUTOFFS, counts how many of the 21 genes each
compound hits, filters to n_targets >= MIN_TARGETS, and saves one CSV per cutoff (compound_id,
smiles, n_targets, targets_hit, score_<gene> for all 21 genes) - written even when empty, since a
zero-hit result at a stringent cutoff is a real finding.

compute_protein_hit_counts(): for each cutoff, counts how many compounds (out of all ~99,105) hit
each of the 21 genes - a library-wide reference stat, not itself plotted by figure_3_plot.py.

compute_selected_set_protein_hits(): the same per-gene x per-cutoff count, but restricted to the
already-selected SELECTED_SET_CUTOFF multi-target hit set (the >= MIN_TARGETS compounds saved by
compute_multi_target_hits) rather than the full library, and across PANEL_A_CUTOFFS (5 cutoffs,
loosest to strictest). Feeds figure_3_plot.py's panel a bars - this matches
notebooks/46_docking_exploration_deliverables.ipynb's precedent exactly, where the plotted `df` was
likewise a pre-filtered "selected compounds" set, not the full screened library.

Usage:
    python figure_3_calculations.py
"""
import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

import json
import tarfile

import pandas as pd

from docking_utils import LIBRARIES, build_matrix, load_scores, lookup_smiles

ROOT = os.path.join(root, "..", "..")
plots_dir = os.path.join(ROOT, "output", "plots", "figure_3")
os.makedirs(plots_dir, exist_ok=True)

CUTOFFS = [-10, -11, -12]
MIN_TARGETS = 4

# Loosest -> strictest; the population figure_3_plot.py's panel a breaks down by protein is the
# already-selected multi-target hit set at this same cutoff (see compute_selected_set_protein_hits).
PANEL_A_CUTOFFS = [-8, -9, -10, -11, -12]
SELECTED_SET_CUTOFF = -11

# Rank 1 (n_targets, then best score) of figure_3_multi_target_hits_cutoff11.csv - user-confirmed
# showcase compound for figure_3_plot.py's panel h.
SHOWCASE_COMPOUND_ID = "m_270196____8770926____26701300"


def banner(title):
    line = "=" * (len(title) + 10)
    print(line)
    print(f"==== {title} ====")
    print(line)


def pocket_to_gene(pocket, uniprot_to_gene):
    uniprot_ac = pocket.split("_model_")[0].split("_")[-1]
    return uniprot_to_gene.get(uniprot_ac)


def compute_gene_min_scores():
    banner("Loading gene mapping from figure 1's color mapping")
    with open(os.path.join(ROOT, "output", "plots", "figure_1", "color_mapping.json")) as f:
        uniprot_to_gene = json.load(f)["uniprot_to_gene"]

    banner("Loading REAL 10B round-2 docking scores (276 pockets)")
    pockets = sorted(os.listdir(LIBRARIES["REAL"]))
    pocket_map = {}
    for pocket in pockets:
        gene = pocket_to_gene(pocket, uniprot_to_gene)
        if gene is None:
            print(f"  Warning: no gene mapping for pocket '{pocket}', skipping.")
            continue
        pocket_map[pocket] = pocket
    pocket_scores = build_matrix(pocket_map, LIBRARIES["REAL"], label="Loading pocket scores")
    print(f"Pocket score matrix: {pocket_scores.shape[0]:,} compounds x {pocket_scores.shape[1]} pockets")

    banner("Reducing pockets -> genes (min score per gene, i.e. best of that gene's own pockets)")
    gene_of_pocket = {p: pocket_to_gene(p, uniprot_to_gene) for p in pocket_scores.columns}
    gene_scores = pocket_scores.T.groupby(gene_of_pocket).min().T
    print(f"Gene score matrix: {gene_scores.shape[0]:,} compounds x {gene_scores.shape[1]} genes")

    banner("Looking up SMILES")
    smiles_map = lookup_smiles(gene_scores.index.tolist(), "REAL")
    gene_scores.insert(0, "smiles", gene_scores.index.map(smiles_map))
    gene_scores.index.name = "compound_id"
    return gene_scores


def compute_multi_target_hits(gene_scores):
    gene_cols = [c for c in gene_scores.columns if c != "smiles"]

    for cutoff in CUTOFFS:
        banner(f"Cutoff <= {cutoff}")
        hit = gene_scores[gene_cols] <= cutoff
        n_targets = hit.sum(axis=1)

        tally = n_targets.value_counts().sort_index()
        for k, count in tally.items():
            print(f"  {k} target(s): {count:,} compound(s)")
        print(f"  >= {MIN_TARGETS} target(s): {(n_targets >= MIN_TARGETS).sum():,} compound(s)")

        selected = n_targets[n_targets >= MIN_TARGETS].index
        out = pd.DataFrame(index=selected)
        out["smiles"] = gene_scores.loc[selected, "smiles"]
        out["n_targets"] = n_targets.loc[selected]
        out["targets_hit"] = hit.loc[selected].apply(lambda row: ",".join(row.index[row]), axis=1)
        for gene in gene_cols:
            out[f"score_{gene}"] = gene_scores.loc[selected, gene]
        out = out.sort_values(["n_targets"] + [f"score_{g}" for g in gene_cols],
                               ascending=[False] + [True] * len(gene_cols))
        out = out.reset_index(names="compound_id")

        out_path = os.path.join(plots_dir, f"figure_3_multi_target_hits_cutoff{abs(cutoff)}.csv")
        out.to_csv(out_path, index=False)
        print(f"  Saved {len(out):,} row(s) to {out_path}")


def compute_protein_hit_counts(gene_scores):
    gene_cols = [c for c in gene_scores.columns if c != "smiles"]

    banner("Per-protein hit counts across cutoffs")
    counts = pd.DataFrame({
        f"count_cutoff{abs(cutoff)}": (gene_scores[gene_cols] <= cutoff).sum(axis=0)
        for cutoff in CUTOFFS
    })
    counts.index.name = "gene"
    counts = counts.reset_index()

    out_path = os.path.join(plots_dir, "figure_3_protein_hit_counts.csv")
    counts.to_csv(out_path, index=False)
    print(counts.to_string(index=False))
    print(f"Saved to {out_path}")


def compute_pocket_scores():
    banner("Loading P2Rank pocket probabilities (276 pockets)")
    pockets = pd.read_csv(os.path.join(ROOT, "output", "pocket_detection_data.csv"))
    pockets = pockets.sort_values("Pocket probability", ascending=False).reset_index(drop=True)

    out = pockets[["Uniprot AC", "File name", "Pocket number", "Pocket probability"]].copy()
    out.insert(0, "pocket_rank", range(1, len(out) + 1))
    out = out.rename(columns={"Pocket probability": "pocket_probability"})

    out_path = os.path.join(plots_dir, "figure_3_pocket_scores.csv")
    out.to_csv(out_path, index=False)
    print(f"Saved {len(out):,} row(s) to {out_path}")


def compute_showcase_compound_pockets():
    banner(f"Locating {SHOWCASE_COMPOUND_ID}'s actual best-scoring pocket per hit gene")

    hits_path = os.path.join(plots_dir, f"figure_3_multi_target_hits_cutoff{abs(SELECTED_SET_CUTOFF)}.csv")
    hits = pd.read_csv(hits_path)
    row = hits[hits["compound_id"] == SHOWCASE_COMPOUND_ID]
    if row.empty:
        print(f"  Error: {SHOWCASE_COMPOUND_ID} not found in {hits_path}.")
        return
    genes = row.iloc[0]["targets_hit"].split(",")

    with open(os.path.join(ROOT, "output", "plots", "figure_1", "color_mapping.json")) as f:
        uniprot_to_gene = json.load(f)["uniprot_to_gene"]
    gene_to_uniprot = {g: u for u, g in uniprot_to_gene.items()}

    pockets = sorted(os.listdir(LIBRARIES["REAL"]))
    pockets_by_gene = {}
    for pocket in pockets:
        gene = pocket_to_gene(pocket, uniprot_to_gene)
        pockets_by_gene.setdefault(gene, []).append(pocket)

    # "Catalytic Domain (ATP Binding Site)" -> CAT badge (panel h), every other curated annotation
    # -> Other badge; a pocket with no row here at all -> NA. Real InterPro-derived classification,
    # not a fabricated category - see output/pocket_detection_data_interpro.tsv's own values.
    interpro = pd.read_csv(os.path.join(ROOT, "output", "pocket_detection_data_interpro.tsv"), sep="\t")

    rows = []
    for gene in genes:
        best_pocket, best_score = None, None
        for pocket in pockets_by_gene.get(gene, []):
            scores = load_scores(os.path.join(LIBRARIES["REAL"], pocket, "report.csv"))
            if SHOWCASE_COMPOUND_ID not in scores.index:
                continue
            score = scores[SHOWCASE_COMPOUND_ID]
            if best_score is None or score < best_score:
                best_pocket, best_score = pocket, score

        structure_name, pocket_number = best_pocket.rsplit("_pocket_", 1)
        has_pose = os.path.isfile(os.path.join(LIBRARIES["REAL"], best_pocket, "docking.tar.gz"))

        uniprot_ac = gene_to_uniprot[gene]
        annotations = interpro[
            (interpro["Uniprot AC"] == uniprot_ac)
            & (interpro["File name"] == f"{structure_name}.pdb")
            & (interpro["Pocket number"] == int(pocket_number))
        ]["Interpro curated annotation"].unique()
        if len(annotations) == 0:
            interpro_categories = "NA"
        else:
            categories = {"CAT" if a == "Catalytic Domain (ATP Binding Site)" else "Other" for a in annotations}
            interpro_categories = ",".join(sorted(categories))

        print(f"  {gene}: {best_pocket} (score={best_score}, has_pose={has_pose}, "
              f"interpro={interpro_categories})")
        rows.append({
            "compound_id": SHOWCASE_COMPOUND_ID,
            "gene": gene,
            "uniprot_ac": uniprot_ac,
            "structure_name": structure_name,
            "pocket_number": int(pocket_number),
            "pocket_name": best_pocket,
            "score": best_score,
            "has_pose": has_pose,
            "interpro_categories": interpro_categories,
        })

    out = pd.DataFrame(rows)
    out_path = os.path.join(plots_dir, "figure_3_showcase_compound_pockets.csv")
    out.to_csv(out_path, index=False)
    print(f"  {out['has_pose'].sum()}/{len(out)} gene(s) have a retained docked pose.")
    print(f"Saved {len(out):,} row(s) to {out_path}")


def compute_selected_set_protein_hits():
    hits_path = os.path.join(plots_dir, f"figure_3_multi_target_hits_cutoff{abs(SELECTED_SET_CUTOFF)}.csv")
    selected = pd.read_csv(hits_path)
    gene_cols = [c.removeprefix("score_") for c in selected.columns if c.startswith("score_")]

    banner(f"Selected set (n={len(selected)}, >= {MIN_TARGETS} targets @ cutoff {SELECTED_SET_CUTOFF}) "
           "- per-protein hit counts")
    rows = []
    for gene in gene_cols:
        row = {"gene": gene}
        for cutoff in PANEL_A_CUTOFFS:
            row[f"count_cutoff{abs(cutoff)}"] = int((selected[f"score_{gene}"] <= cutoff).sum())
        rows.append(row)
    counts = pd.DataFrame(rows)

    out_path = os.path.join(plots_dir, "figure_3_selected_set_protein_hit_counts.csv")
    counts.to_csv(out_path, index=False)
    print(counts.to_string(index=False))
    print(f"Saved to {out_path}")


def main():
    gene_scores = compute_gene_min_scores()
    compute_multi_target_hits(gene_scores)
    compute_protein_hit_counts(gene_scores)
    compute_pocket_scores()
    compute_selected_set_protein_hits()
    compute_showcase_compound_pockets()


if __name__ == "__main__":
    main()
