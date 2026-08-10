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

# User-confirmed per-gene overrides for compute_gene_summary_stats' otherwise-NaN novelty /
# experimental_tractability placeholders. Genes not listed here stay NaN.
NOVELTY_OVERRIDES = {"ileS": False}
EXPERIMENTAL_TRACTABILITY_OVERRIDES = {"glyS": False}


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


def compute_gene_summary_stats(gene_scores):
    banner("Per-gene summary stats (21 proteins): max P2Rank prob + best HL Lib / REAL 10B docking scores")
    with open(os.path.join(ROOT, "output", "plots", "figure_1", "color_mapping.json")) as f:
        uniprot_to_gene = json.load(f)["uniprot_to_gene"]

    pockets = pd.read_csv(os.path.join(ROOT, "output", "pocket_detection_data.csv"))
    pockets["gene"] = pockets["Uniprot AC"].map(uniprot_to_gene)
    max_prob = pockets.groupby("gene")["Pocket probability"].max()

    # HL Lib = Enamine Hit Locator Library 100k, LIBRARIES["DL"] - see scripts/56, 61, 68,
    # figure_2_calculations.py's own "HL" naming for the same library.
    dl_pockets = sorted(os.listdir(LIBRARIES["DL"]))
    pocket_map = {p: p for p in dl_pockets if pocket_to_gene(p, uniprot_to_gene) is not None}
    dl_scores = build_matrix(pocket_map, LIBRARIES["DL"], label="Loading HL Lib pocket scores")
    gene_of_pocket = {p: pocket_to_gene(p, uniprot_to_gene) for p in dl_scores.columns}
    best_hl_score = dl_scores.T.groupby(gene_of_pocket).min().T.min(axis=0)

    # REAL 10B (round-2) best score per gene, reusing the gene_scores matrix compute_gene_min_scores
    # already built (compound x gene, min-per-gene-group over that gene's own pockets) rather than
    # rebuilding the 99,105-compound x 276-pocket matrix from scratch.
    real_gene_cols = [c for c in gene_scores.columns if c != "smiles"]
    best_real_score = gene_scores[real_gene_cols].min(axis=0)

    # Cross-polypharmacology score: raw pairwise hit-overlap count with every OTHER gene, at
    # SELECTED_SET_CUTOFF - i.e. figure_3_plot.py's plot_circos_overlap own node_strength (row sum
    # + column sum of the overlap matrix, user-confirmed via AskUserQuestion as "raw overlap with
    # others, same as the circos plot"), persisted here instead of being discarded after rendering.
    # SELECTED_SET_CUTOFF must stay equal to figure_3_plot.py's CIRCOS_CUTOFF for these to match.
    # Requires compute_multi_target_hits to have already written this cutoff's hits CSV.
    hits_path = os.path.join(plots_dir, f"figure_3_multi_target_hits_cutoff{abs(SELECTED_SET_CUTOFF)}.csv")
    selected = pd.read_csv(hits_path)
    circos_gene_cols = [c.removeprefix("score_") for c in selected.columns if c.startswith("score_")]
    hit_sets = {g: set(selected.loc[selected[f"score_{g}"] <= SELECTED_SET_CUTOFF, "compound_id"])
                for g in circos_gene_cols}
    overlap = pd.DataFrame(0, index=circos_gene_cols, columns=circos_gene_cols, dtype=int)
    for g1 in circos_gene_cols:
        for g2 in circos_gene_cols:
            if g1 != g2:
                overlap.loc[g1, g2] = len(hit_sets[g1] & hit_sets[g2])
    cross_polypharm_score = overlap.sum(axis=0) + overlap.sum(axis=1)

    result = {
        gene: {
            "max_p2rank_prob": float(max_prob[gene]),
            "best_hl_docking_score": float(best_hl_score[gene]),
            "best_real10b_docking_score": float(best_real_score[gene]),
            "cross_polypharmacology_score": int(cross_polypharm_score[gene]),
            # Not computed yet - NaN unless overridden above.
            "novelty": NOVELTY_OVERRIDES.get(gene, float("nan")),
            "experimental_tractability": EXPERIMENTAL_TRACTABILITY_OVERRIDES.get(gene, float("nan")),
        }
        for gene in max_prob.index
    }
    print(json.dumps(result, indent=2))

    out_path = os.path.join(plots_dir, "figure_3_gene_summary_stats.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Saved {len(result)} gene(s) to {out_path}")


# A pocket counts as "catalytic" for panel f's Catalytic row at catalytic_confidence >= 3 (either
# strong direct-PDB or strong AlphaFill ligand evidence for the Catalytic Domain (ATP/ligase) label
# - see scripts/77_pocket_annotation/10_assemble_final_table.py's confidence scale, 0-4), vs. the
# rest (weak/no evidence or no catalytic label at all, confidence 0-2) - a threshold the user
# specified directly, not chosen here.
CATALYTIC_CONFIDENCE_MIN = 3


# The other 3 domain rows (tRNA binding, Editing, Anticodon binding) have no ligand-evidence
# confidence score like catalytic_confidence - only a binary "this InterPro domain label is present
# for this pocket" (curated_labels). Red there means label present AND catalytic_confidence <
# CATALYTIC_CONFIDENCE_MIN - mutually exclusive with the Catalytic row, since a pocket can carry
# both a catalytic and a non-catalytic label at once (e.g. "Catalytic Domain (ATP/ligase)|Editing
# Domain") and 3 pockets do exactly that with catalytic_confidence >= 3, per the user's explicit
# request not to double-count a pocket already flagged red in the Catalytic row.
CURATED_LABEL_COLUMNS = {
    "is_trna_binding": "tRNA Binding Domain",
    "is_editing": "Editing Domain",
    "is_anticodon_binding": "Anticodon Binding Domain",
}


def compute_pocket_scores():
    banner("Loading P2Rank pocket probabilities (276 pockets)")
    pockets = pd.read_csv(os.path.join(ROOT, "output", "pocket_detection_data.csv"))
    pockets = pockets.sort_values("Pocket probability", ascending=False).reset_index(drop=True)

    interpro = pd.read_csv(os.path.join(ROOT, "output", "77_pocket_annotation", "pocket_detection_interpro_updated.csv"),
                            keep_default_na=False)
    pockets = pockets.merge(
        interpro[["Uniprot AC", "File name", "Pocket number", "catalytic_confidence", "curated_labels"]],
        on=["Uniprot AC", "File name", "Pocket number"], how="left",
    )

    out = pockets[["Uniprot AC", "File name", "Pocket number", "Pocket probability", "catalytic_confidence"]].copy()
    out.insert(0, "pocket_rank", range(1, len(out) + 1))
    out = out.rename(columns={"Pocket probability": "pocket_probability"})
    out["is_catalytic"] = out["catalytic_confidence"] >= CATALYTIC_CONFIDENCE_MIN
    has_label = {
        col: pockets["curated_labels"].apply(lambda labels, label=label: label in labels.split("|"))
        for col, label in CURATED_LABEL_COLUMNS.items()
    }
    for col, mask in has_label.items():
        out[col] = mask & ~out["is_catalytic"]

    out_path = os.path.join(plots_dir, "figure_3_pocket_scores.csv")
    out.to_csv(out_path, index=False)
    print(f"Saved {len(out):,} row(s) to {out_path}")


def _compute_pockets_for_compound(compound_id, genes, gene_to_uniprot, pockets_by_gene, interpro):
    """Per-gene best-VERIFIED-retained-pose pocket walk, used by
    compute_top_avg_score_compounds_pockets. For each gene: rank that gene's own pockets by
    this compound's score (best/most negative first - docked poses were only ever retained for a
    curated subset of pockets, so the true best-scoring one often isn't among them), then walk down
    until a pocket with an actually-verified retained pose is found (tar exists AND contains this
    compound's SDF - a tar's mere presence doesn't guarantee this compound is in it), else fall back
    to the true best pocket with has_pose=False. The score reported is THIS pocket's own score
    (matching what's actually depicted), not the gene's true best - showing a fallback pocket's
    ligand pose next to a different pocket's score would be inconsistent with what the image shows.
    Returns a list of row dicts."""
    rows = []
    for gene in genes:
        ranked = []
        for pocket in pockets_by_gene.get(gene, []):
            scores = load_scores(os.path.join(LIBRARIES["REAL"], pocket, "report.csv"))
            if compound_id not in scores.index:
                continue
            ranked.append((scores[compound_id], pocket))
        ranked.sort()

        chosen_pocket, chosen_score, has_pose = ranked[0][1], ranked[0][0], False
        for score, pocket in ranked:
            tar_path = os.path.join(LIBRARIES["REAL"], pocket, "docking.tar.gz")
            if not os.path.isfile(tar_path):
                continue
            with tarfile.open(tar_path, "r|gz") as tf:
                if any(m.name == f"docking/{compound_id}_out.sdf" for m in tf):
                    chosen_pocket, chosen_score, has_pose = pocket, score, True
                    break

        best_pocket, best_score = chosen_pocket, chosen_score
        is_true_best_pocket = best_pocket == ranked[0][1]
        structure_name, pocket_number = best_pocket.rsplit("_pocket_", 1)

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
              f"is_true_best_pocket={is_true_best_pocket}, interpro={interpro_categories})")
        rows.append({
            "compound_id": compound_id,
            "gene": gene,
            "uniprot_ac": uniprot_ac,
            "structure_name": structure_name,
            "pocket_number": int(pocket_number),
            "pocket_name": best_pocket,
            "score": best_score,
            "has_pose": has_pose,
            "is_true_best_pocket": is_true_best_pocket,
            "interpro_categories": interpro_categories,
        })
    return rows


def _load_gene_uniprot_maps():
    with open(os.path.join(ROOT, "output", "plots", "figure_1", "color_mapping.json")) as f:
        uniprot_to_gene = json.load(f)["uniprot_to_gene"]
    gene_to_uniprot = {g: u for u, g in uniprot_to_gene.items()}
    return uniprot_to_gene, gene_to_uniprot


def _load_pockets_by_gene(uniprot_to_gene):
    pockets = sorted(os.listdir(LIBRARIES["REAL"]))
    pockets_by_gene = {}
    for pocket in pockets:
        gene = pocket_to_gene(pocket, uniprot_to_gene)
        pockets_by_gene.setdefault(gene, []).append(pocket)
    return pockets_by_gene


# Best-average-docking-score pair from the cutoff-12 multi-target hit tier (each hits exactly 4
# genes there) - user-confirmed via AskUserQuestion, for figure_3_plot.py's panels c/d (one row per
# compound: 2D structure + 4 docking-pose renders).
TOP_AVG_SCORE_COMPOUND_IDS = ["s_271570____28264988____28567424", "s_51____13974142____77337"]


def compute_top_avg_score_compounds_pockets():
    banner("Locating TOP_AVG_SCORE_COMPOUND_IDS's actual best-scoring pocket per hit gene")

    hits_path = os.path.join(plots_dir, "figure_3_multi_target_hits_cutoff12.csv")
    hits = pd.read_csv(hits_path)

    uniprot_to_gene, gene_to_uniprot = _load_gene_uniprot_maps()
    pockets_by_gene = _load_pockets_by_gene(uniprot_to_gene)
    interpro = pd.read_csv(os.path.join(ROOT, "output", "pocket_detection_data_interpro.tsv"), sep="\t")

    all_rows = []
    for compound_rank, compound_id in enumerate(TOP_AVG_SCORE_COMPOUND_IDS, start=1):
        match = hits[hits["compound_id"] == compound_id]
        if match.empty:
            print(f"  Error: {compound_id} not found in {hits_path}.")
            continue
        genes = match.iloc[0]["targets_hit"].split(",")
        print(f"  Compound {compound_rank} ({compound_id}): {len(genes)} hit gene(s) - {genes}")
        rows = _compute_pockets_for_compound(compound_id, genes, gene_to_uniprot, pockets_by_gene, interpro)
        for r in rows:
            r["compound_rank"] = compound_rank
        all_rows.extend(rows)

    out = pd.DataFrame(all_rows)
    out_path = os.path.join(plots_dir, "figure_3_top_avg_score_compounds_pockets.csv")
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
    compute_gene_summary_stats(gene_scores)
    compute_pocket_scores()
    compute_selected_set_protein_hits()
    compute_top_avg_score_compounds_pockets()


if __name__ == "__main__":
    main()
