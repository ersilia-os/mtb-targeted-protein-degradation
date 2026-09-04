"""
Computes all data feeding figure_3_plot.py: physicochemical properties for HL/REAL 10M/REAL 10B
(panel B), per-pocket docking-score percentiles for all three libraries (panel C), and one
best-compound PyMOL docking snapshot per gene (panel D).

HL and REAL 10M (compute_hl_and_real10m): random 100k samples, physicochemical properties via
docking_utils.compute_properties(). REAL 10B (compute_real10b): the Enamine REAL 10B library is
too large to hold locally, so it's sampled 1000 compounds per chunk (994 chunks), merged and
re-sampled down to 100k, then run through the same property calculation. Resumable: one CSV per
chunk (processed/figure_3b_REAL10B_samples/{chunk}.csv), skips chunks already sampled - each
chunk's ~90MB SMILES mapping is deleted after use to keep tmp/ bounded across 994 chunks (herbert's
root disk has little headroom).

Docking percentiles (compute_docking_percentiles): median, p1/p0.1 (top-1,000/top-100 out of
each pocket's screened compound set) per pocket, for all three libraries. Stores precomputed
summary stats rather than raw score arrays - 276 pockets x ~100k compounds per library would make
for an unnecessarily large file. REAL 10M stats use each pocket's prioritized top-100k
"active" set only (docking_utils.load_real_positive_scores), not the shared ~12,958-compound
background sample also present in that library's docking output - matching HL's and REAL 10B's
fully-screened sets.

Docking snapshots (compute_docking_snapshots): a PyMOL render (cartoon protein + stick ligand +
pocket-residue lines + H-bond dashes, following notebooks/46_docking_exploration_IIa.ipynb's
pymol_screenshot()) of the single best-scoring compound per gene. Docked 3D poses (not just
scores) are only archived (as a per-pocket docking.tar.gz) for a small hand-curated subset of
pockets - 6 of 276 for HL, 14 of 276 for REAL 10B, 0 for REAL 10M - covering 6 genes total
(alaS, aspS, ileS, lysS, pheS, pheT). Candidates are further restricted to each gene's
canonical-pocket winners (best REAL 10B p1 per spatial_cluster_id, matching figure_3_plot.py's
panel c dedup) so panel d's stars always land on one of panel c's columns - a pose-archived
structure that isn't its cluster's winner is skipped even if its own single best compound
scores well. So "best compound" here is the best score among only that gene's pose-archived,
canonical-winner (library, pocket) candidates, restricted to HL/REAL 10B (REAL 10M can never
contribute a renderable snapshot - no poses are archived for it anywhere) - not the true best
across all of a gene's pockets/libraries, most of which only ever had their score kept.

Usage:
    python figure_3_calculations.py [--max-chunks N]
"""
import argparse
import glob
import json
import os
import sys
import tarfile
import tempfile
from collections import defaultdict

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

import numpy as np
import pandas as pd

from default import RANDOM_SEED
from docking_utils import LIBRARIES, compute_properties, load_real_positive_scores, load_scores, sample_prescreened_smiles
from screening_10b_utils import download_file, list_smiles_chunks

ROOT = os.path.join(root, "..", "..")
TMP_DIR = os.path.join(ROOT, "tmp")
SAMPLES_DIR = os.path.join(ROOT, "processed", "figure_3b_REAL10B_samples")
plots_dir = os.path.join(ROOT, "output", "plots", "figure_3")
os.makedirs(TMP_DIR, exist_ok=True)
os.makedirs(SAMPLES_DIR, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)

HL_SAMPLE_N = 100_000
REAL10M_SAMPLE_N = 100_000
N_PER_CHUNK = 1000
N_FINAL = 100_000


def banner(title):
    line = "=" * (len(title) + 10)
    print(line)
    print(f"==== {title} ====")
    print(line)


def compute_hl_and_real10m():
    banner("HL (Enamine Hit Locator 100K) - random 100k sample")
    hl_smiles = sample_prescreened_smiles("DL", exclude_ids=set(), n=HL_SAMPLE_N, seed=RANDOM_SEED)
    print(f"Sampled {len(hl_smiles):,} compounds")
    hl_props = compute_properties(hl_smiles)
    output_path = os.path.join(plots_dir, "figure_3a_HL.csv")
    hl_props.to_csv(output_path)
    print(f"Saved to {output_path}")

    banner("REAL 10M (Enamine REAL 9.56M) - random 100k sample")
    real10m_smiles = sample_prescreened_smiles("REAL_ROUND1", exclude_ids=set(), n=REAL10M_SAMPLE_N, seed=RANDOM_SEED)
    print(f"Sampled {len(real10m_smiles):,} compounds")
    real10m_props = compute_properties(real10m_smiles)
    output_path = os.path.join(plots_dir, "figure_3a_REAL10M.csv")
    real10m_props.to_csv(output_path)
    print(f"Saved to {output_path}")


def sample_real10b_chunk(chunk, chunk_index):
    out_path = os.path.join(SAMPLES_DIR, f"{chunk}.csv")
    if os.path.isfile(out_path):
        return out_path

    ids_path = os.path.join(TMP_DIR, f"{chunk}_SMILES_IDs.tsv.zip")
    download_file(ids_path)
    try:
        df = pd.read_csv(ids_path, sep="\t")
        sample = df.sample(n=min(N_PER_CHUNK, len(df)), random_state=RANDOM_SEED + chunk_index)
        sample[["id", "smiles"]].to_csv(out_path, index=False)
    finally:
        os.remove(ids_path)
    return out_path


def compute_real10b(max_chunks=None):
    banner("Discovering REAL 10B chunks")
    chunks = list_smiles_chunks()
    if max_chunks is not None:
        chunks = chunks[:max_chunks]
    print(f"{len(chunks)} chunk(s) to process")

    banner("Sampling 1000 compounds per chunk")
    n_processed, n_skipped = 0, 0
    for i, chunk in enumerate(chunks):
        out_path = os.path.join(SAMPLES_DIR, f"{chunk}.csv")
        if os.path.isfile(out_path):
            n_skipped += 1
            continue
        print(f"[{i + 1}/{len(chunks)}] {chunk}")
        sample_real10b_chunk(chunk, i)
        n_processed += 1
    print(f"Processed {n_processed}, skipped (already done) {n_skipped}")

    banner("Merging per-chunk samples")
    sample_files = sorted(glob.glob(os.path.join(SAMPLES_DIR, "*.csv")))
    merged = pd.concat([pd.read_csv(f) for f in sample_files], ignore_index=True)
    print(f"Merged pool: {len(merged):,} compounds from {len(sample_files)} chunk file(s)")

    final = merged.sample(n=min(N_FINAL, len(merged)), random_state=RANDOM_SEED)
    print(f"Final random sample: {len(final):,} compounds")

    banner("Computing properties")
    smiles_dict = dict(zip(final["id"], final["smiles"]))
    props = compute_properties(smiles_dict)
    output_path = os.path.join(plots_dir, "figure_3b_REAL10B.csv")
    props.to_csv(output_path)
    print(f"Saved to {output_path}")


def pocket_percentiles(pocket, scores, library, uniprot_to_gene):
    uniprot_ac = pocket.split("_model_")[0].split("_")[-1]
    return {
        "pocket": pocket,
        "uniprot_ac": uniprot_ac,
        "gene": uniprot_to_gene.get(uniprot_ac, "unknown"),
        "library": library,
        "n": len(scores),
        "median": np.median(scores),
        "p0_1": np.percentile(scores, 0.1),
        "p1": np.percentile(scores, 1),
    }


def compute_docking_percentiles():
    banner("Loading gene mapping from figure 1's color mapping")
    with open(os.path.join(ROOT, "output", "plots", "figure_1", "color_mapping.json")) as f:
        uniprot_to_gene = json.load(f)["uniprot_to_gene"]

    rows = []

    banner("HL (output/unidock_docking) - percentiles per pocket")
    for pocket in sorted(os.listdir(LIBRARIES["DL"])):
        report_path = os.path.join(LIBRARIES["DL"], pocket, "report.csv")
        if not os.path.isfile(report_path):
            print(f"  Warning: report not found for pocket '{pocket}', skipping.")
            continue
        scores = load_scores(report_path).values
        rows.append(pocket_percentiles(pocket, scores, "HL", uniprot_to_gene))

    banner("REAL 10M (output/unidock_REAL_docking, prioritized 100k only) - percentiles per pocket")
    for pocket in sorted(os.listdir(LIBRARIES["DL"])):
        scores = load_real_positive_scores(pocket).values
        if len(scores) == 0:
            print(f"  Warning: no REAL 10M positive scores for pocket '{pocket}', skipping.")
            continue
        rows.append(pocket_percentiles(pocket, scores, "REAL 10M", uniprot_to_gene))

    banner("REAL 10B (output/unidock_REAL_docking_2) - percentiles per pocket")
    for pocket in sorted(os.listdir(LIBRARIES["REAL"])):
        report_path = os.path.join(LIBRARIES["REAL"], pocket, "report.csv")
        if not os.path.isfile(report_path):
            print(f"  Warning: report not found for pocket '{pocket}', skipping.")
            continue
        scores = load_scores(report_path).values
        rows.append(pocket_percentiles(pocket, scores, "REAL 10B", uniprot_to_gene))

    df = pd.DataFrame(rows)
    output_path = os.path.join(plots_dir, "figure_3c_docking_percentiles.csv")
    df.to_csv(output_path, index=False)
    print(f"Saved {len(df)} row(s) to {output_path}")
    for library, group in df.groupby("library"):
        print(f"  {library}: {len(group)} pocket(s)")


# Only these two libraries ever archive docked 3D poses (a per-pocket docking.tar.gz) - REAL 10M's
# docking_results/ tree has report.csv (scores) everywhere but no docking.tar.gz anywhere, so it
# can never supply a renderable pose and is excluded here.
POSE_LIBRARIES = {"HL": LIBRARIES["DL"], "REAL 10B": LIBRARIES["REAL"]}


def pockets_with_poses(library_dir):
    return sorted(
        pocket for pocket in os.listdir(library_dir)
        if os.path.isfile(os.path.join(library_dir, pocket, "docking.tar.gz"))
    )


def load_pocket_clusters():
    """structure-level pocket id -> spatial_cluster_id, from the same 6.14 A greedy
    centroid dedup figure_1_calculations.py uses for gene_to_unique_pocket_count
    (scripts/77_pocket_annotation/09_cluster_pockets.py's persisted assignments) -
    lets snapshot candidates be restricted to each canonical pocket's own winner,
    matching figure_3_plot.py's panel c dedup."""
    path = os.path.join(ROOT, "output", "77_pocket_annotation", "pocket_clusters.csv")
    df = pd.read_csv(path)
    df["pocket_id"] = df["File name"].str.replace(".pdb", "", regex=False) + "_pocket_" + df["Pocket number"].astype(str)
    return df.set_index("pocket_id")["spatial_cluster_id"].to_dict()


def canonical_pocket_winners(pocket_to_cluster):
    """{(gene, pocket)} of every canonical pocket's REAL 10B p1 winner - the same
    per-(gene, spatial_cluster_id) selection figure_3_plot.py's gene_docking_stats()
    applies to panel c, recomputed here so panel d only draws snapshots from structures
    that actually get a column in panel c."""
    docking_percentiles = pd.read_csv(os.path.join(plots_dir, "figure_3c_docking_percentiles.csv"))
    docking_percentiles["spatial_cluster_id"] = docking_percentiles["pocket"].map(pocket_to_cluster)
    real10b = docking_percentiles[docking_percentiles["library"] == "REAL 10B"]
    winners = real10b.sort_values("p1").drop_duplicates(["gene", "spatial_cluster_id"])
    return set(zip(winners["gene"], winners["pocket"]))


def best_compound_candidates(uniprot_to_gene, winner_pockets):
    """{gene: [{library, pocket, compound, score}, ...]} - one entry per (library, pose-archived
    pocket) belonging to that gene, restricted to that gene's canonical-pocket winners
    (winner_pockets), each already reduced to that pocket's own best-scoring compound.
    compute_docking_snapshots() flattens these across genes and takes the global top-N by
    score, so a gene can supply more than one of the N snapshots. Also requires the
    structure's own .pdbqt to be resolvable (find_structure_pdbqt) - a winner can have
    docked poses (docking.tar.gz) archived without the protein .pdbqt ever being copied
    into this library's tree (e.g. alphafold3_P9WFV3_model_0, ileS), which would otherwise
    only surface as a pymol_snapshot() crash deep in the rendering step."""
    candidates = defaultdict(list)
    for library, library_dir in POSE_LIBRARIES.items():
        for pocket in pockets_with_poses(library_dir):
            uniprot_ac = pocket.split("_model_")[0].split("_")[-1]
            gene = uniprot_to_gene.get(uniprot_ac)
            if gene is None or (gene, pocket) not in winner_pockets:
                continue
            try:
                find_structure_pdbqt(library_dir, pocket)
            except FileNotFoundError:
                print(f"  Warning: no .pdbqt found for pocket '{pocket}' ({library}), skipping as a snapshot candidate.")
                continue
            scores = load_scores(os.path.join(library_dir, pocket, "report.csv"))
            candidates[gene].append({
                "gene": gene, "library": library, "pocket": pocket,
                "compound": scores.idxmin(), "score": scores.min(),
            })
    return candidates


def find_structure_pdbqt(library_dir, pocket):
    """Usually a pocket's own folder holds a copy of its structure's .pdbqt, but some sibling
    pockets sharing the same structure only keep ONE shared copy, stored under a different
    sibling's folder (e.g. alphafold3_P9WFU9_model_1's .pdbqt lives only under
    .../model_1_pocket_2/, not .../model_1_pocket_1/, even though pocket_1 also uses it) - so
    search every sibling pocket of the same structure, not just this one's own folder."""
    structure = pocket.rsplit("_pocket_", 1)[0]
    direct = os.path.join(library_dir, pocket, structure + ".pdbqt")
    if os.path.isfile(direct):
        return direct
    for sibling in os.listdir(library_dir):
        if sibling.startswith(structure + "_pocket_"):
            candidate = os.path.join(library_dir, sibling, structure + ".pdbqt")
            if os.path.isfile(candidate):
                return candidate
    raise FileNotFoundError(f"{structure}.pdbqt not found in any sibling pocket folder under {library_dir}")


def pymol_snapshot(pocket, compound, library_dir, out_path):
    """One PyMOL render (cartoon protein + stick ligand + pocket-residue lines + H-bond dashes)
    of a single docked compound - following notebooks/46_docking_exploration_IIa.ipynb's
    pymol_screenshot(), trimmed to just the PNG (no pandamap 2D-interaction diagram, no separate
    .pdb dump)."""
    import pymol
    from pymol import cmd

    pocket_dir = os.path.join(library_dir, pocket)
    structure_path = find_structure_pdbqt(library_dir, pocket)

    with tarfile.open(os.path.join(pocket_dir, "docking.tar.gz"), "r|gz") as tf:
        for member in tf:
            if member.name == f"docking/{compound}_out.sdf":
                data = tf.extractfile(member).read()
                break
        else:
            raise FileNotFoundError(f"{compound}_out.sdf not found in {pocket_dir}/docking.tar.gz")

    tmp_fd, tmp_sdf = tempfile.mkstemp(suffix=".sdf")
    os.close(tmp_fd)
    try:
        with open(tmp_sdf, "w") as f:
            f.write(data.decode("utf-8", errors="replace"))

        pymol.finish_launching(["pymol", "-cq"])
        cmd.reinitialize()
        cmd.bg_color("white")
        cmd.set("ray_opaque_background", 0)

        cmd.load(structure_path, "structure")
        cmd.util.cbag("structure")
        cmd.set_color("structure_col", [0x8D / 255, 0xC7 / 255, 0xFA / 255])
        cmd.color("structure_col", "structure and elem C")

        lig = "lig"
        cmd.load(tmp_sdf, lig)
        cmd.util.cbag(lig)
        cmd.set_color("ligC_orange", [0xF5 / 255, 0xA6 / 255, 0x3A / 255])
        cmd.color("ligC_orange", f"{lig} and elem C")

        cmd.hide("everything", "structure")
        cmd.show("cartoon", "structure")
        cmd.set("cartoon_transparency", 0.4, "structure")
        cmd.select("pocket_atoms", f"(structure within 6 of {lig}) and not solvent")
        cmd.show("lines", "pocket_atoms")
        cmd.show("sticks", lig)
        cmd.orient(lig)
        cmd.zoom(lig, buffer=6)
        cmd.h_add("structure")
        cmd.h_add(lig)

        pairs = cmd.find_pairs(f"{lig} and donor", "pocket_atoms and acceptor", mode=1, cutoff=3.5, angle=45)
        pairs += cmd.find_pairs(f"{lig} and acceptor", "pocket_atoms and donor", mode=1, cutoff=3.5, angle=45)
        cmd.delete("hbonds")
        for a1, a2 in list(dict.fromkeys(pairs)):
            cmd.distance("hbonds", a1, a2)
        cmd.color("yellow", "hbonds")
        cmd.hide("labels", "hbonds")
        cmd.set("dash_width", 2)
        cmd.hide("everything", "elem H")

        cmd.set("ray_shadows", 0)
        cmd.set("ray_trace_mode", 1)
        cmd.png(out_path, width=1200, height=1200, dpi=300, ray=1)
    finally:
        os.remove(tmp_sdf)


# Matches panel D's 7-column layout in figure_3_plot.py - a fixed count, not one-per-gene, since
# only 6 genes have any pose-archived candidates at all (see POSE_LIBRARIES above) and repeats
# are wanted rather than leaving a slot empty.
N_SNAPSHOTS = 7


def select_minimizing_repeats(candidates_by_gene, n):
    """Picks n entries with the minimum possible number of repeated genes: every gene's own best
    candidate is taken before any gene's second-best, every gene's second-best before any third,
    etc. - a gene only repeats once all other genes (with a candidate left) have already gotten
    that many picks. Within each round, ties across genes are broken by score."""
    remaining = {gene: sorted(entries, key=lambda e: e["score"]) for gene, entries in candidates_by_gene.items()}
    selected = []
    round_index = 0
    while len(selected) < n:
        round_pool = sorted(
            (entries[round_index] for entries in remaining.values() if len(entries) > round_index),
            key=lambda e: e["score"],
        )
        if not round_pool:
            break  # no gene has any candidate left at this depth
        selected.extend(round_pool[:n - len(selected)])
        round_index += 1
    return selected


def compute_docking_snapshots():
    banner("Loading gene mapping from figure 1's color mapping")
    with open(os.path.join(ROOT, "output", "plots", "figure_1", "color_mapping.json")) as f:
        uniprot_to_gene = json.load(f)["uniprot_to_gene"]

    banner("Finding pose-archived compound candidates (HL/REAL 10B only, canonical-pocket winners only)")
    winner_pockets = canonical_pocket_winners(load_pocket_clusters())
    candidates_by_gene = best_compound_candidates(uniprot_to_gene, winner_pockets)

    # Frozen to the original 6-gene pose-archived pool (2026-08-25, user decision): gatA/glyS/
    # trpS/tyrS picked up pose archives on 2026-08-10/08-19, after this panel was last built
    # (2026-08-07) - excluding them here keeps this revision scoped to the requested lysS -> aspS
    # swap only. Drop this filter in a later, separate update to fold the new genes in.
    ORIGINAL_POSE_GENES = {"alaS", "aspS", "ileS", "lysS", "pheS", "pheT"}
    candidates_by_gene = {g: v for g, v in candidates_by_gene.items() if g in ORIGINAL_POSE_GENES}

    all_candidates = [e for entries in candidates_by_gene.values() for e in entries]
    for e in sorted(all_candidates, key=lambda e: e["score"]):
        print(f"  {e['gene']:<8} {e['library']:<8} {e['pocket']:<45} {e['compound']:<15} {e['score']:.3f}")

    banner(f"Rendering top-{N_SNAPSHOTS} PyMOL snapshots (minimum repeated genes, best score first)")
    top = sorted(select_minimizing_repeats(candidates_by_gene, N_SNAPSHOTS), key=lambda e: e["score"])

    output_dir = os.path.join(plots_dir, "docking_snapshots")
    os.makedirs(output_dir, exist_ok=True)
    for old_file in glob.glob(os.path.join(output_dir, "*.png")):
        os.remove(old_file)

    index_rows = []
    for rank, entry in enumerate(top, start=1):
        filename = f"{rank:02d}_{entry['gene']}.png"
        out_path = os.path.join(output_dir, filename)
        print(f"#{rank} {entry['gene']}: {entry['compound']} in {entry['pocket']} "
              f"({entry['library']}, score {entry['score']:.3f}) -> {out_path}")
        pymol_snapshot(entry["pocket"], entry["compound"], POSE_LIBRARIES[entry["library"]], out_path)
        index_rows.append({**entry, "rank": rank, "filename": filename})

    index_path = os.path.join(output_dir, "index.csv")
    pd.DataFrame(index_rows).to_csv(index_path, index=False)
    print(f"Saved {len(top)} snapshot(s) and {index_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-chunks", type=int, default=None,
                         help="Limit REAL 10B sampling to the first N chunks (for a quick dry run).")
    args = parser.parse_args()

    compute_hl_and_real10m()
    compute_real10b(max_chunks=args.max_chunks)
    compute_docking_percentiles()
    compute_docking_snapshots()


if __name__ == "__main__":
    main()
