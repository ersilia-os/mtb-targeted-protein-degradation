#!/usr/bin/env python3
"""
Identify Enamine REAL 10B compounds that are selective for each of the 7 NON-CAT reference pockets
(catalytic-adjacent/interface sites) independently: top 1% predicted probability (ind_1.npz
threshold, from ersilia-os/gcadda4tb-enamine-real-screening) for that pocket, but NOT top 1% for
any pocket belonging to a different protein. Same selectivity logic as the CAT analysis, applied
per-pocket to the NON-CAT set instead - see src/screening_10b_utils.py for the shared helpers
(download_file, get_pocket_to_ac, load_ind1) and their docstrings for why only top-1% SET
MEMBERSHIP is available (continuous probabilities are never persisted by that external repo).

Processes the first N_CHUNKS chunks (of 994) available under PATH_TO_RESULTS. All 994 are local as
of this run (N_CHUNKS=994 below), but the skip-if-already-done logic (see "Output/resumability"
below) means this stays safe to point at a smaller N_CHUNKS or a partially-downloaded
PATH_TO_RESULTS in the future - only newly available chunks get processed.

The 7 targets are read from output/selected_pockets.csv (site_type == "NON-CAT"), excluding the
dimer-interface pocket 7K98_pocket_6 - same filter scripts/56_NONCAT_top100_selection.py applies:
  pheS: alphafold3_P9WFU3_model_2_pocket_2
  pheT: alphafold2_P9WFU1_model_0_pocket_1
  aspS: alphafold3_P9WFW3_model_3_pocket_2, chai1_P9WFW3_model_0_pocket_2
  lysS: alphafold3_P9WFU9_model_1_pocket_2
  alaS: swissmodel_P9WFW7_model_0_pocket_2, alphafold2_P9WFW7_model_0_pocket_3

Background per target = all pockets NOT belonging to the target's own protein or (if applicable)
its heterodimer partner - pheS<->pheT, PARTNER_AC_OF, applied symmetrically in both directions
(same convention as scripts/53_CAT_selective.py, confirmed with the user to apply both ways).
aspS/lysS/alaS have no documented partner gene.

Per chunk, the tar is opened once and all 276 pockets' ind_1/thr are read into memory up front, then
each of the 7 targets is resolved via in-memory set operations - no need to reopen/reread the tar
per target since all 7 share the same chunk.

Compound ids/SMILES are looked up via each chunk's {chunk}_SMILES_IDs.tsv.zip mapping file,
downloaded on-demand into this repo's tmp/ via screening_10b_utils.download_file. Each mapping is
~90MB - keeping all of them around would balloon to >100GB across ~994 chunks, so once a chunk's
mapping has been used to write all 7 pockets' output for that chunk, the cached .tsv.zip is deleted
immediately (tmp/ only ever holds the one chunk currently being processed, not an accumulating
per-chunk cache).

Output/resumability: one CSV per pocket PER CHUNK
(output/57_NONCAT_selective_10B/{gene}_{pocket}/{chunk}.csv), not one combined file - the
~86MB-per-chunk SMILES-mapping download dominates runtime, so a chunk whose 7 output files already
all exist is skipped entirely (no re-download, no re-read of its tar). This makes the job resumable
and incrementally scalable: bump N_CHUNKS and re-run, only the newly-added chunks get processed.
The end-of-run "Totals" summary is therefore cumulative - it reads every {chunk}.csv present on
disk per pocket, not just chunks processed in this particular invocation.

Usage:
    python 57_NONCAT_selective_10B.py
"""
import glob
import os
import sys
import tarfile
import time

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from screening_10b_utils import download_file, get_pocket_to_ac, load_ind1

PATH_TO_RESULTS = os.path.expanduser("~/github/gcadda4tb-enamine-real-screening/output")
SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")
TMP_DIR = os.path.join(ROOT, "tmp")
OUTPUT_DIR = os.path.join(ROOT, "output", "57_NONCAT_selective_10B")
os.makedirs(TMP_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

N_CHUNKS = 994  # all chunks currently available locally
DIMER_POCKET = "7K98_pocket_6"
PARTNER_AC_OF = {"P9WFU3": "P9WFU1", "P9WFU1": "P9WFU3"}  # pheS <-> pheT heterodimer, symmetric


def load_noncat_targets():
    """[(gene_name, pocket_name), ...] for the 7 NON-CAT pockets, excluding the dimer pocket."""
    df = pd.read_csv(SELECTED_POCKETS_CSV)
    noncat = df[df["site_type"] == "NON-CAT"]
    noncat = noncat[noncat["pocket_name"] != DIMER_POCKET]
    return list(zip(noncat["gene_name"], noncat["pocket_name"]))


def get_smiles_ids(chunk):
    ids_path = os.path.join(TMP_DIR, f"{chunk}_SMILES_IDs.tsv.zip")
    if not os.path.exists(ids_path):
        download_file(ids_path)
    return pd.read_csv(ids_path, sep="\t"), ids_path


def pocket_dir(gene, pocket):
    d = os.path.join(OUTPUT_DIR, f"{gene}_{pocket}")
    os.makedirs(d, exist_ok=True)
    return d


def main():
    t_start = time.perf_counter()

    targets = load_noncat_targets()
    print(f"NON-CAT targets (excl. dimer pocket {DIMER_POCKET}): {len(targets)}")
    for gene, pocket in targets:
        print(f"  {gene}: {pocket}")

    pocket_to_ac = get_pocket_to_ac()

    all_chunks = sorted(f[:-len(".tar")] for f in os.listdir(PATH_TO_RESULTS) if f.endswith(".tar"))
    chunks = all_chunks[:N_CHUNKS]
    print(f"\nChunks available locally: {len(all_chunks)} - considering first {len(chunks)}")

    n_processed, n_skipped = 0, 0

    for chunk in chunks:
        out_paths = {pocket: os.path.join(pocket_dir(gene, pocket), f"{chunk}.csv") for gene, pocket in targets}
        if all(os.path.isfile(p) for p in out_paths.values()):
            print(f"\nChunk {chunk}: already done, skipping")
            n_skipped += 1
            continue

        print(f"\n=== Chunk: {chunk} ===")
        n_processed += 1
        tar_path = os.path.join(PATH_TO_RESULTS, f"{chunk}.tar")
        pocket_data = {}
        with tarfile.open(tar_path) as tar:
            for pocket in pocket_to_ac:
                idx, thr = load_ind1(tar, chunk, pocket)
                pocket_data[pocket] = (set(idx.tolist()), thr)

        smiles_ids, ids_path = get_smiles_ids(chunk)

        for gene, pocket in targets:
            target_ac = pocket_to_ac[pocket]
            partner_ac = PARTNER_AC_OF.get(target_ac)
            excluded_acs = {target_ac} | ({partner_ac} if partner_ac else set())
            background_pockets = [p for p, ac in pocket_to_ac.items() if ac not in excluded_acs]

            target_idx, target_thr = pocket_data[pocket]
            background_idx = set()
            for bg_pocket in background_pockets:
                background_idx.update(pocket_data[bg_pocket][0])

            selective_idx = sorted(target_idx - background_idx)

            print(f"  --- {gene}: {pocket} ---")
            print(f"    Background pockets: {len(background_pockets)}")
            print(f"    Target top-1% hits: {len(target_idx)} (thr={target_thr:.6f})")
            print(f"    Background union: {len(background_idx)}")
            print(f"    Selective hits: {len(selective_idx)}")

            idx_array = pd.Index(selective_idx).to_numpy()
            sub = smiles_ids.iloc[idx_array].reset_index()
            assert (sub["index"].to_numpy() == idx_array).all()
            sub["chunk"] = chunk
            sub["ind1_threshold"] = target_thr
            sub["gene_name"] = gene
            sub["pocket_name"] = pocket
            sub = sub.rename(columns={"id": "compound_id"})
            sub = sub[["chunk", "index", "compound_id", "smiles", "ind1_threshold", "gene_name", "pocket_name"]]
            sub.to_csv(out_paths[pocket], index=False)

        # All 7 pockets are now written for this chunk - the SMILES mapping has served its
        # purpose and would otherwise accumulate to >100GB across all chunks (each ~90MB).
        freed_mb = os.path.getsize(ids_path) / 1e6
        os.remove(ids_path)
        print(f"  Removed cached SMILES mapping for {chunk} ({freed_mb:.0f} MB freed)")

    print(f"\n=== Cumulative totals (all chunk files on disk per pocket) ===")
    for gene, pocket in targets:
        chunk_files = sorted(glob.glob(os.path.join(pocket_dir(gene, pocket), "*.csv")))
        dfs = [pd.read_csv(f) for f in chunk_files]
        n_selective = sum(len(d) for d in dfs)
        print(f"  {gene} ({pocket}): {len(chunk_files)} chunk file(s), {n_selective:,} cumulative selective hits")

    elapsed = time.perf_counter() - t_start
    print(f"\nThis run: {n_processed} chunk(s) processed, {n_skipped} skipped (already done)")
    if n_processed:
        print(f"Total time: {elapsed:.1f} s ({elapsed / n_processed:.1f} s/processed chunk)")
    else:
        print(f"Total time: {elapsed:.1f} s")


if __name__ == "__main__":
    main()
