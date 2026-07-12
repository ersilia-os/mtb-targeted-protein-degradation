#!/usr/bin/env python3
"""
Merge scripts/57_NONCAT_selective_10B.py's per-chunk-per-pocket selective-hit CSVs into one
combined table, then generate 3D conformations for the selected compounds.

Script 57's data has no continuous score to rank by - only top-1% SET MEMBERSHIP survives from the
external screening repo (ersilia-os/gcadda4tb-enamine-real-screening never persists per-compound
probabilities, see src/screening_10b_utils.py) - so any per-pocket downsampling here is random, not
"top-N by score" (unlike the old docking-score-based scripts/XX_merge_scores_select_hits.py, which
this script otherwise mirrors the structure of).

Steps:
  1. Gather - for each of the 7 NON-CAT pockets, concatenate every {chunk}.csv under
     output/57_NONCAT_selective_10B/{gene}_{pocket}/, then dedup by compound_id WITHIN that
     pocket's own set - script 57's 994 chunks are 3 overlapping families (LeadLike/
     NaturalProducts/Sample), so the same compound can appear more than once for one pocket.
  2. Per-pocket cap - random sample down to MAX_PER_POCKET if exceeded (seeded via RANDOM_SEED).
  3. Merge - concatenate the 7 capped tables, then drop any compound selective for more than one
     distinct pocket entirely (not considered for any pocket - possible via the pheS<->pheT
     partner exemption, or aspS/alaS's own 2 pockets each). Saved to merged_selective_hits.csv.
  4. Dedup + conformers - drop_duplicates(subset="compound_id") (a no-op after step 3's exclusion,
     kept as a safety net), then one 3D conformer per unique compound (RDKit ETKDGv3 + UFF,
     matching scripts/44_generate_conformations.py), parallelized.

Usage:
    python 58_generate_conformations_noncat_selective_10B.py
"""
import glob
import os
import sys
from multiprocessing import Pool

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from default import RANDOM_SEED

SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")
SCRIPT57_OUTPUT_DIR = os.path.join(ROOT, "output", "57_NONCAT_selective_10B")
OUTPUT_DIR = os.path.join(ROOT, "output", "58_generate_conformations_noncat_selective_10B")
CONFORMATIONS_DIR = os.path.join(OUTPUT_DIR, "conformations")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(CONFORMATIONS_DIR, exist_ok=True)

DIMER_POCKET = "7K98_pocket_6"
MAX_PER_POCKET = 100_000
N_WORKERS = 16


def load_noncat_targets():
    """[(gene_name, pocket_name), ...] for the 7 NON-CAT pockets, excluding the dimer pocket."""
    df = pd.read_csv(SELECTED_POCKETS_CSV)
    noncat = df[df["site_type"] == "NON-CAT"]
    noncat = noncat[noncat["pocket_name"] != DIMER_POCKET]
    return list(zip(noncat["gene_name"], noncat["pocket_name"]))


def gather_pocket_hits(gene, pocket):
    chunk_files = sorted(glob.glob(os.path.join(SCRIPT57_OUTPUT_DIR, f"{gene}_{pocket}", "*.csv")))
    dfs = [pd.read_csv(f) for f in chunk_files]
    if not dfs:
        return pd.DataFrame()
    df = pd.concat(dfs, ignore_index=True)
    return df.drop_duplicates(subset="compound_id")


def generate_3d_conformer(input_list):
    ID, smi = input_list
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = RANDOM_SEED
    status = AllChem.EmbedMolecule(mol, params=params)
    if status != 0:
        return None
    AllChem.UFFOptimizeMolecule(mol)
    return [ID, mol]


def main():
    targets = load_noncat_targets()
    print(f"NON-CAT targets (excl. dimer pocket {DIMER_POCKET}): {len(targets)}")

    capped_dfs = []
    for gene, pocket in targets:
        df = gather_pocket_hits(gene, pocket)
        n_gathered = len(df)
        if n_gathered > MAX_PER_POCKET:
            df = df.sample(n=MAX_PER_POCKET, random_state=RANDOM_SEED)
            print(f"  {gene} ({pocket}): gathered {n_gathered:,} -> capped to {len(df):,}")
        else:
            print(f"  {gene} ({pocket}): gathered {n_gathered:,} (under cap, kept all)")
        capped_dfs.append(df)

    merged = pd.concat(capped_dfs, ignore_index=True)

    pockets_per_compound = merged.groupby("compound_id")["pocket_name"].nunique()
    multi_pocket_ids = pockets_per_compound[pockets_per_compound > 1].index
    n_excluded_rows = merged["compound_id"].isin(multi_pocket_ids).sum()
    print(f"\nExcluding {len(multi_pocket_ids):,} compounds selective for more than one distinct "
          f"pocket ({n_excluded_rows:,} rows removed) - not considered for any pocket")
    merged = merged[~merged["compound_id"].isin(multi_pocket_ids)].reset_index(drop=True)

    merged_path = os.path.join(OUTPUT_DIR, "merged_selective_hits.csv")
    merged.to_csv(merged_path, index=False)
    print(f"\nMerged (after exclusion): {len(merged):,} rows across {len(targets)} pockets")
    print(f"Unique compounds: {merged['compound_id'].nunique():,} (one pocket each, by construction)")
    print(f"Saved: {merged_path}")

    unique_df = merged.drop_duplicates(subset="compound_id")
    IDs = unique_df["compound_id"].tolist()
    ID_TO_SMILES = dict(zip(unique_df["compound_id"], unique_df["smiles"]))

    # Skip 0-byte files too - a killed run can leave a truncated file at the final path.
    already_done = {
        f[:-len(".sdf")] for f in os.listdir(CONFORMATIONS_DIR)
        if f.endswith(".sdf") and os.path.getsize(os.path.join(CONFORMATIONS_DIR, f)) > 0
    }
    todo_ids = [ID for ID in IDs if ID not in already_done]
    print(f"\n{len(already_done):,} / {len(IDs):,} already have a conformation on disk "
          f"(from a prior run) - generating the remaining {len(todo_ids):,}...")

    input_list = [[ID, ID_TO_SMILES[ID]] for ID in todo_ids]
    processed_molecules = 0
    with Pool(N_WORKERS) as pool:
        for result in tqdm(pool.imap_unordered(generate_3d_conformer, input_list), total=len(input_list)):
            if result is None:
                continue
            ID, mol = result
            mol.SetProp("_Name", ID)
            mol.SetProp("_ID", ID)
            molblock = Chem.MolToV2KMolBlock(mol)
            out_path = os.path.join(CONFORMATIONS_DIR, ID + ".sdf")
            part_path = out_path + ".part"
            with open(part_path, "w") as f:
                f.write(molblock)
                f.write(f">  <_ID>\n{ID}\n\n")
                f.write("$$$$\n")
            os.replace(part_path, out_path)
            processed_molecules += 1

    print(f"\nInitial number of unique molecules: {len(IDs):,}")
    print(f"Already done (skipped): {len(already_done):,}")
    print(f"Newly processed this run: {processed_molecules:,}")
    print(f"Total now on disk: {len(already_done) + processed_molecules:,}")


if __name__ == "__main__":
    main()
