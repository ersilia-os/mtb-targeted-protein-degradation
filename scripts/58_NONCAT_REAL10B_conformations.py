#!/usr/bin/env python3
"""
Merge script 57's per-chunk-per-pocket selective-hit CSVs, cap each pocket at MAX_PER_POCKET
(random sample - no score to rank by), drop compounds selective for more than one pocket, then
generate a 3D conformer per unique compound (RDKit ETKDGv3 + UFF, matching script 44).

Script 57's 994 chunks are 3 overlapping Enamine families, so dedup by compound_id happens
per-pocket before capping, not just at the end.

Usage:
    python 58_NONCAT_REAL10B_conformations.py
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
