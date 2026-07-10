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
  1. Gather - for each of the 7 NON-CAT pockets (output/selected_pockets.csv, excl. the dimer
     pocket 7K98_pocket_6, same filter as script 57's load_noncat_targets), concatenate every
     {chunk}.csv under output/57_NONCAT_selective_10B/{gene}_{pocket}/.
  2. Per-pocket cap - if a pocket's gathered row count exceeds MAX_PER_POCKET, randomly sample
     down to MAX_PER_POCKET (seeded via RANDOM_SEED); otherwise keep all rows.
  3. Merge - concatenate the 7 (possibly capped) per-pocket tables as-is, WITHOUT deduplicating
     (a compound selective for >1 pocket - possible via the pheS<->pheT partner exemption, or
     aspS/alaS each having 2 of their own NON-CAT pockets - appears once per pocket, same
     convention as the old merge script). Saved to merged_selective_hits.csv.
  4. Dedup + conformers - drop_duplicates(subset="compound_id") (dedup happens here, not at merge
     time), then one 3D conformer per unique compound: RDKit ETKDGv3 + UFF, same approach as
     scripts/44_generate_conformations.py / scripts/XX_generate_conformations_noncat_top10k.py,
     parallelized via multiprocessing.Pool. One .sdf per successfully-embedded compound.

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
CONFORMATIONS_DIR = os.path.join(ROOT, "processed", "58_generate_conformations_noncat_selective_10B", "conformations")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(CONFORMATIONS_DIR, exist_ok=True)

DIMER_POCKET = "7K98_pocket_6"
MAX_PER_POCKET = 10_000
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
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


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
    merged_path = os.path.join(OUTPUT_DIR, "merged_selective_hits.csv")
    merged.to_csv(merged_path, index=False)
    print(f"\nMerged: {len(merged):,} rows across {len(targets)} pockets")
    print(f"Saved: {merged_path}")

    n_unique = merged["compound_id"].nunique()
    pockets_per_compound = merged["compound_id"].value_counts()
    n_duplicated = (pockets_per_compound > 1).sum()
    print(f"\nUnique compounds: {n_unique:,} / {len(merged):,} rows "
          f"({n_duplicated:,} compounds hit more than one pocket)")
    print("Compounds by number of pockets hit:")
    print(pockets_per_compound.value_counts().sort_index().rename_axis("n_pockets").rename("n_compounds"))

    unique_df = merged.drop_duplicates(subset="compound_id")
    IDs = unique_df["compound_id"].tolist()
    ID_TO_SMILES = dict(zip(unique_df["compound_id"], unique_df["smiles"]))
    print(f"\nGenerating 3D conformations for {len(IDs):,} unique compounds...")

    input_list = [[ID, ID_TO_SMILES[ID]] for ID in IDs]
    with Pool(N_WORKERS) as pool:
        mols = list(tqdm(pool.imap(generate_3d_conformer, input_list), total=len(input_list)))

    print("Creating individual files...")
    processed_molecules = 0
    for mol in mols:
        if mol is not None:
            ID, mol = mol
            mol.SetProp("_Name", ID)
            mol.SetProp("_ID", ID)
            molblock = Chem.MolToV2KMolBlock(mol)
            out_path = os.path.join(CONFORMATIONS_DIR, ID + ".sdf")
            with open(out_path, "w") as f:
                f.write(molblock)
                f.write(f">  <_ID>\n{ID}\n\n")
                f.write("$$$$\n")
                processed_molecules += 1

    print(f"\nInitial number of molecules: {len(IDs)}")
    print(f"Processed number of molecules: {processed_molecules}")


if __name__ == "__main__":
    main()
