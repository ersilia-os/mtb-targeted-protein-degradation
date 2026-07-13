#!/usr/bin/env python3
"""
Generate a 3D conformer per compound in script 62's aggregated_hits.csv (RDKit ETKDGv3 + UFF,
matching scripts 44/58).

Usage:
    python 63_aggregated_conformations.py
"""
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

AGGREGATED_HITS_CSV = os.path.join(ROOT, "output", "62_aggregate_hits", "aggregated_hits.csv")
OUTPUT_DIR = os.path.join(ROOT, "output", "63_aggregated_conformations")
CONFORMATIONS_DIR = os.path.join(OUTPUT_DIR, "conformations")
os.makedirs(CONFORMATIONS_DIR, exist_ok=True)

N_WORKERS = 16


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
    df = pd.read_csv(AGGREGATED_HITS_CSV, usecols=["compound_id", "smiles"])
    IDs = df["compound_id"].tolist()
    ID_TO_SMILES = dict(zip(df["compound_id"], df["smiles"]))
    print(f"Aggregated compounds: {len(IDs):,}")

    # Skip 0-byte files too - a killed run can leave a truncated file at the final path.
    already_done = {
        f[:-len(".sdf")] for f in os.listdir(CONFORMATIONS_DIR)
        if f.endswith(".sdf") and os.path.getsize(os.path.join(CONFORMATIONS_DIR, f)) > 0
    }
    todo_ids = [ID for ID in IDs if ID not in already_done]
    print(f"{len(already_done):,} / {len(IDs):,} already have a conformation on disk "
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
    print(f"Embedding failures this run: {len(todo_ids) - processed_molecules:,}")
    print(f"Total now on disk: {len(already_done) + processed_molecules:,}")


if __name__ == "__main__":
    main()
