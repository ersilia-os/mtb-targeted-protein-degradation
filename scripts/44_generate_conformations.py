from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.rdmolfiles import MolToV2KMolBlock
from multiprocessing import Pool
from rdkit.Chem import AllChem
from rdkit import Chem
from tqdm import tqdm
import pandas as pd
import numpy as np
import pickle
import shutil
import string
import gzip
import os

root = "/home/acomajuncosa/Documents_GPU/mtb-targeted-protein-degradation"

# Define path to last docking round
PATH_TO_OUTPUT = os.path.join(root, "processed", "unidock_REAL_docking_2")
os.makedirs(PATH_TO_OUTPUT, exist_ok=True)

# Load data from 10B screening
df = pd.read_csv(os.path.join(root, "processed", "unidock_REAL_docking", "inference_10B", "clustered_compounds.csv"))

# Load compound info
IDs = df['id'].tolist()
ID_TO_SMILES = {i: j for i,j in zip(df['id'], df['smiles'])}

# Given (NAME, SMILES) create 3D conformer
def generate_3d_conformer(input_list):
    ID, smi = input_list
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42 
    status = AllChem.EmbedMolecule(mol, params=params)
    if status != 0:
        return None
    AllChem.UFFOptimizeMolecule(mol)
    return [ID, mol]

# Prepare input list
input_list = [[id, ID_TO_SMILES[id]]for id in IDs]

# Create output directories 
OUTPATH = os.path.join(PATH_TO_OUTPUT, "conformations")
os.makedirs(os.path.join(OUTPATH), exist_ok=True)

# Run conformer generation - parallelized
N_WORKERS = 16
with Pool(N_WORKERS) as pool:
    mols = list(tqdm(pool.imap(generate_3d_conformer, input_list), total=len(input_list)))


print("Creating individual files...")

# Write results
processed_molecules = 0
for mol in mols:
    if mol is not None:
        ID, mol = mol
        mol.SetProp("_Name", ID)
        mol.SetProp("_ID", ID)
        molblock = Chem.MolToV2KMolBlock(mol)
        out_path = os.path.join(OUTPATH, ID + ".sdf")
        with open(out_path, "w") as f:
            f.write(molblock)
            f.write(f">  <_ID>\n{ID}\n\n")
            f.write("$$$$\n")
            processed_molecules += 1

# Print some info
print(f"Initial number of molecules: {len(IDs)}")
print(f"Processed number of molecules: {processed_molecules}")