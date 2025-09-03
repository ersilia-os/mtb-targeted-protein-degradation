from multiprocessing import Pool
from rdkit.Chem import AllChem
from rdkit import Chem
import pandas as pd
import numpy as np
from tqdm import tqdm
import pickle
import gzip
import os

# Given (NAME, SMILES) create 3D conformer
def generate_3d_conformer(input_list):
    ID, smi = input_list
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42 
    status = AllChem.EmbedMolecule(mol, params=params)
    if status != 0:
        return [ID, None]
    AllChem.UFFOptimizeMolecule(mol)
    return [ID, mol]

# Load compound info
data = pd.read_csv(os.path.join("..", "processed", "enamine_REAL_characterization", "enamine_REAL.tsv"), sep='\t')
ID_TO_SMILES = {i: j for i,j in zip(data['id'], data['smiles'])}
IDs = sorted(ID_TO_SMILES)

# Prepare input list
input_list = [[id, ID_TO_SMILES[id]]for id in IDs]

# Create output directories 
OUTPATH = "../processed/enamine_REAL_characterization"
os.makedirs(os.path.join(OUTPATH), exist_ok=True)

# Run conformer generation - parallelized
N_WORKERS = 32
with Pool(N_WORKERS) as pool:
    mols = list(tqdm(pool.imap(generate_3d_conformer, input_list), total=len(input_list)))

# Convert to dict
mols = dict(mols)

# Dump results
with gzip.open(os.path.join(OUTPATH, "mols.pkl.gz"), "wb") as f:
    pickle.dump(mols, f, protocol=pickle.HIGHEST_PROTOCOL)