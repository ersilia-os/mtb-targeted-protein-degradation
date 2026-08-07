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
data = pd.read_csv(os.path.join("..", "output", "enamine_REAL_characterization", "enamine_REAL.tsv"), sep='\t')
ID_TO_SMILES = {i: j for i,j in zip(data['id'], data['smiles'])}
IDs = sorted(ID_TO_SMILES)

# Prepare input list
input_list = [[id, ID_TO_SMILES[id]]for id in IDs]

# Splitting
SPLIT_SIZE = 100000
chunks = [input_list[i:i+SPLIT_SIZE] for i in range(0, len(input_list), SPLIT_SIZE)]
id_to_split = {}
for si, chunk in enumerate(chunks):
    for cid, _ in chunk:
        id_to_split[cid] = si

# Create output directories 
OUTPATH = "../output/enamine_REAL_characterization/conformations"
os.makedirs(os.path.join(OUTPATH), exist_ok=True)

N_WORKERS = 32  # TO BE RUN IN A 32 CPU MACHINE
for si, chunk in enumerate(chunks):
    with Pool(N_WORKERS) as pool:
        part = dict(tqdm(pool.imap(generate_3d_conformer, chunk), total=len(chunk)))
    with gzip.open(os.path.join(OUTPATH, f"mols_{si:02d}.pkl.gz"), "wb") as f:
        pickle.dump(part, f, protocol=pickle.HIGHEST_PROTOCOL)

with gzip.open(os.path.join(OUTPATH, "id_to_split.tsv.gz"), "wt", encoding="utf-8") as f:
    for cid, si in id_to_split.items():
        f.write(f"{cid}\t{si}\n")