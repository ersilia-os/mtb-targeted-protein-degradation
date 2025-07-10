from rdkit.Chem import rdFingerprintGenerator
from multiprocessing import Pool
from rdkit.Chem import AllChem
from rdkit import Chem
from tqdm import tqdm
import pandas as pd
import numpy as np
import pickle
import gzip
import os

# Load compound info
IKs = open("../processed/enamine_characterization/IKs.txt", "r").readlines()
IKs = [i.strip() for i in IKs]
IK_TO_SMILES = pickle.load(open("../processed/enamine_characterization/IK_TO_SMI.pkl", "rb"))

# Given (NAME, SMILES) create 3D conformer
def generate_3d_conformer(input_list):
    IK, smi = input_list
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42 
    status = AllChem.EmbedMolecule(mol, params=params)
    if status != 0:
        return None
    AllChem.UFFOptimizeMolecule(mol)
    mol.SetProp("_Name", IK)
    return mol


# Prepare input list
input_list = [[ik, IK_TO_SMILES[ik]]for ik in IKs]

# Run conformer generation - parallelized
N_WORKERS = 12
with Pool(N_WORKERS) as pool:
    mols = list(tqdm(pool.imap(generate_3d_conformer, input_list), total=len(input_list)))

# Write results
processed_molecules = 0
with gzip.open("../processed/enamine_characterization/conformations.sdf.gz", "wt") as gzfile:
    writer = Chem.SDWriter(gzfile)
    for mol in mols:
        if mol is not None:
            writer.write(mol)
            processed_molecules += 1
    writer.close()

print(f"Initial number of molecules: {len(IKs)}")
print(f"Processed number of molecules: {processed_molecules}")