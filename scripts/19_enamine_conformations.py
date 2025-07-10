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
    return [IK, mol]


# Prepare input list
input_list = [[ik, IK_TO_SMILES[ik]]for ik in IKs]

# Get alphabet
# alphabet_upper = string.ascii_uppercase

# Create output directories 
OUTPATH = "../processed/enamine_characterization/conformations"
os.makedirs(os.path.join(OUTPATH), exist_ok=True)
# for i in alphabet_upper:
#     for j in alphabet_upper:
#         os.makedirs(os.path.join(OUTPATH, i, j), exist_ok=True)


# Run conformer generation - parallelized
N_WORKERS = 12
with Pool(N_WORKERS) as pool:
    mols = list(tqdm(pool.imap(generate_3d_conformer, input_list), total=len(input_list)))

print("Creating individual files...")

# Write results
processed_molecules = 0
for mol in mols:
    if mol is not None:
        IK, mol = mol
        mol.SetProp("_Name", IK)
        molblock = Chem.MolToV2KMolBlock(mol)
        out_path = os.path.join(OUTPATH, IK + ".sdf")
        with open(out_path, "w") as f:
            f.write(molblock)
            f.write("$$$$\n")
            processed_molecules += 1

print("Creating a single SDF file...")

# Write results in a single SDF file
processed_molecules = 0
with gzip.open("../processed/enamine_characterization/conformations.sdf.gz", "wt") as gzfile:
    for mol in mols:
        if mol is not None:
            IK, mol = mol
            mol.SetProp("_Name", IK)
            molblock = Chem.MolToV2KMolBlock(mol)
            gzfile.write(molblock)
            gzfile.write("$$$$\n")
            processed_molecules += 1


# Create zip
shutil.make_archive(OUTPATH, 'gztar', OUTPATH)

# Remove directory and keep only the compressed file
shutil.rmtree(OUTPATH)

# Print some info
print(f"Initial number of molecules: {len(IKs)}")
print(f"Processed number of molecules: {processed_molecules}")