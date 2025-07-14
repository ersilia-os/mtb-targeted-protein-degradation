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
IDs = open("../processed/enamine_characterization/IDs.txt", "r").readlines()
IDs = [i.strip() for i in IDs]
ID_TO_SMILES = pickle.load(open("../processed/enamine_characterization/ID_TO_SMI.pkl", "rb"))
IK_TO_SMILES = pickle.load(open("../processed/enamine_characterization/IK_TO_SMI.pkl", "rb"))
ID_TO_IK = pickle.load(open("../processed/enamine_characterization/ID_TO_IK.pkl", "rb"))


# Given (NAME, SMILES) create 3D conformer
def generate_3d_conformer(input_list):
    ID, IK, smi = input_list
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42 
    status = AllChem.EmbedMolecule(mol, params=params)
    if status != 0:
        return None
    AllChem.UFFOptimizeMolecule(mol)
    return [ID, IK, mol]


# Prepare input list
input_list = [[id, ID_TO_IK[id], ID_TO_SMILES[id]]for id in IDs]

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
        ID, IK, mol = mol
        mol.SetProp("_Name", ID)
        mol.SetProp("_ID", ID)
        mol.SetProp("_IK", IK)
        molblock = Chem.MolToV2KMolBlock(mol)
        out_path = os.path.join(OUTPATH, ID + ".sdf")
        with open(out_path, "w") as f:
            f.write(molblock)
            f.write(f">  <_ID>\n{ID}\n\n")
            f.write(f">  <_IK>\n{IK}\n\n")
            f.write("$$$$\n")
            processed_molecules += 1

print("Creating a single SDF file...")

# Write results in a single SDF file
processed_molecules = 0
with gzip.open("../processed/enamine_characterization/conformations.sdf.gz", "wt") as gzfile:
    for mol in mols:
        if mol is not None:
            ID, IK, mol = mol
            mol.SetProp("_Name", ID)
            mol.SetProp("_ID", ID)
            mol.SetProp("_IK", IK)
            molblock = Chem.MolToV2KMolBlock(mol)
            gzfile.write(molblock)
            gzfile.write(f">  <_ID>\n{ID}\n\n")
            gzfile.write(f">  <_IK>\n{IK}\n\n")
            gzfile.write("$$$$\n")
            processed_molecules += 1


# Create zip
shutil.make_archive(OUTPATH, 'gztar', OUTPATH)

# Remove directory and keep only the compressed file
shutil.rmtree(OUTPATH)

# Print some info
print(f"Initial number of molecules: {len(IDs)}")
print(f"Processed number of molecules: {processed_molecules}")