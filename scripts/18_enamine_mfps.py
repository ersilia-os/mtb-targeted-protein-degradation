from rdkit.Chem import rdFingerprintGenerator
from rdkit import Chem
from tqdm import tqdm
import pandas as pd
import numpy as np
import pickle
import os

# Read SMILES - This library corresponds to ENAMINE Diversity Library HLL-100 (Sublibrary of HLL-460)
ENAMINE_SMILES = pd.read_csv("../data/Enamine/Enamine_Hit_Locator_Library_100K_Set_plated_100160cmpds_20250629.smiles", low_memory=False, sep='\t')['SMILES'].tolist()

# Create OUTDIR
OUTDIR = '../processed/enamine_characterization'
os.makedirs(OUTDIR, exist_ok=True)

print("SMILES to InChIKeys...")

# Create and store dict mapping IK to SMILES
IK_TO_SMILES = {}
for smi in tqdm(ENAMINE_SMILES):
    mol = Chem.MolFromSmiles(smi)
    ik = Chem.MolToInchiKey(mol)
    IK_TO_SMILES[ik] = smi

pickle.dump(IK_TO_SMILES, open(os.path.join(OUTDIR, "IK_TO_SMI.pkl"), "wb"))

print("Calculating Morgan Count Fingerprints...")

# Get all features
IKs, X = [], []
mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
for IK in tqdm(sorted(IK_TO_SMILES)):
    try:
        smi = IK_TO_SMILES[IK]
        mol = Chem.MolFromSmiles(smi)
        mfp = mfpgen.GetCountFingerprint(mol)
        X.append(mfp.ToList())
        IKs.append(IK)
    except:
        pass


# Convert to numpy array
X = np.array(X, dtype=np.int16)

# Assert that the number of IKs and MFPs are the same
assert len(IKs) == X.shape[0], f"ERROR: Number of IKs ({len(IKs)}) does not match number of MFPs ({X.shape[0]})"

# Save results
np.savez_compressed(os.path.join(OUTDIR, "X.npz"), X=X)
with open(os.path.join(OUTDIR, "IKs.txt"), "w") as f:
    f.write("\n".join(IKs))