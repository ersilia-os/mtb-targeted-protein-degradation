from rdkit.Chem import rdFingerprintGenerator
from rdkit import Chem
from tqdm import tqdm
import pandas as pd
import numpy as np
import pickle
import os

# Read SMILES - This library corresponds to ENAMINE Diversity Library HLL-100 (Sublibrary of HLL-460)
ENAMINE_SMILES = pd.read_csv("../data/Enamine/Enamine_Hit_Locator_Library_100K_Set_plated_100160cmpds_20250629.smiles", low_memory=False, sep='\t')['SMILES'].tolist()
ENAMINE_IDs = pd.read_csv("../data/Enamine/Enamine_Hit_Locator_Library_100K_Set_plated_100160cmpds_20250629.smiles", low_memory=False, sep='\t')['Catalog ID'].tolist()

# Create OUTDIR
OUTDIR = '../output/enamine_characterization'
os.makedirs(OUTDIR, exist_ok=True)

print("SMILES to InChIKeys/IDs...")

# Create and store dict mapping IK to SMILES
IK_TO_SMILES, ID_TO_SMILES, ID_TO_IK = {}, {}, {}
for smi, id in tqdm(zip(ENAMINE_SMILES, ENAMINE_IDs)):
    mol = Chem.MolFromSmiles(smi)
    ik = Chem.MolToInchiKey(mol)
    IK_TO_SMILES[ik] = smi
    ID_TO_SMILES[id] = smi
    ID_TO_IK[id] = ik

pickle.dump(IK_TO_SMILES, open(os.path.join(OUTDIR, "IK_TO_SMI.pkl"), "wb"))
pickle.dump(ID_TO_SMILES, open(os.path.join(OUTDIR, "ID_TO_SMI.pkl"), "wb"))
pickle.dump(ID_TO_IK, open(os.path.join(OUTDIR, "ID_TO_IK.pkl"), "wb"))


print("Calculating Morgan Count Fingerprints...")

# Get all features
IDs, X = [], []
mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
for ID in tqdm(sorted(ID_TO_SMILES)):
    try:
        smi = ID_TO_SMILES[ID]
        mol = Chem.MolFromSmiles(smi)
        mfp = mfpgen.GetCountFingerprint(mol)
        X.append(mfp.ToList())
        IDs.append(ID)
    except:
        pass


# Convert to numpy array
X = np.array(X, dtype=np.int8)

# Assert that the number of IDs and MFPs are the same
assert len(IDs) == X.shape[0], f"ERROR: Number of IDs ({len(IDs)}) does not match number of MFPs ({X.shape[0]})"

# Save results
np.savez_compressed(os.path.join(OUTDIR, "X.npz"), X=X)
with open(os.path.join(OUTDIR, "IDs.txt"), "w") as f:
    f.write("\n".join(IDs))