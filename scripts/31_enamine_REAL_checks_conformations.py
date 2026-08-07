from tqdm import tqdm
import pandas as pd
import random
import pickle
import gzip
import os

random.seed(42)
root = os.path.dirname(os.path.abspath(__file__))

### FIRST STEP - CHECK WHICH CONFORMATIONS FAILED

# Define paths
PATH_TO_CONFORMATIONS = os.path.join(root, "..", "..", "..", "Documents_GPU", "mtb-targeted-protein-degradation", "output", "enamine_REAL_characterization", 'conformations')
MOLS_FAILED, MOLS_SUCCESS = [], []

# For each split of 100k compounds
for split in tqdm(range(0, 96)):
    with gzip.open(os.path.join(PATH_TO_CONFORMATIONS, f"mols_{split:02d}.pkl.gz"), "rb") as f:
        mols = pickle.load(f)
        failed = [i for i in mols if mols[i] == None]
        success = [i for i in mols if mols[i] != None]
        MOLS_FAILED.extend(failed)
        MOLS_SUCCESS.extend(success)

print(f"Total number of molecules: {len(MOLS_SUCCESS)}")
# assert len(mols) == 9557695
print(f"Total number of failed molecules: {len(MOLS_FAILED)}")

# Save failed molecules
PATH_TO_FAILED = os.path.join(root, "..", "output", "enamine_REAL_characterization", "failed_REAL.csv")
with open(PATH_TO_FAILED, "w") as out:
    out.write("\n".join(sorted(MOLS_FAILED)))

### SECOND STEP - GENERATE A RANDOM SAMPLE OF 10K REAL COMPOUNDS

# Get IDs
IDs = sorted(MOLS_SUCCESS)

# Sample 25k compounds
sample = random.sample(IDs, 25000)
print(f"Total number of sampled compounds: {len(sample)}")

# Save results
PATH_TO_SAMPLE = os.path.join(root, "..", "output", "enamine_REAL_characterization", "random_sample_REAL.csv")
with open(PATH_TO_SAMPLE, "w") as out:
    out.write("\n".join(sorted(sample)))