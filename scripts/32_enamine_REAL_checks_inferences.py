from tqdm import tqdm
import pandas as pd
import random
import pickle
import numpy as np
import gzip
import os

random.seed(42)
root = os.path.dirname(os.path.abspath(__file__))

PATH_TO_INFERENCE_OUTPUTS = os.path.join("/aloy/home/acomajuncosa/Ersilia/mtb/output/unidock_docking/inferences_outputs")
print("Checking inference logs...")
for out in tqdm(sorted(os.listdir(PATH_TO_INFERENCE_OUTPUTS))):
    out = open(os.path.join(PATH_TO_INFERENCE_OUTPUTS, out), "r").readlines()
    out = [i for i in out if "warning" in i.lower() or "error" in i.lower()]
    if len(out) != 0:
        print(f"ERROR: CHECK FILE {out}")
print("If no message appeared, all is good")
print("\n\n")

# Loading failed compounds
FAILED = set(pd.read_csv(os.path.join(root, "..", "output", "enamine_REAL_characterization", "failed_REAL.csv"), header=None)[0])

# Load compound info
data = pd.read_csv(os.path.join("..", "output", "enamine_REAL_characterization", "enamine_REAL.tsv"), sep='\t')

# Loading successful compounds
SUCCESS = set([i for i in data['id'] if i not in FAILED])
SUCCESS = sorted(SUCCESS)

print(f"Loaded {len(FAILED)} failed compounds")
print(f"Loaded {len(SUCCESS)} successful compounds")

# Dump list of successful compounds
PATH_TO_RESULTS = os.path.join("..", "output", "unidock_docking", "inference_probs")
pickle.dump(SUCCESS, open(os.path.join(PATH_TO_RESULTS, "success_mols.pkl"), "wb"))

print("Reading inference results")
PATH_TO_INFERENCES = "/aloy/home/acomajuncosa/Ersilia/mtb/output/unidock_docking/inferences"
os.makedirs(PATH_TO_RESULTS, exist_ok=True)

# For each pocket
for pocket in sorted(os.listdir(PATH_TO_INFERENCES)):

    print(f"Processing {pocket} ...")

    # Load results
    inf = pd.read_csv(os.path.join(PATH_TO_INFERENCES, pocket))
    inf = {i: j for i,j in zip(inf['ids'], inf['probs'])}

    # Get ordered probabilities
    probs = [inf[i] for i in SUCCESS]

    # Store
    np.savez_compressed(os.path.join(PATH_TO_RESULTS, pocket.replace(".csv.gz", '')), np.array(probs).astype('float32'))