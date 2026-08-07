import lazyqsar
import pandas as pd
import numpy as np
import pickle
import joblib
import h5py
import sys
import os

alpha = int(sys.argv[1])

# Set root directory
root = os.path.dirname(os.path.abspath(__file__))

# Load pickle
st, perc = pickle.load(open(os.path.join(root, "..", "..", "output", 'unidock_docking', 'pickle_bins_01.pkl'), "rb"))[alpha]
assert perc == "bin_01"

# Define paths
PATH_TO_MODELS = "/aloy/home/acomajuncosa/Ersilia/mtb/output/unidock_docking/models"
PATH_TO_CHUNKS = "/aloy/home/acomajuncosa/Ersilia/mtb/output/enamine_REAL_characterization/embeddings"
PATH_TO_OUTPUT = "/aloy/home/acomajuncosa/Ersilia/mtb/output/unidock_docking/inferences"

# Load model
model = joblib.load(os.path.join(PATH_TO_MODELS, st, perc, 'LQ_RF.joblib'))

# Set n_jobs to 1
for mm in model.model.models:
    mm.model_.n_jobs = 1

IDS, PROBS = [], []

# For each chunk
for chunk in range(0, 96):

    # Load chunk
    h5_path = os.path.join(PATH_TO_CHUNKS, f"enamine_REAL_chemeleon_chunk_{chunk}.h5")
    with h5py.File(h5_path, "r") as h5f:
        X, ids = h5f['X'][:], h5f['ids'][:]

    # Calculate probabilities
    probs = model.predict_proba(X)

    # Store results
    PROBS.extend(probs[:, 1])
    IDS.extend(ids.astype(str))

assert len(PROBS) == len(IDS)

# Save results
results = pd.DataFrame([])
results['ids'] = np.array(IDS, dtype=str)
results['probs'] = PROBS
results.to_csv(os.path.join(PATH_TO_OUTPUT, f"{st}_bin_01.csv.gz"), index=False, compression={'method': 'gzip','compresslevel': 9})