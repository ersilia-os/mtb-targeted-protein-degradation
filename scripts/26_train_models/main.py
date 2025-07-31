from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
import pandas as pd
import lazyqsar
import numpy as np
import pickle
import joblib
import sys
import os

alpha = int(sys.argv[1])

# Set root directory
root = os.path.dirname(os.path.abspath(__file__))

# Load pickle
st, perc = pickle.load(open(os.path.join(root, "..", "..", "processed", 'unidock_docking', 'pickle.pkl'), "rb"))[alpha]

# Get input data
PATH_TO_REPORTS = os.path.join(root, "..", "..", "processed", "unidock_docking", "binarized_reports")
PATH_TO_OUTPUT = os.path.join(root, "..", "..", "processed", "unidock_docking", "models", st, perc)
PATH_TO_EMBEDDINGS = os.path.join(root, "..", "..", "processed", "enamine_characterization")
os.makedirs(PATH_TO_OUTPUT, exist_ok=True)

# Load compounds and activities
report = pd.read_csv(os.path.join(PATH_TO_REPORTS, f"report_bin_{st}.csv"))
compounds = report["compound"].tolist()
Y = np.array(report[perc].tolist())

# Load ids and embeddings
ids = open(os.path.join(PATH_TO_EMBEDDINGS, "IDs_CheMeleon.txt")).read().splitlines()
embeddings = np.load(os.path.join(PATH_TO_EMBEDDINGS, "X_CheMeleon.npz"))['X']

# Mapping id to embedding
id_to_embedding = {i: j for i,j in zip(ids, embeddings)}

# Creating matrix
X = np.array([id_to_embedding[i] for i in compounds])

# Stratified 3-fold CV
kf = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)
aucs = []

sys.stderr.write(str(st) + " -- " + str(perc) + "\n\n")
sys.stderr.flush()

# For each CV
for train_idx, test_idx in kf.split(X, Y):

    # Train test split
    X_train, X_test = X[train_idx], X[test_idx]
    Y_train, Y_test = Y[train_idx], Y[test_idx]

    # Train model only on training set
    model = lazyqsar.LazyBinaryClassifier(model_type="random_forest", pca=False, min_seen_across_partitions=1, 
                                          num_trials=20, base_num_splits=1, max_samples=10000)

    # Fit and predict
    model.fit(X_train, Y_train)
    probs = model.predict_proba(X_test)[:, 1]
    
    # Evaluate
    auc = roc_auc_score(Y_test, probs)
    aucs.append(auc)

# Save results
with open(os.path.join(PATH_TO_OUTPUT, f"AUROCs.csv"), "w") as f:
    f.write(",".join([str(round(i, 3)) for i in aucs]))

# Train model on all data
model = lazyqsar.LazyBinaryClassifier(model_type="random_forest", pca=False, min_seen_across_partitions=1, 
                                        num_trials=20, base_num_splits=1, max_samples=10000)
model.fit(X, Y)

# Save model
joblib.dump(model, os.path.join(PATH_TO_OUTPUT, "LQ_RF.joblib"))