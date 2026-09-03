#!/usr/bin/env python3
"""
Held-out 3-fold stratified CV for the Gen 2 (script 37) Naive Bayes surrogate models, using the
same StratifiedKFold(n_splits=3, shuffle=True, random_state=RANDOM_SEED) protocol as script 26's
Gen 1 CV, so the two generations' AUROCs are comparable. Script 37 itself only reports an in-sample
AUROC (fit and scored on the same data); this script adds a genuinely held-out estimate without
touching script 37's production models or its existing aurocs/ output.

Fingerprints (enamine_REAL_ECFP6.h5, 9.56M x 2048 int8, ~18GB uncompressed) are read in
chunk-aligned sequential blocks, batching pockets so at most POCKET_BATCH_SIZE pockets' fingerprint
matrices are held in memory at once -- a full random-access read (per pocket, or of the whole file
into a dict as script 37 does) is either far too slow or too memory-heavy on this machine.

Usage:
    python 37b_surrogate_model_cv.py
"""
import os
import sys

import h5py
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.naive_bayes import ComplementNB
from sklearn.preprocessing import normalize
from tqdm import tqdm

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))
from default import RANDOM_SEED

ENAMINE_REAL_TSV = os.path.join(ROOT, "output", "enamine_REAL_characterization", "enamine_REAL.tsv")
ECFP6_H5 = os.path.join(ROOT, "output", "enamine_REAL_characterization", "enamine_REAL_ECFP6.h5")
DOCKING_RESULTS_DIR = os.path.join(ROOT, "output", "unidock_REAL_docking", "docking_results")

OUTPUT_DIR = os.path.join(ROOT, "output", "unidock_REAL_docking", "training_results", "aurocs_cv")
os.makedirs(OUTPUT_DIR, exist_ok=True)

NUMBER_ACTIVES = 1130  # same 1% docking-score threshold used by script 37
POCKET_BATCH_SIZE = 30  # bounds peak memory during the chunked fingerprint read below


def load_row_lookup():
    """Maps compound id (str) -> row position in the ECFP6 h5 file's fps/ids datasets."""
    data = pd.read_csv(ENAMINE_REAL_TSV, sep="\t")
    ind_to_id = data["id"].to_numpy()
    with h5py.File(ECFP6_H5, "r") as h5:
        ids_ind = h5["ids"][:]
    id_str = ind_to_id[ids_ind]
    return pd.Series(np.arange(len(id_str), dtype=np.int64), index=id_str)


def load_pocket_labels(pocket):
    """compound ids + binary labels, ranked by docking score (top NUMBER_ACTIVES = 1)."""
    scores = pd.read_csv(os.path.join(DOCKING_RESULTS_DIR, pocket, "report.csv"))
    scores = scores.sort_values("score").reset_index(drop=True)
    Y = np.zeros(len(scores), dtype=np.int64)
    Y[:NUMBER_ACTIVES] = 1
    return scores["compound"].to_numpy(), Y


def run_cv(X, Y):
    X = normalize(X, norm="l2")
    skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=RANDOM_SEED)
    aucs = []
    for train_idx, test_idx in skf.split(X, Y):
        clf = ComplementNB()
        clf.fit(X[train_idx], Y[train_idx])
        probs = clf.predict_proba(X[test_idx])[:, 1]
        aucs.append(roc_auc_score(Y[test_idx], probs))
    return aucs


def main():
    pockets = sorted(os.listdir(DOCKING_RESULTS_DIR))
    row_lookup = load_row_lookup()

    pocket_data = {}
    for pocket in pockets:
        compounds, Y = load_pocket_labels(pocket)
        row_positions = row_lookup.loc[compounds].to_numpy()
        pocket_data[pocket] = {"Y": Y, "row_positions": row_positions}

    with h5py.File(ECFP6_H5, "r") as h5:
        n_total, n_bits = h5["fps"].shape
        chunk_rows = h5["fps"].chunks[0]

        summary = []
        for batch_start in tqdm(range(0, len(pockets), POCKET_BATCH_SIZE)):
            batch = pockets[batch_start:batch_start + POCKET_BATCH_SIZE]
            X_batch = {p: np.empty((len(pocket_data[p]["row_positions"]), n_bits), dtype=np.int8)
                       for p in batch}

            for chunk_start in range(0, n_total, chunk_rows):
                chunk_end = min(chunk_start + chunk_rows, n_total)
                block = h5["fps"][chunk_start:chunk_end]
                for p in batch:
                    rows = pocket_data[p]["row_positions"]
                    in_chunk = (rows >= chunk_start) & (rows < chunk_end)
                    X_batch[p][in_chunk] = block[rows[in_chunk] - chunk_start]

            for p in batch:
                aucs = run_cv(X_batch[p], pocket_data[p]["Y"])
                with open(os.path.join(OUTPUT_DIR, f"{p}_auroc_cv.csv"), "w") as f:
                    f.write(",".join(str(round(a, 3)) for a in aucs))
                summary.append([p] + aucs + [np.mean(aucs)])
            del X_batch

    summary = pd.DataFrame(summary, columns=["pocket", "fold_1", "fold_2", "fold_3", "mean_auroc"])
    summary.to_csv(os.path.join(OUTPUT_DIR, "summary.csv"), index=False)
    print("Overall mean AUROC (3-fold CV):", round(summary["mean_auroc"].mean(), 4))


if __name__ == "__main__":
    main()
