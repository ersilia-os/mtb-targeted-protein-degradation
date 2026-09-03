"""
Shared pipeline for turning a per-pocket Uni-Dock docking-results directory (one
report.csv of raw scores per pocket, same compound set/order across pockets) into a 2D
t-SNE embedding: rank-transform each pocket's scores (scipy.stats.rankdata; a more negative
score = a better pose = a lower rank), stack into a (n_pockets x n_compounds) rank matrix,
reduce via PCA, then t-SNE on the PCA embedding. Used by scripts/plots/figure_1_5_plot.py's
panel a, for both the 100,154-compound Hit Locator Library and the ~99k-compound Enamine
REAL 10B selection (script 46).
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

PCA_N_COMPONENTS_DEFAULT = 256


def load_docking_rank_matrix(docking_results_dir):
    keys = sorted(os.listdir(docking_results_dir))
    compound_order = None
    ranks = []
    for k in keys:
        report = pd.read_csv(os.path.join(docking_results_dir, k, "report.csv"))
        if compound_order is None:
            compound_order = report["compound"].tolist()
        scores = report.set_index("compound").reindex(compound_order)["score"].values
        assert not np.isnan(scores).any(), f"{k}: missing compound(s) vs. reference order"
        ranks.append(rankdata(scores))
    return keys, np.array(ranks, dtype=np.float32)


def compute_docking_tsne_embedding(docking_results_dir, random_seed, pca_n_components=PCA_N_COMPONENTS_DEFAULT):
    keys, X = load_docking_rank_matrix(docking_results_dir)
    X_pca = PCA(n_components=pca_n_components, random_state=random_seed).fit_transform(X)
    coords = TSNE(n_components=2, perplexity=30, random_state=random_seed, init="random").fit_transform(X_pca)
    return keys, coords
