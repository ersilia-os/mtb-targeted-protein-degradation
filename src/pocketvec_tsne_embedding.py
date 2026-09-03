"""
PocketVec t-SNE embedding (128-dim rDock rank fingerprints, script 16's fps_rank.pkl),
used by scripts/plots/figure_1_5_plot.py's panel a.
"""
import os

import numpy as np
import pandas as pd
from sklearn.manifold import TSNE

OUTLIER_MAX_PENALTY_RANKS = 80  # same threshold as notebooks/16_PocketVec_analyses.ipynb


def compute_pocketvec_embedding(output_dir, random_seed):
    fps = pd.read_pickle(os.path.join(output_dir, "pocketvec_RUN", "fps_rank.pkl"))
    n_penalty_ranks = {p: int(np.sum(fps[p] > len(fps[p]))) for p in fps}
    keys = sorted(p for p in fps if n_penalty_ranks[p] < OUTLIER_MAX_PENALTY_RANKS)
    X = np.array([fps[p] for p in keys])
    coords = TSNE(n_components=2, perplexity=30, random_state=random_seed, init="random").fit_transform(X)
    return keys, coords
