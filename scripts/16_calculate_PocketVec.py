import pickle
import os
import numpy as np
import scipy.stats as ss

path = "../processed/pocketvec_RUN/pocketvec_POST"
order = pickle.load(open("../processed/pocketvec_RUN/TOP128_rDock_LLM/ALL/all.pkl", "rb"))

def rank_fp(fp):
    return ss.rankdata(fp, method='max')


fps_raw, fps_rank = {}, {}

for st in sorted(os.listdir(os.path.join(path))):
    
    raw_scores = {}
    # Take scores
    for scores in sorted(os.listdir(os.path.join(path, st, "rDock_results_LLM"))):
        if "scores" in scores:
            with open(os.path.join(path, st, "rDock_results_LLM", scores), "r") as f:
                for l in [i[:-1] for i in f.readlines()]:
                    name = l.split("\t")[0]
                    aff = float(l.split("\t")[1])
                    raw_scores[name] = aff

    # Store in dict
    fp = np.array([raw_scores[mol] for mol in order])
    fps_raw[st] = fp                

    # Additional fps
    for st in sorted(fps_raw):
        fps_rank[st] = rank_fp(fps_raw[st])

    for st in sorted(fps_raw):
        for c in range(len(fps_raw[st])):
            if fps_raw[st][c] > 0 and fps_raw[st][c] < 50:
                fps_rank[st][c] = len(fps_rank[st]) + 1
            elif fps_raw[st][c] > 50 and fps_raw[st][c] < 100:
                fps_rank[st][c] = len(fps_rank[st]) + 2
            elif fps_raw[st][c] > 100:
                fps_rank[st][c] = len(fps_rank[st]) + 3

pickle.dump(fps_rank, open("../processed/pocketvec_RUN/fps_rank.pkl", "wb"))
