import os
import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.join(root, "..", "..")
output_dir = os.path.join(repo_root, "output", "77_pocket_annotation")
os.makedirs(output_dir, exist_ok=True)

# Same empirical same-pocket distance cutoff used by scripts/plots/figure_1_calculations.py
# (derived in notebooks/08_coherence_detected_pockets.ipynb). This script previously cut
# hierarchical clustering at the single biggest gap in the merge-distance sequence, but that
# heuristic is not robust: the final merge is almost always the largest distance overall, so
# the biggest gap tends to sit near the top of the dendrogram and inflate the threshold,
# merging spatially distinct pockets together (verified on pheT: auto-threshold ~47 A, merging
# pocket-instances up to 28 A apart). Replaced with figure 1's fixed-threshold greedy dedup.
POCKET_DEDUP_DISTANCE_THRESHOLD = 6.14

if __name__ == "__main__":
    df = pd.read_csv(os.path.join(repo_root, "output", "pocket_detection_data.csv"))

    rows = []
    for ac, sub in df.groupby("Uniprot AC"):
        sub = sub.sort_values("Pocket score", ascending=False).reset_index(drop=True)
        coords = np.array([[float(x) for x in c.split()] for c in sub["Pocket centroid coordinate (x y z)"]])

        # greedy dedup: accept a pocket as a new distinct site only if farther than the
        # threshold from every already-accepted centroid (same algorithm/order as figure 1)
        accepted_centroids = []
        for c in coords:
            if all(np.linalg.norm(c - a) > POCKET_DEDUP_DISTANCE_THRESHOLD for a in accepted_centroids):
                accepted_centroids.append(c)
        accepted_centroids = np.array(accepted_centroids)

        # assign every pocket-instance (accepted or not) to its nearest accepted site
        assignments = [int(np.argmin(np.linalg.norm(accepted_centroids - c, axis=1))) for c in coords]

        # renumber sites by size (largest = 1) for a stable, readable id
        sizes = pd.Series(assignments).value_counts()
        remap = {old: new + 1 for new, old in enumerate(sizes.index)}
        cluster_ids = [remap[a] for a in assignments]
        cluster_sizes = [sizes[a] for a in assignments]

        for i, (_, row) in enumerate(sub.iterrows()):
            rows.append({
                "Uniprot AC": ac, "File name": row["File name"], "Pocket number": row["Pocket number"],
                "spatial_cluster_id": cluster_ids[i], "n_models_in_cluster": cluster_sizes[i],
            })

    out_df = pd.DataFrame(rows)
    out_path = os.path.join(output_dir, "pocket_clusters.csv")
    out_df.to_csv(out_path, index=False)
    print("Saved {} pocket cluster assignments -> {}".format(len(out_df), out_path))

    print("\nDistinct pocket counts by protein (should match figure_1_calculations.py's gene_to_unique_pocket_count):")
    print(out_df.groupby("Uniprot AC")["spatial_cluster_id"].nunique().sort_index())

    pheS = out_df[out_df["Uniprot AC"] == "P9WFU3"]
    print("\npheS cluster assignments (sanity check):")
    print(pheS.to_string(index=False))
