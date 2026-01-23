# to be run in herbert
from scipy.spatial import distance
import pandas as pd
import numpy as np
from tqdm import tqdm
import os

root = "/home/acomajuncosa/Documents_GPU/mtb-targeted-protein-degradation/notebooks"

perc = "1"

# Defining path to Enamine screening results and path to output
PATH_TO_RESULTS = "/aloy/home/acomajuncosa/Ersilia/gcadda4tb-enamine-real-screening/results"
PATH_TO_OUTPUT = os.path.join(root, "..", "processed", 'unidock_REAL_docking', 'inference_10B')
os.makedirs(os.path.join(PATH_TO_OUTPUT, "shared_compounds"), exist_ok=True)

# Load chunk info
CHUNKS = pd.read_csv("/aloy/home/acomajuncosa/Ersilia/gcadda4tb-enamine-real-screening/data/chunks/chunks.csv", header=None)[0].tolist()

# Load pocket info
POCKETS = sorted(set(sorted([i.split("_ind_")[0] for i in sorted(os.listdir(os.path.join(PATH_TO_RESULTS, "Enamine_REAL_LeadLike_000")))])))

# Load pocket detection data
pocket_detection_data = pd.read_csv(os.path.join(root, "..", "processed", "pocket_detection_data.csv"))
pocket_to_coords = {}

for file, pocket_numb, coord in zip(pocket_detection_data['File name'], pocket_detection_data['Pocket number'], pocket_detection_data['Pocket centroid coordinate (x y z)']):
    
    # Prepare labels
    label = file.replace(".pdb", "") + "_pocket_" + str(pocket_numb)

    # Create mappings
    pocket_to_coords[label] = np.array(coord.split(), dtype=float)

# Identify pairs
diff_protein, same_protein, same_pocket = set(), set(), set()

for c, pocket1 in tqdm(enumerate(POCKETS)):
    for pocket2 in POCKETS[c+1:]:

        # Define pair
        pair = (pocket1, pocket2)
        
        # If the protein is different
        if pocket1.split("_")[1] != pocket2.split("_")[1]:
                diff_protein.add(pair)

        else:

            # Calculate distance
            d = distance.euclidean(pocket_to_coords[pocket1], pocket_to_coords[pocket2])

            # If the pocket is the same
            if d < 6.14:
                same_pocket.add(pair)

            # If the pocket is not the same
            else:
                same_protein.add(pair)

print(f"Same pocket: {len(same_pocket)}")
print(f"Same protein: {len(same_protein)}")
print(f"Different protein: {len(diff_protein)}")


# For each chunk
for chunk in CHUNKS:

    TOPs = {}
    SHARED_HITS = []

    # For each pocket
    for pocket in tqdm(POCKETS):

        # Get indices at specified percentile
        inds = np.load(os.path.join(PATH_TO_RESULTS, chunk, f"{pocket}_ind_{perc}.npz"))[f"ind_{perc}"]
        cut = float(np.load(os.path.join(PATH_TO_RESULTS, chunk, f"{pocket}_ind_{perc}.npz"))["thr"])

        # Store indices
        TOPs[pocket] = set(inds.tolist())

    # Evaluate shared hits among pairs
    SHARED_HITS_SAME_POCKET = []
    SHARED_HITS_SAME_PROTEIN = []
    SHARED_HITS_DIFF_PROTEIN = []

    # Calculate intersections
    for c1, pocket1 in tqdm(enumerate(sorted(TOPs))):
        for pocket2 in sorted(TOPs)[c1+1:]:
            intersection = len(TOPs[pocket1].intersection(TOPs[pocket2]))
            if (pocket1, pocket2) in same_pocket:
                SHARED_HITS_SAME_POCKET.append(intersection)
            if (pocket1, pocket2) in same_protein:
                SHARED_HITS_SAME_PROTEIN.append(intersection)
            if (pocket1, pocket2) in diff_protein:
                SHARED_HITS_DIFF_PROTEIN.append(intersection)

    # Save chunk results
    SHARED_HITS.append([chunk, perc, "SAME_POCKET", len(SHARED_HITS_SAME_POCKET), ";".join([str(k) for k in sorted(SHARED_HITS_SAME_POCKET)])])
    SHARED_HITS.append([chunk, perc, "SAME_PROTEIN", len(SHARED_HITS_SAME_PROTEIN), ";".join([str(k) for k in sorted(SHARED_HITS_SAME_PROTEIN)])])
    SHARED_HITS.append([chunk, perc, "DIFF_PROTEIN", len(SHARED_HITS_DIFF_PROTEIN), ";".join([str(k) for k in sorted(SHARED_HITS_DIFF_PROTEIN)])])
    SHARED_HITS = pd.DataFrame(SHARED_HITS, columns=['chunk', 'perc', 'set', 'count', 'indices'])
    SHARED_HITS.to_csv(os.path.join(PATH_TO_OUTPUT, "shared_compounds", f"{chunk}.csv"), index=False)