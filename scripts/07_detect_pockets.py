import os
import subprocess
import pandas as pd

# Define paths
root = os.path.dirname(os.path.abspath(__file__))
aligned_dir = os.path.abspath(os.path.join(root, "..", "processed", "aligned_structures"))
detected_pockets_dir = os.path.abspath(os.path.join(root, "..", "processed", "detected_pockets"))

# Define path to P2RANK binary package - downloaded from https://github.com/rdk/p2rank/releases
path_to_p2rank = "/aloy/home/acomajuncosa/programs/p2rank_2.5/prank"  # Change as needed 

# Load alignment RMSD data
alignment_df = pd.read_csv(os.path.abspath(os.path.join(root, "..", "processed", "alignment_rmsd_data.csv")))

# Iterate over Uniprot ACs
for uniprot_ac, file_name in zip(alignment_df["uniprot_ac"], alignment_df["file_name"]):
    print(uniprot_ac, file_name)


    break
