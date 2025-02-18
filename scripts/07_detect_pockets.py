import os
import subprocess
import pandas as pd

# Define paths
root = os.path.dirname(os.path.abspath(__file__))
aligned_dir = os.path.abspath(os.path.join(root, "..", "processed", "aligned_structures"))
detected_pockets_dir = os.path.abspath(os.path.join(root, "..", "processed", "detected_pockets"))

# Load alignment RMSD data
alignment_df = pd.read_csv(alignment_rmsd_file)

# Iterate over Uniprot ACs
for uniprot_ac, file_name in zip(df["uniprot_ac"], df["file_name"]):
    print(uniprot_ac, file_name)


    break
