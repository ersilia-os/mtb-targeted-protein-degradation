import os
import subprocess
import pandas as pd
from pymol import cmd

# # Suppress executive details and warnings
# cmd.feedback("disable", "executive", "details")
# cmd.feedback("disable", "all", "warnings")


# Define paths
root = os.path.dirname(os.path.abspath(__file__))
aligned_dir = os.path.abspath(os.path.join(root, "..", "processed", "aligned_structures"))
aligned_relaxed_dir = os.path.abspath(os.path.join(root, "..", "processed", "aligned_relaxed_structures"))
detected_pockets = os.path.abspath(os.path.join(root, "..", "processed", "detected_pockets"))

# Load data
alignment_df = pd.read_csv(os.path.abspath(os.path.join(root, "..", "processed", "alignment_relaxed_rmsd_data.csv")))
uniprots = sorted(set(alignment_df["uniprot_ac"]))

# For each uniprot
for uni in uniprots:

    # Get reference


    # Create pymol session




    break

# Initialize PyMOL
cmd.reinitialize()

# Load PDB structure
# cmd.load(input_file, "structure")