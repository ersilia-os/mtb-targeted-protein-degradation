import os
import subprocess
import pandas as pd
from pymol import cmd
cmd.feedback("disable", "all", "everything")  # Disables all output, comment this if something doesn't go as expected

def prepare_pymol_session(uni, directory, reference, structures, pymol_sessions):

    # Define some colors
    COLOR_REFERENCE = 'blue'
    COLOR_ALIGNED = 'grey'

    # Initialize PyMOL
    cmd.reinitialize()

    # Prettify session 
    cmd.do("set orthoscopic, on")
    cmd.do("set ray_trace_fog, 0")
    cmd.do("set depth_cue, 0")
    cmd.do("set antialias, 4")
    cmd.do("set ray_trace_mode, 1")
    cmd.do("set ray_trace_gain, 0.005")
    cmd.do("bg_color white")
    cmd.do("set spec_reflect, 0")

    # Load all structures
    for st in reference + structures:

        # Load structure
        cmd.load(os.path.join(directory, f"{st}.pdb"), st)

        # Color reference structure differently
        if reference and st == reference[0]:  
            cmd.color(COLOR_REFERENCE, st)  
        else:
            cmd.color(COLOR_ALIGNED, st)  

    # Save PyMOL session
    cmd.save(os.path.join(pymol_sessions, uni + ".pse.gz"))


# Define paths
root = os.path.dirname(os.path.abspath(__file__))
aligned_dir = os.path.abspath(os.path.join(root, "..", "processed", "aligned_structures"))
aligned_relaxed_dir = os.path.abspath(os.path.join(root, "..", "processed", "aligned_relaxed_structures"))
detected_pockets = os.path.abspath(os.path.join(root, "..", "processed", "detected_pockets"))
pymol_sessions = os.path.abspath(os.path.join(root, "..", "processed", "pymol_sessions"))

# Create dir if necessary
os.makedirs(pymol_sessions, exist_ok=True)

# Load data
alignment_df = pd.read_csv(os.path.abspath(os.path.join(root, "..", "processed", "alignment_relaxed_rmsd_data.csv")))
uniprots = sorted(set(alignment_df["uniprot_ac"]))

# For each uniprot
for uni in uniprots:

    # Specify directory
    directory = os.path.join(aligned_relaxed_dir, uni)

    # Get reference
    reference = [st.replace(".pdb", "") for st in sorted(os.listdir(directory)) if "alphafold2" in st]

    # Get all the other structures, but not the reference
    structures = [st.replace(".pdb", "") for st in sorted(os.listdir(directory)) if "alphafold2" not in st]

    print("\n\n\n\n")
    print(f" -------------- Creating PyMOL session for {uni}  --------------  ")

    # Create pymol session
    prepare_pymol_session(uni, directory, reference, structures, pymol_sessions)

    print(f" -------------- PyMOL session for {uni} created  ---------------  ")
    print("\n\n\n\n")


    break

