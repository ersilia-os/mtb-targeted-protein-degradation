import os
import subprocess
import pandas as pd
from pymol import cmd
cmd.feedback("disable", "all", "everything")  # Disables all output, comment this if something doesn't go as expected

def prepare_pymol_session(uni, directory, reference, structures, pymol_sessions, alignment_df, detected_pockets):

    # Define some colors
    COLOR_REFERENCE = 'blue'
    COLOR_ALIGNED = 'grey'
    COLOR_POCKETS = 'orange'

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
    cmd.do("set transparency, 0.1")  
    cmd.do("set sphere_scale, 2")

    # Load all structures
    for st in reference + structures:

        # Load structure
        cmd.load(os.path.join(directory, f"{st}.pdb"), st)

        # Color reference structure differently
        if st == reference[0]:  
            cmd.color(COLOR_REFERENCE, st)
            cmd.show("cartoon", st)
            cmd.show("surface", st)  
        else:
            cmd.color(COLOR_ALIGNED, st)
            cmd.show("cartoon", st)
            cmd.disable(st)


    # Load all detected pockets:
    for st in reference + structures:

        # Get pockets
        pockets_to_consider = alignment_df[(alignment_df['Uniprot AC'] == uni) & (alignment_df['File name'] == f"{st}.pdb")]['Pocket number'].tolist()

        # Load pockets
        for ptc in sorted(pockets_to_consider):

            # Load pocket
            cmd.load(os.path.join(detected_pockets, uni, st, "pockets", f"pocket_{ptc}.pdb"), f"pocket_{ptc}_{st}")

            # Color reference structure differently
            if st == reference[0]:  
                cmd.color(COLOR_POCKETS, f"pocket_{ptc}_{st}")
                cmd.show("spheres", f"pocket_{ptc}_{st}")  
            else:
                cmd.color(COLOR_ALIGNED, f"pocket_{ptc}_{st}")
                # cmd.show("spheres", f"pocket_{ptc}_{st}")
                cmd.disable(f"pocket_{ptc}_{st}")


    # Save PyMOL session
    cmd.reset()
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
alignment_df = pd.read_csv(os.path.abspath(os.path.join(root, "..", "processed", "pocket_detection_data.csv")))
uniprots = sorted(set(alignment_df["Uniprot AC"]))

# For each uniprot
for uni in uniprots:

    # Specify directory
    directory = os.path.join(aligned_relaxed_dir, uni)

    # Get reference
    reference = [st.replace(".pdb", "") for st in sorted(os.listdir(directory)) if "alphafold2" in st]

    # Get all the other structures, but not the reference
    structures = [st.replace(".pdb", "") for st in sorted(os.listdir(directory)) if "alphafold2" not in st]

    print(f" -------------- Creating PyMOL session for {uni}  --------------  ")

    # Create pymol session
    prepare_pymol_session(uni, directory, reference, structures, pymol_sessions, alignment_df, detected_pockets)

    print(f" -------------- PyMOL session for {uni} created  ---------------  ")


    break

