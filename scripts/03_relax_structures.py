import os
import pandas as pd
import os
import sys
import pyrosetta
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory

root = os.path.dirname(os.path.abspath(__file__))
processed_dir = os.path.abspath(os.path.join(root, "..", "processed"))

sys.path.append(os.path.join(root, "..", "src"))
from utils.conda import SimpleConda

CONDA_ENV = "adda4tb"
N_STRUCTURES_FOR_RELAX = 3

# Initialize PyRosetta
pyrosetta.init()


def calculate_protonation_states(input_file, output_file):
    cwd = os.getcwd()
    file_name = os.path.basename(input_file).split(".")[0]
    input_file = os.path.abspath(input_file)
    output_file = os.path.abspath(output_file)
    print(f"Protonation states calculated and saved to '{output_file}'.")
    sc = SimpleConda()
    sc.run_commandlines(CONDA_ENV, ["pdb2pqr --ff=AMBER --with-ph=7.0 --pdb-output {0} {1} {2}.pqr".format(output_file, input_file, file_name)])
    os.remove(os.path.join(cwd, file_name + ".pqr"))
    os.remove(os.path.join(cwd, file_name + ".log"))


def relax_structure(input_file, output_file):

    # Load the PDB structure
    pose = pyrosetta.pose_from_file(input_file)

    # Set up the scoring function
    scorefxn = ScoreFunctionFactory.create_score_function("ref2015")

    # Set up the Relax protocol
    relax = FastRelax()
    relax.set_scorefxn(scorefxn)

    # Number of structures to generate
    n_structures = N_STRUCTURES_FOR_RELAX
    relaxed_structures = []

    # Relax the structure multiple times and score them
    for i in range(n_structures):
        print(f"Relaxing structure {i + 1}...")
        
        # Create a copy of the original pose to avoid overwriting
        relaxed_pose = pose.clone()
        
        # Apply relaxation
        relax.apply(relaxed_pose)
        
        # Get the total score for the relaxed structure
        score = scorefxn(relaxed_pose)
        
        # Store the pose and its score
        relaxed_structures.append((relaxed_pose, score))
        print(f"Score for structure {i + 1}: {score}")

    # Sort the relaxed structures by score (lower is better)
    relaxed_structures.sort(key=lambda x: x[1])

    # Save the best structure
    best_pose, best_score = relaxed_structures[0]
    print(f"Best structure score: {best_score}")
    print(f"Saving best structure to {output_file}")
    best_pose.dump_pdb(output_file)

df = pd.read_csv(os.path.join(processed_dir, "trna_synthetases_data.csv"))
file_names = df["file_name"].tolist()
for file_name in file_names:
    print("----------------------------------------")
    print(f"Relaxing structure {file_name}...")
    uniprot_ac = file_name.split("_")[1]
    input_file = os.path.join(processed_dir, "structures", uniprot_ac, file_name)
    tmp_file = os.path.join(root, "..", "tmp", file_name)
    calculate_protonation_states(input_file, tmp_file)
    print("Protonation states calculated.")
    if not os.path.exists(os.path.join(processed_dir, "relaxed_structures", uniprot_ac)):
        os.makedirs(os.path.join(processed_dir, "relaxed_structures", uniprot_ac))
    output_file = os.path.join(processed_dir, "relaxed_structures", uniprot_ac, file_name)
    relax_structure(tmp_file, output_file)
    print("Removing temporary file {0}.".format(tmp_file))
    os.remove(tmp_file)
    print("Done!")
    print("----------------------------------------")