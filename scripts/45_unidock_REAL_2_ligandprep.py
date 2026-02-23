### CAUTION: THIS SCRIPT NEEDS TO BE RUN 
### WITHIN THE unidock_env CONDA ENVIRONMENT
### using the in-site unidocktools

import tarfile
import shutil
import os

root = "/home/acomajuncosa/Documents_GPU/mtb-targeted-protein-degradation"

# Define and create directories
PATH_TO_CONFORMATIONS = os.path.join(root, "processed", "unidock_REAL_docking_2", "conformations")
PATH_TO_CONFORMATIONS_PREPARED = os.path.join(root, "processed", "unidock_REAL_docking_2", "conformations_prepared")
os.makedirs(PATH_TO_CONFORMATIONS_PREPARED, exist_ok=True)

# Create file with ligand paths
print("Creating file with all ligand paths...")
PATH_TO_INPUT_LIGANDS = os.path.join(root, "processed", "unidock_REAL_docking_2", "input_ligands.txt")
with open(PATH_TO_INPUT_LIGANDS, "w") as outfile:
    for sdf in sorted(os.listdir(PATH_TO_CONFORMATIONS)):
            outfile.write(os.path.join(".", "..", "processed", "unidock_REAL_docking_2", "conformations", sdf))
            outfile.write("\n")

# Prepare sdfs for docking
# CAUTION: unidocktools NEEDS TO BE INSTALLED TO RUN THIS CODE, please visit https://github.com/dptech-corp/Uni-Dock/tree/main/unidock_tools
print("Preparing compounds for docking using unidocktools...")
COMMAND = f"unidocktools ligandprep --ligand_index {PATH_TO_INPUT_LIGANDS} --savedir {os.path.join(root, "processed", "unidock_REAL_docking_2", "conformations_prepared")} --batch_size 100 --use_file_name"
os.system(COMMAND)

# Create file with ligand paths
print("Creating file with all ligand paths again, discarding those that failed...")
with open(PATH_TO_INPUT_LIGANDS, "w") as outfile:
    for sdf in sorted(os.listdir(PATH_TO_CONFORMATIONS_PREPARED)):
            outfile.write(os.path.join(".", "..", "processed", "unidock_REAL_docking_2", "conformations_prepared", sdf))
            outfile.write("\n")