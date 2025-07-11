### CAUTION: THIS SCRIPT NEEDS TO BE RUN 
### WITHIN THE unidock_env CONDA ENVIRONMENT
### using the in-site unidocktools

import tarfile
import shutil
import os

root = os.path.dirname(os.path.abspath(__file__))

# Define and create directories
PATH_TO_CONFORMATIONS_ENAMINE = "../processed/unidock_docking/conformations"
PATH_TO_CONFORMATIONS_PREPARED = "../processed/unidock_docking/conformations_prepared"
os.makedirs(PATH_TO_CONFORMATIONS_ENAMINE, exist_ok=True)
os.makedirs(PATH_TO_CONFORMATIONS_PREPARED, exist_ok=True)

# Extract conformations
print("Exctracting conformations...")
with tarfile.open("../processed/enamine_characterization/conformations.tar.gz", 'r') as tar:
        tar.extractall(path=PATH_TO_CONFORMATIONS_ENAMINE, filter="data")

# Create file with ligand paths
print("Creating file with all ligand paths...")
with open("../processed/unidock_docking/input_ligands.txt", "w") as outfile:
    for sdf in sorted(os.listdir(PATH_TO_CONFORMATIONS_ENAMINE)):
            outfile.write(os.path.join("..", "processed", "unidock_docking", "conformations", sdf))
            outfile.write("\n")

# Prepare sdfs for docking
# CAUTION: unidocktools NEEDS TO BE INSTALLED TO RUN THIS CODE, please visit https://github.com/dptech-corp/Uni-Dock/tree/main/unidock_tools
print("Preparing compounds for docking using unidocktools...")
COMMAND = f"unidocktools ligandprep --ligand_index ../processed/unidock_docking/input_ligands.txt --savedir ../processed/unidock_docking/conformations_prepared --batch_size 100 --use_file_name"
os.system(COMMAND)

# Create file with ligand paths
print("Creating file with all ligand paths again, discarding those that failed...")
with open("../processed/unidock_docking/input_ligands.txt", "w") as outfile:
    for sdf in sorted(os.listdir(PATH_TO_CONFORMATIONS_PREPARED)):
            outfile.write(os.path.join("..", "processed", "unidock_docking", "conformations", sdf))
            outfile.write("\n")

print("Organizing directory...")

# Tar prepared sdfs
shutil.make_archive(PATH_TO_CONFORMATIONS_PREPARED, 'gztar', PATH_TO_CONFORMATIONS_PREPARED)

# Remove directories
shutil.rmtree(PATH_TO_CONFORMATIONS_ENAMINE)
shutil.rmtree(PATH_TO_CONFORMATIONS_PREPARED)