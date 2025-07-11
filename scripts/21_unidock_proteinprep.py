### CAUTION: THIS SCRIPT NEEDS TO BE RUN 
### WITHIN THE unidock_tools CONDA ENVIRONMENT
### ambertools and reduce need to be
### installed, and np version needs to be 1.26
### in short: pip install git+https://github.com/dptech-corp/Uni-Dock.git#subdirectory=unidock_tools

import pandas as pd
import shutil
import sys
import os

# Define some paths
pocket_detection_data = pd.read_csv("../processed/pocket_detection_data.csv")
PATH_TO_STRUCTURES = "../processed/aligned_relaxed_structures"
PATH_TO_OUTPUT = "../processed/unidock_docking/structures_prepared"
os.makedirs(PATH_TO_OUTPUT, exist_ok=True)

# Structures to prepare
sts = sorted(set([tuple([i,j]) for i,j in zip(pocket_detection_data['Uniprot AC'], 
                                  pocket_detection_data['File name'])]))

print(f"Totalling {len(sts)} protein structures...")

# For each structure
for st in sts:

    # Copy structure - only if not copied already
    if os.path.exists(os.path.join(PATH_TO_OUTPUT, st[1])) is False:
        shutil.copy(os.path.join(PATH_TO_STRUCTURES, st[0], st[1]), os.path.join(PATH_TO_OUTPUT, st[1]))

    # Prepare protein structure
    print(f"Preparing structure {st}...")
    COMMAND = f"unidocktools proteinprep -r {os.path.join(PATH_TO_OUTPUT, st[1])} -o {os.path.join(PATH_TO_OUTPUT, st[1].replace('.pdb', '.pdbqt'))}"
    ret = os.system(COMMAND)
    if ret != 0:
        print(f"Error: Command failed for structure {st[1]}. Please check manually. Exiting.")
        sys.exit(1)

# Tar file
shutil.make_archive(PATH_TO_OUTPUT, 'gztar', PATH_TO_OUTPUT)

# Remove directory
shutil.rmtree(PATH_TO_OUTPUT)