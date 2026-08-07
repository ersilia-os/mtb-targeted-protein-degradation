from tqdm import tqdm
from rdkit import Chem
import tarfile
import shutil
import pickle
import gzip
import os

root = "."
# root = os.path.dirname(os.path.abspath(__file__))

# Define and create directories
PATH_TO_CONFORMATIONS_ENAMINE_REAL = os.path.join(root, "..", "output", "enamine_REAL_characterization", "conformations")
PATH_TO_CONFORMATIONS_ENAMINE_REAL_TMP = os.path.join(root, "..", "output", "enamine_REAL_characterization", "conformations_tmp")
PATH_TO_CONFORMATIONS_ENAMINE_REAL_PREPARED = os.path.join(root, "..", "output", "unidock_REAL_docking", "conformations_prepared")
os.makedirs(PATH_TO_CONFORMATIONS_ENAMINE_REAL_TMP, exist_ok=True)
os.makedirs(PATH_TO_CONFORMATIONS_ENAMINE_REAL_PREPARED, exist_ok=True)

# Load failed molecules
FAILED = set([i.strip() for i in open(os.path.join(root, "..", "output", "enamine_REAL_characterization", "failed_REAL.csv")).readlines()])

for chunk in range(0, 96, 1):

    print(f"Loading chunk {chunk}")

    # Load pickle
    filepath = os.path.join(PATH_TO_CONFORMATIONS_ENAMINE_REAL, f"mols_{chunk:02d}.pkl.gz")
    with gzip.open(filepath, "rb") as f:
        part = pickle.load(f)

    print("Creating individual SDF files...")

    for ID in tqdm(sorted(part)):
        if ID not in FAILED:
            # ID to rdkit mol
            mol = part[ID]
            # Create SDF file
            if mol is not None:
                mol.SetProp("_Name", ID)
                mol.SetProp("_ID", ID)
                molblock = Chem.MolToV2KMolBlock(mol)
                out_path = os.path.join(PATH_TO_CONFORMATIONS_ENAMINE_REAL_TMP, ID + ".sdf")
                with open(out_path, "w") as f:
                    f.write(molblock)
                    f.write(f">  <_ID>\n{ID}\n\n")
                    f.write("$$$$\n")
            else:
                raise Exception(f"Mol with ID={ID} is None and was not included in failed_REAL.csv. Please check.")
            
    # Create file with ligand paths
    print("Creating file with all ligand paths...")
    with open(os.path.join(root, "..", "output", "enamine_REAL_characterization", "input_ligands_tmp.txt"), "w") as outfile:
        for sdf in sorted(os.listdir(PATH_TO_CONFORMATIONS_ENAMINE_REAL_TMP)):
                outfile.write(os.path.join(".", "..", "output", "enamine_REAL_characterization", "conformations_tmp", sdf))
                outfile.write("\n")

    # Create subdirectory
    os.makedirs(os.path.join(PATH_TO_CONFORMATIONS_ENAMINE_REAL_PREPARED, f"{chunk:02d}"), exist_ok=True)

    # Prepare sdfs for docking
    # CAUTION: unidocktools NEEDS TO BE INSTALLED TO RUN THIS CODE, please visit https://github.com/dptech-corp/Uni-Dock/tree/main/unidock_tools
    print("Preparing compounds for docking using unidocktools...")
    COMMAND = f"unidocktools ligandprep --ligand_index {os.path.join(root, "..", "output", "enamine_REAL_characterization", "input_ligands_tmp.txt")} \
        --savedir {os.path.join(PATH_TO_CONFORMATIONS_ENAMINE_REAL_PREPARED, f"{chunk:02d}")} --batch_size 100 --use_file_name"
    os.system(COMMAND)

    print("Removing tmp files and folders")
    shutil.rmtree(PATH_TO_CONFORMATIONS_ENAMINE_REAL_TMP)
    os.makedirs(PATH_TO_CONFORMATIONS_ENAMINE_REAL_TMP, exist_ok=True)
    os.remove(os.path.join(root, "..", "output", "enamine_REAL_characterization", "input_ligands_tmp.txt"))  