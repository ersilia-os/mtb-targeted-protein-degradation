### CAUTION: THIS SCRIPT NEEDS TO BE RUN 
### WITHIN THE unidock_env CONDA ENVIRONMENT
### IN A GPU MACHINE FOR MASSIVE PARALLELIZATION

import pandas as pd
import subprocess
import tarfile
import shutil
import os
import csv

def run_unidock(
    receptor: str,
    ligand_index: str,
    center_x: float,
    center_y: float,
    center_z: float,
    size_x: float = 22.5,
    size_y: float = 22.5,
    size_z: float = 22.5,
    seed: int = 42,
    max_gpu_memory: int = 22000,
    num_modes: int = 1,
    verbosity: int = 1,
    cpu: int = 32,
    search_mode: str = "detail",
    scoring: str = "vina",
    output_dir: str = "output",
    log_file: str = "log.txt"
):
    
    """
    Runs the Uni-Dock docking command with the provided parameters. For further information, please check https://github.com/dptech-corp/Uni-Dock
    """ 


    cmd = [
        "unidock",
        "--receptor", receptor,
        "--ligand_index", ligand_index,
        "--center_x", str(center_x),
        "--center_y", str(center_y),
        "--center_z", str(center_z),
        "--size_x", str(size_x),
        "--size_y", str(size_y),
        "--size_z", str(size_z),
        "--seed", str(seed),
        "--max_gpu_memory", str(max_gpu_memory),
        "--num_modes", str(num_modes),
        "--verbosity", str(verbosity),
        "--cpu", str(cpu),
        "--search_mode", search_mode,
        "--scoring", scoring,
        "--dir", output_dir
    ]

    print(" ".join(cmd))

    with open(log_file, "w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)

def extract_score_from_sdf(filepath: str) -> tuple[str, str | None]:
    """
    Extracts the ENERGY value from the <Uni-Dock RESULT> block in an SDF file.
    """
    try:
        with open(filepath, 'r') as file:
            for line in file:
                if line.strip() == "> <Uni-Dock RESULT>":
                    result_line = next(file).strip()
                    if "ENERGY=" in result_line:
                        return os.path.basename(filepath), result_line.split("ENERGY=")[1].split()[0]
    except Exception:
        pass
    return os.path.basename(filepath), None

def generate_report(directory: str, output_csv: str = "scores.csv", score_field: str = "score"):
    """
    Processes all .sdf files in a directory, extracts score values,
    and writes the results as a CSV: filename,score.
    """

    # Enumerate SDF files
    sdf_files = [os.path.join(directory, f) for f in sorted(os.listdir(directory)) if f.endswith(".sdf")]

    # Get results
    results = [extract_score_from_sdf(file) for file in sdf_files]

    # Write results
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["compound", "score"])
        for compound, score in results:
            writer.writerow([compound.replace("_out.sdf", ""), score if score is not None else ""])


# Define paths
# root = '.'
root = os.path.dirname(os.path.abspath(__file__))
UNIDOCK_PATH = os.path.join(root, "..", "processed", "unidock_docking")
OUTPATH = os.path.join(UNIDOCK_PATH, "docking_results")
os.makedirs(OUTPATH, exist_ok=True)

# Load pocket detection data
pocket_detection_data = pd.read_csv(os.path.join(root, "..", "processed", "pocket_detection_data.csv"))

# Shuffle with fixed seed
df = pocket_detection_data.sample(frac=1, random_state=42).reset_index(drop=True)

# Untar conformations prepared
shutil.unpack_archive(os.path.join(UNIDOCK_PATH, "conformations_prepared.tar.gz"), 
                      os.path.join(UNIDOCK_PATH, "conformations_prepared"), 'gztar')


# For each pocket
for file, pocket_number, centroid in zip(df['File name'], df['Pocket number'], df['Pocket centroid coordinate (x y z)']):

    # Prints
    print("\n\n--- RUNNING DOCKING ---\n\n")
    print(f"\n\n--- File name: {file}; Pocket number: {pocket_number} ---\n\n")

    # Create directories
    label = file.replace(".pdb", "") + "_pocket_" + str(pocket_number)
    outpath = os.path.join(OUTPATH, label)
    os.makedirs(os.path.join(outpath, "docking"), exist_ok=True)

    # Extract structure
    tarfile.open(os.path.join(UNIDOCK_PATH, "structures_prepared.tar.gz")).extract("./" + file.replace(".pdb", ".pdbqt"), 
                                                                                   path=outpath, filter='data')

    # Copy pocket SD file
    shutil.copyfile(os.path.join(root, "..", "output", "pocketvec_RUN", "pocketvec_PRE", label, f"pocket_{label.split('_')[-1]}.sd"), 
                    os.path.join(outpath, f"pocket_{label.split('_')[-1]}.sd"))
    

    # Prepare docking variables
    receptor = os.path.join(outpath, file.replace(".pdb", ".pdbqt"))
    center_x, center_y, center_z = centroid.split()
    ligand_index = os.path.join(UNIDOCK_PATH, "input_ligands.txt")
    search_mode = 'fast'
    output_dir = os.path.join(outpath, "docking")
    log_file = os.path.join(outpath, "logs.log")

    # Run docking
    run_unidock(receptor=receptor, ligand_index=ligand_index, center_x=center_x, center_y=center_y, center_z=center_z,
                search_mode=search_mode, output_dir=output_dir, log_file=log_file)
    
    print("Generating report!")
    
    # Generate report
    generate_report(os.path.join(outpath, "docking"), os.path.join(outpath, "report.csv"))

    print("Compressing results!")

    # Tar results
    with tarfile.open(os.path.join(outpath, "docking.tar.gz"), "w:gz", compresslevel=9) as tar:
        tar.add(os.path.join(outpath, "docking"), arcname=os.path.basename(os.path.join(outpath, "docking")))

    # Tar logs
    with tarfile.open(os.path.join(outpath, "logs.tar.gz"), "w:gz") as tar:
        tar.add(os.path.join(outpath, "logs.log"), arcname="logs.log")

    # Remove file and directory
    os.remove(os.path.join(outpath, "logs.log"))
    shutil.rmtree(os.path.join(outpath, "docking"))



    