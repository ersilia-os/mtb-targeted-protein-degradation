import os
import subprocess
import pandas as pd


# Function to run P2Rank on a PDB file
def detect_pockets(pdb_file, output_dir):
    """
    Runs P2Rank to detect pockets in a given PDB file.

    This function calls the P2Rank binary to predict binding pockets 
    on a protein structure. The results are saved in the specified output directory.

    :param pdb_file: Path to the input PDB file.
    :param output_dir: Directory where the detected pocket results will be stored.
    :return: None
    """


    # Define path to P2RANK binary package - downloaded from https://github.com/rdk/p2rank/releases
    path_to_p2rank = "/aloy/home/acomajuncosa/programs/p2rank_2.5/prank"  # Change as needed 
    
    # Construct command
    command = [path_to_p2rank, "predict", "-f", pdb_file, "-c", "alphafold", "-o", output_dir, "-visualizations", "0"]
    
    try:
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"Error detecting pockets in {pdb_file}: {e}")


def extract_pocket_centers(csv_file):
    """
    Reads a P2Rank output CSV file and maps pocket number (rank) to its center coordinates (x, y, z).

    :param csv_file: Path to the P2Rank output CSV file
    :return: Dictionary mapping pocket number to (x, y, z) center coordinates
    """
    # Load the CSV file
    df = pd.read_csv(csv_file)

    # Strip any leading/trailing spaces from column names
    df.columns = df.columns.str.strip()

    # Create dictionary mapping pocket number to (x, y, z) coordinates
    pocket_centers = {
        int(row["rank"]): (row["center_x"], row["center_y"], row["center_z"])
        for _, row in df.iterrows()
    }

    return pocket_centers


def write_pocket_pdbs(pocket_dict, output_dir):
    """
    Creates a separate PDB file for each pocket, strictly following PDB format specifications.

    :param pocket_dict: Dictionary mapping pocket number to (x, y, z) coordinates.
    :param output_dir: Directory where the output PDB files will be stored.
    """
    os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists

    for pocket_number, (x, y, z) in pocket_dict.items():
        pdb_filename = os.path.join(output_dir, f"pocket_{pocket_number}.pdb")
        with open(pdb_filename, "w") as pdb_file:
            pdb_file.write(f"HEADER    Pocket {pocket_number} Centroid\n")
            pdb_file.write(
                f"HETATM{pocket_number:5d} C   LIG A{pocket_number:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          C\n"
            )
            pdb_file.write("END\n")


# Define paths
root = os.path.dirname(os.path.abspath(__file__))
aligned_dir = os.path.abspath(os.path.join(root, "..", "processed", "aligned_structures"))
detected_pockets_dir = os.path.abspath(os.path.join(root, "..", "processed", "detected_pockets"))

# Load alignment RMSD data
alignment_df = pd.read_csv(os.path.abspath(os.path.join(root, "..", "processed", "alignment_rmsd_data.csv")))

# Iterate over Uniprot ACs
for uniprot_ac, file_name in zip(alignment_df["uniprot_ac"], alignment_df["file_name"]):

    print("---------------  Detecting pockets in: " + uniprot_ac + ", " + file_name + "  -------------")
    
    # Create folder in detected_pockets_dir
    os.makedirs(os.path.join(detected_pockets_dir, uniprot_ac), exist_ok=True)

    # Run P2RANK
    detect_pockets(os.path.join(aligned_dir, uniprot_ac, file_name), 
                   os.path.join(detected_pockets_dir, uniprot_ac, file_name.replace(".pdb", "")))

    # Get pocket centroids
    pocket_centroids = extract_pocket_centers(os.path.join(detected_pockets_dir, uniprot_ac, file_name.replace(".pdb", ""), file_name + "_predictions.csv"))

    # Create a single PDB file per detected pocket
    write_pocket_pdbs(pocket_centroids, os.path.join(detected_pockets_dir, uniprot_ac, file_name.replace(".pdb", ""), "pockets"))

    


    break
