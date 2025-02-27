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

def extract_pocket_scores(csv_file):
    """
    Reads a P2Rank output CSV file and maps pocket number (rank) to its score.

    :param csv_file: Path to the P2Rank output CSV file
    :return: Dictionary mapping pocket number to score
    """
    df = pd.read_csv(csv_file)
    df.columns = df.columns.str.strip()

    return {int(row["rank"]): row["score"] for _, row in df.iterrows()}


def extract_pocket_probabilities(csv_file):
    """
    Reads a P2Rank output CSV file and maps pocket number (rank) to its probability.

    :param csv_file: Path to the P2Rank output CSV file
    :return: Dictionary mapping pocket number to probability
    """
    df = pd.read_csv(csv_file)
    df.columns = df.columns.str.strip()

    return {int(row["rank"]): row["probability"] for _, row in df.iterrows()}

def extract_pocket_residues(csv_file):
    """
    Reads a P2Rank output CSV file and maps pocket number (rank) to its associated residue IDs.

    :param csv_file: Path to the P2Rank output CSV file
    :return: Dictionary mapping pocket number to a list of residue IDs
    """
    df = pd.read_csv(csv_file)
    df.columns = df.columns.str.strip()

    return {int(row["rank"]): row["residue_ids"].split() for _, row in df.iterrows()}



def write_pocket_pdbs(pockets_to_consider, pocket_dict, output_dir):
    """
    Creates a separate PDB file for each pocket, strictly following PDB format specifications.

    :param pocket_dict: Dictionary mapping pocket number to (x, y, z) coordinates.
    :param output_dir: Directory where the output PDB files will be stored.
    """
    os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists

    for pocket_number, (x, y, z) in pocket_dict.items():
        if pocket_number in pockets_to_consider:
            pdb_filename = os.path.join(output_dir, f"pocket_{pocket_number}.pdb")
            with open(pdb_filename, "w") as pdb_file:
                pdb_file.write(f"HEADER    Pocket {pocket_number} Centroid\n")
                pdb_file.write(
                    f"HETATM{pocket_number:5d} C   LIG A{pocket_number:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          C\n"
                )
                pdb_file.write("END\n")

def get_protein_prediction_type(file_name):
    """
    Determines the protein structure prediction type based on the file name.

    :param file_name: Name of the PDB file
    :return: String indicating the prediction type ('alphafold2', 'alphafold3', 'chai1', 'swissmodel', or 'unknown')
    """
    file_name = file_name.lower()  # Convert to lowercase for case insensitivity

    if "alphafold2" in file_name or "af2" in file_name:
        return "alphafold2"
    elif "alphafold3" in file_name or "af3" in file_name:
        return "alphafold3"
    elif "chai1" in file_name:  # Strictly check for "chai1"
        return "chai1"
    elif "swissmodel" in file_name:
        return "swissmodel"
    else:
        return "unknown"  # If no match is found

def extract_residue_confidence_mapping(pdb_file, prediction_type):
    """
    Extracts a mapping of residue numbers to confidence values from a PDB file.
    
    Assumes confidence values are stored in the B-factor column.

    :param pdb_file: Path to the PDB file.
    :param prediction_type: Prediction type ('alphafold2', 'alphafold3', 'chai1', 'swissmodel'). At this moment, we're not using this parameter.
    :return: Dictionary mapping residue numbers (int) to confidence values (float).
    
    alphafold2: pLDDT per residue
    alphafold3: pLDDT per atom (we consider the minimum per residue)
    chai-1: equivalent to alphafold3
    swissmodel: local confidence scores - https://swissmodel.expasy.org/docs/help#qmean

    """

    res_num_to_conf = {}

    with open(pdb_file, "r") as file:
        for line in file:
            if line.startswith("ATOM"):  # Only process atomic coordinate lines
                residue_number = int(line[22:26].strip())  # Extract residue ID
                confidence_value = float(line[60:66].strip())  # Extract B-factor (confidence score)

                # Store the highest confidence value per residue (in case of multiple atoms per residue)
                if residue_number not in res_num_to_conf:
                    res_num_to_conf[residue_number] = confidence_value
                else:
                    res_num_to_conf[residue_number] = min(res_num_to_conf[residue_number], confidence_value)

    return res_num_to_conf  # { residue_number: confidence_value }


# Define paths
root = os.path.dirname(os.path.abspath(__file__))
original_structures_dir = os.path.abspath(os.path.join(root, "..", "processed", "structures"))  # We need original structures to get the confidence values
aligned_dir = os.path.abspath(os.path.join(root, "..", "processed", "aligned_relaxed_structures"))
detected_pockets_dir = os.path.abspath(os.path.join(root, "..", "processed", "detected_pockets"))

# Load alignment RMSD data
alignment_df = pd.read_csv(os.path.abspath(os.path.join(root, "..", "processed", "alignment_relaxed_rmsd_data.csv")))

# Create report
report = []

# Iterate over Uniprot ACs
for uniprot_ac, file_name in zip(alignment_df["uniprot_ac"], alignment_df["file_name"]):

    print("---------------  Detecting pockets in: " + uniprot_ac + ", " + file_name + "  -------------")
    
    # Create folder in detected_pockets_dir
    os.makedirs(os.path.join(detected_pockets_dir, uniprot_ac), exist_ok=True)

    # Run P2RANK
    detect_pockets(os.path.join(aligned_dir, uniprot_ac, file_name), 
                   os.path.join(detected_pockets_dir, uniprot_ac, file_name.replace(".pdb", "")))

    # Get pocket centroids, scores, probabilities and residues involved
    pocket_centroids = extract_pocket_centers(os.path.join(detected_pockets_dir, uniprot_ac, file_name.replace(".pdb", ""), file_name + "_predictions.csv"))
    pocket_scores = extract_pocket_scores(os.path.join(detected_pockets_dir, uniprot_ac, file_name.replace(".pdb", ""), file_name + "_predictions.csv"))
    pocket_probabilities = extract_pocket_probabilities(os.path.join(detected_pockets_dir, uniprot_ac, file_name.replace(".pdb", ""), file_name + "_predictions.csv"))
    pocket_residues = extract_pocket_residues(os.path.join(detected_pockets_dir, uniprot_ac, file_name.replace(".pdb", ""), file_name + "_predictions.csv"))

    # Select only a limited number of pockets - check https://github.com/rdk/p2rank/issues/76
    P = 0.2  # probability
    K = 3  # TOP-K
    pockets_to_consider = set([i for i in sorted(pocket_centroids) if i < K or pocket_probabilities[i] >= P])

    # Get protein structure prediction type
    prediction_type = get_protein_prediction_type(file_name)

    # Get pocket residue confidence scores - assumes bfactors are confidence scores
    residues_confidence = extract_residue_confidence_mapping(os.path.join(original_structures_dir, uniprot_ac, file_name), prediction_type)

    # Define minimum confidence value for a pocket to be considered
    if prediction_type in ['alphafold2', 'alphafold3', 'chai1']:
        MIN = 65
    elif prediction_type in ['swissmodel']:
        MIN = 0.65

    # Keep only those pockets with a MIN confidence value above threshold
    filtered_pockets = set()
    for pocket in sorted(pockets_to_consider):
        residues = [int(res.split("_")[1]) for res in pocket_residues[pocket]]
        confidence = [residues_confidence[i] for i in residues]
        if min(confidence) >= MIN:
            filtered_pockets.add(pocket)

    pockets_to_consider = filtered_pockets


    # Create a single PDB file per detected pocket
    write_pocket_pdbs(pockets_to_consider, pocket_centroids, os.path.join(detected_pockets_dir, uniprot_ac, file_name.replace(".pdb", ""), "pockets"))

    # Save results
    for ptc in sorted(pockets_to_consider):
        report.append([uniprot_ac, file_name, prediction_type, os.path.join("/".join(aligned_dir.split("/")[-2:]), uniprot_ac, file_name), ptc, 
            pocket_scores[ptc], pocket_probabilities[ptc], " ".join([str(cent) for cent in pocket_centroids[ptc]]), " ".join(pocket_residues[ptc]),
            " ".join([str(residues_confidence[int(j.split("_")[1])]) for j in pocket_residues[ptc]])])

    print("---------------  Pockets detected in: " + uniprot_ac + ", " + file_name + "   -------------")

report = pd.DataFrame(report, columns=['Uniprot AC', 'File name', 'Prediction type', 'Full path', 'Pocket number', 'Pocket score', 'Pocket probability', 'Pocket centroid coordinate (x y z)', 'Pocket residues (chain_resn)', 'B-factors'])
report.to_csv(os.path.abspath(os.path.join(root, "..", "processed", "pocket_detection_data.csv")), index=False)