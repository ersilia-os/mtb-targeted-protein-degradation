import pandas as pd
import subprocess
import os

def pdb_to_mol2_obabel(input_pdb, output_mol2=None):
    """
    Converts a PDB file to MOL2 format using Open Babel.

    Parameters:
    - input_pdb (str): Path to the input .pdb file.
    - output_mol2 (str, optional): Path to the output .mol2 file. 
                                   If not provided, uses the same name as input.

    Returns:
    - output_mol2 (str): Path to the generated .mol2 file.
    """

    if not os.path.isfile(input_pdb):
        raise FileNotFoundError(f"Input file '{input_pdb}' does not exist.")

    if output_mol2 is None:
        output_mol2 = os.path.splitext(input_pdb)[0] + ".mol2"

    try:
        subprocess.run(["obabel", "-ipdb", input_pdb, "-omol2", "-O", output_mol2],
                       check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"Converted {input_pdb} to {output_mol2}")
        return output_mol2
    except subprocess.CalledProcessError as e:
        print(f"Error during conversion: {e.stderr.decode().strip()}")
        raise RuntimeError("Open Babel conversion failed.")
    
def pdb_to_sd_obabel(input_pdb, output_sd=None):
    """
    Converts a PDB file to SD (.sd) format using Open Babel.

    Parameters:
    - input_pdb (str): Path to the input .pdb file.
    - output_sd (str, optional): Path to the output .sd file. 
                                 If not provided, uses the same name with .sd extension.

    Returns:
    - output_sd (str): Path to the generated .sd file.
    """

    if not os.path.isfile(input_pdb):
        raise FileNotFoundError(f"Input file '{input_pdb}' does not exist.")

    if output_sd is None:
        output_sd = os.path.splitext(input_pdb)[0] + ".sd"

    try:
        subprocess.run(["obabel", "-ipdb", input_pdb, "-osdf", "-O", output_sd],
                       check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"Converted {input_pdb} to {output_sd}")
        return output_sd
    except subprocess.CalledProcessError as e:
        print(f"Error during conversion: {e.stderr.decode().strip()}")
        raise RuntimeError("Open Babel conversion failed.")

pocket_detection_data = pd.read_csv("../processed/pocket_detection_data.csv")
path_to_pocketvec = "../processed/pocketvec"
path_to_processed = "../processed/aligned_relaxed_structures"
path_to_pockets = "../processed/detected_pockets"

for uni, st, pocket in zip(pocket_detection_data['Uniprot AC'], pocket_detection_data['File name'], pocket_detection_data['Pocket number']):

    name_dir = st.strip(".pdb") + "_pocket_" + str(pocket)

    # Create directory
    os.makedirs(os.path.join(path_to_pocketvec, name_dir), exist_ok=True)

    # Copy/convert structure
    pdb_to_mol2_obabel(os.path.join(path_to_processed, uni, st),
                       os.path.join(path_to_pocketvec, name_dir, st.strip(".pdb") + ".mol2"))
    
    # Copy/convert pocket
    pdb_to_mol2_obabel(os.path.join(path_to_pockets, uni, st.strip(".pdb"), "pockets", f"pocket_{pocket}.pdb"),
                       os.path.join(path_to_pocketvec, name_dir, f"pocket_{pocket}.sd"))