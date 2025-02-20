import os
import pandas as pd

from pymol import cmd
from Bio.PDB import PDBParser, Superimposer, PDBIO

# Align structues
root = os.path.dirname(os.path.abspath(__file__))

# We will align relaxed structures with their unrelaxed counterparts ('aligned')
# to get the aligned_relaxed_structures
print("Aligning relaxed structures to their unrelaxed counterparts...")

def get_sequence_length_in_pdb(pdb_file):
    # Initialize PDB parser
    parser = PDBParser(QUIET=True)

    # Parse the PDB file
    structure = parser.get_structure("structure", pdb_file)

    # Initialize variables to store chain information
    chain_lengths = {}

    # Iterate over chains in the first model
    for chain in structure[0]:
        chain_id = chain.id
        # Count residues in the chain
        residues = [res for res in chain if res.id[0] == " "]  # Exclude heteroatoms
        if len(residues) > 0:
            chain_lengths[chain_id] = len(residues)

    # Print the results
    print(f"Number of chains: {len(chain_lengths)}")
    for k,v in chain_lengths.items():
        print(k,v)
    if len(chain_lengths) > 1:
        raise ValueError("Multiple chains detected. It should only be one...")
    
    for chain_id, length in chain_lengths.items():
        return length
    
def get_sequence_length_from_pdb_from_pymol(pdb_file, chain_id=None):
    """
    Reads a PDB file, loads it into PyMOL, and calculates the sequence length.
    
    Parameters:
        pdb_file (str): Path to the PDB file.
        chain_id (str, optional): Specific chain to calculate the sequence length for.
                                  If None, calculates for all chains.
    
    Returns:
        dict: Dictionary of chain IDs and their sequence lengths.
    """
    print("PDB file for sequence length", pdb_file)
    # Load the PDB file into PyMOL
    cmd.load(pdb_file, "structure")
    
    # Initialize a dictionary to store sequence lengths
    sequence_lengths = {}
    
    # Get chains
    chains = cmd.get_chains("structure")
    
    # Iterate through chains to calculate sequence lengths
    for chain in chains:
        if chain_id is None or chain == chain_id:
            # Get unique residues in the chain
            residues = cmd.get_model(f"structure and chain {chain} and name CA").atom
            unique_residues = {atom.resi for atom in residues}
            
            # Store the length in the dictionary
            sequence_lengths[chain] = len(unique_residues)
    
    # Clean up PyMOL session
    cmd.delete("structure")
    
    return sequence_lengths


def align_structure(reference_file, mobile_file):
    print(mobile_file)
    reference_name = os.path.basename(reference_file).split(".pdb")[0]
    mobile_name = os.path.basename(mobile_file).split(".pdb")[0]

    # Define file paths and parameters
    full_structure_file = reference_file
    subset_structure_file = mobile_file
    input_dir = os.path.join(root, "..", "processed", "relaxed_structures")
    output_dir = os.path.join(root, "..", "processed", "aligned_relaxed_structures")
    aligned_structure_file = mobile_file.replace(input_dir, output_dir)

    # Sequence lengths...
    n_res_ref = get_sequence_length_in_pdb(full_structure_file)
    n_res_mob = get_sequence_length_in_pdb(subset_structure_file)
    n_res_mob2 = get_sequence_length_from_pdb_from_pymol(subset_structure_file)
    print("Reference sequence length: ", n_res_ref)
    print("Mobile sequence length: ", n_res_mob)
    print("Mobile sequence length from PyMOL: ", n_res_mob2)

    # Define the residue range in the full structure that corresponds to the subset
    df = pd.read_csv(os.path.join(root, "..", "processed", "trna_synthetases_data.csv"))
    mobile_file_name = os.path.basename(mobile_file)
    start_residue = df[df["file_name"] == mobile_file_name]["start_resid"].values[0]
    end_residue = df[df["file_name"] == mobile_file_name]["end_resid"].values[0]
    print(start_residue, end_residue)

    # Load the structures
    parser = PDBParser(QUIET=True)
    full_structure = parser.get_structure(reference_name, full_structure_file)
    subset_structure = parser.get_structure(mobile_name, subset_structure_file)

    # Select Cα atoms from the full structure (corresponding to the subset region)
    ref_atoms = [
        residue["CA"] for residue in full_structure[0]["A"]
        if start_residue <= residue.id[1] <= end_residue
    ]

    # Select Cα atoms from the subset structure
    chain_id = list(subset_structure[0])[0].get_id()
    try:
        sample_atoms = [
            residue["CA"] for residue in subset_structure[0][chain_id]
            if start_residue <= residue.id[1] <= end_residue
        ]
    except:
        for residue in subset_structure[0][chain_id]:
            print(residue)

    print(len(ref_atoms), len(sample_atoms))

    # Perform the superimposition
    super_imposer = Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(subset_structure.get_atoms())

    # Print the RMSD value
    print(f"RMSD: {super_imposer.rms:.4f} Å")

    # Save the aligned subset structure
    io = PDBIO()
    io.set_structure(subset_structure)
    io.save(aligned_structure_file)

    print(f"Aligned structure saved to {aligned_structure_file}")

    return round(super_imposer.rms, 2)


def align_all_structures(uniprot_ac):
    input_dir = os.path.join(root, "..", "processed", "relaxed_structures", uniprot_ac)
    output_dir = os.path.join(root, "..", "processed", "aligned_relaxed_structures", uniprot_ac)
    aligned_dir = os.path.join(root, "..", "processed", "aligned_structures", uniprot_ac)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    file_names = []
    for file_name in os.listdir(input_dir):
        if file_name.endswith(".pdb"):
            file_names.append(file_name)
    file_names = sorted(file_names)
    # for file_name in file_names:
    #     if file_name.startswith("alphafold3"):
    #         reference_file = os.path.join(input_dir, file_name)
    #         break
    # print("Reference file is {}".format(reference_file))
    data = {}
    for file_name in file_names:
        if file_name.startswith("alphafill"):
            print("Skipping AlphaFill structure since it is identical to AlphaFold2...")
            continue
        mobile_file = os.path.join(input_dir, file_name)
        reference_file = os.path.join(aligned_dir, file_name)
        print("Aligning {} to the reference...".format(file_name))
        rmsd = align_structure(reference_file, mobile_file)
        print("Structure {} aligned.".format(file_name))
        data[file_name] = rmsd
    return data


if __name__ == "__main__":
    df = pd.read_csv(os.path.join(root, "..", "processed", "trna_synthetases_data.csv"))
    uniprot_acs = df["uniprot_ac"].unique()
    R = []
    for uniprot_ac in uniprot_acs:
        print("----------------------------------------")
        print(f"Aligning structures for {uniprot_ac}...")
        data = align_all_structures(uniprot_ac)
        print(f"All structures aligned for {uniprot_ac}.")
        for k,v in data.items():
            R.append([uniprot_ac, k, v])
        break
    df = pd.DataFrame(R, columns=["uniprot_ac", "file_name", "rmsd"])
    to_remove = df[df["rmsd"] > 10]["file_name"].tolist()
    df = df[df["rmsd"] <= 10]
    df.to_csv(os.path.join(root, "..", "processed", "alignment_relaxed_rmsd_data.csv"), index=False)
    # for fn in to_remove:
    #     uniprot_ac = fn.split("_")[1]
    #     os.remove(os.path.join(root, "..", "processed", "structures", uniprot_ac, fn))
    #     os.remove(os.path.join(root, "..", "processed", "aligned_structures", uniprot_ac, fn))
    #     da = pd.read_csv(os.path.join(root, "..", "processed", "trna_synthetases_data.csv"))
    #     da = da[da["file_name"] != fn]
    #     da.to_csv(os.path.join(root, "..", "processed", "trna_synthetases_data.csv"), index=False)
    #     break