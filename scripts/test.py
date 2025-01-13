import os
import pandas as pd
from pymol import cmd

root = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(root, "..", "data"))
processed_dir = os.path.abspath(os.path.join(root, "..", "processed"))

def get_chain_ids(input_file):
    """
    Extract all chain IDs from a CIF file.

    Args:
        input_file (str): Path to the input CIF file.

    Returns:
        list: A list of chain IDs present in the file.
    """
    # Load the CIF file into PyMOL
    cmd.load(input_file, "structure")
    
    # Get all chain IDs
    chain_ids = cmd.get_chains("structure")
    
    # Clean up PyMOL session
    cmd.delete("all")
    
    return chain_ids


def extract_sequence_from_cif_by_chain(input_file, chain_id):
    """
    Extract the sequence of a specific chain from a CIF file and save it as a FASTA file.

    Args:
        input_file (str): Path to the input CIF file.
        chain_id (str): Chain ID to extract the sequence from (e.g., 'A').
    """
    # Load the CIF file into PyMOL
    cmd.load(input_file, "structure")
    
    # Generate the selection for the specified chain
    selection_name = f"chain_{chain_id}"
    cmd.select(selection_name, f"chain {chain_id}")
    
    # Extract the sequence in FASTA format
    fasta_sequence = cmd.get_fastastr(selection_name)
    
    # Clean up PyMOL session
    cmd.delete("all")

    # Return the sequence
    sequence = ''.join(fasta_sequence.split('\n')[1:])
    return sequence


def extract_sequences_from_cif(file_name):
    chain_ids = get_chain_ids(file_name)
    sequences = {}
    for chain_id in chain_ids:
        sequence = extract_sequence_from_cif_by_chain(file_name, chain_id)
        if chain_id == "_": # In SwissModel, the chain ID for ligands seems to be "_"
            continue
        if len(sequence) == 0:
            continue
        if sequence == "?":
            continue
        sequences[chain_id] = sequence
    return sequences


def find_subset_indices(seq1, seq2):
    """
    Find the 1-based indices of the first and last residues of seq2 in seq1.

    Parameters:
    - seq1 (str): The full protein sequence.
    - seq2 (str): The subset protein sequence.

    Returns:
    - tuple: (start_index, end_index), where both indices are 1-based.
    """
    # Find the starting index of seq2 in seq1
    start_index = seq1.find(seq2)
    
    if start_index == -1:
        raise ValueError("seq2 is not a subset of seq1.")
    
    # Calculate 1-based indices
    start_index_1based = start_index + 1
    end_index_1based = start_index + len(seq2)
    
    return (start_index_1based, end_index_1based)

def find_subset_indices_from_file(file_name):
    uniprot_ac = file_name.split("_")[1]
    with open(os.path.join(data_dir, "sequences", "fasta", "{0}.fasta".format(uniprot_ac)), "r") as f:
        fasta = f.readlines()
        seq1 = ''.join(fasta[1:]).replace("\n", "")
    sequences = extract_sequences_from_cif(os.path.join(processed_dir, "structures", uniprot_ac, file_name))
    seq2 = list(sequences.values())[0]
    start_resid, end_resid = find_subset_indices(seq1, seq2)
    data = {
        "file_name": file_name.replace(".cif", ".pdb"),
        "chain_id": get_chain_ids(file_name)[0],
        "uniprot_ac": uniprot_ac,
        "n_residues": len(seq1),
        "start_resid": start_resid,
        "end_resid": end_resid,
        "coverage": len(seq2)/len(seq1)*100,
        "structure_sequence_length": len(seq2),
        "full_sequence_length": len(seq1),
        "structure_sequence": seq2,
        "full_sequence": seq1,
    }
    return data

if __name__ == "__main__":
    file_name = "/Users/mduranfrigola/Documents/GitHub/mtb-targeted-protein-degradation/scripts/../processed/structures/P9WFS9/swissmodel_P9WFS9_model_1.pdb"
    data = find_subset_indices_from_file(file_name)
    print(data)