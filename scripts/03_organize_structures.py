import os
import pandas as pd
import shutil
import gzip
import collections
from pymol import cmd

root = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(root, "..", "data"))
processed_dir = os.path.abspath(os.path.join(root, "..", "processed"))

df = pd.read_csv(os.path.join(data_dir, "mtb_trna_synthetases_bosch_2021_fig5_annotated.csv"))

uniprot_acs = list(df["uniprot_ac"])
names = list(df["gene_name_in_bosch_2021"])

prot2name = dict(zip(uniprot_acs, names))

print("Getting AlphaFill data")
for uniprot_ac in uniprot_acs:
    dirname = os.path.join(data_dir, "structures", "alphafill_database", uniprot_ac)
    cif_file = os.path.join(dirname, uniprot_ac+".cif")
    output_dir = os.path.join(processed_dir, "structures", uniprot_ac)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    else:
        shutil.rmtree(output_dir)
        os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "alphafill_{0}_model_0.cif".format(uniprot_ac))
    shutil.copy(cif_file, output_file)

print("Getting AlphaFold2 data")
for uniprot_ac in uniprot_acs:
    dirname = os.path.join(data_dir, "structures", "alphafold2_database", uniprot_ac)
    i = 0
    for l in os.listdir(dirname):
        if l.endswith(".cif"):
            cif_file = os.path.join(dirname, l)
            output_dir = os.path.join(processed_dir, "structures", uniprot_ac)
            output_file = os.path.join(output_dir, "alphafold2_{0}_model_{1}.cif".format(uniprot_ac, i))
            shutil.copy(cif_file, output_file)
            i += 1

print("Getting AlphaFold3 data")
for uniprot_ac in uniprot_acs:
    dirname = os.path.join(data_dir, "structures", "alphafold3_webserver", "fold_{0}_mtb_trna_synthetase".format(prot2name[uniprot_ac]))
    print(dirname)
    for l in os.listdir(dirname):
        if l.endswith(".cif"):
            i = int(l.split("_model_")[1].split(".")[0])
            cif_file = os.path.join(dirname, l)
            output_dir = os.path.join(processed_dir, "structures", uniprot_ac)
            output_file = os.path.join(output_dir, "alphafold3_{0}_model_{1}.cif".format(uniprot_ac, i))
            shutil.copy(cif_file, output_file)

print("Getting Chai-1 data")
for uniprot_ac in uniprot_acs:
    dirname = os.path.join(data_dir, "structures", "chai1_server", "{0}_Mtb_tRNA_Synthetase".format(prot2name[uniprot_ac]))
    for l in os.listdir(dirname):
        if l.endswith(".cif"):
            i = int(l.split(".rank_")[1].split(".")[0])
            cif_file = os.path.join(dirname, l)
            output_dir = os.path.join(processed_dir, "structures", uniprot_ac)
            output_file = os.path.join(output_dir, "chai1_{0}_model_{1}.cif".format(uniprot_ac, i))
            shutil.copy(cif_file, output_file)

print("Getting SwissModel data")
for uniprot_ac in uniprot_acs:
    output_dir = os.path.join(processed_dir, "structures", uniprot_ac)
    chosen_name = None
    for name in os.listdir(os.path.join(data_dir, "structures", "swissmodel")):
        if "_{0}_".format(uniprot_ac) in name:
            chosen_name = name
            break
    dirname = os.path.join(data_dir, "structures", "swissmodel", chosen_name, "models")
    for l in os.listdir(dirname):
        if l.startswith("."):
            continue
        i = int(l)-1
        for file_name in os.listdir(os.path.join(dirname, l)):
            if file_name.endswith(".cif.gz"):
                with gzip.open(os.path.join(dirname, l, file_name), 'rb') as f_in:
                    with open(os.path.join(output_dir, "swissmodel_{0}_model_{1}.cif".format(uniprot_ac, i)), 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)


# Some functions

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


def extract_chain_in_file_inplace(file_name, chain_id):
    """
    Extract a single chain from a CIF file and save it to a new file.
    
    Args:
        input_file (str): Path to the input CIF file.
        chain_id (str): Chain ID to extract (e.g., 'A').
        output_file (str): Path to the output CIF file.
    """
    # Load the CIF file into PyMOL
    cmd.load(file_name, "structure")
    
    # Select the chain
    selection_name = f"chain_{chain_id}"
    cmd.select(selection_name, f"chain {chain_id}")
    
    # Save the selected chain to a new CIF file
    cmd.save(file_name, selection_name)
    
    # Clean up PyMOL session
    cmd.delete("all")
    print(f"Chain {chain_id} extracted and saved to {file_name}")


print("Checking that all structures correspond to the same sequence")
uniprot_ac_seqs = collections.defaultdict(list)
uniprot_ac_seqs_with_origin = collections.defaultdict(list)
for uniprot_ac in uniprot_acs:
    output_dir = os.path.join(processed_dir, "structures", uniprot_ac)
    cif_files = [f for f in os.listdir(output_dir) if f.endswith(".cif")]
    for cif_file in cif_files:
        sequences = extract_sequences_from_cif(os.path.join(output_dir, cif_file))
        if len(sequences) > 1:
            lengths = set([len(x) for x in sequences.values()])
            print(lengths)
            if len(lengths) > 1:
                raise Exception("Different chain sequences found for {0}".format(cif_file))
            extract_chain_in_file_inplace(os.path.join(output_dir, cif_file), list(sequences.keys())[0])
            sequences = extract_sequences_from_cif(os.path.join(output_dir, cif_file))
        if len(sequences) > 1:
            raise Exception("Different chain sequences found for {0}".format(cif_file))
        sequence = list(sequences.values())[0]
        uniprot_ac_seqs[uniprot_ac] += [sequence]
        uniprot_ac_seqs_with_origin[(cif_file, uniprot_ac)] = sequence


print("Checking that all structures correspond to the same sequence")
for k,v in uniprot_ac_seqs.items():
    if len(set(v)) > 1:
        print("Different sequences found for {0}".format(k))
        print([len(x) for x in v])
    else:
        print("All sequences for {0} are the same".format(k))

print("Debugging sequences that have a length different than what is expected")
modes = {}
for k, v in uniprot_ac_seqs.items():
    counts = collections.defaultdict(int)
    for x in v:
        counts[len(x)] += 1
    max_count = sorted(counts.items(), key=lambda x: x[1], reverse=True)[0][0]
    modes[k] = max_count

print("Expected lengths")
print(modes)

def check_sequence_continuity(file_name):
    uniprot_ac = file_name.split("_")[1]
    sequences = extract_sequences_from_cif(file_name)
    sequence = list(sequences.values())[0]
    with open(os.path.join(data_dir, "sequences", "fasta", "{0}.fasta".format(uniprot_ac)), "r") as f:
        fasta = f.readlines()
        reference_sequence = ''.join(fasta[1:]).replace("\n", "")
    if sequence in reference_sequence:
        return True
    else:
        return False

print("Finding proteins with a length different than the expected one")
file_names = []
for k, v in uniprot_ac_seqs_with_origin.items():
    name = k[0]
    uniprot_ac = k[1]
    if len(v) != modes[uniprot_ac]:
        if len(v) < modes[uniprot_ac]*0.8:
            print("Sequence for {0} is too short. Removing it".format(k))
            os.remove(os.path.join(processed_dir, "structures", uniprot_ac, name))
            continue
        if not check_sequence_continuity(os.path.join(processed_dir, "structures", uniprot_ac, name)):
            print("Sequence for {0} is not continuous with respect to the reference. Removing it".format(k))
            os.remove(os.path.join(processed_dir, "structures", uniprot_ac, name))
            continue
        print("Different sequences found for {0}. Expected length: {1}. Actual length: {2}".format(k, modes[uniprot_ac], len(v)))
    file_names += [name]


print("Converting CIF files to PDB files")

def convert_cif_to_pdb(input_cif, output_pdb):
    """
    Converts a CIF file to a PDB file using PyMOL.

    Parameters:
    - input_cif: Path to the input CIF file.
    - output_pdb: Path to the output PDB file.
    """
    # Initialize PyMOL
    cmd.reinitialize()
    
    # Load the CIF file
    cmd.load(input_cif, "structure")
    
    # Save as PDB
    cmd.save(output_pdb, "structure")
    
    # Clean up
    cmd.delete("all")

for file_name in file_names:
    uniprot_ac = file_name.split("_")[1]
    file_name = os.path.join(processed_dir, "structures", uniprot_ac, file_name)
    convert_cif_to_pdb(file_name, file_name.replace(".cif", ".pdb"))


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
        "uniprot_ac": uniprot_ac,
        "n_residues": len(seq1),
        "start_resid": start_resid,
        "end_resid": end_resid,
        "coverage": len(seq2)/len(seq1)*100,
        "structure_sequence": seq2,
        "full_sequence": seq1
    }
    return data

print("Finding the subset indices for the structures")
R = []
for file_name in file_names:
    uniprot_ac = file_name.split("_")[1]
    data = find_subset_indices_from_file(file_name)
    R += [[data[k] for k in ["file_name", "uniprot_ac", "n_residues", "start_resid", "end_resid", "coverage", "structure_sequence", "full_sequence"]]]
df = pd.DataFrame(R, columns=["file_name", "uniprot_ac", "n_residues", "start_resid", "end_resid", "coverage", "sequence_structure", "full_sequence"])

print("Saving the data lookup table")
df.to_csv(os.path.join(processed_dir, "trna_synthetases_data.csv"), index=False)