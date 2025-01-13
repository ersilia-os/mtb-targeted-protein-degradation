import os
import numpy as np
import pandas as pd
import shutil
import gzip
import collections
from pymol import cmd

# Suppress executive details and warnings
cmd.feedback("disable", "executive", "details")
cmd.feedback("disable", "all", "warnings")

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
        for fn in os.listdir(os.path.join(dirname, l)):
            if fn.endswith(".cif.gz"):
                with gzip.open(os.path.join(dirname, l, fn), 'rb') as f_in:
                    with open(os.path.join(output_dir, "swissmodel_{0}_model_{1}.cif".format(uniprot_ac, i)), 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)


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
    cmd.load(input_cif)
    
    # Save as PDB
    cmd.save(output_pdb, "all")
    
    # Clean up
    cmd.delete("all")


print("Converting CIF to PDB")
for uniprot_ac in os.listdir(os.path.join(processed_dir, "structures")):
    for file_name in os.listdir(os.path.join(processed_dir, "structures", uniprot_ac)):
        if file_name.endswith(".cif"):
            print("Converting {0}".format(file_name))
            convert_cif_to_pdb(os.path.join(processed_dir, "structures", uniprot_ac, file_name), os.path.join(processed_dir, "structures", uniprot_ac, file_name.replace(".cif", ".pdb")))

print("Removing CIF files")
for uniprot_ac in os.listdir(os.path.join(processed_dir, "structures")):
    for file_name in os.listdir(os.path.join(processed_dir, "structures", uniprot_ac)):
        if file_name.endswith(".cif"):
            os.remove(os.path.join(processed_dir, "structures", uniprot_ac, file_name))

# Some functions

def get_sequence_from_uniprot(uniprot_ac):
    with open(os.path.join(data_dir, "sequences", "fasta", "{0}.fasta".format(uniprot_ac)), "r") as f:
        fasta = f.readlines()
        seq = ''.join(fasta[1:]).replace("\n", "")
        return seq


def get_chain_ids(input_file):
    """
    Extract all chain IDs from a PDB file.

    Args:
        input_file (str): Path to the input PDB file.

    Returns:
        list: A list of chain IDs present in the file.
    """
    # Load the PDB file into PyMOL
    cmd.load(input_file, "structure")
    
    # Get all chain IDs
    chain_ids = cmd.get_chains("structure")
    
    # Clean up PyMOL session
    cmd.delete("all")
    
    return chain_ids


def select_chain_id_by_coverage(input_file, coverage_threshold=0.8):
    chain_ids = get_chain_ids(input_file)
    uniprot_ac = input_file.split("/")[-2]
    full_seq = get_sequence_from_uniprot(uniprot_ac)
    coverages = {}
    for chain_id in chain_ids:
        seq = extract_sequence_from_pdb_by_chain(input_file, chain_id)
        seq = seq.replace("?", "")
        if seq not in full_seq:
            continue
        coverage = float(len(seq))/len(full_seq)
        if coverage < coverage_threshold:
            continue
        coverages[chain_id] = coverage
    if len(coverages) == 0:
        return None
    max_coverage = np.max([v for v in coverages.values()])
    if max_coverage < 0.95:
        raise Exception("Coverage is less than 1", input_file)

    possible_chains = [k for k,v in coverages.items() if v == max_coverage]
    possible_chains = sorted(possible_chains)
    return possible_chains[0]


def extract_sequence_from_pdb_by_chain(input_file, chain_id):
    """
    Extract the sequence of a specific chain from PDB file and save it as a FASTA file.

    Args:
        input_file (str): Path to the input CIF file.
        chain_id (str): Chain ID to extract the sequence from (e.g., 'A').
    """
    # Load the PDB file into PyMOL
    cmd.load(input_file, "structure")
    
    # Generate the selection for the specified chain
    selection_name = f"chain_{chain_id}"
    cmd.select(selection_name, f"chain {chain_id}")
    
    # Extract the sequence in FASTA format
    fasta_sequence = cmd.get_fastastr(selection_name)
    fasta_sequence = fasta_sequence.replace("?", "")
    
    # Clean up PyMOL session
    cmd.delete("all")

    # Return the sequence
    sequence = ''.join(fasta_sequence.split('\n')[1:])
    return sequence


def extract_sequences_from_pdb(input_file):
    chain_ids = get_chain_ids(input_file)
    sequences = {}
    for chain_id in chain_ids:
        sequence = extract_sequence_from_pdb_by_chain(input_file, chain_id)
        if chain_id == "_": # In SwissModel, the chain ID for ligands seems to be "_"
            continue
        if len(sequence) == 0:
            continue
        if sequence == "?":
            continue
        sequences[chain_id] = sequence
    return sequences


def extract_chain_in_file_inplace(input_file, chain_id):
    """
    Extract a single chain from a PDB file and save it to a new file.
    
    Args:
        input_file (str): Path to the input PDB file.
        chain_id (str): Chain ID to extract (e.g., 'A').
        output_file (str): Path to the output PDB file.
    """
    # Load the CIF file into PyMOL
    cmd.load(input_file, "structure")
    
    # Select the chain
    selection_name = f"chain_{chain_id}"
    cmd.select(selection_name, f"chain {chain_id}")
    
    # Save the selected chain to a new CIF file
    cmd.save(input_file, selection_name)
    
    # Clean up PyMOL session
    cmd.delete("all")
    print(f"Chain {chain_id} extracted and saved to {input_file}")


def get_maximum_coverage(file_name):
    uniprot_ac = file_name.split("_")[1]
    full_seq = get_sequence_from_uniprot(uniprot_ac)
    sequences = extract_sequences_from_pdb(os.path.join(processed_dir, "structures", uniprot_ac, file_name))
    max_coverage = 0
    for _, seq in sequences.items():
        if seq not in full_seq:
            print(file_name, "Sequence does not match the reference (it is not a subset)", seq)
            continue
        coverage = float(len(seq))/len(full_seq)
        if coverage > max_coverage:
            max_coverage = coverage
    return max_coverage


coverages = {}
for uniprot_ac in os.listdir(os.path.join(processed_dir, "structures")):
    for fn in os.listdir(os.path.join(processed_dir, "structures", uniprot_ac)):
        if fn.endswith(".pdb"):
            cov = get_maximum_coverage(fn)
            coverages[fn] = cov

print("Checking that all kept structures have a coverage of at least 95%")
for k, v in coverages.items():
    if v < 0.95:
        print("Coverage for {0} is less than 95%: {1}. Removing it".format(k, v))
        os.remove(os.path.join(processed_dir, "structures", k.split("_")[1], k)) 

print("Checking that all structures correspond to the same sequence")
uniprot_ac_seqs = collections.defaultdict(list)
uniprot_ac_seqs_with_origin = collections.defaultdict(list)
for uniprot_ac in uniprot_acs:
    output_dir = os.path.join(processed_dir, "structures", uniprot_ac)
    pdb_files = [f for f in os.listdir(output_dir) if f.endswith(".pdb")]
    for pdb_file in pdb_files:
        sequences = extract_sequences_from_pdb(os.path.join(output_dir, pdb_file))
        if len(sequences) > 1:
            lengths = set([len(x) for x in sequences.values()])
            max_length = np.max(list(lengths))
            chains = sorted(list(sequences.keys()))
            for chain_id in chains:
                if len(sequences[chain_id]) == max_length:
                    break
            extract_chain_in_file_inplace(os.path.join(output_dir, pdb_file), chain_id)
            sequences = extract_sequences_from_pdb(os.path.join(output_dir, pdb_file))
        if len(sequences) > 1:
            raise Exception("Different chain sequences found for {0}".format(pdb_file))
        sequence = list(sequences.values())[0]
        uniprot_ac_seqs[uniprot_ac] += [sequence]
        uniprot_ac_seqs_with_origin[(pdb_file, uniprot_ac)] = sequence


print("Checking that all structures correspond to the same sequence")
for k,v in uniprot_ac_seqs.items():
    if len(set(v)) > 1:
        print("Different sequences found for {0}".format(k))
        print([len(x) for x in v])
    else:
        print("All sequences for {0} are the same".format(k))


print("Getting sequence lengths for all proteins from UniProt")
fullseq_lengths = {}
for k, _ in uniprot_ac_seqs.items():
    seq = get_sequence_from_uniprot(k)
    fullseq_lengths[k] = len(seq)


def check_sequence_continuity(file_name):
    uniprot_ac = file_name.split("_")[1]
    sequences = extract_sequences_from_pdb(file_name)
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
    if len(v) != fullseq_lengths[uniprot_ac]:
        if len(v) < fullseq_lengths[uniprot_ac]*0.8:
            print("Sequence for {0} is too short. Removing it".format(k))
            os.remove(os.path.join(processed_dir, "structures", uniprot_ac, name))
            continue
        if not check_sequence_continuity(os.path.join(processed_dir, "structures", uniprot_ac, name)):
            print("Sequence for {0} is not continuous with respect to the reference. Removing it".format(k))
            os.remove(os.path.join(processed_dir, "structures", uniprot_ac, name))
            continue
        print("Different sequences found for {0}. Expected length: {1}. Actual length: {2}".format(k, fullseq_lengths[uniprot_ac], len(v)))
    file_names += [name]

print("These are the file names that will be kept")
print(file_names)


print("Selecting the chain with the highest coverage for each structure")
selected_chains = {}
for fn in file_names:
    uniprot_ac = fn.split("_")[1]
    selected_chains[fn] = select_chain_id_by_coverage(os.path.join(processed_dir, "structures", uniprot_ac, fn))

print(selected_chains)

print("Extracting the selected chains")
for fn in file_names:
    extract_chain_in_file_inplace(os.path.join(processed_dir, "structures", fn.split("_")[1], fn), selected_chains[fn])


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
    sequences = extract_sequences_from_pdb(os.path.join(processed_dir, "structures", uniprot_ac, file_name))
    seq2 = list(sequences.values())[0]
    start_resid, end_resid = find_subset_indices(seq1, seq2)
    data = {
        "file_name": file_name.split("/")[-1],
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

print("Finding the subset indices for the structures")
R = []
for fn in file_names:
    uniprot_ac = fn.split("_")[1]
    print("Finding subset indices for {0}".format(fn))
    data = find_subset_indices_from_file(os.path.join(root, "..", "processed", "structures", uniprot_ac, fn))
    R += [[data[k] for k in ["file_name", "chain_id", "uniprot_ac", "n_residues", "start_resid", "end_resid", "coverage", "structure_sequence_length", "full_sequence_length", "structure_sequence", "full_sequence"]]]
df = pd.DataFrame(R, columns=["file_name", "chain_id", "uniprot_ac", "n_residues", "start_resid", "end_resid", "coverage", "structure_sequence_length", "full_sequence_length", "sequence_structure", "full_sequence"])
df.sort_values(by=["uniprot_ac", "file_name"], inplace=True)

print("Saving the data lookup table")
df.to_csv(os.path.join(processed_dir, "trna_synthetases_data.csv"), index=False)

print("Removing non-sequence elements from the structures")
def remove_non_sequence_elements(pdb_file, chain_id, output_file):
    """
    Remove non-sequence elements (e.g., water, ligands, heteroatoms) from a specific chain in a PDB file.
    
    Parameters:
        pdb_file (str): Path to the input PDB file.
        chain_id (str): The chain ID to process.
        output_file (str): Path to save the cleaned PDB file.
    """
    # Load the PDB file
    cmd.load(pdb_file, "structure")
    
    # Select non-sequence elements in the specified chain
    cmd.select("non_sequence", f"chain {chain_id} and not polymer.protein")
    
    # Remove the non-sequence elements
    cmd.remove("non_sequence")
    
    # Save the cleaned structure
    cmd.save(output_file, f"chain {chain_id}")
    
    # Clean up PyMOL session
    cmd.delete("all")
    print(f"Cleaned structure saved to {output_file}")

for v in df[["file_name","chain_id"]].values:
    file_name = v[0]
    chain_id = v[1]
    uniprot_ac = file_name.split("_")[1]
    input_file = os.path.join(processed_dir, "structures", uniprot_ac, file_name)
    remove_non_sequence_elements(input_file, chain_id, input_file)
