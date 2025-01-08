import os
import pandas as pd
import shutil
import gzip
import collections

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
        i = int(l)-1
        for file_name in os.listdir(os.path.join(dirname, l)):
            if file_name.endswith(".cif.gz"):
                with gzip.open(os.path.join(dirname, l, file_name), 'rb') as f_in:
                    with open(os.path.join(output_dir, "swissmodel_{0}_model_{1}.cif".format(uniprot_ac, i)), 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)


from Bio.PDB.MMCIFParser import MMCIFParser

def extract_sequence_from_cif(cif_file_path):
    """
    Extracts the amino acid sequence from a CIF file.
    
    Args:
        cif_file_path (str): Path to the CIF file.
    
    Returns:
        dict: A dictionary where keys are chain IDs and values are sequences.
    """
    if "swissprot" in cif_file_path:
        return []
    parser = MMCIFParser(QUIET=False)
    print("Extracting sequence from {0}".format(cif_file_path))
    try:
        structure = parser.get_structure("structure_id", cif_file_path)
        sequences = []
        for model in structure:
            for chain in model:
                residues = [
                    residue.resname for residue in chain.get_residues()
                    if residue.id[0] == " "  # Exclude heteroatoms and water
                ]
                sequence = "".join(residues)
                print("Chain {0}: {1}".format(chain.id, sequence))
                sequences += [sequence]
    except:
        sequence = ""
        with open(cif_file_path, "r") as file:
            for line in file:
                if "_entity_poly.pdbx_seq_one_letter_code" in line:
                    sequence += line.split()[-1]
        sequences = [sequence.strip()]
    sequences = [x for x in sequences if x]
    if len(sequences) > 1:
        raise Exception("Multiple sequences found in {0}".format(cif_file_path))
    return sequences


print("Checking that all structures correspond to the same sequence")
uniprot_ac_seqs = collections.defaultdict(list)
for uniprot_ac in uniprot_acs:
    output_dir = os.path.join(processed_dir, "structures", uniprot_ac)
    cif_files = [f for f in os.listdir(output_dir) if f.endswith(".cif")]
    for cif_file in cif_files:
        try:
            sequences = extract_sequence_from_cif(os.path.join(output_dir, cif_file))
            uniprot_ac_seqs[uniprot_ac] += sequences
            if len(set(sequences)) > 1:
                print("Different sequences found in {0}".format(cif_file))
                print(sequences)
            else:
                print("All sequences in {0} are the same".format(cif_file))
        except Exception as e:
            pass


print("Checking that all structures correspond to the same sequence")
for k,v in uniprot_ac_seqs.items():
    if len(set(v)) > 1:
        print("Different sequences found for {0}".format(k))
        print([len(x) for x in v])
    else:
        print("All sequences for {0} are the same".format(k))
