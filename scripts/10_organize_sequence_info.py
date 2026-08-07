import numpy as np
import pandas as pd
import os

# Define paths
root = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.abspath(os.path.join(root, "..", "data", "sequences", "interpro"))
outdir = os.path.abspath(os.path.join(root, "..", "output", "sequences"))


# Define uniprots
uniprots = [i.split(".")[0].split("-")[2] for i in sorted(os.listdir(data_path))]
uniprots_to_number = {i: c for c,i in enumerate(sorted(uniprots))}

dfs = []

for uni in uniprots:

    # Load TSV file
    df = pd.read_csv(os.path.join(data_path, f"entry-matching-{uni}.tsv"), sep="\t")
    
    # Select only interpro data
    df = df[df['Source Database'] == 'interpro']

    # Select only relevant columns
    df = df[['Accession', 'Name', "Source Database", 'Type', 'Protein Accession', 'Protein Length', 'Matches']]
    df['Protein Accession'] = df['Protein Accession'].str.upper()

    # Append to dfs
    dfs.append(df)

# Merge everything
dfs = pd.concat(dfs, ignore_index=True)

# Create dict - from accession to name and type
acc_to_name_type = {i: [j, k] for i,j,k in zip(dfs['Accession'], dfs['Name'], dfs['Type'])}

# Create dict - from accession to protein length and matches
# acc_to_plength_matches = {i: {} for i in sorted(set(dfs['Accession']))}

interpro_data = []
interpro_data_collapsed = []  # Do we need it?

for interpro_id in sorted(acc_to_name_type):

    # Store data
    name, type_ = acc_to_name_type[interpro_id]
    proteins = ",".join(dfs[dfs['Accession'] == interpro_id]['Protein Accession'].tolist())

    # For each protein, get coverage
    coverages = []
    for protein in proteins.split(","):

        # If len data_protein > 1; error
        data_protein = dfs[(dfs['Accession'] == interpro_id) & (dfs['Protein Accession'] == protein)]
        if len(data_protein) != 1:
            raise Exception("Something did not go as expected. Please check")

        # Get protein length and matches
        protein_length = int(data_protein["Protein Length"].tolist()[0])
        matches = data_protein["Matches"].tolist()[0]
        matches = [i.split("..") for i in matches.split(",")]
        matches = sum([int(i[1]) - int(i[0]) + 1 for i in matches])  # Why +1? If [1,2] --> 2 residues
        coverages.append(round(matches / protein_length, 3))

    # Create binary vector
    vec = np.zeros(len(uniprots), dtype=int)
    ind = [uniprots_to_number[i] for i in sorted(uniprots_to_number) if i in proteins]
    vec[ind] = 1

    # Append to final df
    cols = [interpro_id, name, type_, len(proteins.split(",")), round(np.mean(coverages), 3), proteins, ",".join([str(i) for i in coverages])]
    cols.extend(vec)
    interpro_data.append(cols)

# Create dataframe
interpro_data = pd.DataFrame(interpro_data, columns=['Accession', 'Name', 'Type', 'Number of proteins', 'Average Coverage', 'Proteins', 'Protein Coverages'] + sorted(uniprots))
interpro_data = interpro_data.sort_values('Number of proteins', ascending=False).reset_index(drop=True)

# Save results
interpro_data.to_csv(os.path.join(outdir, 'interpro_summary.tsv'), sep='\t', index=False)


### CURATE INTERPRO ANNOTATIONS ###

def curate_InterPro(name):
    """
    Classifies an InterPro domain name into functional categories with refined classification.

    Parameters:
    name (str): The name of the InterPro domain to be classified.

    Returns:
    str: The category of the domain (e.g., "Catalytic Domain", "Editing Domain", "Other too broad/unspecified functional entities").
    """
    
    # Convert name to lowercase for case-insensitive matching
    name = name.lower()

    # Categorization logic with refined classification
    if "rossmann" in name or "atp" in name or "catalytic" in name or "core domain" in name:
        return "Catalytic Domain (ATP Binding Site)"
    
    elif "anticodon-binding" in name or "anticodon binding" in name or "anti-codon-binding" in name:
        return "Anticodon Binding Domain"
    
    elif "editing" in name:
        return "Editing Domain"
    
    elif "trna-binding" in name or "ob-fold" in name:
        return "tRNA Binding Domain"
    
    elif "rna-binding" in name or "s4 domain" in name:
        return "RNA Binding Domain"
    
    elif "synthetase" in name or "ligase" in name or "amidotransferase" in name or "synthetase-associated" in name or "amidase" in name or "b5-domain" in name or "beta-barrel" in name or "conserved site" in name:
        return "Other too broad/unspecified functional entities"
    
    else:
        return "Other"
    

# root = os.path.dirname(os.path.abspath(__file__))
root = "/home/acomajuncosa/Documents/mtb-targeted-protein-degradation"
df = pd.read_csv(os.path.join(root, "output", "sequences", "interpro_summary.tsv"), sep='\t')

# Remove those with coverage > 60
df = df[df['Average Coverage'] < 0.60].reset_index(drop=True)

# Curate names
curated_annotation = df['Name'].tolist()
curated_annotation = [curate_InterPro(i) for i in curated_annotation]
df.insert(2, column='Curated annotation', value=curated_annotation)

# Save results
df.to_csv(os.path.join(root, "output", "sequences", "interpro_summary_curated.tsv"), sep='\t', index=False)

labels = sorted(set(df['Curated annotation']))
for label in labels:
    names = df[df['Curated annotation'] == label]['Name']
    print(label + " (" + str(len(names)) + ")")
    for name in names:
        print("   -" + name)
    print("\n")