import numpy as np
import pandas as pd
import os

# Define paths
root = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.abspath(os.path.join(root, "..", "data", "sequences", "interpro"))
outdir = os.path.abspath(os.path.join(root, "..", "processed", "sequences"))

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

interpro_data = []
interpro_data_collapsed = []  # Do we need it?

for interpro_id in sorted(acc_to_name_type):

    # Store data
    name, type_ = acc_to_name_type[interpro_id]
    proteins = ",".join(dfs[dfs['Accession'] == interpro_id]['Protein Accession'].tolist())

    # Create binary vector
    vec = np.zeros(len(uniprots), dtype=int)
    ind = [uniprots_to_number[i] for i in sorted(uniprots_to_number) if i in proteins]
    vec[ind] = 1

    # Append to final df
    cols = [interpro_id, name, type_, len(proteins.split(",")), proteins]
    cols.extend(vec)
    interpro_data.append(cols)

# Create dataframe
interpro_data = pd.DataFrame(interpro_data, columns=['Accession', 'Name', 'Type', 'Number of proteins', 'Proteins'] + sorted(uniprots))
interpro_data = interpro_data.sort_values('Number of proteins', ascending=False).reset_index(drop=True)

# Save results
interpro_data.to_csv(os.path.join(outdir, 'interpro_summary.tsv'), sep='\t', index=False)