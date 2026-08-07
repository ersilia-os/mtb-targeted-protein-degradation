import pandas as pd
import os

# Define paths
root = os.path.dirname(os.path.abspath(__file__))
outdir = os.path.abspath(os.path.join(root, "..", "output"))
path_to_sequence_info = os.path.abspath(os.path.join(root, "..", "data", "sequences", "interpro"))
path_to_sequence_summary = os.path.abspath(os.path.join(root, "..", "output", "sequences", "interpro_summary_curated.tsv"))

# Load pocket detection data
pocket_detection_df = pd.read_csv(os.path.join(outdir, "pocket_detection_data.csv"))
uniprots = sorted(set(pocket_detection_df['Uniprot AC']))

# Load all interpro IDs
df = pd.read_csv(path_to_sequence_summary, sep='\t')
interpro_ids = set(df['Accession'])
interpro_id_to_name = {i: j for i,j in zip(df['Accession'], df['Name'])}
interpro_id_to_curated = {i: j for i,j in zip(df['Accession'], df['Curated annotation'])}

uniprot_to_interpro_matches = {}
uniprot_to_interpro_names = {}
uniprot_to_plength = {}

# Load interpro data for each protein
for uni in uniprots:

    # Load TSV file
    df = pd.read_csv(os.path.join(path_to_sequence_info, f"entry-matching-{uni}.tsv"), sep="\t")
    
    # Select only interpro data
    df = df[df['Source Database'] == 'interpro']

    # Select only relevant columns
    df = df[['Accession', 'Name', "Source Database", 'Type', 'Protein Accession', 'Protein Length', 'Matches']]
    df['Protein Accession'] = df['Protein Accession'].str.upper()

    # Get id to matches
    # Select only those interpro that are relevant i.e. coverage < 0.60 in interpro_summary.tsv
    interpro_matches = {i: [k.split("..") for k in j.split(",")] for i,j in zip(df["Accession"], df["Matches"]) if i in interpro_ids}
    interpro_names = {i: j for i,j in zip(df["Accession"], df["Name"]) if i in interpro_ids}

    # Save
    uniprot_to_interpro_matches[uni] = interpro_matches
    uniprot_to_interpro_names[uni] = interpro_names
    uniprot_to_plength[uni] = int(df["Protein Length"].tolist()[0])

    del interpro_matches, interpro_names

interpro_ids = []
interpro_names = []
interpro_curated = []
interpro_matches = []
len_pocket = []
len_interpro = []
len_overlap = []
coverage_pocket = []
coverage_domain = []
overall_coverage_domain = []

pocket_detection_interpro_df = pd.DataFrame(columns=pocket_detection_df.columns).astype(pocket_detection_df.dtypes.to_dict())

# For each detected pocket
for pocket in pocket_detection_df.iterrows():

    # Get data from pocket
    uni = pocket[1]['Uniprot AC']
    pocket_residues = pocket[1]['Pocket residues (chain_resn)'].split()
    pocket_residues = set([int(i.split("_")[1]) for i in pocket_residues])

    # Append as many rows as needed in the new df
    rows = pd.DataFrame([pocket[1]] * len(uniprot_to_interpro_matches[uni]))
    pocket_detection_interpro_df = pd.concat([pocket_detection_interpro_df, rows], ignore_index=True)

    # For each interpro within the protein
    for interpro in sorted(uniprot_to_interpro_matches[uni]):

        # Store data from interpro domains
        interpro_ids.append(interpro)
        interpro_names.append(uniprot_to_interpro_names[uni][interpro])
        interpro_curated.append(interpro_id_to_curated[interpro])
        interpro_matches.append(uniprot_to_interpro_matches[uni][interpro])

        # Number of residues in the pocket
        len_pocket.append(len(pocket_residues))

        # Nomber of residues in the interpro domain
        dom_residues = set([int(res) for start, end in uniprot_to_interpro_matches[uni][interpro] 
                            for res in range(int(start), int(end) + 1)])
        len_interpro.append(len(dom_residues))

        # Check overlap
        len_overlap.append(len(dom_residues.intersection(pocket_residues)))
        dom_in_pocket = set([i for i in dom_residues if i in pocket_residues]) 
        pocket_in_dom = set([i for i in pocket_residues if i in dom_residues])

        # Save
        coverage_pocket.append(round(len(pocket_in_dom) / len(pocket_residues), 3))
        coverage_domain.append(round(len(dom_in_pocket) / len(dom_residues), 3))
        overall_coverage_domain.append(round(len(dom_residues) / uniprot_to_plength[uni], 3))

# Store data
pocket_detection_interpro_df["Interpro ID"] = interpro_ids
pocket_detection_interpro_df["Interpro name"] = interpro_names
pocket_detection_interpro_df["Interpro curated annotation"] = interpro_curated
pocket_detection_interpro_df["Interpro Matches"] = interpro_matches
pocket_detection_interpro_df["Residues in pocket"] = len_pocket
pocket_detection_interpro_df["Residues in interpro"] = len_interpro
pocket_detection_interpro_df["Residues overlap"] = len_overlap
pocket_detection_interpro_df["Coverage pocket"] = coverage_pocket
pocket_detection_interpro_df["Coverage domain"] = coverage_domain
pocket_detection_interpro_df["Overall coverage domain"] = overall_coverage_domain

# Omit those pairs having no overlapping residues 
pocket_detection_interpro_df = pocket_detection_interpro_df[pocket_detection_interpro_df["Residues overlap"] > 0].reset_index(drop=True)

# Save file
pocket_detection_interpro_df.to_csv(os.path.join(outdir, "pocket_detection_data_interpro.tsv"), sep='\t', index=False)