import os
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))

df = pd.read_csv(os.path.join("..", "processed/trna_synthetases_data.csv"))

uniprot_acs = sorted(set(df["uniprot_ac"].tolist()))

R = []
for uniprot_ac in uniprot_acs:
    di = pd.read_csv(os.path.join(root, "..", "data", "sequences", "interpro", "entry-matching-{}.tsv".format(uniprot_ac)), sep="\t")
    print(di.columns)
    di = di[["Accession", "Name", "Source Database", "Type", "Integrated Into", "Integrated Signatures", "GO Terms", "Protein Length", "Matches"]]
    columns = ["uniprot_ac", "interpro_ac", "name", "source_database", "type", "integrated_into", "integrated_signatures", "go_terms", "protein_length", "matches"]
    for i, row in di.iterrows():
        R += [[uniprot_ac] + row.tolist()]

df = pd.DataFrame(R, columns=columns)
df.to_csv(os.path.join(root, "..", "processed", "sequences", "interpro_data.csv"), index=False)