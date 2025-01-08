import os
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(root, "..", "data"))

df = pd.read_csv(os.path.join(data_dir, "mtb_trna_synthetases_bosch_2021_fig5_annotated.csv"))

uniprot_acs = list(df["uniprot_ac"])
names = list(df["gene_name_in_bosch_2021"])

name2prot = dict(zip(names, uniprot_acs))

print(name2prot)
print(len(name2prot))