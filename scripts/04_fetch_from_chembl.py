import os
import csv
import json
import pandas as pd
from chembl_webresource_client.new_client import new_client

root = os.path.dirname(os.path.abspath(__file__))

activity = new_client.activity

df = pd.read_csv(os.path.join(root, "..", "data", "mtb_trna_synthetases_bosch_2021_fig5_annotated.csv"))
uniprot_acs = list(df["uniprot_ac"])

available_uniprot_acs = {}
with open(os.path.join(root, "..", "data", "chembl_uniprot_mapping.txt"), "r") as f:
    reader = csv.reader(f, delimiter="\t")
    for r in reader:
        if r[0] in uniprot_acs:
            available_uniprot_acs[r[0]] = r[1]


for uniprot_ac, chembl_id in available_uniprot_acs.items():
    print(f"Fetching data for {uniprot_ac} ({chembl_id})")
    activities = activity.filter(target_chembl_id=chembl_id).filter(standard_type="IC50")
    activities = list(activities)
    file_name = os.path.join(root, "..", "data", "ligands", "chembl", f"{uniprot_ac}.json")
    with open(file_name, "w") as f:
        json.dump(activities, f, indent=4)