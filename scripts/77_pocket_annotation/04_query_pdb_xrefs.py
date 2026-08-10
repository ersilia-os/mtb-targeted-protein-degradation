import os
import json
import time
import requests
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.join(root, "..", "..")
output_dir = os.path.join(repo_root, "output", "77_pocket_annotation")
os.makedirs(output_dir, exist_ok=True)


if __name__ == "__main__":
    bosch = pd.read_csv(os.path.join(repo_root, "data", "mtb_trna_synthetases_bosch_2021_fig5_annotated.csv"))

    results = {}
    for _, row in bosch.iterrows():
        ac = row["uniprot_ac"]
        gene = row["gene_name_in_bosch_2021"]
        r = requests.get("https://rest.uniprot.org/uniprotkb/{}.json".format(ac), timeout=30)
        r.raise_for_status()
        d = r.json()
        pdb_ids = [x["id"] for x in d.get("uniProtKBCrossReferences", []) if x["database"] == "PDB"]
        results[ac] = {"gene": gene, "pdb_ids": pdb_ids}
        print("{:8s} {:10s} n_pdb={:2d}  {}".format(gene, ac, len(pdb_ids), pdb_ids))
        time.sleep(0.2)

    out_path = os.path.join(output_dir, "pdb_xrefs.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print("\nSaved PDB cross-references for {} proteins -> {}".format(len(results), out_path))
