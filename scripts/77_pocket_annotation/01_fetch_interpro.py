import os
import sys
import json
import requests

root = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.join(root, "..", "..")
raw_dir = os.path.join(repo_root, "processed", "pocket_annotation", "interpro_raw")
os.makedirs(raw_dir, exist_ok=True)

UNIPROT_AC = sys.argv[1] if len(sys.argv) > 1 else "P9WFU3"  # default: pheS
BASE_URL = "https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{}/?page_size=200"


def fetch_all(uniprot_ac):
    url = BASE_URL.format(uniprot_ac)
    results = []
    while url:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        results.extend(data["results"])
        url = data.get("next")
    return results


if __name__ == "__main__":
    results = fetch_all(UNIPROT_AC)
    out_path = os.path.join(raw_dir, "{}_interpro.json".format(UNIPROT_AC))
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print("Fetched {} entries for {} -> {}".format(len(results), UNIPROT_AC, out_path))
