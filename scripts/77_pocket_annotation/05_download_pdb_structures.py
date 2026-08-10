import os
import json
import requests

root = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.join(root, "..", "..")
output_dir = os.path.join(repo_root, "output", "77_pocket_annotation")
pdb_dir = os.path.join(repo_root, "processed", "pocket_annotation", "pdb_structures")
os.makedirs(pdb_dir, exist_ok=True)

# A legacy-PDB-format 404 comes back as a small HTML error page, not a real
# structure file -- anything under this size is treated as "not available",
# triggering the mmCIF fallback (RCSB retires legacy PDB format for very
# large/complex depositions).
MIN_VALID_PDB_SIZE = 2000


def download(pdb_id, fmt):
    url = "https://files.rcsb.org/download/{}.{}".format(pdb_id, fmt)
    out_path = os.path.join(pdb_dir, "{}.{}".format(pdb_id, fmt))
    resp = requests.get(url, timeout=30)
    with open(out_path, "wb") as f:
        f.write(resp.content)
    return out_path, len(resp.content)


if __name__ == "__main__":
    xrefs = json.load(open(os.path.join(output_dir, "pdb_xrefs.json")))
    all_ids = sorted({pid for info in xrefs.values() for pid in info["pdb_ids"]})
    print("{} distinct PDB IDs to download".format(len(all_ids)))

    for pdb_id in all_ids:
        pdb_path, size = download(pdb_id, "pdb")
        if size < MIN_VALID_PDB_SIZE:
            os.remove(pdb_path)
            cif_path, cif_size = download(pdb_id, "cif")
            print("{}: legacy PDB unavailable, using mmCIF ({} bytes)".format(pdb_id, cif_size))
        else:
            print("{}: pdb OK ({} bytes)".format(pdb_id, size))
