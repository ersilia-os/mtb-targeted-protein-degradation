import os
import sys
import json
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.join(root, "..", "..")
raw_dir = os.path.join(repo_root, "processed", "pocket_annotation", "interpro_raw")
output_dir = os.path.join(repo_root, "output", "77_pocket_annotation")
os.makedirs(output_dir, exist_ok=True)

UNIPROT_AC = sys.argv[1] if len(sys.argv) > 1 else "P9WFU3"  # default: pheS


def locations_str(entry_protein_locations):
    ranges = []
    for loc in entry_protein_locations:
        for frag in loc["fragments"]:
            ranges.append((frag["start"], frag["end"]))
    ranges = sorted(set(ranges))
    return ranges


def range_str(ranges):
    return ",".join("{}-{}".format(s, e) for s, e in ranges)


def total_coverage(ranges, protein_length):
    covered = sum(e - s + 1 for s, e in ranges)
    return covered / protein_length


def go_terms_str(go_terms):
    if not go_terms:
        return ""
    return ";".join("{}|{}".format(g["identifier"], g["name"]) for g in go_terms)


if __name__ == "__main__":
    raw_path = os.path.join(raw_dir, "{}_interpro.json".format(UNIPROT_AC))
    with open(raw_path) as f:
        results = json.load(f)

    parsed = {}
    for r in results:
        m = r["metadata"]
        prot = r["proteins"][0]
        ranges = locations_str(prot["entry_protein_locations"])
        parsed[m["accession"]] = {
            "accession": m["accession"],
            "name": m["name"],
            "source_database": m["source_database"],
            "type": m["type"],
            "integrated_into": m["integrated"],
            "go_terms": go_terms_str(m["go_terms"]),
            "start": ranges[0][0],
            "end": ranges[-1][1],
            "ranges": range_str(ranges),
            "protein_length": prot["protein_length"],
            "coverage": round(total_coverage(ranges, prot["protein_length"]), 3),
            "member_db_accession": [],
            "member_db_name": [],
        }

    rows = []
    for acc, entry in parsed.items():
        if entry["integrated_into"]:
            parent_acc = entry["integrated_into"]
            if parent_acc in parsed:
                parsed[parent_acc]["member_db_accession"].append(acc)
                parsed[parent_acc]["member_db_name"].append(entry["name"])
                continue
        rows.append(entry)

    for entry in rows:
        entry["member_db_accession"] = ";".join(entry["member_db_accession"])
        entry["member_db_name"] = ";".join(n or "" for n in entry["member_db_name"])
        entry["name"] = entry["name"] or ""
        entry["specificity"] = "specific" if entry["type"] in {"domain", "family", "site"} else "structural/broad"

    columns = [
        "accession", "name", "source_database", "type", "specificity",
        "start", "end", "ranges", "coverage", "protein_length",
        "member_db_accession", "member_db_name", "go_terms",
    ]
    df = pd.DataFrame(rows)[columns].sort_values("start").reset_index(drop=True)

    out_path = os.path.join(output_dir, "{}_annotation_table.csv".format(UNIPROT_AC))
    df.to_csv(out_path, index=False)
    print("Built annotation table with {} rows (from {} raw InterPro entries) -> {}".format(
        len(df), len(results), out_path))
    print(df[["accession", "name", "source_database", "type", "specificity", "start", "end", "coverage"]].to_string(index=False))
