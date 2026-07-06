import argparse
import os
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument("--raw",       required=True, help="Script 48 multi-target CSV")
parser.add_argument("--selective", required=True, help="Script 49 selective hits CSV")
args = parser.parse_args()

eos_dir    = os.path.join(root, "..", "output", "50_ersilia_eos42ez")
output_dir = os.path.join(root, "..", "output", "51_filtered_hits")
os.makedirs(output_dir, exist_ok=True)

CYTOTOX_COLS = ["cytotoxicity_hepg2", "cytotoxicity_hskmc", "cytotoxicity_imr90"]


def _read_eos42ez(path):
    """Read eos42ez CSV robustly: SMILES (col 1) may contain commas (extended SMILES notation).
    Structure is always: key, smiles, cytotox_hepg2, cytotox_hskmc, cytotox_imr90."""
    rows = []
    with open(path) as f:
        header = next(f).strip().split(",")
        for line in f:
            parts = line.strip().split(",")
            key = parts[0]
            smiles = ",".join(parts[1:-3])
            rows.append([key, smiles] + parts[-3:])
    df = pd.DataFrame(rows, columns=header)
    df[header[2:]] = df[header[2:]].astype(float)
    return df


# --- Load QED and is_pains from source CSVs ---
df48 = pd.read_csv(args.raw)[["smiles", "QED", "is_pains"]]

df49_raw = pd.read_csv(args.selective)
df49 = df49_raw[df49_raw["selected"].notna()][["smiles", "QED", "is_pains"]]

props = (
    pd.concat([df48, df49], ignore_index=True)
    .drop_duplicates(subset="smiles", keep="first")
    .set_index("smiles")
)

# --- Load cytotox predictions ---
eos = _read_eos42ez(os.path.join(eos_dir, "eos42ez.csv"))
eos = eos.rename(columns={"input": "smiles"})

# --- Merge ---
df = eos.join(props, on="smiles", how="left")

# --- Filter with reporting ---
print(f"Starting compounds:             {len(df)}")

df1 = df[df["QED"] > 0.5]
print(f"After QED > 0.5:                {len(df1)}")

df2 = df1[~df1["is_pains"].astype(bool)]
print(f"After PAINS removal:            {len(df2)}")

df3 = df2[(df2[CYTOTOX_COLS] < 0.3).all(axis=1)]
print(f"After cytotoxicity (all < 0.3): {len(df3)}")

out_cols = ["key", "smiles"] + CYTOTOX_COLS + ["QED", "is_pains"]
df3[out_cols].to_csv(os.path.join(output_dir, "filtered_hits.csv"), index=False)
print(f"Written to output/51_filtered_hits/filtered_hits.csv")
