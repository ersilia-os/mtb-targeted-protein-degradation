"""Shared constants and helpers for docking hit analysis scripts (48_docking_hits_raw, 49_docking_hits_selective)."""
import os
import sys

import pandas as pd
from tqdm import tqdm

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

LIBRARIES = {
    "DL":   os.path.join(ROOT, "output", "unidock_docking",        "docking_results"),
    "REAL": os.path.join(ROOT, "output", "unidock_REAL_docking_2", "docking_results"),
}

REFERENCE_POCKET_CSV = os.path.join(ROOT, "output", "reference_pocket.csv")

SMILES_PATHS = {
    "DL":   os.path.join(ROOT, "data", "Enamine", "Enamine_Hit_Locator_Library_100K_Set_plated_100160cmpds_20250629.smiles"),
    "REAL": os.path.join(ROOT, "output", "unidock_REAL_docking", "inference_10B", "clustered_compounds.csv"),
}

SMILES_SEPS = {
    "DL":   "\t",
    "REAL": ",",
}

SMILES_COLS = {
    "DL":   ("Catalog ID", "SMILES"),
    "REAL": ("id",         "smiles"),
}


def lookup_smiles(compound_ids, lib):
    """Stream the library SMILES file in chunks, returning {id: smiles} for the requested ids."""
    path = SMILES_PATHS[lib]
    if not os.path.isfile(path):
        print(f"  Warning: SMILES file not found at {path}, omitting SMILES column.")
        return {}
    id_col, smi_col = SMILES_COLS[lib]
    sep = SMILES_SEPS[lib]
    needed = set(compound_ids)
    found = {}
    for chunk in pd.read_csv(path, sep=sep, usecols=[id_col, smi_col], chunksize=200_000):
        match = chunk[chunk[id_col].isin(needed)]
        found.update(zip(match[id_col], match[smi_col]))
        if len(found) == len(needed):
            break
    return found


def load_gene_map():
    path = os.path.join(ROOT, "data", "mtb_trna_synthetases_bosch_2021_fig5.csv")
    df = pd.read_csv(path)
    return {row["gene_name_in_bosch_2021"]: row["uniprot_ac"] for _, row in df.iterrows()}


def load_reference_pockets():
    if not os.path.isfile(REFERENCE_POCKET_CSV):
        print(f"Error: {REFERENCE_POCKET_CSV} not found.")
        print("Run script 47 and manually populate output/reference_pocket.csv")
        print("  (columns: gene_name, pocket_name)")
        sys.exit(1)
    df = pd.read_csv(REFERENCE_POCKET_CSV)
    return dict(zip(df["gene_name"], df["pocket_name"]))


def load_scores(pocket_path):
    """Return a pd.Series of scores indexed by compound id."""
    df = pd.read_csv(pocket_path)
    return df.set_index("compound")["score"]


def build_matrix(pocket_map, results_dir, label=""):
    """
    Build a raw-score DataFrame from a dict {column_label: pocket_name}.
    Rows = compounds, columns = column_labels. Missing cells → NaN.
    """
    series = {}
    for col_label, pocket_name in tqdm(pocket_map.items(), desc=label, unit="pocket"):
        report = os.path.join(results_dir, pocket_name, "report.csv")
        if not os.path.isfile(report):
            print(f"  Warning: report not found for pocket '{pocket_name}', skipping.")
            continue
        series[col_label] = load_scores(report)
    if not series:
        return pd.DataFrame()
    return pd.concat(series, axis=1)
