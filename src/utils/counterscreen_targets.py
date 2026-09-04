"""Shared gene-list loader for the structure-based counter-screen pipeline (scripts 90-97),
which runs the exact same AF2-monomer-only, no-curation recipe against two organisms' aaRS
targets: the 38 human genes (the original counter-screen) and, in parallel, all 21 Mtb genes
(scripts 90-97's own --organism flag)."""
import os

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")

HUMAN_AARS_CSV = os.path.join(ROOT, "data", "human_trna_synthetases_uniprot.csv")
MTB_AARS_CSV = os.path.join(ROOT, "data", "mtb_trna_synthetases_bosch_2021_fig5_annotated.csv")

TARGET_COLUMNS = ["gene_name", "uniprot_ac", "sequence_length"]


def load_targets(organism):
    """Returns a DataFrame with columns gene_name, uniprot_ac, sequence_length for the given
    organism ("human" or "mtb")."""
    if organism == "human":
        df = pd.read_csv(HUMAN_AARS_CSV)
    elif organism == "mtb":
        df = pd.read_csv(MTB_AARS_CSV)
        df = df.rename(columns={"gene_name_in_bosch_2021": "gene_name"})
        df["sequence_length"] = df["sequence"].str.len()
    else:
        raise ValueError(f"Unknown organism: {organism!r} (expected 'human' or 'mtb')")
    return df[TARGET_COLUMNS]
