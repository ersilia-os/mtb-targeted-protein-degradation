"""
Extended Data Table 1: one row per target protein (21 tRNA synthetases + gatA/gatB), with gene
name, UniProt AC, aaRS Class I/II designation, and full sequence.

Source: data/mtb_trna_synthetases_bosch_2021_fig5_annotated.csv (gene_name_in_bosch_2021,
uniprot_ac, sequence columns) - the same 21-target list used throughout the project. Class comes
from src/default.py's AARS_CLASS_LABELS (same mapping figure_1_plot.py/ExtendedDataFigure2.py use);
gatA/gatB are transamidases, not aaRS ligases, so they carry no Class I/II designation and are
labeled "N/A" here (user-confirmed) rather than being dropped, since all 21 targets must appear.

Rows are sorted alphabetically by gene name for a stable, citable order.

Usage:
    python ExtendedDataTable1.py
"""
import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "..", "src"))

import pandas as pd

from default import AARS_CLASS_LABELS

output_dir = os.path.join(root, "..", "..", "..", "output")
table_dir = os.path.join(output_dir, "plots", "FigSupp", "ExtendedDataTable1")
os.makedirs(table_dir, exist_ok=True)

SOURCE_CSV = os.path.join(root, "..", "..", "..", "data",
                           "mtb_trna_synthetases_bosch_2021_fig5_annotated.csv")


def main():
    source = pd.read_csv(SOURCE_CSV)
    assert len(source) == 21, f"Expected 21 target proteins in {SOURCE_CSV}, got {len(source)}."

    table = pd.DataFrame({
        "Gene name": source["gene_name_in_bosch_2021"],
        "Uniprot AC": source["uniprot_ac"],
        "Class": source["gene_name_in_bosch_2021"].map(AARS_CLASS_LABELS).fillna("N/A"),
        "Sequence": source["sequence"],
    })
    table = table.sort_values("Gene name").reset_index(drop=True)

    output_path = os.path.join(table_dir, "ExtendedDataTable1.csv")
    table.to_csv(output_path, index=False)
    print(f"Saved {len(table)} row(s) to {output_path}")
    print(table[["Gene name", "Uniprot AC", "Class"]].to_string(index=False))


if __name__ == "__main__":
    main()
