#!/usr/bin/env python3
"""
Prepares Nesso-1 inputs for the human tRNA synthetase counter-screen: reads all 38 reviewed
human aminoacyl-tRNA synthetases (data/human_trna_synthetases_uniprot.csv, fetched from UniProt
this session -- query `ec:6.1.1.* AND organism_id:9606 AND reviewed:true`, verified one-by-one
against known human aaRS biology) and writes them out in the same shape script 78 uses for the
Mtb targets, so scripts 85-89 can be near-identical copies of 79-83.

No filtering -- all 38 genes, per the user's explicit scope decision, not just the direct
orthologs of the 5 Mtb targets. This is a separate, parallel pipeline (scripts 84-89, own output
paths) so it never touches the Mtb screen's (scripts 78-83) inputs/outputs, which may still be
running.

Usage:
    python 84_nesso1_human_prepare_inputs.py
"""
import os

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

HUMAN_AARS_CSV = os.path.join(ROOT, "data", "human_trna_synthetases_uniprot.csv")

OUTPUT_DIR = os.path.join(ROOT, "output", "84_nesso1_human_prepare_inputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def main():
    df = pd.read_csv(HUMAN_AARS_CSV)
    df = df[["gene_name", "uniprot_ac", "protein_name", "ec_number", "sequence", "sequence_length"]]
    df = df.sort_values("gene_name").reset_index(drop=True)

    out_path = os.path.join(OUTPUT_DIR, "protein_sequences.csv")
    df.to_csv(out_path, index=False)

    for _, row in df.iterrows():
        print(f"  {row['gene_name']} ({row['uniprot_ac']}): seq_len={row['sequence_length']}")
    print(f"\nSaved {len(df)} protein sequences to {out_path}")
    print("Compounds: reused directly from output/71_boltz2_prepare_inputs/compounds.csv, "
          "not copied.")


if __name__ == "__main__":
    main()
