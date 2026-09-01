#!/usr/bin/env python3
"""
Prepares Nesso-1 inputs at the PROTEIN level, not the pocket/structure level (revised from this
script's original pocket-based design). Nesso-1 has no pocket-conditioning (see
docs/prediction.md's "Coming Soon: Pocket Conditioning" in the recursionpharma/nesso repo) and
takes only a sequence + ligand as input -- no 3D coordinates at all. This was empirically
confirmed to matter: an earlier pocket/structure-level run showed several differently-sourced
structures for the same gene (e.g. alaS's AlphaFold2/AlphaFold3/SwissModel models) producing
byte-identical Nesso-1 predictions, because Nesso-1 only ever sees whatever sequence was
extracted from each structure file, and multiple structure-prediction methods for the same
full-length, gap-free protein extract to the exact same sequence. Given that, pocket- or
structure-level granularity adds nothing Nesso-1 can act on -- the right unit is the protein
itself: all 21 Mtb tRNA synthetases (the project's full CRISPR-fitness-screen target set, see
CLAUDE.md), each tested against the full compound set. This also drops the 7K98 pheS+pheT dimer
entirely -- it isn't a named protein in its own right, and Nesso-1 can't be pointed at its
dimerization-interface pocket any more than it could at any other pocket.

Sequences come from data/mtb_trna_synthetases_bosch_2021_fig5_annotated.csv's own `sequence`
column (the canonical UniProt sequence used at the very start of this pipeline, script 00) rather
than from any structure file -- this sidesteps the exact ambiguity above (which structure-derived
variant to trust) by using the one sequence that was never subject to a homology model's
truncations or a particular structure-prediction method's residue coverage.

Usage:
    python 78_nesso1_prepare_inputs.py
"""
import os

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

PROTEINS_CSV = os.path.join(ROOT, "data", "mtb_trna_synthetases_bosch_2021_fig5_annotated.csv")

OUTPUT_DIR = os.path.join(ROOT, "output", "78_nesso1_prepare_inputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def main():
    df = pd.read_csv(PROTEINS_CSV)
    sub = df[["gene_name_in_bosch_2021", "uniprot_ac", "sequence"]].rename(
        columns={"gene_name_in_bosch_2021": "gene_name"}
    )

    sub["sequence_length"] = sub["sequence"].str.len()
    sub = sub.sort_values("gene_name").reset_index(drop=True)

    out_path = os.path.join(OUTPUT_DIR, "protein_sequences.csv")
    sub.to_csv(out_path, index=False)

    for _, row in sub.iterrows():
        print(f"  {row['gene_name']} ({row['uniprot_ac']}): seq_len={row['sequence_length']}")
    print(f"\nSaved {len(sub)} protein sequences to {out_path}")
    print("Compounds: reused directly from output/71_boltz2_prepare_inputs/compounds.csv, "
          "not copied.")


if __name__ == "__main__":
    main()
