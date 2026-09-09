"""
Extended Data Table 1: one row per target protein (21 tRNA synthetases + gatA/gatB) - identity
(gene name, UniProt AC, protein name, EC number, aaRS Class I/II), structural coverage (number of
structures, number of raw pocket detections, number of deduplicated canonical pockets), the same
4-domain presence flags and AlphaFold2 pLDDT confidence breakdown as ExtendedDataFigure1.py's
panels a/b, and full sequence.

Sources, all already computed elsewhere in the project - this script joins them, it doesn't
recompute anything:
- data/mtb_trna_synthetases_bosch_2021_fig5_annotated.csv: gene_name_in_bosch_2021, uniprot_ac,
  protein_names, ec_number, sequence - the same 21-target list used throughout the project.
- src/default.py's AARS_CLASS_LABELS: Class I/II (gatA/gatB are transamidases, not aaRS ligases,
  so they carry no Class I/II designation - labeled "N/A" here, user-confirmed, rather than being
  dropped, since all 21 targets must appear).
- output/plots/FigSupp/ExtendedDataFigure1/ExtendedDataFigure1_domain_presence.csv: Catalytic/
  tRNA binding/Editing/Anticodon presence (transposed here into 4 bool columns per gene).
- .../ExtendedDataFigure1_plddt_fractions.csv: per-residue AlphaFold2 pLDDT class fractions,
  reported here as percentages (rounded to 1 decimal).
- .../ExtendedDataFigure1_structure_counts.csv: number of structures per source (AlphaFold2/
  AlphaFold3/Chai-1/SwissModel), summed here into one "Number of structures" total. AlphaFill is
  deliberately NOT counted (matches Extended Data Figure 1, user-confirmed) - every protein has
  exactly one AlphaFill structure in output/trna_synthetases_data.csv's fuller per-structure table,
  but it's the same coordinates as that protein's own AlphaFold2 model with a homology-transplanted
  ligand grafted on (03_align_structures.py skips aligning it for exactly this redundancy), and it
  never gets its own independent pocket-detection/docking run - so "Number of structures" here
  means structurally-independent models actually carried through the pipeline, one less per
  protein than the full organized-structure count.
- output/pocket_detection_data.csv: one row per detected pocket per structure - grouped by
  UniProt AC for "Number of pocket structures" (raw, pre-dedup detections - figure_1_plot.py's own
  term for this exact count, e.g. "the total pocket structures for lysS").
- output/plots/figure_1/color_mapping.json's gene_to_unique_pocket_count: the same 6.14 A
  greedy-centroid-dedup canonical pocket count figure_1_calculations.py computes and every other
  figure in this project reuses (not recomputed here).

Rows are sorted alphabetically by gene name for a stable, citable order.

Usage:
    python ExtendedDataTable1.py
"""
import json
import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "..", "src"))

import pandas as pd

from default import AARS_CLASS_LABELS
from xlsx_utils import save_table_with_legend

output_dir = os.path.join(root, "..", "..", "..", "output")
table_dir = os.path.join(output_dir, "plots", "FigSupp", "ExtendedDataTable1")
os.makedirs(table_dir, exist_ok=True)

SOURCE_CSV = os.path.join(root, "..", "..", "..", "data",
                           "mtb_trna_synthetases_bosch_2021_fig5_annotated.csv")
COLOR_MAPPING_JSON = os.path.join(output_dir, "plots", "figure_1", "color_mapping.json")
FIG1_SUPP_DIR = os.path.join(output_dir, "plots", "FigSupp", "ExtendedDataFigure1")
DOMAIN_PRESENCE_CSV = os.path.join(FIG1_SUPP_DIR, "ExtendedDataFigure1_domain_presence.csv")
PLDDT_FRACTIONS_CSV = os.path.join(FIG1_SUPP_DIR, "ExtendedDataFigure1_plddt_fractions.csv")
STRUCTURE_COUNTS_CSV = os.path.join(FIG1_SUPP_DIR, "ExtendedDataFigure1_structure_counts.csv")
POCKET_DETECTION_CSV = os.path.join(output_dir, "pocket_detection_data.csv")

STRUCTURE_SOURCES = ["alphafold2", "alphafold3", "chai1", "swissmodel"]
DOMAIN_COLUMNS = ["Catalytic", "tRNA binding", "Editing", "Anticodon"]
PLDDT_COLUMNS = ["Very high (>90)", "Confident (70-90)", "Low (50-70)", "Very low (<50)"]


def build_legend():
    explanations = {
        "Gene name": "Gene symbol used throughout this project.",
        "Uniprot AC": "UniProt accession code for the M. tuberculosis protein.",
        "Protein name": "Full protein name(s), from UniProt (data/mtb_trna_synthetases_bosch_2021_fig5_annotated.csv).",
        "EC number": "Enzyme Commission number.",
        "Class": "aaRS Class I or Class II (fold/active-site architecture). N/A for gatA/gatB, "
                 "which are transamidases, not aaRS ligases.",
        "Number of structures": "Number of structurally-independent models used downstream (pocket "
                                 "detection, docking) - AlphaFold2, AlphaFold3, Chai-1, SwissModel. "
                                 "Excludes the redundant AlphaFill structure every protein also has "
                                 "(same coordinates as its own AlphaFold2 model, not independently "
                                 "carried through the pipeline).",
        "Number of pocket structures": "Total raw P2Rank pocket detections across all of this "
                                        "protein's structures (one count per pocket per structure, "
                                        "before deduplication).",
        "Number of canonical pockets": "Number of physically distinct pockets after deduplicating "
                                        "raw detections (6.14 A greedy centroid clustering across "
                                        "this protein's aligned structures).",
        "Catalytic": "Whether any of this protein's pockets carry the curated Catalytic Domain "
                     "(ATP/ligase) label.",
        "tRNA binding": "Whether any of this protein's pockets carry the curated tRNA Binding "
                        "Domain label.",
        "Editing": "Whether any of this protein's pockets carry the curated Editing Domain label.",
        "Anticodon": "Whether any of this protein's pockets carry the curated Anticodon Binding "
                     "Domain label.",
        "pLDDT Very high (>90) (%)": "Percent of residues with AlphaFold2 pLDDT > 90 (very high "
                                      "confidence).",
        "pLDDT Confident (70-90) (%)": "Percent of residues with AlphaFold2 pLDDT in [70, 90) "
                                        "(confident).",
        "pLDDT Low (50-70) (%)": "Percent of residues with AlphaFold2 pLDDT in [50, 70) (low "
                                  "confidence).",
        "pLDDT Very low (<50) (%)": "Percent of residues with AlphaFold2 pLDDT < 50 (very low "
                                     "confidence).",
        "Sequence": "Full amino acid sequence.",
    }
    return explanations


def main():
    source = pd.read_csv(SOURCE_CSV)
    assert len(source) == 21, f"Expected 21 target proteins in {SOURCE_CSV}, got {len(source)}."

    with open(COLOR_MAPPING_JSON) as f:
        color_mapping = json.load(f)
    uniprot_to_gene = color_mapping["uniprot_to_gene"]
    gene_to_unique_pocket_count = color_mapping["gene_to_unique_pocket_count"]

    domain_presence = pd.read_csv(DOMAIN_PRESENCE_CSV, index_col=0).T.astype(bool)
    plddt_fractions = pd.read_csv(PLDDT_FRACTIONS_CSV, index_col=0)
    structure_counts = pd.read_csv(STRUCTURE_COUNTS_CSV, index_col=0)
    n_structures = structure_counts[STRUCTURE_SOURCES].sum(axis=1)

    pocket_detection = pd.read_csv(POCKET_DETECTION_CSV)
    pocket_detection["gene"] = pocket_detection["Uniprot AC"].map(uniprot_to_gene)
    n_pocket_structures = pocket_detection.groupby("gene").size()

    genes = source["gene_name_in_bosch_2021"]
    table = pd.DataFrame({
        "Gene name": genes,
        "Uniprot AC": source["uniprot_ac"],
        "Protein name": source["protein_names"],
        "EC number": source["ec_number"],
        "Class": genes.map(AARS_CLASS_LABELS).fillna("N/A"),
        "Number of structures": genes.map(n_structures).astype(int),
        "Number of pocket structures": genes.map(n_pocket_structures).astype(int),
        "Number of canonical pockets": genes.map(gene_to_unique_pocket_count).astype(int),
    })
    for domain in DOMAIN_COLUMNS:
        table[domain] = genes.map(domain_presence[domain])
    for interval in PLDDT_COLUMNS:
        table[f"pLDDT {interval} (%)"] = (genes.map(plddt_fractions[interval]) * 100).round(1)
    table["Sequence"] = source["sequence"]

    table = table.sort_values("Gene name").reset_index(drop=True)

    output_path = os.path.join(table_dir, "ExtendedDataTable1.xlsx")
    save_table_with_legend(table, build_legend(), output_path)
    print(f"Saved {len(table)} row(s) to {output_path}")
    print(table.drop(columns=["Sequence"]).to_string(index=False))


if __name__ == "__main__":
    main()
