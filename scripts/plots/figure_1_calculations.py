import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

import json
import pandas as pd
import stylia
from matplotlib.colors import to_hex

data_dir = os.path.join(root, "..", "..", "data")
plots_dir = os.path.join(root, "..", "..", "plots", "figure_1")
os.makedirs(plots_dir, exist_ok=True)


def banner(title):
    line = "=" * (len(title) + 10)
    print(line)
    print(f"==== {title} ====")
    print(line)


uniprot_to_gene_df = pd.read_csv(os.path.join(data_dir, "mtb_trna_synthetases_bosch_2021_fig5.csv"))
uniprot_to_gene = {i: j for i, j in zip(uniprot_to_gene_df['uniprot_ac'], uniprot_to_gene_df['gene_name_in_bosch_2021'])}

uniprot_ids = ["P9WFW5", "P9WFW7", "P9WFW3", "P9WQA1", "P9WN61", "P9WFT5", "P9WFV3", "P9WFS9", "P9WFV1", "P9WFV9", "P9WFT9", "P9WFV7", "P9WFT7",
               "P9WFW1", "P9WFU5", "P9WFU9", "P9WFV5", "P9WFT3", "P9WFU3", "P9WFU1", "P9WFT1"]

# stylia SpectralColormap ("npg" preset), one color per protein, assigned in
# alphabetical gene-name order so the gradient lines up with PROTEINS in figure_1_plot.py.
genes_sorted = sorted(uniprot_to_gene[uid] for uid in uniprot_ids)
palette = stylia.SpectralColormap("npg").sample(len(genes_sorted))
gene_to_color = {gene: to_hex(palette[i]) for i, gene in enumerate(genes_sorted)}

banner("SAVING MAPPINGS")
mappings = {
    "uniprot_to_gene": uniprot_to_gene,
    "gene_to_color": gene_to_color,
}
output_path = os.path.join(plots_dir, "color_mapping.json")
with open(output_path, "w") as f:
    json.dump(mappings, f, indent=2)
print(f"Saved to {output_path}")
