import pandas as pd
import os

root = os.path.dirname(os.path.abspath(__file__))


def create_fasta(uniprot_id, sequence, output_file, description="protein|empty"):
    """
    Create a FASTA file with the given UniProt ID and sequence.

    :param uniprot_id: The UniProt identifier (e.g., P69905)
    :param sequence: The protein sequence
    :param output_file: Output file name
    :param description: Optional description for the FASTA header
    """
    # Format the header
    header = f">{uniprot_id}|{description}"
    
    # Wrap sequence to 80 characters per line
    wrapped_sequence = "\n".join([sequence[i:i+80] for i in range(0, len(sequence), 80)])
    
    # Write to file
    with open(output_file, "w") as fasta_file:
        fasta_file.write(f"{header}\n{wrapped_sequence}\n")
    
    print(f"FASTA file created: {output_file}")

df = pd.read_csv(os.path.join(root, "..", "data", "mtb_trna_synthetases_bosch_2021_fig5_annotated.csv"))
for v in df[["gene_name_in_bosch_2021", "uniprot_ac", "sequence"]].values:
    gene_name, uniprot_id, sequence = v
    output_file = os.path.join(root, "..", "data", "sequences", f"{uniprot_id}.fasta")
    create_fasta(uniprot_id, sequence, output_file)