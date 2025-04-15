from pymol import cmd
import pandas as pd
import os

def align_structure_with_pymol(reference_file, mobile_file):

    # Extract base names for the reference and mobile structure files
    reference_name = os.path.basename(reference_file).split(".pdb")[0]
    mobile_name = os.path.basename(mobile_file).split(".pdb")[0]

    # Load both structures into PyMOL
    cmd.reinitialize()
    cmd.load(reference_file, reference_name)
    cmd.load(mobile_file, mobile_name)

    # Perform the alignment using PyMOL's super command
    alignment_result = cmd.super(mobile_name, reference_name)

    # Extract RMSD value (first element in the result tuple)
    rmsd = round(alignment_result[0], 2)

    return rmsd


# Dict mapping uniprot to filename
df = pd.read_csv("../processed/alignment_relaxed_rmsd_data.csv")
uniprots = sorted(set(df["uniprot_ac"]))
uniprot_to_filename = {i: sorted(df[df['uniprot_ac'] == i]['file_name']) for i in uniprots}

# Get all possible uniprot pairs
uni_pairs = [(i,j) for i in uniprots for j in uniprots]

PATH_TO_STRUCTURES = "../processed/aligned_relaxed_structures"
PATH_TO_OUTPUT = "../processed/structural_comparisons"

os.makedirs(PATH_TO_OUTPUT, exist_ok=True)

for uni1, uni2 in uni_pairs:

    if uni1 != uni2:

        print(uni1, uni2)
        results = dict()

        # Get associated structures
        sts1 = uniprot_to_filename[uni1]
        sts2 = uniprot_to_filename[uni2]

        # Get paths to sts
        sts1 = [os.path.join(PATH_TO_STRUCTURES, uni1, i) for i in sts1]
        sts2 = [os.path.join(PATH_TO_STRUCTURES, uni2, i) for i in sts2]

        # For each structure pair
        for st1 in sts1:
            for st2 in sts2:

                # Calculate RMSD
                rmsd = align_structure_with_pymol(st1, st2)
                results[(os.path.basename(st1), os.path.basename(st2))] = rmsd

        # Save results to CSV
        results_df = pd.DataFrame([[i[0], i[1], results[i]] for i in sorted(results)], columns=["file_name_1", "file_name_2", "rmsd"])
        results_df.to_csv(os.path.join(PATH_TO_OUTPUT, f"{uni1}_{uni2}_rmsd.csv"), index=False)
        print(f"Results saved to {os.path.join(PATH_TO_OUTPUT, f'{uni1}_{uni2}_rmsd.csv')}")