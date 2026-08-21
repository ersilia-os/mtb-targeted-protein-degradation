from pymol import cmd
import pandas as pd
import os

def align_structure_with_cealign(reference_file, mobile_file):

    # Extract base names for the reference and mobile structure files
    reference_name = os.path.basename(reference_file).split(".pdb")[0]
    mobile_name = os.path.basename(mobile_file).split(".pdb")[0]

    # Load both structures into PyMOL
    cmd.reinitialize()
    cmd.load(reference_file, reference_name)
    cmd.load(mobile_file, mobile_name)

    # Perform the alignment using PyMOL's cealign command (Combinatorial Extension:
    # a local, topology-independent structural alignment - unlike super/align it
    # does not assume a global 1:1 residue correspondence across the whole chain)
    result = cmd.cealign(target=reference_name, mobile=mobile_name)

    rmsd = round(result["RMSD"], 2)
    alignment_length = result["alignment_length"]

    return rmsd, alignment_length


# Dict mapping uniprot to filename
df = pd.read_csv("../output/alignment_relaxed_rmsd_data.csv")
uniprots = sorted(set(df["uniprot_ac"]))

PATH_TO_STRUCTURES = "../output/aligned_relaxed_structures"
PATH_TO_OUTPUT = "../output/structural_comparisons_local"

os.makedirs(PATH_TO_OUTPUT, exist_ok=True)

results = []

# Unlike scripts/15_calculate_StSim.py (which compares every predicted structure of
# protein A against every predicted structure of protein B), this uses only the AF2
# model per protein, and each unordered pair is computed once.
for i in range(len(uniprots)):
    for j in range(i + 1, len(uniprots)):

        uni1, uni2 = uniprots[i], uniprots[j]
        print(uni1, uni2)

        st1 = os.path.join(PATH_TO_STRUCTURES, uni1, f"alphafold2_{uni1}_model_0.pdb")
        st2 = os.path.join(PATH_TO_STRUCTURES, uni2, f"alphafold2_{uni2}_model_0.pdb")

        rmsd, alignment_length = align_structure_with_cealign(st1, st2)
        results.append((uni1, uni2, rmsd, alignment_length))

results_df = pd.DataFrame(results, columns=["uniprot_ac_1", "uniprot_ac_2", "rmsd", "alignment_length"])
out_path = os.path.join(PATH_TO_OUTPUT, "StSim_local_pairwise.csv")
results_df.to_csv(out_path, index=False)
print(f"Results saved to {out_path}")
