from pymol import cmd
import pandas as pd
import os
import time

def compare_structures(ac1, file1, ac2, file2, struct_dir):

    f1 = os.path.join(struct_dir, ac1, file1)
    f2 = os.path.join(struct_dir, ac2, file2)
    name1 = file1.replace(".pdb", "")
    name2 = file2.replace(".pdb", "")

    # Global: PyMOL's sequence-independent, whole-chain cmd.super (script 15's method)
    cmd.reinitialize()
    cmd.load(f1, name1)
    cmd.load(f2, name2)

    res_1 = cmd.count_atoms(f"{name1} and name CA")
    res_2 = cmd.count_atoms(f"{name2} and name CA")

    global_result = cmd.super(name2, name1)
    global_rmsd = round(global_result[0], 2)
    global_aligned_res = global_result[6]
    global_coverage = round(global_aligned_res / min(res_1, res_2), 4)

    # Local: cealign (Combinatorial Extension, script 15b's method) - reload fresh so
    # cealign's own optimization isn't seeded by cmd.super's prior rotation/translation
    cmd.reinitialize()
    cmd.load(f1, name1)
    cmd.load(f2, name2)

    local_result = cmd.cealign(target=name1, mobile=name2)
    local_rmsd = round(local_result["RMSD"], 2)
    local_aligned_res = local_result["alignment_length"]
    local_coverage = round(local_aligned_res / min(res_1, res_2), 4)

    return res_1, res_2, global_rmsd, global_aligned_res, global_coverage, local_rmsd, local_aligned_res, local_coverage


# Dict mapping uniprot to filename (same source as scripts/15_calculate_StSim.py)
df = pd.read_csv("../output/alignment_relaxed_rmsd_data.csv")
uniprots = sorted(set(df["uniprot_ac"]))
uniprot_to_filename = {i: sorted(df[df["uniprot_ac"] == i]["file_name"]) for i in uniprots}

STRUCT_DIR = "../output/aligned_relaxed_structures"
OUTPUT_DIR = "../output/structural_comparisons_full"
ERROR_LOG = os.path.join(OUTPUT_DIR, "errors.log")

os.makedirs(OUTPUT_DIR, exist_ok=True)

COLUMNS = [
    "uniprot_ac_1", "uniprot_ac_2", "file_name_1", "file_name_2", "res_1", "res_2",
    "global_rmsd", "global_aligned_res", "global_coverage",
    "local_rmsd", "local_aligned_res", "local_coverage",
]

# Unordered pairs INCLUDING self-comparisons (i <= j) - self-pairs compare a protein's
# own structures (AF2/AF3/Chai-1/SwissModel/...) against each other, filling the
# within-protein consistency gap noted (but never computed) in figure_1_calculations.py.
pairs = [(uniprots[i], uniprots[j]) for i in range(len(uniprots)) for j in range(i, len(uniprots))]

start_time = time.time()
for ac1, ac2 in pairs:

    out_path = os.path.join(OUTPUT_DIR, f"{ac1}_{ac2}_global_local.csv")
    if os.path.exists(out_path):
        print(f"Skipping {ac1}/{ac2} (already done)")
        continue

    files1 = uniprot_to_filename[ac1]
    files2 = uniprot_to_filename[ac2]

    if ac1 == ac2:
        combos = [(files1[k], files1[l]) for k in range(len(files1)) for l in range(k + 1, len(files1))]
    else:
        combos = [(f1, f2) for f1 in files1 for f2 in files2]

    rows = []
    for file1, file2 in combos:
        try:
            res_1, res_2, g_rmsd, g_ares, g_cov, l_rmsd, l_ares, l_cov = compare_structures(
                ac1, file1, ac2, file2, STRUCT_DIR
            )
            rows.append((ac1, ac2, file1, file2, res_1, res_2, g_rmsd, g_ares, g_cov, l_rmsd, l_ares, l_cov))
        except Exception as e:
            with open(ERROR_LOG, "a") as f:
                f.write(f"{ac1},{ac2},{file1},{file2},{e}\n")
            continue

    pd.DataFrame(rows, columns=COLUMNS).to_csv(out_path, index=False)
    elapsed_min = (time.time() - start_time) / 60
    print(f"{ac1}/{ac2}: {len(combos)} combinations -> {out_path} (elapsed {elapsed_min:.1f} min)")

print("Done.")
