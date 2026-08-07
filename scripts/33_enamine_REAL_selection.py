from collections import Counter
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
import numpy as np
import pickle
import os

# Load labels and pockets
root = '.'
PATH_TO_PROBS = os.path.join(root, "..", "output", "unidock_docking", "inference_probs")
labels = np.array(pickle.load(open(os.path.join(PATH_TO_PROBS, "success_mols.pkl"), "rb")))
pockets = [i.replace("_bin_01.npz", "") for i in sorted(os.listdir(PATH_TO_PROBS)) if i != "success_mols.pkl"]
success_mols = np.array(pickle.load(open(os.path.join(PATH_TO_PROBS, "success_mols.pkl"), "rb")))
assert (success_mols == labels).all()

print(f"Number of pockets: {len(pockets)}")
print(f"Number of unique molecules: {len(set(labels))}")

# Load probabilities
print("Loading probabilities...")
pocket_to_all_probs = {}
for c, pocket in tqdm(enumerate(pockets)):
    probs = np.load(os.path.join(PATH_TO_PROBS, f"{pocket}_bin_01.npz"))["arr_0"]
    pocket_to_all_probs[pocket] = probs

# Prepare matrix
M = np.column_stack([pocket_to_all_probs[p] for p in pockets])
del pocket_to_all_probs
print(f"Matrix shape: {M.shape}")
print("Normalizing by columns...")
means_1 = np.mean(M, axis=0, keepdims=True)
stds_1 = np.std(M, axis=0, keepdims=True)
M -= means_1
M /= stds_1

print("Normalizing by rows...")
means_2 = np.mean(M, axis=1, keepdims=True)
stds_2 = np.std(M, axis=1, keepdims=True)
M -= means_2
M /= stds_2

print(f"Matrix normalized, shape: {M.shape}")

# Load smiles
df = pd.read_csv(os.path.join(root, "..", "output", "enamine_REAL_characterization", "enamine_REAL.tsv"), sep='\t')

# From compound ID to SMILES
synthons = {i: j for i, j in zip(df['id'], df['smiles'])}

# Dictionary with counts of each synthon
syns_counter = dict(Counter([j for i in synthons for j in i.split("____")[1:]]))

print(f"Number of unique 'synthons': {len(syns_counter)}")
print(f"1 occ: {len([i for i in syns_counter if syns_counter[i] == 1])}")
print(f">100 occ: {len([i for i in syns_counter if syns_counter[i] > 100])}")
print(f">1000 occ: {len([i for i in syns_counter if syns_counter[i] > 1000])}")
print(f"MAX occ: {max(syns_counter.values())}")

del df, synthons

print("Getting indices that sort every column...")

# Get indices that sort every column
inds = {}
for col in tqdm(range(M.shape[1])):
    inds[col] = np.argsort(M[:, col])[::-1].astype(np.int32)

# Remove huge matrix M from memory
del M

MAX_MOL = 100000
MAX_SYN = 5

ALL_CONSIDERED_MOLECULES = {}
ALL_CONSIDERED_SYNTHONS = {}
EVALUATED_MOLECULES = []

for c1, pocket in enumerate(pockets):

    # Molecules sorted
    sorted_molecules = labels[inds[c1]]

    # Counting synthons in the sorted molecules
    syns_counter_tmp = {i: 0 for i in syns_counter}
    considered_molecules = []

    # For each molecule
    for c2, mol in enumerate(sorted_molecules):

        syns = mol.split("____")[1:]
        include = [syns_counter_tmp[syn] < MAX_SYN for syn in syns]

        # If we can include the molecule
        if all(include) == True:
            considered_molecules.append(mol)
            for syn in syns:
                syns_counter_tmp[syn] += 1

        # If we can not include the molecule
        else:
            continue

        # If we have reached MAX_MOL
        if len(considered_molecules) >= MAX_MOL:
            break
    
    # Sanity check
    syns = [syn for mol in considered_molecules for syn in mol.split("____")[1:]]
    syns = dict(Counter(syns))
    assert all([syns[i] <= MAX_SYN for i in syns])

    print("\n")
    print(f"Pocket model: {pocket}")
    print(f"Number of considered molecules: {len(considered_molecules)}")
    print(f"Number of evaluated molecules: {c2+1}")
    print(f"Number of considered synthons: {len(syns)}")
    EVALUATED_MOLECULES.append(c2+1)

    ALL_CONSIDERED_MOLECULES[pocket] = set(considered_molecules)
    ALL_CONSIDERED_SYNTHONS[pocket] = set(syns)

active_molecules = set([j for pocket in ALL_CONSIDERED_MOLECULES for j in ALL_CONSIDERED_MOLECULES[pocket]])
active_synthons = set([j for pocket in ALL_CONSIDERED_SYNTHONS for j in ALL_CONSIDERED_SYNTHONS[pocket]])

print(f"Molecules being active at least once: {len(active_molecules)}")
print(f"Synthons being active at least once: {len(active_synthons)}")

print("Creating set of inactive compounds")
INACTIVES = set()
for molecule in tqdm(labels):
    # If molecule is not active
    if molecule not in active_molecules:
        # If None of the synthons is in active_synthons
        if any([syn in active_synthons for syn in molecule.split("____")[1:]]) == False:
            INACTIVES.add(molecule)

print(f"{len(INACTIVES)} compounds for which none of their synthons were found to be active in any of the 276 pocket models")

del active_synthons, active_molecules, ALL_CONSIDERED_SYNTHONS

PATH_TO_INPUT_LIGANDS_FILES = os.path.join(root, "..", "output", "unidock_REAL_docking", "input_ligands")
os.makedirs(PATH_TO_INPUT_LIGANDS_FILES, exist_ok=True)
print("Loading mol ID to split dict")

# Loading mol ID to split dict
df = pd.read_csv(os.path.join(root, "..", "..", "..", "Documents_GPU", "mtb-targeted-protein-degradation",
                              "output", "enamine_REAL_characterization", "conformations", "id_to_split.tsv.gz"), sep='\t', header=None)

id_to_split = {i: j for i,j in zip(df[0], df[1])}
del df

print("Creating input_ligands.txt files for each pocket")

for pocket in tqdm(pockets):

    # A couple of checks
    assert len(ALL_CONSIDERED_MOLECULES[pocket]) == 100000
    assert len(ALL_CONSIDERED_MOLECULES[pocket].intersection(INACTIVES)) == 0

    # Create file
    with open(os.path.join(PATH_TO_INPUT_LIGANDS_FILES, f"input_ligands_{pocket}.txt"), "w") as outfile:
        # Active molecules
        for mol_id in sorted(ALL_CONSIDERED_MOLECULES[pocket]):
            outfile.write(os.path.join(".", "..", "output", "unidock_REAL_docking", "conformations_prepared", f"{id_to_split[mol_id]:02d}", mol_id + ".sdf"))
            outfile.write("\n")

        # Inactive molecules
        for mol_id in sorted(INACTIVES):
            outfile.write(os.path.join(".", "..", "output", "unidock_REAL_docking", "conformations_prepared", f"{id_to_split[mol_id]:02d}", mol_id + ".sdf"))
            outfile.write("\n")

print("Finished!")