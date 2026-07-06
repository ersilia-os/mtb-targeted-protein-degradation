from collections import Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import os

root_gpu = "/home/acomajuncosa/Documents_GPU/mtb-targeted-protein-degradation"
PATH_TO_OUTPUT = os.path.join(root_gpu, "processed", "unidock_REAL_docking", "inference_10B")

def get_splits(path):
    splits = os.listdir(path)
    return [i.replace(".csv", "") for i in sorted(splits)]

def get_pockets():
    df = pd.read_csv("/home/acomajuncosa/Documents_GPU/mtb-targeted-protein-degradation/processed/pocket_detection_data.csv")
    pockets = [f"{i.replace('.pdb', '')}_pocket_{j}" for i, j in zip(df['File name'], df['Pocket number'])]
    return sorted(pockets)

# Get splits (994)
SPLITS = get_splits(os.path.join(PATH_TO_OUTPUT, 'A_proteins'))

# Get pockets (276)
POCKETS = get_pockets()

# Get proteins (21)
PROTEINS = sorted(set([i.split("_")[1] for i in POCKETS]))

print(f"Getting TOP 250k most promiscuous compounds (by pocket structure)")
df = []
for split in tqdm(SPLITS):
    df_ = pd.read_csv(os.path.join(os.path.join(PATH_TO_OUTPUT, "A_pockets", f"{split}.csv")))
    df.append(df_)
df = pd.concat(df, ignore_index=True)
rng = np.random.default_rng(42)
df["_r"] = rng.random(len(df))
df = (df.sort_values(["n_targets", "_r"], ascending=[False, True]).drop(columns="_r"))
print(f"Number of considered chunks in the TOP250k: {len(set(df[:250_000]['split']))}")
print(f"Number of considered chunks overall (sanity check): {len(set(df['split']))}")
print(f"Expected number of compounds per split: {round(250_000 / len(SPLITS), 1)}")
print("Minimum number of compounds per split")
print(sorted(Counter(df[:250_000]['split']).values())[:5])
print("Maximum number of compounds per split")
print(sorted(Counter(df[:250_000]['split']).values())[::-1][:5])
A_pockets = df[:250_000]
print(f"Number of unique compounds: {len(set([tuple(i) for i in A_pockets[['split', 'index']].values]))}")
print(f"Number of targets: {round(np.mean(A_pockets['n_targets'].tolist()), 1)} ± {round(np.std(A_pockets['n_targets'].tolist()), 1)}")
A_pockets.to_csv(os.path.join(PATH_TO_OUTPUT, "A_pockets.csv"), index=False)
del df, A_pockets




print(f"Getting TOP 276k most selective* compounds (by pocket structure)")
df = []
for split in tqdm(SPLITS):
    df_ = pd.read_csv(os.path.join(os.path.join(PATH_TO_OUTPUT, "B_pockets", f"{split}.csv")))
    df.append(df_)
df = pd.concat(df, ignore_index=True)
print(f"Getting only 1,000 compounds per pocket")
numb_splits = []
B_pockets = []
for c, pocket in tqdm(enumerate(POCKETS)):
    rng = np.random.default_rng(c)
    df_pocket = df[df['pocket'] == pocket].copy()
    df_pocket["_r"] = rng.random(len(df_pocket))
    df_pocket = df_pocket.sort_values(["n_targets", "_r"], ascending=True).drop(columns="_r")[:1_000]
    numb_splits.append(len(set(df_pocket['split'])))
    B_pockets.append(df_pocket)
B_pockets = pd.concat(B_pockets, ignore_index=True)
print(f"Number of unique compounds: {len(set([tuple(i) for i in B_pockets[['split', 'index']].values]))}")
print(f"Number of targets (overall): {round(np.mean(B_pockets['n_targets'].tolist()), 1)} ± {round(np.std(B_pockets['n_targets'].tolist()), 1)}")
print(f"Number of splits per compound: {round(np.mean(numb_splits), 1)} ± {round(np.std(numb_splits), 1)}")
B_pockets.to_csv(os.path.join(PATH_TO_OUTPUT, "B_pockets.csv"), index=False)
del df, B_pockets


print(f"Getting TOP 250k most promiscuous compounds (by protein)")
df = []
for split in tqdm(SPLITS):
    df_ = pd.read_csv(os.path.join(os.path.join(PATH_TO_OUTPUT, "A_proteins", f"{split}.csv")))
    df.append(df_)
df = pd.concat(df, ignore_index=True)
rng = np.random.default_rng(42)
df["_r"] = rng.random(len(df))
df = (df.sort_values(["n_targets", "_r"], ascending=[False, True]).drop(columns="_r"))
print(f"Number of considered chunks in the TOP250k: {len(set(df[:250_000]['split']))}")
print(f"Number of considered chunks overall (sanity check): {len(set(df['split']))}")
print(f"Expected number of compounds per split: {round(250_000 / len(SPLITS), 1)}")
print("Minimum number of compounds per split")
print(sorted(Counter(df[:250_000]['split']).values())[:5])
print("Maximum number of compounds per split")
print(sorted(Counter(df[:250_000]['split']).values())[::-1][:5])
A_proteins = df[:250_000]
print(f"Number of unique compounds: {len(set([tuple(i) for i in A_proteins[['split', 'index']].values]))}")
print(f"Number of targets: {round(np.mean(A_proteins['n_targets'].tolist()), 1)} ± {round(np.std(A_proteins['n_targets'].tolist()), 1)}")
A_proteins.to_csv(os.path.join(PATH_TO_OUTPUT, "A_proteins.csv"), index=False)
del df, A_proteins


print(f"Getting TOP 273k most selective* compounds (by protein)")
df = []
for split in tqdm(SPLITS):
    df_ = pd.read_csv(os.path.join(os.path.join(PATH_TO_OUTPUT, "B_proteins", f"{split}.csv")))
    df.append(df_)
df = pd.concat(df, ignore_index=True)

print(f"Getting only 13,000 compounds per protein")
numb_splits = []
B_proteins = []
for c, protein in tqdm(enumerate(PROTEINS)):
    rng = np.random.default_rng(c)
    df_protein = df[df['pocket'] == protein].copy()
    df_protein["_r"] = rng.random(len(df_protein))
    df_protein = df_protein.sort_values(["n_targets", "_r"], ascending=True).drop(columns="_r")[:13_000]
    numb_splits.append(len(set(df_protein['split'])))
    B_proteins.append(df_protein)
B_proteins = pd.concat(B_proteins, ignore_index=True)
print(f"Number of unique compounds: {len(set([tuple(i) for i in B_proteins[['split', 'index']].values]))}")
print(f"Number of targets (overall): {round(np.mean(B_proteins['n_targets'].tolist()), 1)} ± {round(np.std(B_proteins['n_targets'].tolist()), 1)}")
print(f"Number of splits per compound: {round(np.mean(numb_splits), 1)} ± {round(np.std(numb_splits), 1)}")
B_proteins.to_csv(os.path.join(PATH_TO_OUTPUT, "B_proteins.csv"), index=False)
del df, B_proteins