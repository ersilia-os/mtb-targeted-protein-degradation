from rdkit.Chem import Descriptors, Crippen, QED
from rdkit.Chem.Draw import IPythonConsole
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import Draw
from collections import Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import random
import os

root = "/home/acomajuncosa/Documents_GPU/mtb-targeted-protein-degradation"
PATH_TO_SELECTED_COMPOUNDS = os.path.join(root, "processed", "unidock_REAL_docking", "inference_10B", "selected_compounds")

def get_splits(path):
    splits = os.listdir(path)
    return [i.replace(".csv", "") for i in sorted(splits)]

def get_pockets():
    df = pd.read_csv("/home/acomajuncosa/Documents_GPU/mtb-targeted-protein-degradation/processed/pocket_detection_data.csv")
    pockets = [f"{i.replace('.pdb', '')}_pocket_{j}" for i, j in zip(df['File name'], df['Pocket number'])]
    return sorted(pockets)

def smiles_grid(smiles, per_row=4, size=(200, 200)):
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    return Draw.MolsToGridImage(mols, molsPerRow=per_row, subImgSize=size, maxMols=10_000)

# Get splits (994)
SPLITS = get_splits(os.path.join(root, "processed", "unidock_REAL_docking", "inference_10B", 'A_proteins'))

# Get pockets (276)
POCKETS = get_pockets()

df = []
for split in SPLITS:
    df_ = pd.read_csv(os.path.join(os.path.join(PATH_TO_SELECTED_COMPOUNDS, f"{split}.csv")))
    df.append(df_)
df = pd.concat(df, ignore_index=True)
rng = np.random.default_rng(42)
df["rd"] = rng.random(len(df))
CUSTOM_ORDER = ['A_pockets', 'A_proteins', 'B_proteins', 'B_pockets']
df["label"] = pd.Categorical(df["label"], categories=CUSTOM_ORDER, ordered=True)
df = df.sort_values(["label", "rd"], ascending=[True, True], kind="stable").reset_index(drop=True)
# df[['smiles']].to_csv(os.path.join(root, "processed", "unidock_REAL_docking", "inference_10B", "selected_compounds_smiles.csv"), index=False)  # to run ersilia models

print(Counter(df['index'] == df['index.1']))

# random.seed(42)
# rd_cpds = random.sample(df['smiles'].tolist(), 10 * 25)
# rd_cpds = [i.split()[0] for i in rd_cpds]
# smiles_grid(rd_cpds, per_row=25)

mols = [Chem.MolFromSmiles(smi) for smi in tqdm(df['smiles'])]
df['MW'] = [Descriptors.MolWt(m) for m in tqdm(mols)]
df['logp'] = [Crippen.MolLogP(m) for m in tqdm(mols)]
df['QED'] = [QED.qed(m) for m in tqdm(mols)]
eos12x7 = pd.read_csv(os.path.join(root, "processed", "unidock_REAL_docking", "inference_10B", "ersilia", "eos12x7.csv"))[["input", "nsps_score"]]


assert (eos12x7['input'] == df['smiles']).all
results = pd.concat([df, eos12x7], axis=1)
assert (results['input'] == results['smiles']).all
del results['index.1'], results['input']


# Filter
results = results[(results['MW'] > 250) & (results['MW'] < 450) & 
                  (results['logp'] > -1) & (results['logp'] < 5) & 
                  (results['QED'] > 0.4) & 
                  (results['nsps_score'] > 10) & (results['nsps_score'] < 40)].reset_index(drop=True)


print(len(results))
results.to_csv(os.path.join(root, "processed", "unidock_REAL_docking", "inference_10B", "filtered_compounds.csv"), index=False)