from scipy.sparse import coo_matrix, csr_matrix, hstack
from tqdm import tqdm
import pandas as pd
import numpy as np
import os

root_gpu = "/home/acomajuncosa/Documents_GPU/mtb-targeted-protein-degradation"
screening_aloy = "/aloy/home/acomajuncosa/Ersilia/gcadda4tb-enamine-real-screening"

PATH_TO_OUTPUT = os.path.join(root_gpu, "processed", "unidock_REAL_docking", "inference_10B")
os.makedirs(os.path.join(PATH_TO_OUTPUT, "A_pockets"), exist_ok=True)
os.makedirs(os.path.join(PATH_TO_OUTPUT, "B_pockets"), exist_ok=True)
os.makedirs(os.path.join(PATH_TO_OUTPUT, "A_proteins"), exist_ok=True)
os.makedirs(os.path.join(PATH_TO_OUTPUT, "B_proteins"), exist_ok=True)

def get_splits(screening_aloy):
    splits = os.listdir(os.path.join(screening_aloy, "results"))
    return sorted(splits)

def get_pockets():
    df = pd.read_csv("/home/acomajuncosa/Documents_GPU/mtb-targeted-protein-degradation/processed/pocket_detection_data.csv")
    pockets = [f"{i.replace('.pdb', '')}_pocket_{j}" for i, j in zip(df['File name'], df['Pocket number'])]
    return sorted(pockets)

def get_matrix(split, pockets, perc=1):
    rows, cols = [], []
    for c, pocket in enumerate(pockets):
        ind = np.load(os.path.join(screening_aloy, "results", split, f"{pocket}_ind_{perc}.npz"))[f'ind_{perc}']
        rows.append(ind.astype(np.int32, copy=False))
        cols.append(np.full(ind.size, c, dtype=np.int32))
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    M = coo_matrix((np.ones(rows.size, dtype=np.uint8), (rows, cols)),
                shape=(10_000_000, len(pockets))).tocsr()
    return M

def get_A_compounds(M, split, top=10_000, seed=42):
    sum_cpds = np.asarray(M.sum(axis=1)).ravel().astype(np.float32)
    rng = np.random.default_rng(seed)
    inds = np.lexsort((rng.random(len(sum_cpds)), sum_cpds))[::-1]
    A_compounds = [[split, i, j] for i, j in zip(inds[:top], sum_cpds[inds][:top])]
    A_compounds = pd.DataFrame(A_compounds, columns=['split', 'index', 'n_targets'])
    return A_compounds

def get_B_compounds(M, split, POCKETS, top_per_pocket=100, min_targets=50, seed=42):
    sum_cpds = np.asarray(M.sum(axis=1)).ravel().astype(np.float32)
    eligible = sum_cpds >= float(min_targets)
    rng = np.random.default_rng(seed)
    B_compounds = []
    for column in range(M.shape[1]):
        r = M[:, column].nonzero()[0]  # Get non-zero indices
        r = r[eligible[r]]  # Only eligible cpds 
        n_targets = sum_cpds[r]
        order = np.lexsort((rng.random(len(n_targets)), n_targets))
        top_inds = r[order][:top_per_pocket]
        top_n_targets  = n_targets[order][:top_per_pocket]
        B_compounds.extend([split, POCKETS[column], i, j] for i,j in zip(top_inds, top_n_targets))
    B_compounds= pd.DataFrame(B_compounds, columns=['split', 'pocket', 'index', 'n_targets'])
    return B_compounds

def collapse_per_protein(M, proteins, PROTEINS, PROTEINS_TO_NUMBER):
    cols = []
    for protein in PROTEINS:
        ind_column = np.where(np.array(proteins) == PROTEINS_TO_NUMBER[protein])[0]
        col = M[:, ind_column].max(axis=1)
        cols.append(csr_matrix(col, dtype=np.uint8))
    M_protein = hstack(cols, format="csr").astype(np.uint8)
    return M_protein


# Get splits (994)
SPLITS = get_splits(screening_aloy)

# Get pockets (276)
POCKETS = get_pockets()

# Get proteins (21)
PROTEINS = sorted(set([i.split("_")[1] for i in POCKETS]))
PROTEINS_TO_NUMBER = {i: c for c,i in enumerate(PROTEINS)}
proteins = [PROTEINS_TO_NUMBER[i.split("_")[1]] for i in POCKETS]

for split in tqdm(SPLITS):

    # Get matrix
    M = get_matrix(split, POCKETS)
  
    # Condition A
    A_pockets = get_A_compounds(M, split, top=10_000)

    # Condition B
    B_pockets = get_B_compounds(M, split, POCKETS, top_per_pocket=100, min_targets=50)

    # Collase matrix
    M_protein = collapse_per_protein(M, proteins, PROTEINS, PROTEINS_TO_NUMBER)

    # Get A compounds protein
    A_proteins = get_A_compounds(M_protein, split, top=10_000)

    # Get B compounds protein
    B_proteins = get_B_compounds(M_protein, split, PROTEINS, top_per_pocket=100, min_targets=2)

    # Save results
    A_pockets.to_csv(os.path.join(os.path.join(PATH_TO_OUTPUT, "A_pockets", f"{split}.csv")), index=False)
    B_pockets.to_csv(os.path.join(os.path.join(PATH_TO_OUTPUT, "B_pockets", f"{split}.csv")), index=False)
    A_proteins.to_csv(os.path.join(os.path.join(PATH_TO_OUTPUT, "A_proteins", f"{split}.csv")), index=False)
    B_proteins.to_csv(os.path.join(os.path.join(PATH_TO_OUTPUT, "B_proteins", f"{split}.csv")), index=False)