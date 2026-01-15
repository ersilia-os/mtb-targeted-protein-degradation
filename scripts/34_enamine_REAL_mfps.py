# This script needs to be run with rdkit version 2025.9.1  <====== (!!!!!)
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import AllChem
from rdkit import Chem
from tqdm import tqdm
import pandas as pd
import numpy as np
import h5py
import os

# Load smiles
root = "."
df = pd.read_csv(os.path.join(root, "..", "processed", "enamine_REAL_characterization", "enamine_REAL.tsv"), sep='\t')
df['index'] = df.index

# Build the generator once
mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)

def clip_sparse(vect, nBits=2048):
    """
    Convert a sparse RDKit fingerprint vector to a dense int8 list.

    Parameters
    ----------
    vect : rdkit.DataStructs.cDataStructs.ExplicitBitVect or similar
        Sparse fingerprint vector with nonzero elements.
    nBits : int, optional, default=2048
        Length of the dense output vector.

    Returns
    -------
    list of int
        Dense list representation of the fingerprint, where values are
        clipped to the maximum representable int8 value.
    """
    l = [0] * nBits
    MAX_I8 = 127
    for i, v in vect.GetNonzeroElements().items():
        l[i] = v if v < MAX_I8 else MAX_I8
    return np.array(l, dtype=np.int8)

# Given a SMILES, return its fingerprint
def calculate_ecfp6(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return None
    ecfp = mfpgen.GetCountFingerprint(mol)
    ecfp = clip_sparse(ecfp)
    return np.array(ecfp, dtype=np.int8)

# ---------- Compute & save to HDF5 (compressed) ----------
H5_OUT = os.path.join(root, "..", "processed", "enamine_REAL_characterization", "enamine_REAL_ECFP6.h5")
CHUNK_ROWS = 1000000

with h5py.File(H5_OUT, "w") as h5:

    # Create datasets
    ids_ds = h5.create_dataset("ids",
        shape=(0,), maxshape=(None,),
        dtype='i8',
        chunks=True, compression="gzip", compression_opts=4
    )
    fps_ds = h5.create_dataset("fps",
        shape=(0, 2048), maxshape=(None, 2048),
        dtype="i1",
        chunks=(CHUNK_ROWS, 2048),
        compression="gzip", compression_opts=4, shuffle=True
    )

    offset, n_total = 0, len(df)
    for start in range(0, n_total, CHUNK_ROWS):

        # Define start and stop
        stop = min(start + CHUNK_ROWS, n_total)
        sub = df.iloc[start:stop]

        # Calculate fingerprints
        fps_results = [calculate_ecfp6(smi) for smi in 
               tqdm(sub["smiles"].astype(str), desc=f"Fingerprinting {start/10**6}M:{stop/10**6}M")]

        # Keep only valid results
        keep = [fp for fp in fps_results if fp is not None]
        if not keep:
            continue

        fps_batch = np.stack(keep, axis=0)
        ids_batch = sub.loc[[fp is not None for fp in fps_results], "index"].astype(np.int64).to_numpy()


        # Append to datasets
        new_n = offset + len(ids_batch)
        ids_ds.resize((new_n,))
        fps_ds.resize((new_n, 2048))
        ids_ds[offset:new_n] = ids_batch
        fps_ds[offset:new_n, :] = fps_batch
        offset = new_n