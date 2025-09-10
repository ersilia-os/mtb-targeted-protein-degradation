from rdkit.Chem import rdFingerprintGenerator
from joblib import Parallel, delayed
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
morganGen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

# Given a SMILES, return its fingerprint as packed bytes
def fp_words_and_popcount(smiles):
    m = Chem.MolFromSmiles(smiles)
    if not m: return None
    fp = morganGen.GetFingerprint(m)
    arr = np.zeros((2048,), dtype=np.uint8)
    Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    popc = int(arr.sum())
    fp_u64 = arr.reshape(32, 64)
    fp_u64 = np.packbits(fp_u64, axis=1, bitorder="little")
    fp_u64 = fp_u64.view(np.uint64).reshape(-1) 
    return fp_u64, popc

# ---------- Compute & save to HDF5 (compressed) ----------
H5_OUT = os.path.join(root, "..", "processed", "enamine_REAL_characterization", "enamine_REAL_ECFP4.h5")
CHUNK_ROWS = 1000000

with h5py.File(H5_OUT, "w") as h5:

    # Create datasets
    ids_ds = h5.create_dataset("ids",
        shape=(0,), maxshape=(None,),
        dtype='i8',
        chunks=True, compression="gzip", compression_opts=4
    )
    fps_ds = h5.create_dataset("fps",
        shape=(0, 32), maxshape=(None, 32),
        dtype="u8",
        chunks=(CHUNK_ROWS, 32),
        compression="gzip", compression_opts=4, shuffle=True
    )
    popc_ds = h5.create_dataset("popc",
    shape=(0,), maxshape=(None,),
    dtype="u2",
    chunks=True, compression="gzip", compression_opts=4
    )

    offset, n_total = 0, len(df)
    for start in range(0, n_total, CHUNK_ROWS):

        # Define start and stop
        stop = min(start + CHUNK_ROWS, n_total)
        sub = df.iloc[start:stop]

        # Parallel FP computation
        fps_results = [fp_words_and_popcount(smi) for smi in 
               tqdm(sub["smiles"].astype(str).tolist(), desc=f"Fingerprinting {start/10**6}M:{stop/10**6}M")]

        # Keep only valid results
        keep = [(fp, pc) for fp, pc in fps_results if fp is not None]
        if not keep:
            continue

        fps_batch, popc_batch = zip(*keep)
        fps_batch = np.stack(fps_batch, axis=0).astype(np.uint64)
        popc_batch = np.array(popc_batch, dtype=np.uint16)
        ids_batch = sub.loc[[fp is not None for fp, _ in fps_results], "index"].astype(np.int64).to_numpy()

        # Append to datasets
        new_n = offset + len(ids_batch)
        ids_ds.resize((new_n,))
        fps_ds.resize((new_n, 32))
        popc_ds.resize((new_n,))
        ids_ds[offset:new_n] = ids_batch
        fps_ds[offset:new_n, :] = fps_batch
        popc_ds[offset:new_n] = popc_batch
        offset = new_n