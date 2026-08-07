# FULL ENAMINE REAL CHARACTERIZATION
# CAUTION: THIS SCRIPT NEEDS TO BE RUN WITH A CONDA ENVIRONMENT HAVING CHEMPROP >= 2.2.0 ON A GPU MACHINE FOR MASSIVE PARALLELIZATION

from rdkit.Chem import MolFromSmiles, Mol
from urllib.request import urlretrieve
from chemprop import featurizers, nn
from chemprop.data import BatchMolGraph
from chemprop.nn import RegressionFFN
from chemprop.models import MPNN
from pathlib import Path
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
import torch
import h5py
import os

root = os.path.dirname(os.path.abspath(__file__))

class CheMeleonFingerprint:
    def __init__(self, device: str | torch.device | None = None):
        self.featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
        agg = nn.MeanAggregation()
        root = os.path.dirname(os.path.abspath(__file__))
        # ckpt_dir = Path().home() / ".chemprop"
        # ckpt_dir.mkdir(exist_ok=True)
        # mp_path = ckpt_dir / "chemeleon_mp.pt"
        mp_path = Path(root) / ".." / "data" / "CheMeleon" / "chemeleon_mp.pt"
        mp_path = mp_path.resolve()
        if not mp_path.exists():
            urlretrieve(
                r"https://zenodo.org/records/15460715/files/chemeleon_mp.pt",
                mp_path,
            )
        chemeleon_mp = torch.load(mp_path, weights_only=True)
        mp = nn.BondMessagePassing(**chemeleon_mp['hyper_parameters'])
        mp.load_state_dict(chemeleon_mp['state_dict'])
        self.model = MPNN(
            message_passing=mp,
            agg=agg,
            predictor=RegressionFFN(input_dim=mp.output_dim),  # not actually used
        )
        self.model.eval()
        if device is not None:
            self.model.to(device=device)

    def __call__(self, molecules: list[str | Mol]) -> np.ndarray:
        bmg = BatchMolGraph([self.featurizer(MolFromSmiles(m) if isinstance(m, str) else m) for m in molecules])
        bmg.to(device=self.model.device)
        return self.model.fingerprint(bmg).numpy(force=True)

# Get input data
PATH_TO_MFPS = os.path.join(root, "..", "output", "enamine_REAL_characterization")
PATH_TO_DATA = os.path.join(root, "..", "data", "Enamine", "2024.07_Enamine_REAL_DB_9.6M.cxsmiles")

# Load data
data = pd.read_csv(PATH_TO_DATA, sep="\t", low_memory=False)

# Create sorted list of records using list comprehension
records = sorted(
    [tuple([ik, id_, smi]) for ik, id_, smi in zip(data["InChiKey"], data["id"], data["smiles"])],
    key=lambda x: x[0]  
)

print(f"Number of compounds in the original dataset: {len(data)}")
print(f"Number of IDs in the final compound set: {len(set([i[1] for i in records]))}")


# Store data in a csv file
os.makedirs(PATH_TO_MFPS, exist_ok=True)
output_df = pd.DataFrame(records, columns=["InChiKey", "id", "smiles"])
output_df.to_csv(os.path.join(PATH_TO_MFPS, "enamine_REAL.tsv"), index=False, sep='\t')


# Instantiate CheMeleon embeddings
chemeleon_fingerprint = CheMeleonFingerprint(device="cuda")  # remove cuda if no GPU is available
ID_TO_SMILES = {rec[1]: rec[2] for rec in records}
ALL_IDs = sorted(ID_TO_SMILES)
BATCH_SIZE = 1000
CHUNK_SIZE = 100000

print("Generating CheMeleon embeddings and writing to chunked HDF5 files...")
outpath = os.path.join(PATH_TO_MFPS, "embeddings")
os.makedirs(outpath, exist_ok=True)

# Create chunks and write to HDF5
for file_idx, chunk_start in enumerate(range(0, len(ALL_IDs), CHUNK_SIZE)):

    # Chunk size
    chunk_end = min(chunk_start + CHUNK_SIZE, len(ALL_IDs))
    chunk_ids = ALL_IDs[chunk_start:chunk_end]

    # Create h5 file
    h5_path = os.path.join(outpath, f"enamine_REAL_chemeleon_chunk_{file_idx}.h5")
    with h5py.File(h5_path, "w") as h5f:

        # Create datasets
        emb_ds = h5f.create_dataset("X", shape=(len(chunk_ids), 2048), dtype="float16", compression="gzip")
        id_ds = h5f.create_dataset("ids", shape=(len(chunk_ids),), dtype=h5py.string_dtype(encoding="utf-8"), compression="gzip")

        # For each batch in the chunk
        for batch_start in tqdm(range(0, len(chunk_ids), BATCH_SIZE), desc=f"Chunk {file_idx}"):
            batch_end = min(batch_start + BATCH_SIZE, len(chunk_ids))
            batch_ids = chunk_ids[batch_start:batch_end]
            batch_smiles = [ID_TO_SMILES[i] for i in batch_ids]

            embeddings = chemeleon_fingerprint(batch_smiles)
            embeddings = np.array(embeddings, dtype="float16")

            # Write batch to HDF5
            emb_ds[batch_start:batch_end, :] = embeddings
            id_ds[batch_start:batch_end] = np.array(batch_ids, dtype=h5py.string_dtype(encoding="utf-8"))