# This file is conceptually analogous to 18_enamine_mfps.py but using CheMeleon embeddings instead
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
import pickle
import torch
import os

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
PATH_TO_MFPS = "../processed/enamine_characterization"
OUTDIR = PATH_TO_MFPS
IK_TO_SMILES = pickle.load(open(os.path.join(PATH_TO_MFPS, "IK_TO_SMI.pkl"), "rb"))
ID_TO_SMILES = pickle.load(open(os.path.join(PATH_TO_MFPS, "ID_TO_SMI.pkl"), "rb"))
ID_TO_IK = pickle.load(open(os.path.join(PATH_TO_MFPS, "ID_TO_IK.pkl"), "rb"))

# Instantiate CheMeleong embeddings
chemeleon_fingerprint = CheMeleonFingerprint(device="cuda")  # remove cuda if no GPU is available

print("Calculating CheMeleon embeddings...")
BATCH_SIZE = 1000

# Get all features
ALL_IDs = sorted(ID_TO_SMILES)
X, IDs = [], []
for batch in tqdm(range(0, len(ID_TO_SMILES), BATCH_SIZE)):
    ids = ALL_IDs[batch: batch + BATCH_SIZE]
    smiles = [ID_TO_SMILES[i] for i in ids]
    chemeleons = chemeleon_fingerprint(smiles)
    X.extend(chemeleons)
    IDs.extend(ids)

# Convert to numpy array
X = np.array(X, dtype="float16")

# Assert that the number of IDs and MFPs are the same
assert len(IDs) == X.shape[0], f"ERROR: Number of IDs ({len(IDs)}) does not match number of CheMeleon embeddings ({X.shape[0]})"

# Save results
np.savez_compressed(os.path.join(OUTDIR, "X_CheMeleon.npz"), X=X)
with open(os.path.join(OUTDIR, "IDs_CheMeleon.txt"), "w") as f:
    f.write("\n".join(IDs))