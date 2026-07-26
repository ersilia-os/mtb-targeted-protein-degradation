import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

from default import RANDOM_SEED
from docking_utils import compute_properties, sample_prescreened_smiles

plots_dir = os.path.join(root, "..", "..", "plots", "figure_2")
os.makedirs(plots_dir, exist_ok=True)

HL_SAMPLE_N = 100_000
REAL10M_SAMPLE_N = 100_000


def banner(title):
    line = "=" * (len(title) + 10)
    print(line)
    print(f"==== {title} ====")
    print(line)


banner("HL (Enamine Hit Locator 100K) - random 100k sample")
hl_smiles = sample_prescreened_smiles("DL", exclude_ids=set(), n=HL_SAMPLE_N, seed=RANDOM_SEED)
print(f"Sampled {len(hl_smiles):,} compounds")
hl_props = compute_properties(hl_smiles)
output_path = os.path.join(plots_dir, "figure_2a_HL.csv")
hl_props.to_csv(output_path)
print(f"Saved to {output_path}")

banner("REAL 10M (Enamine REAL 9.56M) - random 100k sample")
real10m_smiles = sample_prescreened_smiles("REAL_ROUND1", exclude_ids=set(), n=REAL10M_SAMPLE_N, seed=RANDOM_SEED)
print(f"Sampled {len(real10m_smiles):,} compounds")
real10m_props = compute_properties(real10m_smiles)
output_path = os.path.join(plots_dir, "figure_2a_REAL10M.csv")
real10m_props.to_csv(output_path)
print(f"Saved to {output_path}")
