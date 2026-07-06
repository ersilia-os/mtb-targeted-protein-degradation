#!/bin/bash
# Run ersilia model eos42ez on combined SMILES from scripts 48 and 49.
# Usage: ./50_ersilia_eos42ez.sh <file_48.csv> <file_49.csv>
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <file_48.csv> <file_49.csv>"
    exit 1
fi

FILE_48="$1"
FILE_49="$2"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTPUT_DIR="$ROOT/output/50_ersilia_eos42ez"
mkdir -p "$OUTPUT_DIR"

# --- 1. Conda env ---
if conda env list | grep -qE "^ersilia_env[[:space:]]"; then
    echo "Conda env 'ersilia_env' already exists, skipping creation."
else
    conda create -n ersilia_env python=3.12 -y
fi

# --- 2. Install ersilia ---
conda run -n ersilia_env pip install ersilia

# --- 3. Merge and deduplicate SMILES ---
conda run -n ersilia_env python -c "
import pandas as pd

df1 = pd.read_csv('$FILE_48')[['smiles']]

raw2 = pd.read_csv('$FILE_49')
df2 = raw2[raw2['selected'].notna()][['smiles']].reset_index(drop=True)

n48 = len(df1)
n49 = len(df2)
print(f'Compounds from script 48: {n48}')
print(f'Compounds from script 49 (selected only): {n49}')

combined = pd.concat([df1, df2], ignore_index=True)
duplicates = combined[combined.duplicated(subset='smiles', keep=False)]['smiles'].unique()
n_dups = len(duplicates)
if n_dups:
    print(f'Duplicate SMILES ({n_dups}):')
    for smi in duplicates:
        print(f'  {smi}')
else:
    print('No duplicate SMILES found.')

combined = combined.drop_duplicates(subset='smiles').reset_index(drop=True)
print(f'Total unique compounds: {len(combined)}')

combined.to_csv('$OUTPUT_DIR/smiles.csv', index=False)
"

# --- 4. Source conda and activate env for ersilia commands ---
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ersilia_env

# --- 5. Ersilia workflow ---
ersilia fetch eos42ez
ersilia serve eos42ez
ersilia run -i "$OUTPUT_DIR/smiles.csv" -o "$OUTPUT_DIR/eos42ez.csv" --batch_size 100
ersilia close
ersilia delete eos42ez

echo "Done. Predictions written to $OUTPUT_DIR/eos42ez.csv"
