#!/usr/bin/env python3
"""
Renders script 70's 1,095 filtered hits (output/70_filtering/filtered_hits.csv) as 6 GIFs, one
molecule per frame, for quick visual browsing. Compounds are split into 6 sequential,
as-equal-as-possible chunks in file order (no reordering/sorting). Uses chemgifs'
(github.com/ersilia-os/chemical-library-gifs) new --style rdkit renderer (RDKit CoordGen +
MolDraw2D, matching the depiction style of Ersilia's molecule-auditing skill), added in PR
https://github.com/ersilia-os/chemical-library-gifs/pull/1 — requires the "chemgifs" conda env,
editable-installed against a local clone of that branch/PR.

Usage:
    python 70b_compound_gifs.py
"""
import os
import subprocess

import numpy as np
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

FILTERED_HITS_CSV = os.path.join(ROOT, "output", "70_filtering", "filtered_hits.csv")

CHUNKS_DIR = os.path.join(ROOT, "processed", "70b_compound_gifs")
OUTPUT_DIR = os.path.join(ROOT, "output", "70b_compound_gifs")
os.makedirs(CHUNKS_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

N_GROUPS = 6
GIF_KWARGS = ["--style", "rdkit", "--n_rows", "1", "--n_cols", "1", "--size", "512", "--duration", "1000", "--color", "white"]


def main():
    df = pd.read_csv(FILTERED_HITS_CSV)
    print(f"{len(df):,} compounds, split into {N_GROUPS} sequential groups")

    for i, idx in enumerate(np.array_split(np.arange(len(df)), N_GROUPS), start=1):
        group = df.iloc[idx]
        chunk_csv = os.path.join(CHUNKS_DIR, f"group_{i}of{N_GROUPS}.csv")
        group.to_csv(chunk_csv, index=False)

        output_gif = os.path.join(OUTPUT_DIR, f"filtered_hits_group{i}of{N_GROUPS}.gif")
        print(f"Group {i}/{N_GROUPS}: {len(group)} compounds -> {output_gif}")
        subprocess.run(
            ["conda", "run", "-n", "chemgifs", "chemgifs", "-i", chunk_csv, "-o", output_gif, *GIF_KWARGS],
            check=True,
        )


if __name__ == "__main__":
    main()
