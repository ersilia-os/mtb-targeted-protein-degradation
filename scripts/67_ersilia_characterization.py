#!/usr/bin/env python3
### CAUTION: THIS SCRIPT NEEDS TO BE RUN ON herbert
### WITH THE ersilia CONDA ENV ACTIVE (ersilia CLI on PATH)
"""
Characterizes script 62's aggregated_hits.csv compounds using Ersilia Model Hub models, in place of
the RDKit physicochemical properties script 62 itself computes (docking_utils.compute_properties).

Models are passed at runtime via --models (comma-separated Ersilia model IDs), not hardcoded. Models
used so far for this project's compound set:
    eos12x7 - Spatial Score topological indicator (outputs sps_score, nsps_score - nsps_score is the
              same NSPS metric already used for the drug-likeness filter in script 42)
    eos5jv3 - MycoPermeNet (outputs mycomembrane_permeation - predicted M. tuberculosis outer-
              membrane permeability)
    eos2xeq - antibacterial structural-alert / motif flags (has_pains, has_brenk, is_sim_known_ab,
              nitrofuran_motif, fluoroquinolone_motif, carbepenem_motif, betalactam_motif)
    eos42ez - cytotoxicity (HepG2, HSkMC, IMR90)
    eos3ujl - Mtb cell wall permeability (SVC on Mordred descriptors, permeable-vs-not probability)
    eos8d8a - MycPermCheck (Mtb cell membrane permeability, physicochemical-property classifier)
    eos1lb5 - Mycobacterium cell wall penetration (6-output compendium: lung diffusion, caseum
              diffusion, and cell wall permeation from 3 literature sources incl. eos5jv3's own)

For each model, in order: fetch (skipped if already present in `ersilia catalog --local`) -> serve
-> run -i smiles_input.csv -o {model_id}.csv -> close -> delete (always, regardless of whether the
model pre-existed locally). Resumable: a model whose {model_id}.csv already exists in OUTPUT_DIR is
skipped entirely (fetch/serve/run/close/delete not repeated).

Each model's raw output is left as its own {model_id}.csv (key, input, <model outputs>) - no merging
or combining across models is done here.

Usage:
    conda activate ersilia
    python 67_ersilia_characterization.py --models eos12x7,eos5jv3
"""
import argparse
import os
import subprocess
import sys

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(os.path.join(ROOT, "src"))

AGGREGATED_HITS_CSV = os.path.join(ROOT, "output", "62_aggregate_hits", "aggregated_hits.csv")
OUTPUT_DIR = os.path.join(ROOT, "output", "67_ersilia_characterization")
SMILES_INPUT_CSV = os.path.join(OUTPUT_DIR, "smiles_input.csv")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--models", required=True,
        help="Comma-separated Ersilia model IDs to run, e.g. eos12x7,eos5jv3",
    )
    return parser.parse_args()


def is_fetched_locally(model_id):
    result = subprocess.run(["ersilia", "catalog", "--local"], capture_output=True, text=True, check=True)
    return model_id in result.stdout


def run_ersilia_model(model_id, input_csv, output_csv):
    """fetch (if needed) -> serve -> run -> close -> delete (always)."""
    if is_fetched_locally(model_id):
        print(f"  {model_id}: already fetched locally, skipping fetch.")
    else:
        print(f"  {model_id}: fetching...")
        subprocess.run(["ersilia", "fetch", model_id], check=True)

    try:
        print(f"  {model_id}: serving...")
        subprocess.run(["ersilia", "serve", model_id], check=True)

        print(f"  {model_id}: running inference...")
        subprocess.run(["ersilia", "run", "-i", input_csv, "-o", output_csv], check=True)
    finally:
        print(f"  {model_id}: closing session...")
        subprocess.run(["ersilia", "close"], check=True)

        print(f"  {model_id}: deleting from local storage...")
        subprocess.run(["ersilia", "delete", model_id], check=True)


def main():
    args = parse_args()
    model_ids = [m.strip() for m in args.models.split(",") if m.strip()]
    if not model_ids:
        sys.exit("Error: --models must contain at least one Ersilia model ID.")
    print(f"Models to run: {model_ids}")

    hits = pd.read_csv(AGGREGATED_HITS_CSV, usecols=["smiles"])
    print(f"Aggregated compounds: {len(hits):,}")

    hits.to_csv(SMILES_INPUT_CSV, index=False)
    print(f"Saved SMILES input for Ersilia: {SMILES_INPUT_CSV}")

    for model_id in model_ids:
        print(f"\n--- {model_id} ---")
        model_out_csv = os.path.join(OUTPUT_DIR, f"{model_id}.csv")

        if os.path.isfile(model_out_csv):
            print(f"  {model_id}: {model_out_csv} already exists, skipping.")
            continue

        run_ersilia_model(model_id, SMILES_INPUT_CSV, model_out_csv)

        model_df = pd.read_csv(model_out_csv)
        if len(model_df) != len(hits):
            print(f"  Warning: {model_id} output has {len(model_df):,} rows, expected {len(hits):,}.")
        print(f"  Saved: {model_out_csv} ({len(model_df):,} rows, {len(model_df.columns)} columns)")


if __name__ == "__main__":
    main()
