#!/usr/bin/env python3
"""
Prints a per-structure progress table for script 80's run: structure predictions done, affinity
predictions done, out of 1,095 compounds each. Read-only -- inspects output/80_nesso1_docking/
on disk, does not touch anything. Meant to be re-run repeatedly over the course of a run (on
nebula, where nesso_results/ actually lives), not a pipeline stage.

No MSA-cache column here (unlike script 74) -- Nesso-1 needs no MSA bootstrap step.

Usage:
    python 81_nesso1_monitor.py [--out-subdir nesso_results]
"""
import argparse
import os

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

STRUCTURE_SEQUENCES_CSV = os.path.join(ROOT, "output", "78_nesso1_prepare_inputs", "structure_sequences.csv")
DOCKING_DIR = os.path.join(ROOT, "output", "80_nesso1_docking")

N_COMPOUNDS = 1095


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-subdir", default="nesso_results",
                         help="Output subdirectory under output/80_nesso1_docking/ to inspect (default: nesso_results)")
    return parser.parse_args()


def structure_progress(structure_id, out_subdir):
    """(n_structures_done, n_affinities_done) for structure_id, both capped implicitly by what's
    actually on disk under predictions/."""
    predictions_dir = os.path.join(DOCKING_DIR, out_subdir, structure_id, "predictions")
    if not os.path.isdir(predictions_dir):
        return 0, 0, False

    compound_ids = os.listdir(predictions_dir)
    n_structures = len(compound_ids)
    n_affinities = sum(
        os.path.isfile(os.path.join(predictions_dir, cid, "affinity.json"))
        for cid in compound_ids
    )
    used_fallback = os.path.isfile(os.path.join(DOCKING_DIR, out_subdir, structure_id, "USED_TRIMMED_DIMER_FALLBACK"))
    return n_structures, n_affinities, used_fallback


def main():
    args = parse_args()
    structures = pd.read_csv(STRUCTURE_SEQUENCES_CSV)

    rows = []
    for _, row in structures.iterrows():
        structure_id = row["structure_id"]
        n_pred, n_aff, used_fallback = structure_progress(structure_id, args.out_subdir)

        if n_aff >= N_COMPOUNDS:
            status = "DONE"
        elif n_pred > 0:
            status = "in progress"
        else:
            status = "not started"

        total_len = row["sequence_length"] + (row["sequence_length_b"] if row["is_dimer"] else 0)
        rows.append({
            "gene": row["gene_name"],
            "structure_id": structure_id,
            "pocket_names": row["pocket_names"],
            "seq_len": total_len,
            "dimer": "yes" if row["is_dimer"] else "no",
            "predictions": f"{n_pred}/{N_COMPOUNDS}",
            "affinities": f"{n_aff}/{N_COMPOUNDS}",
            "used_trimmed_fallback": "yes" if used_fallback else "no",
            "status": status,
        })

    df = pd.DataFrame(rows).sort_values("seq_len")
    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(df.to_string(index=False))

    total_pred = sum(int(r["predictions"].split("/")[0]) for r in rows)
    total_aff = sum(int(r["affinities"].split("/")[0]) for r in rows)
    total_slots = len(rows) * N_COMPOUNDS
    print(f"\nOverall: {total_pred:,}/{total_slots:,} predictions, {total_aff:,}/{total_slots:,} affinities")


if __name__ == "__main__":
    main()
