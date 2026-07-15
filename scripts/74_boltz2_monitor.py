#!/usr/bin/env python3
"""
Prints a per-pocket progress table for script 73's run: MSA cache status, structure predictions
done, affinity predictions done, out of 1,095 compounds each. Read-only -- inspects
output/73_boltz2_docking/ on disk, does not touch anything. Meant to be re-run repeatedly over the
course of a multi-day run (on nebula, where boltz_results/ actually lives), not a pipeline stage.

Usage:
    python 74_boltz2_monitor.py [--out-subdir boltz_results]
"""
import argparse
import os

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

POCKET_SEQUENCES_CSV = os.path.join(ROOT, "output", "71_boltz2_prepare_inputs", "pocket_sequences.csv")
DOCKING_DIR = os.path.join(ROOT, "output", "73_boltz2_docking")
MSA_CACHE_DIR = os.path.join(DOCKING_DIR, "msa_cache")

N_COMPOUNDS = 1095


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-subdir", default="boltz_results",
                         help="Output subdirectory under output/73_boltz2_docking/ to inspect (default: boltz_results)")
    return parser.parse_args()


def pocket_progress(pocket_name, out_subdir):
    """(n_structures_done, n_affinities_done) for pocket_name, both capped implicitly by what's
    actually on disk under predictions/."""
    predictions_dir = os.path.join(DOCKING_DIR, out_subdir, pocket_name, f"boltz_results_{pocket_name}", "predictions")
    if not os.path.isdir(predictions_dir):
        return 0, 0

    compound_ids = os.listdir(predictions_dir)
    n_structures = len(compound_ids)
    n_affinities = sum(
        os.path.isfile(os.path.join(predictions_dir, cid, f"affinity_{cid}.json"))
        for cid in compound_ids
    )
    return n_structures, n_affinities


def main():
    args = parse_args()
    pockets = pd.read_csv(POCKET_SEQUENCES_CSV)

    rows = []
    for _, row in pockets.iterrows():
        pocket_name = row["pocket_name"]
        msa_cached = os.path.isfile(os.path.join(MSA_CACHE_DIR, f"{pocket_name}.csv"))
        n_struct, n_aff = pocket_progress(pocket_name, args.out_subdir)

        if n_aff >= N_COMPOUNDS:
            status = "DONE"
        elif n_struct > 0 or msa_cached:
            status = "in progress"
        else:
            status = "not started"

        rows.append({
            "gene": row["gene_name"],
            "site_type": row["site_type"],
            "pocket_name": pocket_name,
            "seq_len": row["sequence_length"],
            "msa_cached": "yes" if msa_cached else "no",
            "structures": f"{n_struct}/{N_COMPOUNDS}",
            "affinities": f"{n_aff}/{N_COMPOUNDS}",
            "status": status,
        })

    df = pd.DataFrame(rows).sort_values("seq_len")
    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(df.to_string(index=False))

    total_struct = sum(int(r["structures"].split("/")[0]) for r in rows)
    total_aff = sum(int(r["affinities"].split("/")[0]) for r in rows)
    total_slots = len(rows) * N_COMPOUNDS
    print(f"\nOverall: {total_struct:,}/{total_slots:,} structures, {total_aff:,}/{total_slots:,} affinities")


if __name__ == "__main__":
    main()
