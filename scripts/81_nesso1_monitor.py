#!/usr/bin/env python3
"""
Prints a per-gene progress table for script 80's run: predictions done, affinity predictions
done, out of 1,095 compounds each (revised for the per-gene design, see script 78's docstring --
21 genes, no more per-structure/dimer columns). Read-only -- inspects output/80_nesso1_docking/ on
disk, does not touch anything. Meant to be re-run repeatedly over the course of a run (on the IRB
cluster, where nesso_results/ actually lives, via script 83's array), not a pipeline stage.

No MSA-cache column here (unlike script 74) -- Nesso-1 needs no MSA bootstrap step.

Also prints this user's live SLURM job state (squeue) alongside the on-disk progress, when run on
a machine with squeue available -- silently skipped otherwise (e.g. if run from a machine without
Slurm, such as the primary dev machine, purely to inspect files rsynced back from the cluster).

Usage:
    python 81_nesso1_monitor.py [--out-subdir nesso_results]
"""
import argparse
import os
import shutil
import subprocess

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

PROTEIN_SEQUENCES_CSV = os.path.join(ROOT, "output", "78_nesso1_prepare_inputs", "protein_sequences.csv")
DOCKING_DIR = os.path.join(ROOT, "output", "80_nesso1_docking")

N_COMPOUNDS = 1095


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-subdir", default="nesso_results",
                         help="Output subdirectory under output/80_nesso1_docking/ to inspect (default: nesso_results)")
    return parser.parse_args()


def gene_progress(gene_name, out_subdir):
    """(n_predictions_done, n_affinities_done) for gene_name, both capped implicitly by what's
    actually on disk under predictions/."""
    predictions_dir = os.path.join(DOCKING_DIR, out_subdir, gene_name, "predictions")
    if not os.path.isdir(predictions_dir):
        return 0, 0

    compound_ids = os.listdir(predictions_dir)
    n_predictions = len(compound_ids)
    n_affinities = sum(
        os.path.isfile(os.path.join(predictions_dir, cid, "affinity.json"))
        for cid in compound_ids
    )
    return n_predictions, n_affinities


def print_slurm_jobs():
    """Best-effort: prints `squeue -u $USER` for the nesso1 job (script 83's --job-name) if
    Slurm is available on this machine. No-op elsewhere (e.g. the primary dev machine)."""
    if shutil.which("squeue") is None:
        return
    user = os.environ.get("USER", "")
    result = subprocess.run(
        ["squeue", "-u", user, "--name=nesso1"],
        capture_output=True, text=True,
    )
    print("SLURM jobs (squeue):")
    print(result.stdout.strip() or "  (none running/queued)")
    print()


def main():
    args = parse_args()
    print_slurm_jobs()
    proteins = pd.read_csv(PROTEIN_SEQUENCES_CSV)

    rows = []
    for _, row in proteins.iterrows():
        gene_name = row["gene_name"]
        n_pred, n_aff = gene_progress(gene_name, args.out_subdir)

        if n_aff >= N_COMPOUNDS:
            status = "DONE"
        elif n_pred > 0:
            status = "in progress"
        else:
            status = "not started"

        rows.append({
            "gene": gene_name,
            "uniprot_ac": row["uniprot_ac"],
            "seq_len": row["sequence_length"],
            "predictions": f"{n_pred}/{N_COMPOUNDS}",
            "affinities": f"{n_aff}/{N_COMPOUNDS}",
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
