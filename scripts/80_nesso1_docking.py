#!/usr/bin/env python3
### CAUTION: THIS SCRIPT NEEDS TO BE RUN
### WITHIN THE adda4tb CONDA ENVIRONMENT (envs/adda4tb on the IRB SBNB-IRB cluster)
### IN A GPU JOB (sbnb_gpu_3090 / sbnb_gpu_h200 partitions only -- see script 83)
"""
Runs Nesso-1 on script 79's YAML inputs (`nesso predict`) and aggregates the resulting
affinity.json files into one CSV. Unlike Boltz-2's script 73, there is no MSA-bootstrap phase --
Nesso-1 computes and caches ESM-2 embeddings internally per sequence, no external API call
involved. Revised from an earlier per-structure/dimer design to a per-gene design (21 Mtb tRNA
synthetases, see script 78's docstring) -- there is no more dimer complex or OOM-fallback logic
here. The 3090 partition was only empirically confirmed up to alaS (904 residues); 3 of the 21
genes are larger (ileS 1041, leuS 969, valS 886) -- smoke-test the new largest (ileS) before
submitting the full array (script 83), same convention as below.

Scope-limiting flags (--genes, --max-compounds, --out-subdir) let the exact same code path serve
both a cheap single-gene/single-compound smoke test and the full production run -- run the smoke
test first and inspect its output before ever invoking this with no flags.

On the IRB cluster this runs as a SLURM array (script 83, one task per gene on sbnb_gpu_3090) --
since array tasks run concurrently, each invocation is called with --no-aggregate (skip the
per-gene aggregate-and-save step, which would otherwise race across parallel tasks writing the
same CSV), and --aggregate-only is run once at the end, on the login node, after every task has
finished.

Confirmed by reading the recursionpharma/nesso source (nesso/main.py) rather than assumed:
- `nesso predict` takes exactly one DATA path (a file or a directory, non-recursive); its own
  filter_records() skips a compound whose predictions/<id>/affinity.json already exists unless
  --override is passed, so no manual skip-check is needed here (unlike script 72's YAML-write
  skip-check, which still applies at generation time) -- this is also the resumability mechanism:
  a crashed and re-run job never recomputes an already-finished (gene, compound) pair.
- Directory-mode output is flat: <out_dir>/predictions/<compound_id>/affinity.json, no extra
  nesting like Boltz-2's boltz_results_<dirname>/ level.

Usage:
    # Smoke test first: one gene, one compound, isolated output (e.g. via srun).
    envs/adda4tb/bin/python scripts/80_nesso1_docking.py --genes pheS \\
        --max-compounds 1 --out-subdir smoke_test

    # One SLURM array task's worth of work (called from script 83 -- no shared-CSV writes).
    envs/adda4tb/bin/python scripts/80_nesso1_docking.py --genes <gene_name> \\
        --no-aggregate --out-subdir nesso_results

    # After every array task has finished: aggregate everything on disk once.
    envs/adda4tb/bin/python scripts/80_nesso1_docking.py --aggregate-only --out-subdir nesso_results
"""
import argparse
import json
import os
import subprocess
import sys

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

# Resolved relative to the running interpreter (sys.executable), not looked up on PATH: this
# script is always invoked as `envs/adda4tb/bin/python scripts/80_nesso1_docking.py` (SLURM
# script 83 and the srun smoke test alike), never via `conda activate`, so a bare "nesso" would
# not be found on PATH -- its sibling script in the same env's bin/ directory is.
NESSO_BIN = os.path.join(os.path.dirname(sys.executable), "nesso")

INPUT_YAMLS_DIR = os.path.join(ROOT, "output", "79_nesso1_yaml_generation", "input_yamls")
PROTEIN_SEQUENCES_CSV = os.path.join(ROOT, "output", "78_nesso1_prepare_inputs", "protein_sequences.csv")
COMPOUNDS_CSV = os.path.join(ROOT, "output", "71_boltz2_prepare_inputs", "compounds.csv")

OUTPUT_DIR = os.path.join(ROOT, "output", "80_nesso1_docking")
SCOPE_DIR = os.path.join(OUTPUT_DIR, "_scope")
os.makedirs(SCOPE_DIR, exist_ok=True)

AFFINITY_FIELDS = ["affinity_pred_value", "affinity_pred_value1", "affinity_pred_value2",
                   "affinity_logits_binary", "affinity_probability_binary"]
ENTROPY_FIELDS = ["entropy_pp", "entropy_pl", "entropy_ll",
                   "entropy_crop_pp", "entropy_crop_pl", "entropy_crop_ll"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genes", default=None,
                         help="Comma-separated gene names (default: all in protein_sequences.csv)")
    parser.add_argument("--max-compounds", type=int, default=None,
                         help="Limit to the first N compounds per gene, sorted by compound_id (default: all)")
    parser.add_argument("--out-subdir", default="nesso_results",
                         help="Output subdirectory under output/80_nesso1_docking/ (default: nesso_results)")
    parser.add_argument("--no-aggregate", action="store_true",
                         help="Skip the aggregate-and-save step after each gene. Required when "
                              "called from concurrent SLURM array tasks (script 83) -- otherwise "
                              "parallel tasks race on the same results CSV.")
    parser.add_argument("--aggregate-only", action="store_true",
                         help="Skip running any predictions; just aggregate everything already on "
                              "disk under --out-subdir and save. Meant to run once, after every "
                              "array task has finished.")
    return parser.parse_args()


def prepare_scope_dir(gene_name, max_compounds):
    """Returns the directory `nesso predict` should be pointed at: the real input_yamls dir for
    a full run, or a small directory of symlinks to just the first max_compounds (sorted) for a
    restricted run -- avoids symlinking all 1,095 files in the common full-run case."""
    real_dir = os.path.join(INPUT_YAMLS_DIR, gene_name)
    if max_compounds is None:
        return real_dir

    scope_dir = os.path.join(SCOPE_DIR, gene_name)
    os.makedirs(scope_dir, exist_ok=True)
    selected = sorted(os.listdir(real_dir))[:max_compounds]
    for fname in selected:
        link_path = os.path.join(scope_dir, fname)
        if not os.path.islink(link_path):
            os.symlink(os.path.join(real_dir, fname), link_path)
    return scope_dir


def run_nesso_predict(data_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    subprocess.run([
        NESSO_BIN, "predict", data_dir,
        "--out_dir", out_dir,
        "--no_kernels",
        "--num_workers", "4",
    ], check=True)


def parse_prediction(gene_out_dir, compound_id):
    """One result row for compound_id, or None if no affinity output exists (compound
    failed/skipped during processing -- not fatal, just excluded from the summary)."""
    affinity_path = os.path.join(gene_out_dir, "predictions", compound_id, "affinity.json")
    if not os.path.isfile(affinity_path):
        return None
    with open(affinity_path) as f:
        affinity = json.load(f)
    return {field: affinity.get(field) for field in AFFINITY_FIELDS + ENTROPY_FIELDS}


def aggregate_results(out_subdir):
    """Scans every gene directory actually present under out_subdir -- not just this
    invocation's --genes scope -- so re-running the aggregation after an incremental/partial run
    reflects everything computed so far, not just the latest invocation."""
    out_root = os.path.join(OUTPUT_DIR, out_subdir)
    columns = ["gene_name", "compound_id"] + AFFINITY_FIELDS + ENTROPY_FIELDS
    if not os.path.isdir(out_root):
        return pd.DataFrame(columns=columns)

    rows = []
    n_missing = 0
    for gene_name in sorted(os.listdir(out_root)):
        gene_out_dir = os.path.join(out_root, gene_name)
        predictions_dir = os.path.join(gene_out_dir, "predictions")
        if not os.path.isdir(predictions_dir):
            continue
        for compound_id in sorted(os.listdir(predictions_dir)):
            row = parse_prediction(gene_out_dir, compound_id)
            if row is None:
                n_missing += 1
                print(f"  Warning: no affinity result for {gene_name}/{compound_id}, skipping.")
                continue
            row["gene_name"] = gene_name
            row["compound_id"] = compound_id
            rows.append(row)

    print(f"  Aggregated {len(rows)} results, {n_missing} missing/failed.")
    return pd.DataFrame(rows, columns=columns)


def main():
    args = parse_args()

    if args.aggregate_only:
        print(f"Aggregate-only: scanning everything on disk under --out-subdir {args.out_subdir}...")
        result = aggregate_results(args.out_subdir)
        out_path = os.path.join(OUTPUT_DIR, f"{args.out_subdir}_affinity_results.csv")
        result.to_csv(out_path, index=False)
        print(f"Saved {len(result)} rows to {out_path}")
        return

    proteins = pd.read_csv(PROTEIN_SEQUENCES_CSV)
    if args.genes:
        wanted = set(args.genes.split(","))
        proteins = proteins[proteins["gene_name"].isin(wanted)]
    compounds = pd.read_csv(COMPOUNDS_CSV)
    n_total_compounds = len(compounds)
    print(f"Target genes: {proteins['gene_name'].tolist()}")
    print(f"Max compounds per gene: {args.max_compounds or 'all (' + str(n_total_compounds) + ')'}")
    print(f"Output subdir: {args.out_subdir}")

    for _, prow in proteins.iterrows():
        gene_name = prow["gene_name"]
        print(f"\n--- {gene_name} ---")
        data_dir = prepare_scope_dir(gene_name, args.max_compounds)
        out_dir = os.path.join(OUTPUT_DIR, args.out_subdir, gene_name)

        try:
            run_nesso_predict(data_dir, out_dir)
        except subprocess.CalledProcessError as e:
            print(f"  Warning: {gene_name} failed (exit {e.returncode}), skipping.")
            continue

        if args.no_aggregate:
            continue

        # Re-aggregate after every gene (not just once at the end) so a long multi-gene serial
        # run has an up-to-date summary CSV at any point it's interrupted or checked on. Only
        # safe for a single-process serial run -- concurrent SLURM array tasks must pass
        # --no-aggregate and rely on a separate --aggregate-only call once everything is done.
        print("  Aggregating results so far...")
        result = aggregate_results(args.out_subdir)
        out_path = os.path.join(OUTPUT_DIR, f"{args.out_subdir}_affinity_results.csv")
        result.to_csv(out_path, index=False)
        print(f"  Saved {len(result)} rows to {out_path}")


if __name__ == "__main__":
    main()
