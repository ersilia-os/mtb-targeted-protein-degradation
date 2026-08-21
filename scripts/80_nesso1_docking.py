#!/usr/bin/env python3
### CAUTION: THIS SCRIPT NEEDS TO BE RUN
### WITHIN THE nesso CONDA ENVIRONMENT (on nebula)
### IN A GPU MACHINE
"""
Runs Nesso-1 on script 79's YAML inputs (`nesso predict`) and aggregates the resulting
affinity.json files into one CSV. Unlike Boltz-2's script 73, there is no MSA-bootstrap phase --
Nesso-1 computes and caches ESM-2 embeddings internally per sequence, no external API call
involved.

Scope-limiting flags (--structures, --max-compounds, --out-subdir) let the exact same code path
serve both a cheap single-structure/single-compound smoke test and the full production run --
run the smoke test first and inspect its output before ever invoking this with no flags.

Confirmed by reading the recursionpharma/nesso source (nesso/main.py) rather than assumed:
- `nesso predict` takes exactly one DATA path (a file or a directory, non-recursive); its own
  filter_records() skips a compound whose predictions/<id>/affinity.json already exists unless
  --override is passed, so no manual skip-check is needed here (unlike script 72's YAML-write
  skip-check, which still applies at generation time).
- Directory-mode output is flat: <out_dir>/predictions/<compound_id>/affinity.json, no extra
  nesting like Boltz-2's boltz_results_<dirname>/ level.

The 7K98 dimer (script 78's structure_sequences.csv `is_dimer` row) is tried at its full,
untrimmed size first (user decision: Nesso-1 is a different, reportedly coarser/faster
architecture than Boltz-2, which OOM'd on this complex, so its memory ceiling isn't assumed).
Only on a CalledProcessError does this script fall back once to the pre-trimmed chain windows
already validated for Boltz-2 (script 71's DIMER_TRIM_MARGIN_A/B, carried through in script 78's
dimer_trimmed_sequence/_b columns) -- logged loudly, and marked with a USED_TRIMMED_DIMER_FALLBACK
sentinel file in that structure's output dir so script 82 can flag it rather than silently
presenting it as an untrimmed-complex result.

Usage:
    conda activate nesso

    # Smoke test first: one structure, one compound, isolated output.
    python 80_nesso1_docking.py --structures swissmodel_P9WFU3_model_0 \\
        --max-compounds 1 --out-subdir smoke_test

    # Full production run (11 structures x 1,095 compounds) -- only after the smoke test is
    # reviewed.
    python 80_nesso1_docking.py
"""
import argparse
import json
import os
import subprocess

import pandas as pd
import yaml

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

INPUT_YAMLS_DIR = os.path.join(ROOT, "output", "79_nesso1_yaml_generation", "input_yamls")
STRUCTURE_SEQUENCES_CSV = os.path.join(ROOT, "output", "78_nesso1_prepare_inputs", "structure_sequences.csv")
COMPOUNDS_CSV = os.path.join(ROOT, "output", "71_boltz2_prepare_inputs", "compounds.csv")

OUTPUT_DIR = os.path.join(ROOT, "output", "80_nesso1_docking")
SCOPE_DIR = os.path.join(OUTPUT_DIR, "_scope")
DIMER_FALLBACK_YAMLS_DIR = os.path.join(OUTPUT_DIR, "dimer_trimmed_fallback_yamls")
os.makedirs(SCOPE_DIR, exist_ok=True)

AFFINITY_FIELDS = ["affinity_pred_value", "affinity_pred_value1", "affinity_pred_value2",
                   "affinity_logits_binary", "affinity_probability_binary"]
ENTROPY_FIELDS = ["entropy_pp", "entropy_pl", "entropy_ll",
                   "entropy_crop_pp", "entropy_crop_pl", "entropy_crop_ll"]
FALLBACK_MARKER = "USED_TRIMMED_DIMER_FALLBACK"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--structures", default=None,
                         help="Comma-separated structure_ids (default: all in structure_sequences.csv)")
    parser.add_argument("--max-compounds", type=int, default=None,
                         help="Limit to the first N compounds per structure, sorted by compound_id (default: all)")
    parser.add_argument("--out-subdir", default="nesso_results",
                         help="Output subdirectory under output/80_nesso1_docking/ (default: nesso_results)")
    return parser.parse_args()


def prepare_scope_dir(structure_id, max_compounds):
    """Returns the directory `nesso predict` should be pointed at: the real input_yamls dir for
    a full run, or a small directory of symlinks to just the first max_compounds (sorted) for a
    restricted run -- avoids symlinking all 1,095 files in the common full-run case."""
    real_dir = os.path.join(INPUT_YAMLS_DIR, structure_id)
    if max_compounds is None:
        return real_dir

    scope_dir = os.path.join(SCOPE_DIR, structure_id)
    os.makedirs(scope_dir, exist_ok=True)
    selected = sorted(os.listdir(real_dir))[:max_compounds]
    for fname in selected:
        link_path = os.path.join(scope_dir, fname)
        if not os.path.islink(link_path):
            os.symlink(os.path.join(real_dir, fname), link_path)
    return scope_dir


def build_fallback_dimer_yaml(seq_a, seq_b, smiles):
    return {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": seq_a}},
            {"protein": {"id": "B", "sequence": seq_b}},
            {"ligand": {"id": "C", "smiles": smiles}},
        ],
        "properties": [
            {"affinity": {"binder": "C"}}
        ],
    }


def prepare_fallback_dimer_dir(structure_id, dimer_trimmed_sequence, dimer_trimmed_sequence_b,
                                compounds, max_compounds):
    """Writes trimmed-sequence dimer YAMLs (script 78's OOM-fallback columns) into a separate
    directory. Only ever called from the except branch below, after the untrimmed complex has
    already failed -- never generated pre-emptively."""
    fallback_dir = os.path.join(DIMER_FALLBACK_YAMLS_DIR, structure_id)
    os.makedirs(fallback_dir, exist_ok=True)
    compounds_scoped = compounds if max_compounds is None else compounds.head(max_compounds)
    for _, crow in compounds_scoped.iterrows():
        out_path = os.path.join(fallback_dir, f"{crow['compound_id']}.yaml")
        if os.path.isfile(out_path):
            continue
        yaml_dict = build_fallback_dimer_yaml(dimer_trimmed_sequence, dimer_trimmed_sequence_b, crow["smiles"])
        with open(out_path, "w") as f:
            yaml.safe_dump(yaml_dict, f, sort_keys=False)
    return fallback_dir


def run_nesso_predict(data_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    subprocess.run([
        "nesso", "predict", data_dir,
        "--out_dir", out_dir,
        "--no_kernels",
        "--num_workers", "8",
    ], check=True)


def parse_prediction(structure_out_dir, compound_id):
    """One result row for compound_id, or None if no affinity output exists (compound
    failed/skipped during processing -- not fatal, just excluded from the summary)."""
    affinity_path = os.path.join(structure_out_dir, "predictions", compound_id, "affinity.json")
    if not os.path.isfile(affinity_path):
        return None
    with open(affinity_path) as f:
        affinity = json.load(f)
    return {field: affinity.get(field) for field in AFFINITY_FIELDS + ENTROPY_FIELDS}


def aggregate_results(out_subdir):
    """Scans every structure directory actually present under out_subdir -- not just this
    invocation's --structures scope -- so re-running the aggregation after an
    incremental/partial run reflects everything computed so far, not just the latest
    invocation."""
    out_root = os.path.join(OUTPUT_DIR, out_subdir)
    columns = ["structure_id", "compound_id", "used_trimmed_dimer_fallback"] + AFFINITY_FIELDS + ENTROPY_FIELDS
    if not os.path.isdir(out_root):
        return pd.DataFrame(columns=columns)

    rows = []
    n_missing = 0
    for structure_id in sorted(os.listdir(out_root)):
        structure_out_dir = os.path.join(out_root, structure_id)
        predictions_dir = os.path.join(structure_out_dir, "predictions")
        if not os.path.isdir(predictions_dir):
            continue
        used_fallback = os.path.isfile(os.path.join(structure_out_dir, FALLBACK_MARKER))
        for compound_id in sorted(os.listdir(predictions_dir)):
            row = parse_prediction(structure_out_dir, compound_id)
            if row is None:
                n_missing += 1
                print(f"  Warning: no affinity result for {structure_id}/{compound_id}, skipping.")
                continue
            row["structure_id"] = structure_id
            row["compound_id"] = compound_id
            row["used_trimmed_dimer_fallback"] = used_fallback
            rows.append(row)

    print(f"  Aggregated {len(rows)} results, {n_missing} missing/failed.")
    return pd.DataFrame(rows, columns=columns)


def main():
    args = parse_args()
    structures = pd.read_csv(STRUCTURE_SEQUENCES_CSV)
    if args.structures:
        wanted = set(args.structures.split(","))
        structures = structures[structures["structure_id"].isin(wanted)]
    compounds = pd.read_csv(COMPOUNDS_CSV)
    n_total_compounds = len(compounds)
    print(f"Target structures: {structures['structure_id'].tolist()}")
    print(f"Max compounds per structure: {args.max_compounds or 'all (' + str(n_total_compounds) + ')'}")
    print(f"Output subdir: {args.out_subdir}")

    for _, srow in structures.iterrows():
        structure_id = srow["structure_id"]
        print(f"\n--- {structure_id} ---")
        data_dir = prepare_scope_dir(structure_id, args.max_compounds)
        out_dir = os.path.join(OUTPUT_DIR, args.out_subdir, structure_id)
        used_trimmed_fallback = False

        try:
            run_nesso_predict(data_dir, out_dir)
        except subprocess.CalledProcessError as e:
            if not srow["is_dimer"]:
                print(f"  Warning: {structure_id} failed (exit {e.returncode}), skipping.")
                continue
            print(f"  Warning: {structure_id} (untrimmed dimer, "
                  f"{srow['sequence_length']}+{srow['sequence_length_b']} residues) failed "
                  f"with exit {e.returncode} -- likely OOM. Falling back to the pre-trimmed "
                  f"chain windows already validated for Boltz-2 (script 71's "
                  f"DIMER_TRIM_MARGIN_A/B) and retrying once.")
            fallback_dir = prepare_fallback_dimer_dir(
                structure_id, srow["dimer_trimmed_sequence"], srow["dimer_trimmed_sequence_b"],
                compounds, args.max_compounds)
            try:
                run_nesso_predict(fallback_dir, out_dir)
                used_trimmed_fallback = True
            except subprocess.CalledProcessError as e2:
                print(f"  Warning: {structure_id} failed again with the trimmed fallback "
                      f"(exit {e2.returncode}), skipping.")
                continue

        if used_trimmed_fallback:
            open(os.path.join(out_dir, FALLBACK_MARKER), "w").close()

        # Re-aggregate after every structure (not just once at the end) so a long multi-structure
        # run has an up-to-date summary CSV at any point it's interrupted or checked on.
        print("  Aggregating results so far...")
        result = aggregate_results(args.out_subdir)
        out_path = os.path.join(OUTPUT_DIR, f"{args.out_subdir}_affinity_results.csv")
        result.to_csv(out_path, index=False)
        print(f"  Saved {len(result)} rows to {out_path}")


if __name__ == "__main__":
    main()
