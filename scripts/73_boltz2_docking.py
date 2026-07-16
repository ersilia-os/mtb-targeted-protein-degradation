#!/usr/bin/env python3
### CAUTION: THIS SCRIPT NEEDS TO BE RUN
### WITHIN THE boltz CONDA ENVIRONMENT (on nebula)
### IN A GPU MACHINE
"""
Runs Boltz-2 on script 72's YAML inputs: bootstraps one MSA per pocket, runs `boltz predict`, and
aggregates the resulting affinity predictions into one CSV.

Scope-limiting flags (--pockets, --max-compounds, --out-subdir) let the exact same code path serve
both a cheap single-pocket/single-compound smoke test and the full production run -- run the smoke
test first and inspect its output before ever invoking this with no flags (see Usage below).

Confirmed by reading the installed Boltz-2 source (boltz 2.2.1) rather than assumed:
- `boltz predict` takes exactly one `data` path (a file or a directory); resumability is real at
  the manifest level (a compound whose prediction/affinity output already exists is skipped on a
  re-run unless --override is passed), and a single bad/unembeddable compound is logged and
  skipped, not fatal to the whole batch.
- Output is named after `data`'s own path stem, not per compound: a single-file run produces
  <out_dir>/boltz_results_<yaml_stem>/, with the MSA cache at msa/<yaml_stem>_0.csv; a
  directory-mode run (the main prediction step here) produces one shared
  <out_dir>/boltz_results_<dir_name>/predictions/<compound_id>/ per compound inside it.

Usage:
    conda activate boltz

    # Smoke test first: one pocket, one compound, isolated output.
    python 73_boltz2_docking.py --pockets swissmodel_P9WFU3_model_0_pocket_1 \\
        --max-compounds 1 --out-subdir smoke_test

    # Full production run (11 pockets x 1,095 compounds) -- only after the smoke test is reviewed.
    python 73_boltz2_docking.py
"""
import argparse
import json
import os
import shutil
import subprocess
import time

import pandas as pd
import yaml

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

INPUT_YAMLS_DIR = os.path.join(ROOT, "output", "72_boltz2_yaml_generation", "input_yamls")
POCKET_SEQUENCES_CSV = os.path.join(ROOT, "output", "71_boltz2_prepare_inputs", "pocket_sequences.csv")
COMPOUNDS_CSV = os.path.join(ROOT, "output", "71_boltz2_prepare_inputs", "compounds.csv")

OUTPUT_DIR = os.path.join(ROOT, "output", "73_boltz2_docking")
MSA_CACHE_DIR = os.path.join(OUTPUT_DIR, "msa_cache")
MSA_BOOTSTRAP_DIR = os.path.join(OUTPUT_DIR, "msa_bootstrap")
SCOPE_DIR = os.path.join(OUTPUT_DIR, "_scope")
os.makedirs(MSA_CACHE_DIR, exist_ok=True)
os.makedirs(MSA_BOOTSTRAP_DIR, exist_ok=True)

MAX_BOOTSTRAP_ATTEMPTS = 3
# api.colabfold.com rate-limits to ~1 request per 50s per IP (X-Rate-Limit-Limit: 0.02, seen on a
# 429 response); a fixed 2min sleep before every MSA request keeps us comfortably under that.
MSA_REQUEST_BACKOFF_S = 120
AFFINITY_FIELDS = ["affinity_pred_value", "affinity_probability_binary",
                   "affinity_pred_value1", "affinity_probability_binary1",
                   "affinity_pred_value2", "affinity_probability_binary2"]
CONFIDENCE_FIELDS = ["confidence_score", "ptm", "iptm", "ligand_iptm", "complex_plddt"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pockets", default=None,
                         help="Comma-separated pocket names (default: all pockets in pocket_sequences.csv)")
    parser.add_argument("--max-compounds", type=int, default=None,
                         help="Limit to the first N compounds per pocket, sorted by compound_id (default: all)")
    parser.add_argument("--out-subdir", default="boltz_results",
                         help="Output subdirectory under output/73_boltz2_docking/ (default: boltz_results)")
    return parser.parse_args()


def load_pocket_sequences():
    return pd.read_csv(POCKET_SEQUENCES_CSV).set_index("pocket_name")


def load_compounds_by_smiles_length():
    """Compounds sorted ascending by SMILES length -- shorter SMILES are simpler and more likely
    to embed successfully, used to pick MSA-bootstrap candidates."""
    df = pd.read_csv(COMPOUNDS_CSV)
    return df.assign(smiles_len=df["smiles"].str.len()).sort_values("smiles_len")


def build_bootstrap_yaml(sequence, contacts, smiles):
    """Same schema as script 72's production YAMLs but WITHOUT an `msa:` field. A chain with any
    custom msa: path -- even one that doesn't exist yet -- is treated by Boltz-2 as already having
    a custom MSA and skips auto-generation entirely (confirmed the hard way: pointing
    --use_msa_server at a production YAML, which already has msa: baked in, fails immediately with
    FileNotFoundError instead of generating anything). Omitting the key entirely is required for
    --use_msa_server to actually fire."""
    return {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": sequence}},
            {"ligand": {"id": "B", "smiles": smiles}},
        ],
        "constraints": [
            {"pocket": {"binder": "B", "contacts": [["A", pos] for pos in contacts]}}
        ],
        "properties": [
            {"affinity": {"binder": "B"}}
        ],
    }


def bootstrap_msa(pocket_name, sequence, contacts, candidates):
    """Ensures msa_cache/<pocket_name>.csv exists, generating it via a single-compound
    `boltz predict --use_msa_server` call (on a purpose-built, msa-less YAML) if missing. Tries up
    to MAX_BOOTSTRAP_ATTEMPTS candidate compounds (shortest SMILES first) in case one fails to
    embed."""
    cache_path = os.path.join(MSA_CACHE_DIR, f"{pocket_name}.csv")
    if os.path.isfile(cache_path) and os.path.getsize(cache_path) > 0:
        print(f"  MSA cache already present: {cache_path}")
        return

    bootstrap_out = os.path.join(MSA_BOOTSTRAP_DIR, pocket_name)
    os.makedirs(bootstrap_out, exist_ok=True)
    for compound_id, smiles in zip(candidates["compound_id"].head(MAX_BOOTSTRAP_ATTEMPTS),
                                    candidates["smiles"].head(MAX_BOOTSTRAP_ATTEMPTS)):
        yaml_stem = f"bootstrap_{compound_id}"
        yaml_path = os.path.join(bootstrap_out, f"{yaml_stem}.yaml")
        with open(yaml_path, "w") as f:
            yaml.safe_dump(build_bootstrap_yaml(sequence, contacts, smiles), f, sort_keys=False)

        print(f"  Waiting {MSA_REQUEST_BACKOFF_S}s before MSA request (rate-limit backoff)...")
        time.sleep(MSA_REQUEST_BACKOFF_S)
        print(f"  Bootstrapping MSA with {compound_id}...")
        subprocess.run([
            "boltz", "predict", yaml_path,
            "--out_dir", bootstrap_out,
            "--use_msa_server", "--no_kernels",
        ], check=True)

        # named "<yaml_stem>_0.csv", not a fixed "input_0.csv"
        generated = os.path.join(bootstrap_out, f"boltz_results_{yaml_stem}", "msa", f"{yaml_stem}_0.csv")
        if os.path.isfile(generated) and os.path.getsize(generated) > 0:
            shutil.copyfile(generated, cache_path)
            print(f"  MSA cache written: {cache_path}")
            return
        print(f"  Warning: bootstrap with {compound_id} did not produce a usable MSA, trying next candidate.")

    raise RuntimeError(f"Failed to bootstrap MSA for {pocket_name} after {MAX_BOOTSTRAP_ATTEMPTS} attempts.")


def prepare_scope_dir(pocket_name, max_compounds):
    """Returns the directory `boltz predict` should be pointed at: the real input_yamls dir for a
    full run, or a small directory of symlinks to just the first max_compounds (sorted) for a
    restricted run -- avoids symlinking all 1,095 files in the common full-run case."""
    real_dir = os.path.join(INPUT_YAMLS_DIR, pocket_name)
    if max_compounds is None:
        return real_dir

    scope_dir = os.path.join(SCOPE_DIR, pocket_name)
    os.makedirs(scope_dir, exist_ok=True)
    selected = sorted(os.listdir(real_dir))[:max_compounds]
    for fname in selected:
        link_path = os.path.join(scope_dir, fname)
        if not os.path.islink(link_path):
            os.symlink(os.path.join(real_dir, fname), link_path)
    return scope_dir


def check_yaml_count(pocket_name, expected_n):
    real_dir = os.path.join(INPUT_YAMLS_DIR, pocket_name)
    n_found = len([f for f in os.listdir(real_dir) if f.endswith(".yaml")])
    if n_found != expected_n:
        print(f"  Warning: {pocket_name} has {n_found} YAMLs on disk, expected {expected_n} "
              "(possible incomplete rsync from herbert).")


def run_boltz_predict(data_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    subprocess.run([
        "boltz", "predict", data_dir,
        "--out_dir", out_dir,
        "--no_kernels",
        "--num_workers", "8",
        "--preprocessing-threads", "8",
    ], check=True)


def parse_prediction(boltz_results_dir, compound_id):
    """One result row for compound_id, or None if no affinity output exists (compound
    failed/skipped during processing -- not fatal, just excluded from the summary)."""
    pred_dir = os.path.join(boltz_results_dir, "predictions", compound_id)
    affinity_path = os.path.join(pred_dir, f"affinity_{compound_id}.json")
    if not os.path.isfile(affinity_path):
        return None

    with open(affinity_path) as f:
        affinity = json.load(f)
    row = {field: affinity.get(field) for field in AFFINITY_FIELDS}

    confidence_path = os.path.join(pred_dir, f"confidence_{compound_id}_model_0.json")
    if os.path.isfile(confidence_path):
        with open(confidence_path) as f:
            confidence = json.load(f)
        row.update({field: confidence.get(field) for field in CONFIDENCE_FIELDS})

    return row


def aggregate_results(out_subdir):
    """Scans every pocket directory actually present under out_subdir -- not just this
    invocation's --pockets scope -- so re-running the aggregation after an incremental/partial
    run reflects everything computed so far, not just the latest invocation."""
    out_root = os.path.join(OUTPUT_DIR, out_subdir)
    if not os.path.isdir(out_root):
        return pd.DataFrame(columns=["pocket_name", "compound_id"] + AFFINITY_FIELDS + CONFIDENCE_FIELDS)

    rows = []
    n_missing = 0
    for pocket_name in sorted(os.listdir(out_root)):
        pocket_dir = os.path.join(out_root, pocket_name)
        if not os.path.isdir(pocket_dir):
            continue
        # A directory-mode `boltz predict` call names its output after the input directory, not
        # per compound -- one boltz_results_<pocket_name>/ holds every compound's predictions/.
        boltz_results_dir = os.path.join(pocket_dir, f"boltz_results_{pocket_name}")
        predictions_dir = os.path.join(boltz_results_dir, "predictions")
        if not os.path.isdir(predictions_dir):
            continue
        for compound_id in sorted(os.listdir(predictions_dir)):
            row = parse_prediction(boltz_results_dir, compound_id)
            if row is None:
                n_missing += 1
                print(f"  Warning: no affinity result for {pocket_name}/{compound_id}, skipping.")
                continue
            row["pocket_name"] = pocket_name
            row["compound_id"] = compound_id
            rows.append(row)

    print(f"  Aggregated {len(rows)} results, {n_missing} missing/failed.")
    columns = ["pocket_name", "compound_id"] + AFFINITY_FIELDS + CONFIDENCE_FIELDS
    return pd.DataFrame(rows, columns=columns)


def main():
    args = parse_args()
    pocket_sequences = load_pocket_sequences()
    pockets = args.pockets.split(",") if args.pockets else pocket_sequences.index.tolist()
    compounds_by_len = load_compounds_by_smiles_length()
    n_total_compounds = len(compounds_by_len)
    print(f"Target pockets: {pockets}")
    print(f"Max compounds per pocket: {args.max_compounds or 'all (' + str(n_total_compounds) + ')'}")
    print(f"Output subdir: {args.out_subdir}")

    for pocket_name in pockets:
        print(f"\n--- {pocket_name} ---")
        try:
            prow = pocket_sequences.loc[pocket_name]
            contacts = [int(p) for p in prow["pocket_contacts"].split()]
            bootstrap_msa(pocket_name, prow["sequence"], contacts, compounds_by_len)

            if args.max_compounds is None:
                check_yaml_count(pocket_name, n_total_compounds)
            data_dir = prepare_scope_dir(pocket_name, args.max_compounds)

            out_dir = os.path.join(OUTPUT_DIR, args.out_subdir, pocket_name)
            run_boltz_predict(data_dir, out_dir)
        except (RuntimeError, subprocess.CalledProcessError) as e:
            # A pocket-level failure (e.g. the remote MSA server being transiently unreachable)
            # shouldn't kill the rest of the batch, same as a single bad compound doesn't.
            print(f"  Warning: {pocket_name} failed, skipping. Error: {e}")
            continue

        # Re-aggregate after every pocket (not just once at the end) so a long multi-pocket run
        # has an up-to-date summary CSV at any point it's interrupted or checked on.
        print("  Aggregating results so far...")
        result = aggregate_results(args.out_subdir)
        out_path = os.path.join(OUTPUT_DIR, f"{args.out_subdir}_affinity_results.csv")
        result.to_csv(out_path, index=False)
        print(f"  Saved {len(result)} rows to {out_path}")


if __name__ == "__main__":
    main()
