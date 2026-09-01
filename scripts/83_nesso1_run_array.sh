#!/bin/bash
# SLURM array job: one task per gene (all 21 Mtb tRNA synthetases -- see script 78's docstring
# for why this replaced the earlier 11-structure/dimer design and the earlier 5-gene subset).
# Modeled on the existing cluster precedent at
# chembl-antimicrobial-models/scripts/09_run_models.sh: an --array job indexed by
# $SLURM_ARRAY_TASK_ID, calling a project-local envs/<name>/bin/python directly, no `conda
# activate` needed.
#
# Single partition (no separate H200 job): only confirmed empirically up to alaS (904 residues).
# 3 of the 21 genes are larger (ileS 1041, leuS 969, valS 886) -- smoke-test the new largest
# (ileS) via script 80's --genes/--max-compounds flags before submitting this array; fall back to
# sbnb_gpu_h200 (as the human counter-screen, script 89, already does for its larger sequences)
# if it OOMs. --partition/--nodelist pinned to sbnb_gpu_3090 (node irbgcn07) deliberately: it's
# one of only two GPU partitions free for this lab (the other, sbnb_gpu_h200/irbgcn10, is used by
# the human counter-screen, script 89). spot_gpu/irb_gpu_3090/irb_gpu_h200 now bill the lab per
# PI -- do not switch to those without explicit sign-off.
#
# Each task runs with --no-aggregate (concurrent array tasks would otherwise race on the same
# results CSV); run script 80 with --aggregate-only once, on the login node, after this whole
# array has finished.
#
# Usage (on irb, from the repo root):
#   mkdir -p output/80_nesso1_docking/logs
#   sbatch --chdir=$(pwd) scripts/83_nesso1_run_array.sh

#SBATCH --job-name=nesso1
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --partition=sbnb_gpu_3090
#SBATCH --nodelist=irbgcn07
#SBATCH --array=0-20
#SBATCH --output=output/80_nesso1_docking/logs/%x_%a.out
#SBATCH --error=output/80_nesso1_docking/logs/%x_%a.err
#SBATCH --requeue

export PYTHONDONTWRITEBYTECODE=1
export PYTHONUNBUFFERED=1

GENE=$(envs/adda4tb/bin/python -c "
import pandas as pd
df = pd.read_csv('output/78_nesso1_prepare_inputs/protein_sequences.csv').sort_values('gene_name').reset_index(drop=True)
print(df.loc[$SLURM_ARRAY_TASK_ID, 'gene_name'])
")

echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID -> gene=$GENE"

envs/adda4tb/bin/python -u scripts/80_nesso1_docking.py \
    --genes "$GENE" \
    --no-aggregate \
    --out-subdir nesso_results
