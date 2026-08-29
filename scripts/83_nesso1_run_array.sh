#!/bin/bash
# SLURM array job: one task per gene (5 canonical protein sequences -- pheS, pheT, aspS, lysS,
# alaS -- see script 78's docstring for why this replaced the earlier 11-structure/dimer design).
# Modeled on the existing cluster precedent at
# chembl-antimicrobial-models/scripts/09_run_models.sh: an --array job indexed by
# $SLURM_ARRAY_TASK_ID, calling a project-local envs/<name>/bin/python directly, no `conda
# activate` needed.
#
# Single partition now (no separate H200 job): the largest single gene, alaS (904 residues), was
# already confirmed empirically to run cleanly on a 3090 -- there's no dimer complex left to need
# the H200's extra headroom for. --partition/--nodelist pinned to sbnb_gpu_3090 (node irbgcn07)
# deliberately: it's one of only two GPU partitions free for this lab (the other,
# sbnb_gpu_h200/irbgcn10, isn't needed here but remains available if a future gene ever needs
# it). spot_gpu/irb_gpu_3090/irb_gpu_h200 now bill the lab per PI -- do not switch to those
# without explicit sign-off.
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
#SBATCH --array=0-4
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
