#!/bin/bash
# SLURM array job: one task per human aaRS gene (38 total), each docking the 1,095 filtered Mtb
# hits against that gene's own pockets (avg. ~10/gene, 389 total). Modeled on scripts 83/89's own
# per-gene array convention, adapted for Uni-Dock instead of Nesso-1's ESM-based model: uses the
# "unidock_tools" env (not "adda4tb") and the lab's GPU partitions since Uni-Dock docking itself
# needs a GPU (script 94's receptor prep is CPU-only and already ran on the login node).
#
# Unlike scripts 83/89, no --no-aggregate/--aggregate-only split: script 96 writes one independent
# report.csv per pocket (no shared file for concurrent array tasks to race on), so there's nothing
# to aggregate until script 97 merges everything, once, after the whole array finishes.
#
# Partition: spot_gpu, explicitly authorized by the user for this job (overriding scripts 83/89's
# own caution about spot_gpu billing the lab per PI without sign-off; that sign-off is this line)
# -- but --nodelist is pinned to ONLY the two SBNB-owned nodes within that partition (irbgcn07 =
# sbnb_gpu_3090, irbgcn10 = sbnb_gpu_h200), per explicit user instruction: spot_gpu also spans 3
# nodes owned by other labs (irbgcn[01,06,08]), which must not be used. So this still lets the
# 38-task array run two-at-a-time (rather than queueing behind one pinned node), just not the
# further parallelism a full 5-node partition would allow.
#
# --time is a rough estimate (no prior Uni-Dock run on this cluster to calibrate against) -- adjust
# after the smoke test / first array pass shows real per-pocket timing.
#
# --mem=128G (raised from an initial 32G after the first full-array run): the same 5 largest
# proteins that needed the valence fix in script 94 (EPRS1, IARS1, LARS1, VARS1, plus VARS2) hit
# SLURM's OOM killer within their first pocket's docking call at 32G. Both irbgcn07 (515GB total)
# and irbgcn10 (2TB total) have ample headroom for this even running 2 tasks/node concurrently.
#
# Usage (on irb, from the repo root):
#   mkdir -p output/96_human_docking/logs
#   sbatch --chdir=$(pwd) scripts/96_human_run_array.sh

#SBATCH --job-name=human-docking
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --partition=spot_gpu
#SBATCH --nodelist=irbgcn07,irbgcn10
#SBATCH --array=0-37
#SBATCH --output=output/96_human_docking/logs/%x_%a.out
#SBATCH --error=output/96_human_docking/logs/%x_%a.err
#SBATCH --requeue

export PYTHONDONTWRITEBYTECODE=1
export PYTHONUNBUFFERED=1

GENE=$(envs/adda4tb/bin/python -c "
import pandas as pd
df = pd.read_csv('output/91_human_detect_pockets/pocket_detection_data.csv')
genes = sorted(df['Gene name'].unique())
print(genes[$SLURM_ARRAY_TASK_ID])
")

echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID -> gene=$GENE"

envs/unidock_tools/bin/python -u scripts/96_human_docking.py \
    --genes "$GENE" \
    --out-subdir docking_results
