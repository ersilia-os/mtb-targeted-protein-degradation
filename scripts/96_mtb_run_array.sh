#!/bin/bash
# SLURM array job: Mtb counterpart of 96_human_run_array.sh -- one task per Mtb aaRS gene (21
# total, all CRISPR-screen targets, not just the 5 covered by the 12 hand-curated pockets), each
# docking the 1,095 filtered Mtb hits against that gene's own AF2-monomer pockets (script 91
# --organism mtb). Exact same recipe as the human counter-screen (script 96_human_docking.py
# --organism mtb) so Mtb on-target scores are directly comparable to the human off-target scores --
# see that script's own docstring for the full rationale. See 96_human_run_array.sh's own header
# for the shared design notes (no --no-aggregate/--aggregate-only split, spot_gpu authorization,
# --nodelist pinned to the two SBNB-owned nodes).
#
# --time/--mem are carried over from the human array's own calibrated values (48h, 128G) but were
# calibrated on human protein sizes -- sanity-check against the Mtb smoke test before the full
# array, since Mtb protein sizes differ (the human values already had to be raised once, from an
# initial 32G, after 5 unexpectedly large human proteins OOM'd).
#
# Usage (on irb, from the repo root):
#   mkdir -p output/96_mtb_docking/logs
#   sbatch --chdir=$(pwd) scripts/96_mtb_run_array.sh

#SBATCH --job-name=mtb-docking
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --partition=spot_gpu
#SBATCH --nodelist=irbgcn07,irbgcn10
#SBATCH --array=0-20
#SBATCH --output=output/96_mtb_docking/logs/%x_%a.out
#SBATCH --error=output/96_mtb_docking/logs/%x_%a.err
#SBATCH --requeue

export PYTHONDONTWRITEBYTECODE=1
export PYTHONUNBUFFERED=1

GENE=$(envs/adda4tb/bin/python -c "
import pandas as pd
df = pd.read_csv('output/91_mtb_detect_pockets/pocket_detection_data.csv')
genes = sorted(df['Gene name'].unique())
print(genes[$SLURM_ARRAY_TASK_ID])
")

echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID -> gene=$GENE"

envs/unidock_tools/bin/python -u scripts/96_human_docking.py \
    --organism mtb \
    --genes "$GENE" \
    --out-subdir docking_results
