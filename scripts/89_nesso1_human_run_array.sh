#!/bin/bash
# SLURM array job: one task per human tRNA synthetase gene (38 total). Modeled on script 83's
# Mtb array, but on a separate partition -- sbnb_gpu_h200/irbgcn10, not sbnb_gpu_3090/irbgcn07 --
# so this never competes with the Mtb screen's own array (script 83) for CPU/GPU on the same
# node while both may be running. irbgcn10 also has far more memory headroom, useful given
# several human sequences (up to EPRS1's 1,512 residues) are larger than anything tested in the
# Mtb screen. spot_gpu/irb_gpu_3090/irb_gpu_h200 remain off limits -- billed to the lab.
#
# --cpus-per-task kept modest (4) despite this node's extra headroom, same lesson learned from
# the Mtb array: it's shared with other users too.
#
# Each task runs with --no-aggregate (concurrent array tasks would otherwise race on the same
# results CSV); run script 86 with --aggregate-only once, on the login node, after this whole
# array has finished.
#
# Usage (on irb, from the repo root):
#   mkdir -p output/86_nesso1_human_docking/logs
#   sbatch --chdir=$(pwd) scripts/89_nesso1_human_run_array.sh

#SBATCH --job-name=nesso1-human
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --partition=sbnb_gpu_h200
#SBATCH --nodelist=irbgcn10
#SBATCH --array=0-37
#SBATCH --output=output/86_nesso1_human_docking/logs/%x_%a.out
#SBATCH --error=output/86_nesso1_human_docking/logs/%x_%a.err
#SBATCH --requeue

export PYTHONDONTWRITEBYTECODE=1
export PYTHONUNBUFFERED=1

GENE=$(envs/adda4tb/bin/python -c "
import pandas as pd
df = pd.read_csv('output/84_nesso1_human_prepare_inputs/protein_sequences.csv').sort_values('gene_name').reset_index(drop=True)
print(df.loc[$SLURM_ARRAY_TASK_ID, 'gene_name'])
")

echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID -> gene=$GENE"

envs/adda4tb/bin/python -u scripts/86_nesso1_human_docking.py \
    --genes "$GENE" \
    --no-aggregate \
    --out-subdir nesso_results
