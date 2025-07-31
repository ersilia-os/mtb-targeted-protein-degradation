#!/bin/bash
#

#SBATCH -J mtb_docking
#SBATCH --chdir=/aloy/home/acomajuncosa/Ersilia/mtb
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-827%3
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --output=./processed/unidock_docking/training_outputs/%x_%a.out
#SBATCH -p sbnb_cpu_sphr,sbnb_cpu_zen3

# Loads default environment configuration
export SINGULARITYENV_LD_LIBRARY_PATH=$LD_LIBRARY_PATH #/.singularity.d/libs
export SINGULARITY_BINDPATH="/home/sbnb:/aloy/home,/data/sbnb/data:/aloy/data,/data/sbnb/scratch:/aloy/scratch"

# Load cuda libraries
export LD_LIBRARY_PATH=/apps/manual/software/CUDA/11.6.1/lib64:/apps/manual/software/CUDA/11.6.1/targets/x86_64-linux/lib:/apps/manual/software/CUDA/11.6.1/extras/CUPTI/lib64/:/apps/manual/software/CUDA/11.6.1/nvvm/lib64/:$LD_LIBRARY_PATH

alpha=($(seq 0 827))

singularity exec --cleanenv /apps/singularity/ood_images/docker_irb_intel-optimized-tensorflow-avx512-2.13-pip-conda-jupyter-v6.sif ./scripts/26_train_models/run.sh ${alpha[$SLURM_ARRAY_TASK_ID]}