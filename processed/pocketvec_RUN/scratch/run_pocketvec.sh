#!/bin/bash
#
#PBS -N pocketvec
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -q YOUR_QUEUE_NAME

cd /home/fcarli/pocketvec_RUN/scratch

# Loads default environment configuration
if [[ -f $HOME/.bashrc ]]
then
  source $HOME/.bashrc
fi

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 singularity exec /home/fcarli/pocketvec_RUN/rDock_image_2.simg python /home/fcarli/pocketvec_RUN/gen_fps_rDock_center_HT.py 1 /home/fcarli/pocketvec_RUN/scratch/array.pkl
