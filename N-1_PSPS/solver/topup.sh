#!/bin/bash
#SBATCH -A mth260052
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH -o logs/topup_%j.out
module purge 2>/dev/null; unset PYTHONPATH PYTHONHOME
cd $SLURM_SUBMIT_DIR
bash submit_all.sh all >> logs/autosubmit.log 2>&1
