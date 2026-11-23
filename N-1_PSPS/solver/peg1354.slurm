#!/bin/bash
# One PEG-1354 run on an Anvil highmem node (48 h QOS limit).
# Submitted by submit_all.sh; RUNARGS carries the instance and cut configuration.
#SBATCH -A mth260052
#SBATCH -p highmem
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH -t 48:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

source $HOME/psps-env/bin/activate
source $SLURM_SUBMIT_DIR/lic.sh

# Gurobi node files go to fast scratch, one folder per job
NODEDIR=$SCRATCH/grb_nodefiles/$SLURM_JOB_ID
mkdir -p "$NODEDIR"

cd $SLURM_SUBMIT_DIR
sbatch -d afterany:$SLURM_JOB_ID topup.sh
echo "host $(hostname)  start $(date)  args: $RUNARGS"
python run_peg1354.py $RUNARGS \
    --data ../../matpower/case1354pegase \
    --out results_anvil \
    --threads 32 \
    --timelimit 168000 \
    --nodefilestart 400 \
    --nodefiledir "$NODEDIR"
echo "end $(date)"
rm -rf "$NODEDIR"

