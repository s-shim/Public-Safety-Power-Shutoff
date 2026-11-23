#!/bin/bash
# Run the new PEG-1354 instances one after another on this VM (no Slurm on Jetstream2).
#   bash run_jet.sh nine     # VM 1: nine families, aggressive (L2), instances 11, 12, 13
#   bash run_jet.sh all      # VM 2: all families,  aggressive (L2), instances 11, 12, 13
# Both on ONE r3.xl VM (977 GB; one run peaks near 317 GB on Anvil): run each in its own tmux window with NFS=250.
# Same settings as Anvil (peg1354.slurm): MIPFocus 3, GomoryPasses 100, 32 threads, 168,000 s.
# Run it inside tmux so it keeps going when the browser tab is closed. Finished runs are skipped.
FAM=${1:?usage: bash run_jet.sh nine|all}
LEVEL=${LEVEL:-2}
NFS=${NFS:-400}        # Gurobi NodefileStart (GB); use NFS=250 when two runs share one 1 TB VM
INSTS=${INSTS:-"11 12 13"}
source $HOME/psps-env/bin/activate
export GRB_LICENSE_FILE=$HOME/gurobi.lic
cd "$(dirname "$0")"
mkdir -p results_jet logs
VOLROOT=$(df --output=target . | tail -1)
for i in $INSTS; do
  if [ -f "results_jet/summary_n1354_inst${i}_${FAM}_L${LEVEL}_mf3.json" ]; then
    echo "skip inst$i $FAM L$LEVEL (done)"; continue
  fi
  NODEDIR=$VOLROOT/grb_nodefiles/inst${i}_${FAM}
  mkdir -p "$NODEDIR"
  echo "start inst$i $FAM L$LEVEL  $(date)" | tee -a logs/jet.log
  python run_peg1354.py --inst $i --families $FAM --level $LEVEL --mipfocus 3 \
      --data ../../matpower/case1354pegase --out results_jet \
      --threads 32 --timelimit 168000 --nodefilestart $NFS --nodefiledir "$NODEDIR" \
      > logs/jet_inst${i}_${FAM}_L${LEVEL}.out 2>&1
  echo "end   inst$i $FAM L$LEVEL  $(date)  exit $?" | tee -a logs/jet.log
  rm -rf "$NODEDIR"
done
echo "all done $(date)" | tee -a logs/jet.log
