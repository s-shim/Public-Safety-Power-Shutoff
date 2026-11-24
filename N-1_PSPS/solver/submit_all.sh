#!/bin/bash
# Submit the 2x2 cut experiment (MIPFocus=3 in every cell) on PEG-1354.
#   bash submit_all.sh test     -> one test run (inst 1, nine families, aggressive)
#   bash submit_all.sh all      -> the remaining runs, in priority order
# Run from N-1_PSPS/solver (where run_peg1354.py and myDictionary.py are).
#
# The Gurobi WLS Academic license allows LANES concurrent sessions (baseline 2).
# Jobs are spread over LANES Slurm "singleton" chains, so at most LANES run at once;
# each chain starts its next job only when the previous one in that chain ends.
mkdir -p logs results_anvil
LANES=${LANES:-4}

# Priority: the 2x2 on the two hardest existing instances first, then new instances.
RUNS="
1:nine:2 8:nine:2 1:all:2 8:all:2 1:nine:1 8:nine:1 1:all:1 8:all:1
11:nine:2 12:nine:2 13:nine:2 11:all:2 12:all:2 13:all:2
11:nine:1 12:nine:1 13:nine:1 11:all:1 12:all:1 13:all:1
"

pick_lane () {   # the lane with the fewest queued/running jobs
  best=1; bestn=999999
  for l in $(seq 1 $LANES); do
    n=$(squeue -u $USER -h -n "psps_lane$l" | wc -l)
    if [ "$n" -lt "$bestn" ]; then best=$l; bestn=$n; fi
  done
  echo $best
}

submit () {   # $1 inst, $2 families, $3 level
  if [ -f "results_anvil/grblog_n1354_inst$1_$2_L$3_mf3.txt" ]; then
    echo "skip inst$1 $2 L$3 (done)"; return
  fi
  if squeue -u $USER -h -o "%k" | grep -qx "inst$1_$2_L$3"; then
    echo "skip inst$1 $2 L$3 (already queued)"; return
  fi
  # loop guard: a run that crashes before Gurobi writes its log is never "done";
  # never submit the same run more than twice
  tries=$(grep -cx "inst$1_$2_L$3" logs/submitted.txt 2>/dev/null)
  if [ "${tries:-0}" -ge 2 ]; then
    echo "skip inst$1 $2 L$3 (submitted twice already; check logs/*.err)"; return
  fi
  lane=$(pick_lane)
  if ! sbatch -J "psps_lane$lane" --comment="inst$1_$2_L$3" --dependency=singleton \
         --export=ALL,RUNARGS="--inst $1 --families $2 --level $3 --mipfocus 3" peg1354.slurm; then
    echo "queue limit reached; the rest will be submitted automatically as jobs finish"
    return 1
  fi
  echo "inst$1_$2_L$3" >> logs/submitted.txt
}

case "$1" in
  test)
    submit 1 nine 2
    ;;
  all)
    for r in $RUNS; do
      IFS=: read i f l <<< "$r"
      submit "$i" "$f" "$l" || break
    done
    ;;
  *)
    echo "usage: bash submit_all.sh test|all"; exit 1 ;;
esac
squeue -u $USER -o "%.10i %.14j %.22k %.8T %.10M %R"
