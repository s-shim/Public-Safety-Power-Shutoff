#!/bin/bash
# One-time setup on a Jetstream2 r3.xl VM (Ubuntu), run from N-1_PSPS/solver of the cloned repo:
#   bash setup_jetstream.sh
# Installs Python + gurobipy 11.0.2 (same as Anvil) and generates instances 11-13 (seed 2026).
set -e
cd "$(dirname "$0")"
sudo apt-get update -y
sudo apt-get install -y python3-venv tmux
ENV=$HOME/psps-env
[ -d "$ENV" ] || python3 -m venv "$ENV"
source "$ENV/bin/activate"
pip install --upgrade pip
pip install "gurobipy==11.0.2" pandas networkx
python gen_new_instances.py --data ../../matpower/case1354pegase --first 11 --count 3 --seed 2026
md5sum ../../matpower/case1354pegase/lines_1354_inst1[123].csv
python -c "import gurobipy; print('gurobipy', gurobipy.gurobi.version())"
echo "Next: nano ~/gurobi.lic (Northwestern WLS license), then: tmux new -s psps ; bash run_jet.sh nine (or all)"
