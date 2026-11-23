"""Generate new PEG-1354 risk realizations with a fixed seed (Reviewer 3.10).

Same rule as matpower/case1354pegase/generateinstances.py: every line switchable,
risk drawn uniformly from {1,...,10}. Only the seed is new.

Usage (from the repository root):
    python gen_new_instances.py --data matpower/case1354pegase --first 11 --count 3 --seed 2026
"""
import argparse
import copy
import random

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--data', default='matpower/case1354pegase')
ap.add_argument('--nbus', type=int, default=1354)
ap.add_argument('--first', type=int, default=11)
ap.add_argument('--count', type=int, default=3)
ap.add_argument('--seed', type=int, default=2026)
args = ap.parse_args()

lines = pd.read_csv('%s/lines_%s.csv' % (args.data, args.nbus))
rng = random.Random(args.seed)
for inst in range(args.first, args.first + args.count):
    instance = copy.deepcopy(lines)
    instance['Switchable'] = 1
    instance['Risk'] = [int(rng.random() * 10 + 1) for _ in range(len(lines))]
    out = '%s/lines_%s_inst%s.csv' % (args.data, args.nbus, inst)
    instance.to_csv(out, index=False)
    print(out, 'total risk =', instance['Risk'].sum())
print('seed', args.seed)
