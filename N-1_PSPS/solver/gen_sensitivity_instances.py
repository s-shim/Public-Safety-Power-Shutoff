"""Create the instance files for the sensitivity analyses (Reviewers 1.4, 2.5, 2.11, 3.12).

For each IEEE system and each original instance 1..10:
  * switchable share 50% and 75%: the same risks, but only a random share of lines switchable
      -> lines_<n>_sw50_inst<i>.csv, lines_<n>_sw75_inst<i>.csv
  * spatially clustered risk: three random "fire centres"; lines within a hop radius of a
    centre get high risk U{6..10}, all other lines low risk U{1..3}; the radius is the smallest
    one that covers at least 25% of the lines
      -> lines_<n>_clu_inst<i>.csv
Usage (from N-1_PSPS/solver):  python gen_sensitivity_instances.py
"""
import argparse
import random

import networkx as nx
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--root', default='../../matpower')
ap.add_argument('--systems', default='118,300')
ap.add_argument('--n_inst', type=int, default=10)
ap.add_argument('--seed', type=int, default=2026)
args = ap.parse_args()

for n in [int(x) for x in args.systems.split(',')]:
    folder = '%s/case%d' % (args.root, n)
    nodes = pd.read_csv('%s/nodes_%d.csv' % (folder, n))
    for inst in range(1, args.n_inst + 1):
        base = pd.read_csv('%s/lines_%d_inst%d.csv' % (folder, n, inst))
        rng = random.Random(args.seed * 1000 + n * 100 + inst)

        # switchable share
        for share in (50, 75):
            df = base.copy()
            k = round(len(df) * share / 100)
            sw = set(rng.sample(list(df['Line']), k))
            df['Switchable'] = [1 if l in sw else 0 for l in df['Line']]
            df.to_csv('%s/lines_%d_sw%d_inst%d.csv' % (folder, n, share, inst), index=False)

        # clustered risk
        G = nx.Graph()
        G.add_nodes_from(nodes['Node'])
        G.add_edges_from(zip(base['Source'], base['Target']))
        centres = rng.sample(list(nodes['Node']), 3)
        dist = {}
        for c in centres:
            for v, d in nx.single_source_shortest_path_length(G, c).items():
                dist[v] = min(d, dist.get(v, 10**9))
        line_d = [min(dist.get(a, 10**9), dist.get(b, 10**9)) for a, b in zip(base['Source'], base['Target'])]
        radius = 0
        while sum(d <= radius for d in line_d) < 0.25 * len(base):
            radius += 1
        df = base.copy()
        df['Risk'] = [rng.randint(6, 10) if d <= radius else rng.randint(1, 3) for d in line_d]
        df.to_csv('%s/lines_%d_clu_inst%d.csv' % (folder, n, inst), index=False)
        hot = sum(d <= radius for d in line_d)
        print('case%d inst%d: centres %s radius %d, high-risk lines %d of %d'
              % (n, inst, centres, radius, hot, len(base)))
print('seed', args.seed)
