"""Collect all finished runs into one table for the response letter.

    python summarize_results.py results_anvil   ->  results_anvil/summary_table.csv
Columns: run settings, final status/objective/bound/gap, time, B&B nodes, peak memory,
and incumbent/gap at 1 h, 6 h, 24 h (from the Gurobi log).
"""
import csv
import glob
import json
import os
import sys

from parse_grblog import parse

STATUS = {2: 'optimal', 9: 'time limit'}
folder = sys.argv[1] if len(sys.argv) > 1 else 'results_anvil'
rows = []
for js in sorted(glob.glob(os.path.join(folder, 'summary_*.json'))):
    s = json.load(open(js))
    log = os.path.join(folder, 'grblog_%s.txt' % s['run'])
    g = parse(log) if os.path.exists(log) else {}
    row = {k: s[k] for k in ['inst', 'families', 'level', 'mipfocus', 'threads', 'gurobi', 'host']}
    row.update(status=STATUS.get(s['status'], s['status']), obj=s['obj'], bound=round(s['bound'], 2),
               gap_pct=round(100 * s['gap'], 3), hours=round(s['runtime'] / 3600, 2),
               bb_nodes=int(s['nodes']), peak_mem_gb=s['peak_mem_gb'])
    for k in ['1h', '6h', '24h']:
        v = g.get(k)
        row['inc_%s' % k] = v[0] if v else ''
        row['gap_%s' % k] = v[2] if v else ''
    rows.append(row)

if not rows:
    sys.exit('no summary_*.json files in %s' % folder)
out = os.path.join(folder, 'summary_table.csv')
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
for r in rows:
    print(r)
print('wrote', out)
