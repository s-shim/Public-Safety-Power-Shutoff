"""Tables for the sensitivity analyses: mean and standard error over instances.

    python summarize_sensitivity.py results_sens  ->  results_sens/sensitivity_table.csv
"""
import glob
import json
import math
import os
import re
import sys
from collections import defaultdict

folder = sys.argv[1] if len(sys.argv) > 1 else 'results_sens'
groups = defaultdict(list)
for js in glob.glob(os.path.join(folder, 'summary_*.json')):
    s = json.load(open(js))
    m = re.match(r'n(\d+)_inst\d+_(.*)$', s['run'])
    if not m:
        continue
    groups[(int(m.group(1)), m.group(2))].append(s)

cols = ['risk_pct', 'served_pct', 'off_buses', 'off_lines', 'off_gens', 'runtime', 'gap']
rows = []
for (n, setting), runs in sorted(groups.items()):
    row = {'system': 'IEEE-%d' % n, 'setting': setting, 'runs': len(runs),
           'optimal': sum(r['status'] == 2 for r in runs)}
    for c in cols:
        v = [r[c] for r in runs]
        mean = sum(v) / len(v)
        se = math.sqrt(sum((x - mean) ** 2 for x in v) / (len(v) - 1) / len(v)) if len(v) > 1 else 0.0
        row[c] = round(mean, 3)
        row[c + '_se'] = round(se, 3)
    rows.append(row)

out = os.path.join(folder, 'sensitivity_table.csv')
with open(out, 'w') as f:
    keys = list(rows[0])
    f.write(','.join(keys) + '\n')
    for r in rows:
        f.write(','.join(str(r[k]) for k in keys) + '\n')
for r in rows:
    print(r['system'], r['setting'], 'n=%d opt=%d' % (r['runs'], r['optimal']),
          'risk %.2f%%' % r['risk_pct'], 'served %.2f%%' % r['served_pct'],
          'off buses %.1f lines %.1f gens %.1f' % (r['off_buses'], r['off_lines'], r['off_gens']),
          'time %.1fs' % r['runtime'])
print('wrote', out)
