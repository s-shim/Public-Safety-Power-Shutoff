"""Run dc_check.py on every PEG-1354 solution in a folder and collect one table.

    python dc_check_all.py results_anvil       ->  dc_results/dc_table.csv
Each solution takes about 4 minutes (1,288 contingencies). Finished ones are skipped.
"""
import glob
import json
import os
import re
import subprocess
import sys

folder = sys.argv[1] if len(sys.argv) > 1 else 'results_anvil'
out = 'dc_results'
data = '../../matpower/case1354pegase'
for opt in sorted(glob.glob(os.path.join(folder, 'opt_n1354_inst*.csv'))):
    name = os.path.splitext(os.path.basename(opt))[0].replace('opt_', '')
    if os.path.exists(os.path.join(out, 'dcsummary_%s.json' % name)):
        print('skip', name); continue
    inst = re.search(r'inst(\d+)', name).group(1)
    print('DC check', name, flush=True)
    subprocess.call([sys.executable, 'dc_check.py', '--nbus', '1354', '--data', data,
                     '--lines', '%s/lines_1354_inst%s.csv' % (data, inst), '--opt', opt, '--out', out])

rows = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(out, 'dcsummary_*.json')))]
if rows:
    keys = [k for k in rows[0] if k != 'base_binding_note']
    with open(os.path.join(out, 'dc_table.csv'), 'w') as f:
        f.write(','.join(keys) + '\n')
        for r in rows:
            f.write(','.join(str(r.get(k, '')) for k in keys) + '\n')
    print('wrote', os.path.join(out, 'dc_table.csv'))
