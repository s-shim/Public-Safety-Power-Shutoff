"""Run all sensitivity experiments one after another (skips runs that already finished).

    python run_sensitivity.py --phase ratio      demand requirement 80, 85, 90, 95, 100 %
    python run_sensitivity.py --phase switch     switchable share 50 % and 75 % (at 90 %)
    python run_sensitivity.py --phase cluster    clustered risk (at 90 %)
    python run_sensitivity.py --phase all        everything (IEEE-118 first, then IEEE-300)

Results go to results_sens/ (one summary_*.json, grblog_*.txt and opt_*.csv per run).
"""
import argparse
import os
import subprocess
import sys

ap = argparse.ArgumentParser()
ap.add_argument('--phase', choices=['ratio', 'switch', 'cluster', 'all'], default='all')
ap.add_argument('--systems', default='118,300')
ap.add_argument('--n_inst', type=int, default=10)
ap.add_argument('--root', default='../../matpower')
ap.add_argument('--out', default='results_sens')
ap.add_argument('--timelimit', type=float, default=7200)
args = ap.parse_args()

runs = []
for n in [int(x) for x in args.systems.split(',')]:
    data = '%s/case%d' % (args.root, n)
    for i in range(1, args.n_inst + 1):
        if args.phase in ('ratio', 'all'):
            for r in (80, 85, 90, 95, 100):
                runs.append(('n%d_inst%d_r%d' % (n, i, r), n, data, '%s/lines_%d_inst%d.csv' % (data, n, i), r))
        if args.phase in ('switch', 'all'):
            for sh in (50, 75):
                runs.append(('n%d_inst%d_sw%d_r90' % (n, i, sh), n, data,
                             '%s/lines_%d_sw%d_inst%d.csv' % (data, n, sh, i), 90))
        if args.phase in ('cluster', 'all'):
            runs.append(('n%d_inst%d_clu_r90' % (n, i), n, data, '%s/lines_%d_clu_inst%d.csv' % (data, n, i), 90))

os.makedirs(args.out, exist_ok=True)
for k, (label, n, data, lines, r) in enumerate(runs, 1):
    if os.path.exists(os.path.join(args.out, 'summary_%s.json' % label)):
        print('[%d/%d] skip %s (done)' % (k, len(runs), label))
        continue
    print('[%d/%d] run %s' % (k, len(runs), label), flush=True)
    cmd = [sys.executable, 'run_flow.py', '--nbus', str(n), '--data', data, '--lines', lines,
           '--ratio', str(r / 100), '--label', label, '--out', args.out,
           '--timelimit', str(args.timelimit), '--families', 'default', '--mipfocus', '0']
    with open(os.path.join(args.out, 'stdout_%s.txt' % label), 'w') as f:
        rc = subprocess.call(cmd, stdout=f, stderr=subprocess.STDOUT)
    if rc != 0:
        print('   FAILED (see %s/stdout_%s.txt)' % (args.out, label))
print('done')
