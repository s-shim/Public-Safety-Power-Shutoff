"""DC power-flow validation of an optimal reference network (Reviewers 2.3, 2.7, 3.1).

For the energized buses and lines of a solution (opt_*.csv), solve a DC dispatch LP that
minimizes load shed subject to
    nodal balance   g_i - (d_i - s_i) = sum of flows out of i
    DC flow         f_l = (theta_a - theta_b) / x_l
    line limits     |f_l| <= normal limit (base case) or emergency limit (contingency)
    generation      0 <= g_i <= gmax_i at energized buses
    shed            0 <= s_i <= d_i, where d_i is the demand served in the PSPS solution.
This is done for the base case and for every energized contingency line removed in turn
(corrective re-dispatch allowed). Lines with limit 99999 are treated as unlimited.

    python dc_check.py --nbus 1354 --data ../../matpower/case1354pegase \
        --lines ../../matpower/case1354pegase/lines_1354_inst1.csv \
        --opt results_anvil/opt_n1354_inst1_nine_L2_mf3.csv --out dc_results
Needs numpy, scipy, pandas (no Gurobi license needed).
"""
import argparse
import json
import os
import time

import numpy as np
import pandas as pd
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

ap = argparse.ArgumentParser()
ap.add_argument('--nbus', type=int, required=True)
ap.add_argument('--data', required=True)
ap.add_argument('--lines', required=True)
ap.add_argument('--opt', required=True)
ap.add_argument('--out', default='dc_results')
ap.add_argument('--unlimited', type=float, default=99999)
ap.add_argument('--free', action='store_true',
                help="model's demand definition: 0 <= d_i <= dbar_i free at each energized bus, "
                     "report the max total served against target = frac * sum(dbar)")
ap.add_argument('--frac', type=float, default=0.9)
ap.add_argument('--limits', default=None, help='CSV with Line, Normal/Emergency Flow Limit (MW) overriding the lines file')
args = ap.parse_args()

nodes = pd.read_csv('%s/nodes_%d.csv' % (args.data, args.nbus))
lines = pd.read_csv(args.lines)
if args.limits:
    lim = pd.read_csv(args.limits).set_index('Line')
    for col in ('Normal Flow Limit (MW)', 'Emergency Flow Limit (MW)'):
        lines[col] = lines['Line'].map(lim[col])
sol = pd.read_csv(args.opt)
val = dict(zip(sol['varName'], sol['varVal']))

on_bus = [i for i in nodes['Node'] if val.get('X[%s]' % i, 0) > 0.5]
on_line = [l for l in lines['Line'] if val.get('Z[%s]' % l, 0) > 0.5]
bidx = {b: k for k, b in enumerate(on_bus)}
gmax = dict(zip(nodes['Node'], nodes['max_generation']))
if args.free:
    dbar = dict(zip(nodes['Node'], nodes['Load']))
    dem = {i: max(0.0, dbar[i]) for i in on_bus}
    target = args.frac * nodes['Load'].clip(lower=0).sum()
else:
    dem = {i: max(0.0, val.get('D[%s]' % i, 0.0)) for i in on_bus}
L = lines.set_index('Line').loc[on_line]
gen_bus = [i for i in on_bus if gmax[i] > 0]
ref = 1 if 1 in bidx else on_bus[0]


def solve(active, limit_col):
    """min total shed for the given active lines; returns (shed, n_binding, max_loading)."""
    nb, nl, ng = len(on_bus), len(active), len(gen_bus)
    # variable order: theta (nb), f (nl), g (ng), s (nb)
    nv = nb + nl + ng + nb
    rows, cols, vals, rhs = [], [], [], []
    r = 0
    # DC flow: x_l f_l - theta_a + theta_b = 0
    for k, l in enumerate(active):
        a, b, x = bidx[L.at[l, 'Source']], bidx[L.at[l, 'Target']], L.at[l, 'Reactance']
        rows += [r, r, r]; cols += [nb + k, a, b]; vals += [x, -1.0, 1.0]; rhs.append(0.0); r += 1
    # balance: sum_out f - sum_in f - g_i - s_i = -d_i
    base = r
    for k, l in enumerate(active):
        a, b = bidx[L.at[l, 'Source']], bidx[L.at[l, 'Target']]
        rows += [base + a, base + b]; cols += [nb + k, nb + k]; vals += [1.0, -1.0]
    for k, i in enumerate(gen_bus):
        rows.append(base + bidx[i]); cols.append(nb + nl + k); vals.append(-1.0)
    for i in on_bus:
        rows.append(base + bidx[i]); cols.append(nb + nl + ng + bidx[i]); vals.append(-1.0)
    rhs += [-dem[i] for i in on_bus]
    r = base + nb
    # reference angle
    rows.append(r); cols.append(bidx[ref]); vals.append(1.0); rhs.append(0.0); r += 1
    A = coo_matrix((vals, (rows, cols)), shape=(r, nv)).tocsr()
    lim = [L.at[l, limit_col] for l in active]
    bounds = ([(None, None)] * nb
              + [(None, None) if c >= args.unlimited else (-c, c) for c in lim]
              + [(0.0, gmax[i]) for i in gen_bus]
              + [(0.0, dem[i]) for i in on_bus])
    c = np.zeros(nv); c[nb + nl + ng:] = 1.0
    res = linprog(c, A_eq=A, b_eq=np.array(rhs), bounds=bounds, method='highs')
    if res.status != 0:
        return float('nan'), 0, float('nan')
    f = res.x[nb:nb + nl]
    load = [abs(fv) / cv for fv, cv in zip(f, lim) if cv < args.unlimited]
    binding = sum(v > 1 - 1e-6 for v in load)
    return res.fun, binding, max(load) if load else 0.0


t0 = time.time()
served = sum(dem.values())
base_shed, base_bind, base_max = solve(on_line, 'Normal Flow Limit (MW)')
cont = [l for l in on_line if L.at[l, 'Vulnerable?'] == 1]
records = []
for c in cont:
    act = [l for l in on_line if l != c]
    shed, bind, mx = solve(act, 'Emergency Flow Limit (MW)')
    records.append(dict(contingency=c, shed_mw=shed, binding_lines=bind, max_loading=mx))
df = pd.DataFrame(records)
os.makedirs(args.out, exist_ok=True)
name = os.path.splitext(os.path.basename(args.opt))[0].replace('opt_', '')
df.to_csv(os.path.join(args.out, 'dc_%s.csv' % name), index=False)
summary = dict(run=name, energized_buses=len(on_bus), energized_lines=len(on_line),
               base_binding_note='lines at their limit in the min-shed dispatch',
               served_mw=round(served, 2), limited_lines=int((L['Normal Flow Limit (MW)'] < args.unlimited).sum()),
               base_shed_mw=round(base_shed, 3), base_binding=base_bind, base_max_loading=round(base_max, 3),
               contingencies=len(cont),
               cont_with_shed=int((df['shed_mw'] > 1e-3).sum()) if len(df) else 0,
               cont_max_shed_mw=round(float(df['shed_mw'].max()), 3) if len(df) else 0.0,
               cont_max_shed_pct=round(100 * float(df['shed_mw'].max()) / served, 3) if len(df) and served else 0.0,
               cont_mean_shed_mw=round(float(df['shed_mw'].mean()), 3) if len(df) else 0.0,
               seconds=round(time.time() - t0, 1))
if args.free:
    worst = served - (float(df['shed_mw'].max()) if len(df) else base_shed)
    summary.update(mode='free', target_mw=round(target, 2),
                   base_max_served_mw=round(served - base_shed, 2),
                   cont_min_served_mw=round(worst, 2),
                   base_shortfall_mw=round(max(0.0, target - (served - base_shed)), 2),
                   cont_max_shortfall_mw=round(max(0.0, target - worst), 2),
                   cont_below_target=int((served - df['shed_mw'] < target - 1e-3).sum()) if len(df) else 0)
summary = {k: (v.item() if hasattr(v, 'item') else v) for k, v in summary.items()}
json.dump(summary, open(os.path.join(args.out, 'dcsummary_%s.json' % name), 'w'), indent=1)
print(summary)
