"""ARR (adaptive randomized rounding) repair of a PSPS reference network for DC line limits.

Starts from a (topologically N-1 connected) reference network, e.g. an optimal MILP solution
opt_*.csv, and searches nearby networks with the ARR procedure of Kody, Xavier, Molzahn, Qiu and
Shim (randomized rounding of a moving seed, decelerated exponential smoothing, random reset).

Feasibility of a candidate network G' = (V', L'):
  (a) G' is the connected component of bus 1 after opening the switches rounded to 0;
  (b) topological N-1: no energized contingency line is a bridge of G';
  (c) DC capacity (corrective): for the base case and for every energized contingency line c
      removed, the LP  max sum_i d_i  s.t. nodal balance, x_l f_l = theta_a - theta_b,
      |f_l| <= u_l, 0 <= g_i <= gbar_i, 0 <= d_i <= dbar_i  serves at least Dmin
      (the dispatch may change after each outage). Dcap(G') = min over these states.

Phase 1 maximizes  min(0, Dcap(G') - Dmin)  (ties: lower risk) until a capacity-feasible
network is found. Phase 2 then minimizes the risk R(G') = sum of risks of L' over
capacity-feasible networks. Line ratings: --limits CSV (Line, Normal/Emergency Flow Limit (MW));
99999 = unrated.

    python arr_cap.py --inst 8 --start results_anvil/opt_n1354_inst8_nine_L2_mf3.csv \
        --limits ../../matpower/case1354pegase/limits_1354_adjusted.csv \
        --workers 32 --phase1 21600 --phase2 64800 --out results_arr

Needs numpy, scipy, pandas, networkx (no Gurobi license).
"""
import argparse
import math
import multiprocessing as mp
import os
import time

import networkx as nx
import numpy as np
import pandas as pd
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components

ap = argparse.ArgumentParser()
ap.add_argument('--inst', type=int, required=True)
ap.add_argument('--start', required=True, help="opt_*.csv with Z[l] of the starting network, or 'full'")
ap.add_argument('--data', default='../../matpower/case1354pegase')
ap.add_argument('--nbus', type=int, default=1354)
ap.add_argument('--limits', default=None)
ap.add_argument('--ratio', type=float, default=0.9, help='demand requirement ratio rho')
ap.add_argument('--workers', type=int, default=os.cpu_count())
ap.add_argument('--phase1', type=float, default=6 * 3600, help='time limit of phase 1 (s)')
ap.add_argument('--phase2', type=float, default=18 * 3600, help='time limit of phase 2 (s)')
ap.add_argument('--kappa', type=float, default=0.8, help='initial seed = 0.5 + kappa (z_start - 0.5)')
ap.add_argument('--kappa_reset', type=float, default=0.6, help='reset seed = 0.5 + kappa_reset (z_best - 0.5)')
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--total', type=float, default=0,
                help='overall time limit (s) for phase 1 followed by phase 2; overrides --phase1/--phase2')
ap.add_argument('--independent', action='store_true',
                help='workers do not share the best solution; each saves its own phase-1 solution and stops '
                     'after phase 1 when --phase2 0')
ap.add_argument('--out', default='results_arr')
ap.add_argument('--unlimited', type=float, default=99999)
args = ap.parse_args()
if args.total > 0:
    args.phase1, args.phase2 = args.total, 0.0      # phase 1 may use all of it; phase 2 runs until the end

UNL = args.unlimited

# ------------------------------------------------------------------ data
nodes = pd.read_csv('%s/nodes_%d.csv' % (args.data, args.nbus))
lines = pd.read_csv('%s/lines_%d_inst%d.csv' % (args.data, args.nbus, args.inst))
if args.limits:
    lim = pd.read_csv(args.limits).set_index('Line')
    for col in ('Normal Flow Limit (MW)', 'Emergency Flow Limit (MW)'):
        lines[col] = lines['Line'].map(lim[col])
NB, NL = len(nodes), len(lines)
bidx = {b: k for k, b in enumerate(nodes['Node'])}
SRC = lines['Source'].map(bidx).to_numpy()
TGT = lines['Target'].map(bidx).to_numpy()
X = lines['Reactance'].to_numpy(float)
UN = lines['Normal Flow Limit (MW)'].to_numpy(float)
UE = lines['Emergency Flow Limit (MW)'].to_numpy(float)
CONT = (lines['Vulnerable?'].to_numpy() == 1)
RISK = lines['Risk'].to_numpy(float)
DBAR = nodes['Load'].clip(lower=0).to_numpy(float)
GBAR = nodes['max_generation'].clip(lower=0).to_numpy(float)
REF = bidx[1]
DMIN = args.ratio * min(DBAR.sum(), GBAR.sum())
LINE_ID = lines['Line'].to_numpy()

if args.start == 'full':
    Z0 = np.ones(NL)                                   # start from the full network
else:
    sol = pd.read_csv(args.start)
    val = dict(zip(sol['varName'], sol['varVal']))
    Z0 = np.array([1.0 if val.get('Z[%s]' % l, 0) > 0.5 else 0.0 for l in LINE_ID])


# ------------------------------------------------------------------ evaluation
def component(z):
    """buses and lines of the component of bus 1 in the graph of closed lines."""
    on = np.flatnonzero(z > 0.5)
    A = csr_matrix((np.ones(len(on)), (SRC[on], TGT[on])), shape=(NB, NB))
    _, lab = connected_components(A, directed=False)
    vb = lab == lab[REF]
    on = on[vb[SRC[on]]]
    return np.flatnonzero(vb), on


def has_contingency_bridge(vbus, on):
    G = nx.Graph()
    mult = {}
    for l in on:
        e = (min(SRC[l], TGT[l]), max(SRC[l], TGT[l]))
        mult[e] = mult.get(e, 0) + 1
        G.add_edge(*e)
    br = set(nx.bridges(G))
    for l in on:
        if CONT[l]:
            e = (min(SRC[l], TGT[l]), max(SRC[l], TGT[l]))
            if mult[e] == 1 and (e in br or (e[1], e[0]) in br):
                return True
    return False


def served(vbus, act, U):
    """max total served demand of the network (vbus, act) with ratings U (corrective DC LP)."""
    nb, nl = len(vbus), len(act)
    pos = -np.ones(NB, dtype=int); pos[vbus] = np.arange(nb)
    gen = vbus[GBAR[vbus] > 0]; ng = len(gen)
    a, b = pos[SRC[act]], pos[TGT[act]]
    k = np.arange(nl)
    # rows 0..nl-1: x f - theta_a + theta_b = 0 ; rows nl..nl+nb-1: balance ; last: theta_ref = 0
    rows = np.concatenate([k, k, k, nl + a, nl + b, nl + pos[gen], nl + np.arange(nb), [nl + nb]])
    cols = np.concatenate([nb + k, a, b, nb + k, nb + k, nb + nl + np.arange(ng),
                           nb + nl + ng + np.arange(nb), [pos[REF]]])
    vals = np.concatenate([X[act], -np.ones(nl), np.ones(nl), np.ones(nl), -np.ones(nl),
                           -np.ones(ng), np.ones(nb), [1.0]])
    A = coo_matrix((vals, (rows, cols)), shape=(nl + nb + 1, nb + nl + ng + nb)).tocsr()
    u = U[act]
    lo = np.concatenate([np.full(nb, -np.inf), np.where(u >= UNL, -np.inf, -u), np.zeros(ng), np.zeros(nb)])
    hi = np.concatenate([np.full(nb, np.inf), np.where(u >= UNL, np.inf, u), GBAR[gen], DBAR[vbus]])
    c = np.zeros(nb + nl + ng + nb); c[nb + nl + ng:] = -1.0
    res = linprog(c, A_eq=A, b_eq=np.zeros(nl + nb + 1), bounds=np.column_stack([lo, hi]), method='highs')
    return -res.fun if res.status == 0 else 0.0


def evaluate(z, stop_below, order):
    """returns (status, Dcap, risk, vbus, on). status: 'topo' infeasible, 'cut' stopped early, 'ok'."""
    vbus, on = component(z)
    risk = RISK[on].sum()
    if not np.any(vbus == REF) or has_contingency_bridge(vbus, on):
        return 'topo', -np.inf, risk, vbus, on
    dcap = min(DBAR[vbus].sum(), GBAR[vbus].sum())
    if dcap < stop_below:
        return 'cut', dcap, risk, vbus, on
    s = served(vbus, on, UN)
    dcap = min(dcap, s)
    if dcap < stop_below:
        return 'cut', dcap, risk, vbus, on
    onset = set(on.tolist())
    cl = [l for l in order if l in onset and CONT[l]]
    seen = set(cl)
    cl += [l for l in on if CONT[l] and l not in seen]
    for l in cl:
        s = served(vbus, on[on != l], UE)
        if s < dcap:
            dcap = s
            order.remove(l) if l in order else None
            order.insert(0, l)               # check the hardest contingencies first next time
        if dcap < stop_below:
            return 'cut', dcap, risk, vbus, on
    return 'ok', dcap, risk, vbus, on


def better(a, b):
    """a, b = (phase-1 value min(0, Dcap - Dmin), risk); larger value, then lower risk."""
    return a[0] > b[0] + 1e-6 or (abs(a[0] - b[0]) <= 1e-6 and a[1] < b[1] - 1e-9)


# ------------------------------------------------------------------ ARR worker
def worker(w, shared, lock, t0, name):
    rng = np.random.default_rng(args.seed * 1000 + w)
    order = []
    with lock:
        zbest = np.frombuffer(shared['z'], dtype=float).copy(); best = shared['best']
        lphase = shared['phase']
    seed = 0.5 + args.kappa * (zbest - 0.5)
    nlocal, trial = 0, 0
    while True:
        el = time.time() - t0
        if args.independent:
            phase = lphase
            if el > args.phase1 + args.phase2 or (phase == 2 and args.phase2 <= 0 and args.total <= 0):
                save(zbest, 'w%d' % w, name)
                return
        else:
            with lock:
                phase = shared['phase']
                if el > args.phase1 + args.phase2:
                    return
                if better(shared['best'], best):                  # adopt the global best
                    best = shared['best']; zbest = np.frombuffer(shared['z'], dtype=float).copy()
        trial += 1
        r0, r1 = rng.random(NL), rng.random(NL)
        z = np.where(seed * r0 < (1 - seed) * (1 - r1), 0.0, 1.0)
        vbus, on = component(z)
        zc = np.zeros(NL); zc[on] = 1.0
        same = np.array_equal(zc, zbest)
        if same:
            nlocal += 1
            rmsd = math.sqrt(np.mean((seed - 0.5) ** 2))
            if rng.random() < min(1.0, nlocal / 20) * rmsd:
                seed = 0.5 + args.kappa_reset * (zbest - 0.5); nlocal = 0
                continue
        else:
            nlocal = 0
            if phase == 1:
                stop = DMIN + best[0]                          # cannot improve below the best Dcap
                if RISK[on].sum() >= best[1] and best[0] >= 0:
                    stop = np.inf
            else:
                stop = DMIN
                if RISK[on].sum() >= best[1] - 1e-9:           # phase 2: only lower risk can improve
                    stop = np.inf
            if np.isfinite(stop):
                st, dcap, risk, vb, onl = evaluate(zc, stop - 1e-6, order)
                if st != 'topo':
                    cand = (min(0.0, dcap - DMIN), risk)
                    ok = (st == 'ok') if phase == 2 else (st == 'ok' or better(cand, best))
                    if ok and better(cand, best) and args.independent:
                        best = cand; zbest = zc
                        if cand[0] >= -1e-6 and lphase == 1:
                            lphase = 2
                            save(zbest, 'phase1_w%d' % w, name)
                        with lock:
                            with open(shared['log'], 'a') as f:
                                f.write('%.1f,%d,%d,%d,%.2f,%.1f,%d,%d\n' % (
                                    time.time() - t0, w, trial, lphase, dcap, risk, len(vb), len(onl)))
                    elif ok and better(cand, best):
                        with lock:
                            if better(cand, shared['best']):
                                shared['best'] = cand; shared['z'] = zc.tobytes()
                                if cand[0] >= -1e-6 and shared['phase'] == 1:
                                    shared['phase'] = 2
                                rec = '%.1f,%d,%d,%d,%.2f,%.1f,%d,%d\n' % (
                                    time.time() - t0, w, trial, shared['phase'], dcap, risk, len(vb), len(onl))
                                with open(shared['log'], 'a') as f:
                                    f.write(rec)
                        best = cand; zbest = zc
        rmsd = math.sqrt(np.mean((seed - 0.5) ** 2))
        alpha = 1.0 / (1.0 + math.exp(4 * rmsd))
        seed = (1 - alpha) * seed + alpha * zbest


def save(z, tag, name):
    vbus, on = component(z)
    rows = [('X[%s]' % nodes['Node'][i], 1.0 if i in set(vbus.tolist()) else 0.0) for i in range(NB)]
    rows += [('Z[%s]' % LINE_ID[l], z[l]) for l in range(NL)]
    pd.DataFrame(rows, columns=['varName', 'varVal']).to_csv(
        os.path.join(args.out, 'opt_%s_%s.csv' % (name, tag)), index=False)


if __name__ == '__main__':
    os.makedirs(args.out, exist_ok=True)
    name = 'arr_' + (('n%d_inst%d_full' % (args.nbus, args.inst)) if args.start == 'full'
                     else os.path.splitext(os.path.basename(args.start))[0].replace('opt_', ''))
    t0 = time.time()
    order = []
    st, dcap, risk, vb, on = evaluate(Z0, -np.inf, order)
    print('start: status %s, Dcap %.1f MW (Dmin %.1f), risk %.0f, buses %d, lines %d, %.1f s'
          % (st, dcap, DMIN, risk, len(vb), len(on), time.time() - t0), flush=True)
    zs = np.zeros(NL); zs[on] = 1.0
    mgr = mp.Manager()
    shared = mgr.dict()
    shared['best'] = (min(0.0, dcap - DMIN), risk) if st != 'topo' else (-np.inf, risk)
    shared['z'] = zs.tobytes()
    shared['phase'] = 2 if shared['best'][0] >= -1e-6 else 1
    shared['log'] = os.path.join(args.out, 'log_%s.csv' % name)
    with open(shared['log'], 'w') as f:
        f.write('time,worker,trial,phase,dcap,risk,buses,lines\n')
        f.write('0,-1,0,%d,%.2f,%.1f,%d,%d\n' % (shared['phase'], dcap, risk, len(vb), len(on)))
    lock = mgr.Lock()
    t0 = time.time()
    procs = [mp.Process(target=worker, args=(w, shared, lock, t0, name)) for w in range(args.workers)]
    for p in procs:
        p.start()
    phase1_saved = args.independent
    while any(p.is_alive() for p in procs):
        time.sleep(30)
        if shared['phase'] == 2 and not phase1_saved:
            save(np.frombuffer(shared['z'], dtype=float), 'phase1', name); phase1_saved = True
        if shared['phase'] == 1 and time.time() - t0 > args.phase1 + 600:
            break
    for p in procs:
        p.join(timeout=60)
    if args.independent:
        lg = pd.read_csv(shared['log'])
        last = lg[lg.worker >= 0].groupby('worker').tail(1)
        print('workers that reached phase 2: %d of %d; best Dcap %.1f MW'
              % ((last.phase == 2).sum(), args.workers, lg.dcap.max()), flush=True)
        raise SystemExit
    zb = np.frombuffer(shared['z'], dtype=float)
    save(zb, 'final', name)
    b = shared['best']
    print('final: phase %d, Dcap - Dmin = %.1f MW, risk %.0f (start %.0f), %.0f s'
          % (shared['phase'], b[0], b[1], risk, time.time() - t0), flush=True)
