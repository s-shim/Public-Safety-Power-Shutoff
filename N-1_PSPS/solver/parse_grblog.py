"""Parse Gurobi logs: model size, presolve size, and incumbent/bound/gap at time checkpoints."""
import re, sys, os, glob, csv

NODE_RE = re.compile(r'^\s*[H*]?\s*\d+\s+\d+.*?\s(\S+)\s+(\S+)\s+([\d.]+)%\s+\S+\s+(\d+)s\s*$')
CHECK = [3600, 6 * 3600, 24 * 3600, 48 * 3600]


def parse(path):
    d = {'file': os.path.basename(path)}
    rows = []
    with open(path, errors='replace') as f:
        for line in f:
            m = re.search(r'Optimize a model with (\d+) rows, (\d+) columns and (\d+) nonzeros', line)
            if m and 'rows0' not in d:
                d['rows0'], d['cols0'], d['nnz0'] = map(int, m.groups())
            m = re.search(r'^Presolved: (\d+) rows, (\d+) columns, (\d+) nonzeros', line)
            if m:
                d['rows1'], d['cols1'], d['nnz1'] = map(int, m.groups())
            m = re.search(r'^Presolve time: ([\d.]+)s', line)
            if m:
                d['presolve_s'] = float(m.group(1))
            m = re.search(r'using up to (\d+) threads', line)
            if m:
                d['threads'] = int(m.group(1))
            m = re.search(r'^Root relaxation: objective ([\d.e+-]+)', line)
            if m:
                d['root'] = float(m.group(1))
            m = re.search(r'^Explored (\d+) nodes .* in ([\d.]+) seconds', line)
            if m:
                d['nodes'], d['time_s'] = int(m.group(1)), float(m.group(2))
            m = re.search(r'^Best objective ([\d.e+-]+), best bound ([\d.e+-]+), gap ([\d.]+)%', line)
            if m:
                d['final_obj'], d['final_bd'], d['final_gap'] = map(float, m.groups())
            if line.startswith('Optimal solution found'):
                d['status'] = 'optimal'
            if line.startswith('Time limit reached'):
                d['status'] = 'time limit'
            m = NODE_RE.match(line)
            if m:
                inc, bd, gap, t = m.groups()
                try:
                    rows.append((int(t), float(inc), float(bd), float(gap)))
                except ValueError:
                    pass
    d.setdefault('status', 'unfinished (log ends)')
    if rows:
        d['last_t'] = rows[-1][0]
    for c in CHECK:
        prior = [r for r in rows if r[0] <= c]
        key = f'{c // 3600}h'
        if not prior:
            d[key] = None
        elif 'time_s' in d and d['time_s'] <= c:
            d[key] = (d.get('final_obj'), d.get('final_bd'), d.get('final_gap'), 'solved')
        elif rows[-1][0] < c and 'time_s' not in d:
            d[key] = (rows[-1][1], rows[-1][2], rows[-1][3], 'log ended before')
        else:
            r = prior[-1]
            d[key] = (r[1], r[2], r[3], '')
    return d


if __name__ == '__main__':
    out = []
    for pat in sys.argv[1:]:
        for p in sorted(glob.glob(pat)):
            out.append(parse(p))
    for d in out:
        print(d)
