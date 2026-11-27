"""Original line limits from the UnitCommitment.jl benchmark file (2017-02-01.json).

The lines_1354*.csv files carry 'Normal/Emergency Flow Limit (MW)' columns that were set to 40% of
the benchmark's normal limit (and emergency = normal) for other experiments. For the DC power-flow
validation we use the benchmark's own limits: normal for the base case, emergency (1.1 x normal)
after a contingency; lines without a limit in the benchmark stay unlimited (99999).

    python make_ucjl_limits.py   ->  ../../matpower/case1354pegase/limits_1354_ucjl.csv
"""
import json
import pandas as pd

folder = '../../matpower/case1354pegase'
data = json.load(open(folder + '/2017-02-01.json/2017-02-01.json'))
lines = pd.read_csv(folder + '/lines_1354.csv')
keys = list(data['Transmission lines'])          # same order as lines_1354.csv (checked: reactances match)
rows = []
for i, line in enumerate(lines['Line']):
    v = data['Transmission lines'][keys[i]]
    assert abs(v['Reactance (ohms)'] - lines.at[i, 'Reactance']) < 1e-9
    rows.append((line, v.get('Normal flow limit (MW)', 99999), v.get('Emergency flow limit (MW)', 99999)))
out = pd.DataFrame(rows, columns=['Line', 'Normal Flow Limit (MW)', 'Emergency Flow Limit (MW)'])
out.to_csv(folder + '/limits_1354_ucjl.csv', index=False)
print('limited lines', int((out['Normal Flow Limit (MW)'] < 99999).sum()), 'of', len(out))
