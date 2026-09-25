#!/usr/bin/env bash
# Numerical verification of the explorer.
#  1. JS statistics library vs SciPy/NumPy (24 quantities)
#  2. Independent Python re-implementation of alignment/season/lag/detrend/composites/reclassification
#     (compare its output with the browser: window.__setState(cfg); window.__APP — see README)
#  3. LMR ensemble correlation summary
set -euo pipefail
cd "$(dirname "$0")"
echo "== 1. stats.js vs scipy"
node stats_node.js > js_out.json
python3 - <<'EOF'
import json, numpy as np
from scipy import stats
o = json.load(open('js_out.json')); x, y = np.array(o['x']), np.array(o['y']); n = len(x)
def ar1(v): v = v - v.mean(); return np.sum(v[1:]*v[:-1])/np.sum(v*v)
r1, r2 = ar1(x), ar1(y); ne = min(n, max(3, n*(1-r1*r2)/(1+r1*r2))); r = stats.pearsonr(x, y)[0]
p = 2*stats.t.sf(abs(r)*np.sqrt((ne-2)/(1-r*r)), ne-2)
checks = {'r': (o['r'], r), 'rho': (o['rho'], stats.spearmanr(x, y)[0]), 'neff': (o['neff'], ne), 'p': (o['p'], p),
          'welch_p': (o['welch']['p'], stats.ttest_ind(o['a'], o['b'], equal_var=False).pvalue)}
bad = [k for k, (a, b) in checks.items() if abs(a-b) > 1e-6*max(1, abs(b))]
print('  OK' if not bad else f'  FAIL {bad}', {k: round(v[0], 6) for k, v in checks.items()})
EOF
echo "== 2. independent end-to-end expected values (configs A–D)"
python3 verify_endtoend.py > py_expected.json && python3 -c "import json;d=json.load(open('py_expected.json'));[print(' ',k,'r_ONI',v['ONI']['r'],'r_RONI',v['RONI']['r'],'reclass',v['reclass']) for k,v in d.items()]"
echo "== 3. LMR ensemble (Illimani NH4, 1800–1998)"
python3 verify_lmr.py
