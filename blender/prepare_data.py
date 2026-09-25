"""Compute everything the 3-D infographic needs and write blender/enso3d.json.

Same definitions as the 2-D app: proxy year Y vs DJF index of ENSO year Y (lag 0), no detrending,
Pearson r, N_eff (Bretherton 1999) p-values, Meng–Rosenthal–Rubin test for ΔR.
Full-period stats: official CPC ONI/RONI over each core's overlap with 1950–present.
Running stats: ERSST.v6-recomputed ONI/RONI (validated vs CPC, r=0.999), 21-yr centred windows, 1850–present.
"""
import json, os, numpy as np, pandas as pd
from scipy import stats
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src = open(os.path.join(ROOT, 'app/data/dataset.js')).read()
D = json.loads(src[src.index('=') + 1:].strip().rstrip(';'))

def djf(block, name):
    a, y0 = block[name], block['start']
    return {y0 + i // 12: v for i, v in enumerate(a) if i % 12 == 0 and v is not None}
CPC = {k: djf(D['cpc'], k.lower()) for k in ('ONI', 'RONI')}
ERS = {k: djf(D['ersst'], k.lower()) for k in ('ONI', 'RONI')}

def ar1(v): v = v - v.mean(); return np.sum(v[1:] * v[:-1]) / np.sum(v * v)
def corr(x, y):
    n = len(x); r = np.corrcoef(x, y)[0, 1]
    ne = min(n, max(3, n * (1 - ar1(x) * ar1(y)) / (1 + ar1(x) * ar1(y))))
    p = 2 * stats.t.sf(abs(r) * np.sqrt((ne - 2) / (1 - r * r)), ne - 2)
    return r, p, ne
def mrr(r1, r2, r12, n):
    rm2 = (r1 ** 2 + r2 ** 2) / 2; f = min(1, (1 - r12) / (2 * (1 - rm2))); h = (1 - f * rm2) / (1 - rm2)
    Z = (np.arctanh(r1) - np.arctanh(r2)) * np.sqrt((n - 3) / (2 * (1 - r12) * h)); return 2 * stats.norm.sf(abs(Z))

VAR = {'ILL99': ('ions', 'NH4')}
out = {'cores': [], 'meta': {}}
for c in D['cores']:
    vk = VAR.get(c['id'], ('d18O', None))
    v = c['vars'][vk[0]] if vk[1] is None else c['vars'][vk[0]][vk[1]]
    P = dict(zip(v['years'], v['values']))
    rec = dict(id=c['id'], short={'QUE13': 'Q03', 'QUE18': 'Q18', 'HUACOL': 'HC', 'HUASUM': 'HS', 'HUA93': 'H93', 'ILL99': 'I99', 'ILL17': 'I17'}[c['id']],
               name=c['name'], site=c['site'], lat=c['lat'], lon=c['lon'], elev=c['elev'], ref=c['ref'],
               proxy=f"{v['label']} ({v['units']})", years=[v['years'][0], v['years'][-1]])
    # full period vs official CPC
    Y = [y for y in sorted(P) if P[y] is not None and y in CPC['ONI'] and y in CPC['RONI']]
    p = np.array([P[y] for y in Y]); full = {}
    for k in ('ONI', 'RONI'):
        r, pv, ne = corr(p, np.array([CPC[k][y] for y in Y])); full[k] = dict(r=round(r, 3), p=float('%.3g' % pv), neff=round(ne, 1))
    r12 = np.corrcoef([CPC['ONI'][y] for y in Y], [CPC['RONI'][y] for y in Y])[0, 1]
    full.update(n=len(Y), y0=Y[0], y1=Y[-1], dR=round(full['RONI']['r'] - full['ONI']['r'], 3),
                p_dR=float('%.3g' % mrr(full['RONI']['r'], full['ONI']['r'], r12, min(full['ONI']['neff'], full['RONI']['neff']))))
    ne = min(full['ONI']['neff'], full['RONI']['neff']); tc = stats.t.ppf(0.975, ne - 2)
    full['rcrit'] = round(float(tc / np.sqrt(ne - 2 + tc * tc)), 3)
    rec['full'] = full
    # running 21-yr, ERSST.v6
    W, h = 21, 10; run = {}
    for yc in range(1851 + h, 2026 - h):
        yy = [y for y in range(yc - h, yc + h + 1) if P.get(y) is not None and y in ERS['ONI']]
        if len(yy) < 0.8 * W: continue
        pp = np.array([P[y] for y in yy])
        run[yc] = [round(float(np.corrcoef(pp, [ERS[k][y] for y in yy])[0, 1]), 3) for k in ('ONI', 'RONI')]
    rec['running'] = run
    out['cores'].append(rec)
    print(c['id'], full['n'], full['ONI']['r'], full['RONI']['r'], 'running', (min(run), max(run)) if run else None)
out['meta']['rcrit21'] = round(float(stats.t.ppf(0.975, 19) / np.sqrt(19 + stats.t.ppf(0.975, 19) ** 2)), 3)
out['meta']['years'] = [1861, 2015]
# Huascarán Col time series (for the "time wall")
c = next(c for c in D['cores'] if c['id'] == 'HUACOL'); v = c['vars']['d18O']
P = dict(zip(v['years'], v['values'])); Y = [y for y in sorted(P) if y in CPC['ONI']]
pz = np.array([P[y] for y in Y]); pz = (pz - pz.mean()) / pz.std(ddof=1)
def epi(block, name, thr=0.5):  # CPC episode flag for Jan of each year
    a, y0 = block[name], block['start']; f = {}
    for sg in (1, -1):
        i = 0
        while i < len(a):
            if a[i] is not None and sg * a[i] >= thr:
                j = i
                while j + 1 < len(a) and a[j + 1] is not None and sg * a[j + 1] >= thr: j += 1
                if j - i + 1 >= 5:
                    for k in range(i, j + 1):
                        if k % 12 == 0: f[y0 + k // 12] = sg
                i = j + 1
            else: i += 1
    return f
fo, fr = epi(D['cpc'], 'oni'), epi(D['cpc'], 'roni')
out['wall'] = dict(years=Y, proxy=[round(float(z), 3) for z in pz], ONI=[CPC['ONI'][y] for y in Y], RONI=[CPC['RONI'][y] for y in Y],
                   clsONI=[fo.get(y, 0) for y in Y], clsRONI=[fr.get(y, 0) for y in Y])
out['other_sites'] = D['other_sites']
# terrain (0.2° grid) for Blender
et = pd.read_csv(os.path.join(ROOT, 'data/raw/etopo_andes.csv'), skiprows=[1])
g = et.pivot(index='latitude', columns='longitude', values='altitude').sort_index()
out['terrain'] = dict(lats=[round(float(x), 2) for x in g.index], lons=[round(float(x), 2) for x in g.columns], z=g.values.astype(int).tolist())
json.dump(out, open(os.path.join(ROOT, 'blender/enso3d.json'), 'w'), separators=(',', ':'))
print('wrote enso3d.json', os.path.getsize(os.path.join(ROOT, 'blender/enso3d.json')) // 1024, 'kB; rcrit21', out['meta']['rcrit21'])
