"""Independent re-implementation (from raw NOAA files) of the app's alignment, seasonal
aggregation, detrending, lag, N_eff significance, composites and ENSO-year reclassification.
Prints expected values for several configurations; compared against the browser app."""
import json, numpy as np, pandas as pd
from scipy import stats
R = '../data/raw/'
SEAS = ['DJF','JFM','FMA','MAM','AMJ','MJJ','JJA','JAS','ASO','SON','OND','NDJ']
def cpc(fn, col):
    d = {}
    for l in open(R+fn).read().splitlines()[1:]:
        p = l.split()
        if len(p) >= 3: d[(int(p[1]), SEAS.index(p[0])+1)] = float(p[col])
    return d
CPC = {'ONI': cpc('oni.ascii.txt', 3), 'RONI': cpc('RONI.ascii.txt', 2)}
src = open('../app/data/dataset.js').read(); DS = json.loads(src[src.index('=')+1:].strip().rstrip(';'))
def ersst(name):
    a = DS['ersst'][name.lower()]; y0 = DS['ersst']['start']
    return {(y0 + i//12, i%12+1): v for i, v in enumerate(a) if v is not None}
ERS = {'ONI': ersst('ONI'), 'RONI': ersst('RONI')}
def months(Y, season):
    if season == 'djf': return [(Y, 1)]
    if season == 'ann': return [(Y, m) for m in range(1, 13)]
    if season == 'hydro': return [(Y-1, m) for m in range(8, 13)] + [(Y, m) for m in range(1, 8)]
def seas(d, Y, season):
    ms = months(Y, season)
    return np.mean([d[k] for k in ms]) if all(k in d for k in ms) else None
def epi_flags(d, thr):
    keys = sorted(d); f = {}
    for sg in (1, -1):
        run = []
        for k in keys + [None]:
            if k is not None and sg*d[k] >= thr - 1e-9 and (not run or (k[0]*12+k[1]) - (run[-1][0]*12+run[-1][1]) == 1):
                run.append(k)
            else:
                if len(run) >= 5:
                    for q in run: f[q] = sg
                run = [k] if (k is not None and sg*d[k] >= thr-1e-9) else []
    return f
def rd(fn): return pd.read_csv(R+fn, sep='\t', comment='#', na_values=['NA',''])
def proxy(core, var):
    if core == 'HUACOL': df = rd('weber2023-hs2019.txt'); return dict(zip(df.age_CE, df['col_'+var]))
    if core == 'QUE13':  df = rd('quelccaya2013-noaa.txt'); c = {'d18O':'d18OiceSMOW','dust':'dust'}[var]; df = df[['age_CE', c]].dropna(); return dict(zip(df.age_CE, df[c]))
    if core == 'HUA93':  df = rd('thompson1995-annualt-noaa.txt'); c = {'d18O':'d18O_HScore2','dust':'part>0.63_HScore2'}[var]; df = df[['age_CE', c]].dropna(); return dict(zip(df.age_CE, df[c]))
def ar1(v): v = v - v.mean(); return np.sum(v[1:]*v[:-1])/np.sum(v*v)
def detr(t, y): return y - np.polyval(np.polyfit(t, y, 1), t)
def run(core, var, source, season, lag, detrend, y0, y1, thr=0.5):
    P = proxy(core, var); IDX = CPC if source == 'cpc' else ERS
    out = {}
    for nm in ('ONI', 'RONI'):
        Y, p, x = [], [], []
        for yr in sorted(P):
            if yr < y0 or yr > y1 or pd.isna(P[yr]): continue
            v = seas(IDX[nm], yr - lag, season)
            if v is None: continue
            Y.append(yr); p.append(P[yr]); x.append(v)
        Y, p, x = map(np.array, (Y, p, x))
        if detrend: p = detr(Y, p); x = detr(Y, x)
        p = (p - p.mean())/p.std(ddof=1)
        n = len(Y); r = stats.pearsonr(x, p)[0]; rho = stats.spearmanr(x, p)[0]
        r1, r2 = ar1(p), ar1(x); ne = min(n, max(3, n*(1-r1*r2)/(1+r1*r2)))
        pe = 2*stats.t.sf(abs(r)*np.sqrt((ne-2)/(1-r*r)), ne-2)
        f = epi_flags(IDX[nm], thr)
        cls = np.array([f.get((yy - lag, 1), 0) for yy in Y])
        en, ln = p[cls == 1], p[cls == -1]
        out[nm] = dict(n=n, r=round(r, 4), rho=round(rho, 4), neff=round(ne, 2), p=float('%.4g' % pe),
                       EN=round(en.mean(), 4), nEN=len(en), LN=round(ln.mean(), 4), nLN=len(ln),
                       p_welch=float('%.4g' % stats.ttest_ind(en, ln, equal_var=False).pvalue))
    fo, fr = epi_flags(IDX['ONI'], thr), epi_flags(IDX['RONI'], thr)
    out['reclass'] = [Y for Y in range(y0, y1+1) if (Y, 1) in IDX['ONI'] and (Y,1) in IDX['RONI'] and fo.get((Y,1),0) != fr.get((Y,1),0)]
    return out
cfgs = {
 'A': dict(core='HUACOL', var='d18O', source='cpc', season='djf', lag=0, detrend=False, y0=1960, y1=2019),
 'B': dict(core='QUE13', var='d18O', source='cpc', season='hydro', lag=1, detrend=True, y0=1952, y1=2009),
 'C': dict(core='HUA93', var='dust', source='ersst', season='ann', lag=0, detrend=False, y0=1895, y1=1993),
 'D': dict(core='HUA93', var='d18O', source='ersst', season='djf', lag=-1, detrend=True, y0=1895, y1=1993, thr=1.0),
}
res = {k: run(**v) for k, v in cfgs.items()}
print(json.dumps(res))
