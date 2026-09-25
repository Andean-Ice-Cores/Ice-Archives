"""Check the app's LMR ensemble correlation summary (Illimani NH4, 1800-1998, annual, lag 0,
no detrending): median member r for ONI and RONI must match the browser (0.4158 / 0.2553)."""
import json, numpy as np, pandas as pd
src = open('../app/data/dataset.js').read(); DS = json.loads(src[src.index('=')+1:].strip().rstrip(';'))
il = pd.read_csv('../data/raw/illimani2010nh4ann-noaa.txt', sep='\t', comment='#').dropna()
il = il[(il.age_CE >= 1800) & (il.age_CE <= 1998)]
Y = il.age_CE.values.astype(int); p = il['nh4_ueq/l'].values
for nm in ('oni', 'roni'):
    M = np.array(DS['lmr'][nm]) / 100.0
    rs = [np.corrcoef(M[k, Y], p)[0, 1] for k in range(M.shape[0])]
    print(nm.upper(), 'n =', len(Y), ' median r = %.4f' % np.median(rs), ' 2.5-97.5%% = [%.3f, %.3f]' % tuple(np.percentile(rs, [2.5, 97.5])))
