"""Recompute Nino3.4, tropical-mean SST (20S-20N) and a RONI-style index from
gridded ERSST (v6 from NOAA PSL), to (1) verify the CPC RONI definition and
(2) extend ONI/RONI back to 1854 with the identical algorithm.

RONI (CPC, Feb 2026): 3-month running mean of [Nino3.4 anom - tropical-mean
(20N-20S) anom], rescaled so its variance equals that of Nino3.4; base 1991-2020.
"""
import sys, json, datetime as dt
import numpy as np, h5py

def load(fn):
    f = h5py.File(fn, 'r')
    lat = f['lat'][:].astype(float); lon = f['lon'][:].astype(float)
    t = f['time'][:]
    units = f['time'].attrs['units'].decode()
    base = dt.datetime(1800, 1, 1)
    assert units.startswith('days since 1800-1-1'), units
    dates = [base + dt.timedelta(days=float(x)) for x in t]
    sst = f['sst'][:].astype(float)
    mv = f['sst'].attrs.get('missing_value', [np.nan])[0]
    sst[(sst < -5) | (sst > 50) | (sst == mv)] = np.nan
    return lat, lon, np.array([d.year for d in dates]), np.array([d.month for d in dates]), sst

def box_mean(sst, lat, lon, la0, la1, lo0, lo1):
    mlat = (lat >= la0) & (lat <= la1)
    mlon = (lon >= lo0) & (lon <= lo1)
    sub = sst[:, mlat][:, :, mlon]
    w = np.cos(np.deg2rad(lat[mlat]))[None, :, None] * np.ones_like(sub)
    w[np.isnan(sub)] = 0
    return np.nansum(np.nan_to_num(sub) * w, axis=(1, 2)) / w.sum(axis=(1, 2))

def anom(x, yr, mo, y0=1991, y1=2020):
    out = x.copy()
    for m in range(1, 13):
        sel = mo == m
        clim = x[sel & (yr >= y0) & (yr <= y1)].mean()
        out[sel] = x[sel] - clim
    return out

def run3(x):
    y = np.full_like(x, np.nan)
    y[1:-1] = (x[:-2] + x[1:-1] + x[2:]) / 3
    return y

def compute(fn):
    lat, lon, yr, mo, sst = load(fn)
    n34 = box_mean(sst, lat, lon, -5, 5, 190, 240)      # 170W-120W
    trop = box_mean(sst, lat, lon, -20, 20, 0, 360)
    a34 = anom(n34, yr, mo); atr = anom(trop, yr, mo)
    rel = a34 - atr
    return yr, mo, a34, atr, rel

if __name__ == '__main__':
    fn = sys.argv[1]
    yr, mo, a34, atr, rel = compute(fn)
    # compare with CPC
    seasons = ['DJF','JFM','FMA','MAM','AMJ','MJJ','JJA','JAS','ASO','SON','OND','NDJ']
    def read_season(path, col):
        d = {}
        for l in open(path).read().split('\n')[1:]:
            p = l.split()
            if len(p) < 3: continue
            m = seasons.index(p[0]) + 1  # centre month
            d[(int(p[1]), m)] = float(p[col])
        return d
    roni = read_season('../data/raw/RONI.ascii.txt', 2)
    oni = read_season('../data/raw/oni.ascii.txt', 3)
    key = [(int(a), int(b)) for a, b in zip(yr, mo)]
    r3 = run3(rel); n3 = run3(a34)
    idx = [i for i, k in enumerate(key) if k in roni and k in oni and not np.isnan(r3[i])]
    cpc_r = np.array([roni[key[i]] for i in idx]); cpc_o = np.array([oni[key[i]] for i in idx])
    x = r3[idx]
    for lab, per in [('1950-now', (1950, 2100)), ('1991-2020', (1991, 2020))]:
        m = (yr >= per[0]) & (yr <= per[1]) & ~np.isnan(a34)
        k = np.std(a34[m]) / np.std(rel[m])
        rr = x * k
        print(f'scale over {lab}: k={k:.3f}  corr(myRONI,CPC)={np.corrcoef(rr, cpc_r)[0,1]:.4f}  RMSE={np.sqrt(np.mean((rr-cpc_r)**2)):.3f}  bias={np.mean(rr-cpc_r):+.3f}')
    print(f'least-squares implied k = {np.sum(x*cpc_r)/np.sum(x*x):.3f}')
    print('corr(my 3-mo N34 anom [fixed 1991-2020 base], CPC ONI) =', round(np.corrcoef(n3[idx], cpc_o)[0,1],4),
          ' RMSE', round(float(np.sqrt(np.mean((n3[idx]-cpc_o)**2))),3))
    print('corr(CPC ONI, CPC RONI) =', round(np.corrcoef(cpc_o, cpc_r)[0,1],4))
    print('first year with data:', yr[0], 'n months', len(yr))
