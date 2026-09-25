"""Build app/data/dataset.js from raw NOAA downloads.

Inputs (data/raw/):
  oni.ascii.txt, RONI.ascii.txt            NOAA CPC official indices (ERSST.v6)
  large/ersst.v6.sst.mnmean.nc             NOAA ERSST.v6 grid (PSL) -> extended ONI/RONI 1850-
  large/posterior_climate_indices_...nc    LMRv2.1 Nino3.4 ensemble (20 MC x 100 members)
  large/sst_MCruns_ensemble_{mean,spread}_LMRv2.1.nc  LMRv2.1 SST fields -> tropical mean
  NCEI ice-core text files (see CORES below)
  ne_50m_admin_0_countries.geojson, etopo_andes.csv   base map

Run:  python3 pipeline/build_dataset.py
"""
import os, json, base64, io, datetime as dt
import numpy as np, pandas as pd, h5py

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, 'data', 'raw')
OUT = os.path.join(ROOT, 'app', 'data', 'dataset.js')
SEAS = ['DJF', 'JFM', 'FMA', 'MAM', 'AMJ', 'MJJ', 'JJA', 'JAS', 'ASO', 'SON', 'OND', 'NDJ']
rng = np.random.default_rng(20260924)


def r2(x, nd=3):
    return [None if (v is None or not np.isfinite(v)) else round(float(v), nd) for v in x]


# ---------------------------------------------------------------- CPC indices
def read_cpc_season(fn, col):
    """Return dict {(year, centre_month): value}. DJF 1950 -> Jan 1950."""
    d = {}
    for line in open(os.path.join(RAW, fn)).read().splitlines()[1:]:
        p = line.split()
        if len(p) < 3:
            continue
        d[(int(p[1]), SEAS.index(p[0]) + 1)] = float(p[col])
    return d


oni_cpc = read_cpc_season('oni.ascii.txt', 3)
roni_cpc = read_cpc_season('RONI.ascii.txt', 2)


def to_monthly(d, y0):
    keys = sorted(d)
    y1, m1 = keys[-1]
    n = (y1 - y0) * 12 + m1
    arr = [None] * n
    for (y, m), v in d.items():
        i = (y - y0) * 12 + (m - 1)
        if 0 <= i < n:
            arr[i] = v
    return arr


cpc = dict(start=1950, oni=to_monthly(oni_cpc, 1950), roni=to_monthly(roni_cpc, 1950))
print('CPC months', len(cpc['oni']), 'last', max(oni_cpc))

# ---------------------------------------------------------------- ERSST.v6
from ersst_indices import load, box_mean, anom, run3  # noqa: E402

lat, lon, yr, mo, sst = load(os.path.join(RAW, 'large', 'ersst.v6.sst.mnmean.nc'))
n34 = box_mean(sst, lat, lon, -5, 5, 190, 240)
trop = box_mean(sst, lat, lon, -20, 20, 0, 360)
a34_fix = anom(n34, yr, mo)          # fixed 1991-2020 base
atr = anom(trop, yr, mo)
rel = a34_fix - atr
sel = (yr >= 1950)
k_roni = float(np.std(a34_fix[sel]) / np.std(rel[sel]))
roni_x = run3(rel) * k_roni


def oni_sliding(x, yr, mo):
    """CPC ONI convention: centred 30-yr base periods updated every 5 years
    (e.g. 1950-55 -> 1936-65, 1956-60 -> 1941-70). Base ends are clamped to the
    available record (>=1850 and <= last complete 30-yr period 1991-2020)."""
    out = np.full_like(x, np.nan)
    for y in np.unique(yr):
        b = 5 * ((y - 1) // 5) + 1   # 5-yr block start: 1951, 1956, ...
        if y == 1950:
            b = 1951
        s, e = b - 15, b + 14
        if e > 2020:
            s, e = 1991, 2020
        if s < int(yr.min()):
            s, e = int(yr.min()), int(yr.min()) + 29
        for m in range(1, 13):
            idx = (yr == y) & (mo == m)
            base = (yr >= s) & (yr <= e) & (mo == m)
            out[idx] = x[idx] - np.nanmean(x[base])
    return out


oni_x = run3(oni_sliding(n34, yr, mo))
t20_x = run3(atr)
nmon = len(yr)
# validation against CPC
def cmp(ours, cpcd):
    a, b = [], []
    for i in range(nmon):
        k = (int(yr[i]), int(mo[i]))
        if k in cpcd and np.isfinite(ours[i]):
            a.append(ours[i]); b.append(cpcd[k])
    a, b = np.array(a), np.array(b)
    return dict(r=round(float(np.corrcoef(a, b)[0, 1]), 4),
                rmse=round(float(np.sqrt(np.mean((a - b) ** 2))), 3), n=len(a))


val_oni = cmp(oni_x, oni_cpc)
val_roni = cmp(roni_x, roni_cpc)
print('ERSSTv6 ONI vs CPC', val_oni, ' RONI vs CPC', val_roni, 'k', k_roni)
ersst = dict(start=int(yr[0]), oni=r2(oni_x), roni=r2(roni_x), t20=r2(t20_x), n34rel=r2(run3(rel)),
             k=round(k_roni, 3), val_oni=val_oni, val_roni=val_roni)
# annual (Jan-Dec) Nino3.4 / RONI from ERSST for LMR validation
yrs_e = np.arange(int(yr[0]), int(yr[-1]) + 1)
ann = lambda x: np.array([np.nanmean(x[(yr == y)]) if (yr == y).sum() == 12 else np.nan for y in yrs_e])
e_n34_ann = ann(oni_sliding(n34, yr, mo)); e_rel_ann = ann(rel * k_roni); e_t20_ann = ann(atr)

# ---------------------------------------------------------------- LMRv2.1
lmr_dir = os.path.join(RAW, 'large')
g = h5py.File(os.path.join(lmr_dir, 'posterior_climate_indices_MCruns_ensemble_full_LMRv2.1.nc'), 'r')
nino = g['nino34'][:]                                  # (2001, 20, 100), years 0..2000
fm = h5py.File(os.path.join(lmr_dir, 'sst_MCruns_ensemble_mean_LMRv2.1.nc'), 'r')
fs = h5py.File(os.path.join(lmr_dir, 'sst_MCruns_ensemble_spread_LMRv2.1.nc'), 'r')
llat = fm['lat'][:]; tl = (llat >= -20) & (llat <= 20)
i0, i1 = np.where(tl)[0][[0, -1]]
w = np.cos(np.deg2rad(llat[i0:i1 + 1]))[:, None]


def trop_mean(ds):
    out = np.zeros((2001, 20))
    for t0 in range(0, 2001, 334):
        blk = ds['sst'][t0:t0 + 334, :, i0:i1 + 1, :].astype(float)
        blk[np.abs(blk) > 1e3] = np.nan
        ww = np.broadcast_to(w, blk.shape[2:]).copy()
        for m in range(20):
            for j in range(blk.shape[0]):
                f = blk[j, m]; ok = np.isfinite(f)
                out[t0 + j, m] = np.sum(f[ok] * ww[ok]) / np.sum(ww[ok])
    return out


def sliding_annual(x):
    """Apply the CPC ONI climatology convention to an annual series indexed by year 0..2000:
    each 5-yr block uses the centred 30-yr base (block start b -> b-15..b+14), clamped to the record."""
    out = np.empty_like(x); n = len(x)
    for y in range(n):
        b = 5 * ((y - 1) // 5) + 1
        s0, e0 = b - 15, b + 14
        if s0 < 0: s0, e0 = 0, 29
        if e0 > n - 1: s0, e0 = n - 30, n - 1
        out[y] = x[y] - x[s0:e0 + 1].mean()
    return out


print('LMR tropical mean ...')
T_mean = trop_mean(fm)          # tropical-mean SST anomaly per MC run (ensemble mean of that run)
T_sprd = trop_mean(fs)          # area-mean of grid-point spread (NOT the spread of the area mean)


def box_spread(la0, la1, lo0, lo1, t0=1800, t1=2001):
    il = np.where((llat >= la0) & (llat <= la1))[0]; ll = fs['lon'][:]; jl = np.where((ll >= lo0) & (ll <= lo1))[0]
    b = fs['sst'][t0:t1, :, il[0]:il[-1] + 1, jl[0]:jl[-1] + 1].astype(float); b[np.abs(b) > 1e3] = np.nan
    ww = np.broadcast_to(np.cos(np.deg2rad(llat[il]))[:, None], b.shape[2:]); ok = np.isfinite(b)
    return float((np.nansum(b * ww, axis=(2, 3)) / np.sum(ok * ww, axis=(2, 3))).mean()), float(np.sum(np.isfinite(b[0, 0]) * ww))


# Within-iteration ensemble spread of the tropical-mean is not archived. Estimate it from the
# Nino3.4 box, where both the member spread of the box mean and the grid-point spread are known:
# N_box = (grid spread / box-mean spread)^2 independent error regions; scale by ocean area.
s_member = float(nino[1800:2001].std(axis=2).mean())
sg_box, a_box = box_spread(-5, 5, 190, 240)
sg_tr, a_tr = box_spread(-20, 20, 0, 360)
n_box = (sg_box / s_member) ** 2
sigma_T = 2.0 * sg_tr / np.sqrt(n_box * a_tr / a_box)      # x2 for conservatism
print(f'sigma_T estimate: N_box={n_box:.2f}, area ratio={a_tr / a_box:.1f}, sigma_T(x2)={sigma_T:.3f}')
lmr_years = np.arange(0, 2001)
ND = 100
pick_m = rng.integers(0, 20, ND); pick_k = rng.integers(0, 100, ND)
oni_r, roni_r = [], []
ref = (lmr_years >= 1951) & (lmr_years <= 1980)
for d in range(ND):
    m, k = pick_m[d], pick_k[d]
    nn = nino[:, m, k].astype(float)
    tt = T_mean[:, m] + rng.standard_normal(2001) * sigma_T          # MC-iteration mean + within-iteration draw
    rr = nn - tt
    rr = rr * np.std(nn) / np.std(rr)                  # variance rescaling as in CPC RONI
    nn = sliding_annual(nn); rr = rr - rr[ref].mean()   # ONI: CPC-style centred 30-yr bases; RONI: fixed base
    oni_r.append(np.round(nn * 100).astype(int).tolist())
    roni_r.append(np.round(rr * 100).astype(int).tolist())
oni_arr = np.array(oni_r) / 100; roni_arr = np.array(roni_r) / 100
# skill vs ERSSTv6 annual means over 1880-2000
ov = (yrs_e >= 1880) & (yrs_e <= 2000)
lv = (lmr_years >= 1880) & (lmr_years <= 2000)
sk = dict(
    r_oni=round(float(np.corrcoef(np.median(oni_arr, 0)[lv], e_n34_ann[ov])[0, 1]), 3),
    r_roni=round(float(np.corrcoef(np.median(roni_arr, 0)[lv], e_rel_ann[ov])[0, 1]), 3),
    r_t20=round(float(np.corrcoef(T_mean.mean(1)[lv], e_t20_ann[ov])[0, 1]), 3),
    spread_oni=round(float(oni_arr.std(0)[lv].mean()), 3),
    spread_roni=round(float(roni_arr.std(0)[lv].mean()), 3),
    period='1880-2000', sigma_T=round(sigma_T, 3), n_box=round(n_box, 2), mc_spread_T=round(float(T_mean.std(1)[1800:].mean()), 3))
print('LMR skill', sk)
lmr = dict(start=0, members=ND, oni=oni_r, roni=roni_r, skill=sk,
           t20=r2(T_mean.mean(1) - T_mean.mean(1)[ref].mean()))

# ---------------------------------------------------------------- ice cores
def rd(fn):
    return pd.read_csv(os.path.join(RAW, fn), sep='\t', comment='#', na_values=['NA', '', 'nan'])


def series(df, ycol, vcol, y0=None):
    s = df[[ycol, vcol]].dropna()
    s = s[pd.to_numeric(s[vcol], errors='coerce').notna()]
    s = s.groupby(ycol)[vcol].mean().sort_index()
    if y0:
        s = s[s.index >= y0]
    return dict(years=[int(v) for v in s.index], values=r2(s.values.astype(float), 4))


def hydro_year(df, tcol, vcols, min_frac=0.6):
    """Aggregate decimal-year samples into Aug(Y-1)-Jul(Y) years (Y labels the
    wet season D(Y-1)JF(Y)); keep years with samples spanning >= min_frac of the year."""
    df = df.copy()
    df['hy'] = np.floor(df[tcol] - 0.5833 + 1).astype(int)
    out = {}
    for v in vcols:
        s = df.groupby('hy').agg(val=(v, 'mean'), tmin=(tcol, 'min'), tmax=(tcol, 'max'))
        s = s[(s.tmax - s.tmin) >= min_frac * 0.99]
        out[v] = dict(years=[int(y) for y in s.index], values=r2(s.val.values, 4))
    return out


V = lambda lab, units, d, **kw: dict(label=lab, units=units, **d, **kw)

cores = []
# 1 Quelccaya Summit Dome (Thompson et al. 2013)
q = rd('quelccaya2013-noaa.txt')
qi = {c: series(q, 'age_CE', c) for c in q.columns[1:]}
cores.append(dict(
    id='QUE13', site='Quelccaya', name='Quelccaya Summit Dome core (2003)', lat=-13.93, lon=-70.83, elev=5670,
    country='Peru', ref='Thompson et al. (2013) Science 340:945', doi='10.1126/science.1234210',
    url='https://www.ncei.noaa.gov/pub/data/paleo/icecore/trop/quelccaya/quelccaya2013-noaa.txt',
    dating='Annual layer counting (visible dust and δ18O seasonality); verified by volcanic horizons (e.g. Huaynaputina 1600 CE).',
    year_def='Ice year (wet season D(Y−1)JF(Y) forms most of layer Y)', in_lmr=True,
    vars=dict(
        d18O=V('δ18O', '‰ VSMOW', qi['d18OiceSMOW']),
        accum=V('Accumulation', 'm w.e. yr⁻¹', qi['iceacc_m']),
        dust=V('Dust (0.63–20 µm)', 'particles mL⁻¹', qi['dust']),
        ions={k: V(lab, 'ppb', qi[c]) for k, lab, c in [
            ('NO3', 'NO₃⁻', 'no3_ppb'), ('SO4', 'SO₄²⁻', 'so4_ppb'), ('Cl', 'Cl⁻', 'cl_ppb'),
            ('Na', 'Na⁺', 'na_ppb'), ('NH4', 'NH₄⁺', 'nh4_ppb'), ('K', 'K⁺', 'k_ppb'),
            ('Mg', 'Mg²⁺', 'mg_ppb'), ('Ca', 'Ca²⁺', 'ca_ppb'), ('F', 'F⁻', 'f_ppb')]})))

# 2 Quelccaya 2003-2018 firn core (dos Reis et al. 2022; NCEI 'mayewski2023')
qi2 = rd('mayewski2023-qu-18_isotopes.txt'); qc2 = rd('mayewski2023-qu-18_ic.txt')
qi2['t'] = (qi2.top_year + qi2.bottom_year) / 2; qc2['t'] = (qc2.top_year + qc2.bottom_year) / 2
h1 = hydro_year(qi2, 't', ['d18O']); h2 = hydro_year(qc2, 't', ['Na+', 'K+', 'Mg2+', 'Ca2+', 'Cl-', 'NO3-', 'SO42-'])
cores.append(dict(
    id='QUE18', site='Quelccaya', name='Quelccaya firn core (2018)', lat=-13.9333, lon=-70.8333, elev=5470,
    country='Peru', ref='dos Reis et al. (2022) / NCEI Mayewski 2023 archive', doi='',
    url='https://www.ncei.noaa.gov/pub/data/paleo/icecore/trop/quelccaya/mayewski2023/',
    dating='Seasonal δ18O / ion cycles; sub-annual samples averaged here into Aug–Jul years.',
    year_def='Aug(Y−1)–Jul(Y) mean of sub-annual samples', in_lmr=False,
    vars=dict(d18O=V('δ18O', '‰ VSMOW', h1['d18O']),
              ions={k: V(lab, 'µg L⁻¹', h2[c]) for k, lab, c in [
                  ('NO3', 'NO₃⁻', 'NO3-'), ('SO4', 'SO₄²⁻', 'SO42-'), ('Cl', 'Cl⁻', 'Cl-'), ('Na', 'Na⁺', 'Na+'),
                  ('K', 'K⁺', 'K+'), ('Mg', 'Mg²⁺', 'Mg2+'), ('Ca', 'Ca²⁺', 'Ca2+')]})))

# 3/4 Huascarán Col & Summit composites (Weber et al. 2023) + dust (Weber et al. 2026)
w23 = rd('weber2023-hs2019.txt')
wd = rd('weber2026-dust.txt')
col_ch = [c for c in wd.columns if c.startswith('colB_channel')]
sum_ch = [c for c in wd.columns if c.startswith('summitB_channel')]
wd['col_dust'] = wd[col_ch].sum(axis=1, min_count=len(col_ch))
wd['summit_dust'] = wd[sum_ch].sum(axis=1, min_count=len(sum_ch))
for pre, nm, el in [('col', 'Col', 6050), ('summit', 'Summit', 6768)]:
    cores.append(dict(
        id='HUA' + pre[:3].upper(), site='Huascarán', name=f'Huascarán {nm} composite (2019)', lat=-9.11, lon=-77.61,
        elev=el, country='Peru', ref='Weber et al. (2023); dust: Weber et al. (2026)', doi='',
        url='https://www.ncei.noaa.gov/pub/data/paleo/icecore/trop/huascaran/weber2023/weber2023-hs2019.txt',
        dating='Annual layer counting of δ18O seasonality, 1960–2019 CE.',
        year_def='Annual layer (ice year)', in_lmr=False,
        vars=dict(d18O=V('δ18O', '‰ VSMOW', series(w23, 'age_CE', f'{pre}_d18O')),
                  accum=V('Accumulation', 'm ice eq. yr⁻¹', series(w23, 'age_CE', f'{pre}_accum')),
                  dexcess=V('d-excess', '‰', series(w23, 'age_CE', f'{pre}_dexcess')),
                  dust=V('Dust (0.63–16 µm)', 'particles mL⁻¹', series(wd, 'age_CE', f'{pre}_dust')))))

# 5 Huascarán 1993 core 2 (Thompson et al. 1995), thermal years Aug-Jul
t95 = rd('thompson1995-annualt-noaa.txt')
cores.append(dict(
    id='HUA93', site='Huascarán', name='Huascarán Col core 2 (1993)', lat=-9.11, lon=-77.61, elev=6048,
    country='Peru', ref='Thompson et al. (1995) Science 269:46', doi='10.1126/science.269.5220.46',
    url='https://www.ncei.noaa.gov/pub/data/paleo/icecore/trop/huascaran/thompson1995-annualt-noaa.txt',
    dating='Annual layer counting, 1894–1993 CE; values are thermal-year (Aug–Jul) averages.',
    year_def='Thermal year Aug(Y−1)–Jul(Y)', in_lmr=False,
    vars=dict(d18O=V('δ18O', '‰ VSMOW', series(t95, 'age_CE', 'd18O_HScore2')),
              dust=V('Particles > 0.63 µm', 'particles mL⁻¹', series(t95, 'age_CE', 'part>0.63_HScore2')),
              ions={'NO3': V('NO₃⁻', 'ppb', series(t95, 'age_CE', 'NO3_HScore2'))})))

# 6 Illimani 1999 deep core NH4 (Kellerhals et al. 2010)
il = rd('illimani2010nh4ann-noaa.txt')
cores.append(dict(
    id='ILL99', site='Illimani', name='Illimani deep core (1999)', lat=-16.6167, lon=-67.7667, elev=6300,
    country='Bolivia', ref='Kellerhals et al. (2010) JGR 115:D16123', doi='10.1029/2009JD012603',
    url='https://www.ncei.noaa.gov/pub/data/paleo/icecore/trop/illimani/illimani2010nh4ann-noaa.txt',
    dating='Annual layer counting + 210Pb, tritium and volcanic markers, 1800–1998 CE.',
    year_def='Annual layer', in_lmr=False,
    vars=dict(ions={'NH4': V('NH₄⁺', 'µeq L⁻¹', series(il, 'age_CE', 'nh4_ueq/l'))})))

# 7 Illimani 2017 firn core (Lindau et al. 2021)
li = rd('lindau2021-il2017_isotopes.txt'); lc = rd('lindau2021-il2017_ic.txt')
g1 = hydro_year(li, 'year', ['d18O']); g2 = hydro_year(lc, 'year', ['Na+', 'K+', 'Mg2+', 'Ca2+', 'Cl-', 'NO3-', 'SO42-'])
cores.append(dict(
    id='ILL17', site='Illimani', name='Illimani firn core (2017)', lat=-16.6167, lon=-67.7667, elev=6350,
    country='Bolivia', ref='Lindau et al. (2021) Geophys. Res. Lett.', doi='',
    url='https://www.ncei.noaa.gov/pub/data/paleo/icecore/trop/illimani/lindau2021/',
    dating='Seasonal δ18O / ion cycles; sub-annual samples averaged into Aug–Jul years.',
    year_def='Aug(Y−1)–Jul(Y) mean of sub-annual samples', in_lmr=False,
    vars=dict(d18O=V('δ18O', '‰ VSMOW', g1['d18O']),
              ions={k: V(lab, 'ng g⁻¹', g2[c]) for k, lab, c in [
                  ('NO3', 'NO₃⁻', 'NO3-'), ('SO4', 'SO₄²⁻', 'SO42-'), ('Cl', 'Cl⁻', 'Cl-'), ('Na', 'Na⁺', 'Na+'),
                  ('K', 'K⁺', 'K+'), ('Mg', 'Mg²⁺', 'Mg2+'), ('Ca', 'Ca²⁺', 'Ca2+')]})))

for c in cores:
    for k, v in list(c['vars'].items()):
        if k == 'ions':
            c['vars'][k] = {i: s for i, s in v.items() if len(s['years']) >= 8}
            if not c['vars'][k]:
                del c['vars'][k]
        elif len(v['years']) < 8:
            del c['vars'][k]
    rng_ = [(v['years'][0], v['years'][-1]) for k, v in c['vars'].items() if k != 'ions'] + \
           [(s['years'][0], s['years'][-1]) for s in c['vars'].get('ions', {}).values()]
    print(c['id'], {k: (len(v['years']) if k != 'ions' else list(v)) for k, v in c['vars'].items()}, rng_[:2])

# Andean ice-core sites without annual NCEI data (shown, not analysed)
other_sites = [
    dict(site='Sajama', lat=-18.0, lon=-69.0, elev=6540, country='Bolivia',
         status='In NCEI (Thompson et al. 1998) but only 5-m / 100-yr resolution — no annual 1950–present overlap.'),
    dict(site='Coropuna', lat=-15.54, lon=-72.65, elev=6450, country='Peru',
         status='Cores drilled (2003) but no annual record in the NCEI ice-core archive (checked Sep 2026).'),
    dict(site='Hualcán', lat=-9.23, lon=-77.56, elev=6125, country='Peru',
         status='Cores drilled (2007) but no annual record in the NCEI ice-core archive (checked Sep 2026).'),
    dict(site='Chimborazo', lat=-1.47, lon=-78.82, elev=6268, country='Ecuador',
         status='Cores drilled (2000) but no annual record in the NCEI ice-core archive (checked Sep 2026).'),
]

# ---------------------------------------------------------------- base map
ext = dict(lon0=-84, lon1=-60, lat0=-24, lat1=4)
gj = json.load(open(os.path.join(RAW, 'ne_50m_admin_0_countries.geojson')))


def clip_feat(f):
    def ring_ok(r):
        xs = [p[0] for p in r]; ys = [p[1] for p in r]
        return max(xs) > ext['lon0'] - 5 and min(xs) < ext['lon1'] + 5 and max(ys) > ext['lat0'] - 5 and min(ys) < ext['lat1'] + 5
    geo = f['geometry']
    polys = geo['coordinates'] if geo['type'] == 'MultiPolygon' else [geo['coordinates']]
    keep = [[[[round(x, 2), round(y, 2)] for x, y in ring] for ring in poly] for poly in polys if ring_ok(poly[0])]
    if not keep:
        return None
    return dict(type='Feature', properties=dict(name=f['properties'].get('NAME')),
                geometry=dict(type='MultiPolygon', coordinates=keep))


countries = [c for c in (clip_feat(f) for f in gj['features']) if c]

et = pd.read_csv(os.path.join(RAW, 'etopo_andes.csv'), skiprows=[1])
grid = et.pivot(index='latitude', columns='longitude', values='altitude').sort_index(ascending=False)
Z = grid.values.astype(float)
import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import LightSource
import matplotlib.pyplot as plt
ls = LightSource(azdeg=315, altdeg=40)
hs = ls.hillshade(np.where(Z > 0, Z, 0), vert_exag=0.02, dx=11000, dy=11000)
land = Z > 0
rgb = np.ones(Z.shape + (4,))
elev_norm = np.clip(Z / 6000, 0, 1)
base = 0.97 - 0.35 * elev_norm           # higher = darker grey
shade = base * (0.72 + 0.28 * hs)
for c in range(3):
    rgb[..., c] = np.where(land, shade, 1.0)
rgb[..., 3] = np.where(land, 1.0, 0.0)   # ocean transparent (page colour shows through)
fig = plt.figure(figsize=(Z.shape[1] / 100, Z.shape[0] / 100), dpi=100)
ax = fig.add_axes([0, 0, 1, 1]); ax.axis('off'); ax.imshow(rgb, interpolation='bilinear', aspect='auto')
buf = io.BytesIO(); fig.savefig(buf, format='png', dpi=200, transparent=True); plt.close(fig)
relief_png = 'data:image/png;base64,' + base64.b64encode(buf.getvalue()).decode()
# 3000 m contour as polylines for the SVG
cs = plt.contour(grid.columns.values, grid.index.values, Z, levels=[3000])
contours = []
for seg in cs.allsegs[0]:
    if len(seg) > 12:
        contours.append([[round(float(x), 2), round(float(y), 2)] for x, y in seg[::2]])
plt.close('all')
# Andean crest: max elevation per 0.5° latitude band (for the latitude–elevation inset)
crest = []
for la in np.arange(ext['lat0'], ext['lat1'] + 0.01, 0.5):
    band = (grid.index.values >= la - 0.25) & (grid.index.values < la + 0.25)
    crest.append([round(float(la), 2), int(np.nanmax(Z[band]))])

meta = dict(
    built=dt.date.today().isoformat(),
    sources=dict(
        oni='NOAA CPC oni.ascii.txt (ERSST.v6, centred 30-yr base periods updated every 5 yr)',
        roni='NOAA CPC RONI.ascii.txt (ERSST.v6; Niño3.4 minus 20°S–20°N mean, variance-rescaled; 1991–2020 base)',
        ersst='NOAA ERSST.v6 monthly grid (PSL), indices recomputed here',
        lmr='LMRv2.1 (Tardif et al. 2019, Clim. Past 15:1251) posterior Niño3.4 ensemble + SST fields'),
    cpc_last=f"{SEAS[max(oni_cpc)[1]-1]} {max(oni_cpc)[0]}",
)
data = dict(meta=meta, cpc=cpc, ersst=ersst, lmr=lmr, cores=cores, other_sites=other_sites,
            map=dict(ext=ext, countries=countries, relief=relief_png, contour3000=contours, crest=crest))
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, 'w') as f:
    f.write('// Generated by pipeline/build_dataset.py on ' + meta['built'] + ' — do not edit by hand\n')
    f.write('window.DATA = ' + json.dumps(data, separators=(',', ':'), ensure_ascii=False) + ';\n')
print('wrote', OUT, round(os.path.getsize(OUT) / 1e6, 2), 'MB')
