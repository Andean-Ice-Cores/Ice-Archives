# Andean ice cores × ENSO: ONI vs RONI — interactive figure

An interactive, publication-style figure testing whether NOAA's **Relative Oceanic Niño Index (RONI)**, which removes the
tropical-mean (20°S–20°N) SST anomaly from Niño-3.4, changes the apparent ENSO coupling of Andean ice-core records
relative to the classic **ONI**.

## Website

**https://andean-ice-cores.github.io/Ice-Archives/** — the 3-D explorer (Blender scene), the full analysis figure and
the data/methods notes in one site, plus a **Your data** tab: load your own annual, monthly or sub-annual series
(CSV/TSV/NOAA template, decimal years, ISO dates or years BP; processed only in the browser) and correlate it with
ONI/RONI (CPC 1950–, ERSST.v6 1850–, LMRv2.1 1–2000 CE with ensemble uncertainty), tropical-mean SST and every NCEI
Andean ice-core series — N_eff p-values, ΔR (MRR) test, lag scan, running correlation, composites and a screen of all
historical series with false-discovery-rate control. It is built from this repository into `docs/` (served by GitHub Pages):
`python3 site/build_site.py --pages docs`.

## Open it locally

Double-click **`app/index.html`** (2-D figure) or **`blender/web/index.html`** (3-D viewer), or `docs/index.html` (whole site). It runs entirely offline in any modern browser (no server, no internet; d3 is vendored).
Optionally: `python3 -m http.server 8765 --directory app` and open <http://localhost:8765>.

## What is in the figure

| Panel | Content |
|---|---|
| **a** | Map of the Andes (ETOPO1 relief) with every annually resolved NCEI core, coloured by *r* with the selected index; latitude–elevation inset with the Andean crest; click a core to open it. Split mode shows ONI and RONI maps side by side. |
| **b** | Proxy (z-score) overlaid on ONI and RONI; monthly El Niño/La Niña **episodes** under each definition (CPC rule) plus year shading. |
| **c** | Pearson r, Spearman ρ, 95% CI, naive and **N_eff-corrected** p (Bretherton et al. 1999), ΔR significance (Meng–Rosenthal–Rubin 1992); scatter. |
| **d** | El Niño / neutral / La Niña composites under each index (Welch + permutation tests). |
| **e** | *What changed?* ΔR per core, the ENSO years **reclassified** between ONI and RONI, and a plain-language reason for each (tropical-mean anomaly vs. ONI's base period, post-volcanic cooling, variance rescaling). |
| **f** | **Running correlation** (11–31-yr window) for ONI and RONI; pre-1950 via ERSST.v6 (1850–) or the **LMRv2.1** multi-proxy reconstruction (1–2000 CE) with the ensemble uncertainty band and an instrumental cross-check. |
| **g** | *Own feature 1 — hypothesis test:* regression of the proxy on the two components ONI mixes (ENSO-relative SST and tropical-mean SST), for every core. |
| **h** | *Own feature 2 — age-model uncertainty:* banded age model (Comboul et al. 2014) Monte Carlo of layer-counting errors accumulating from the core top. |
| text | Live interpretation of every parameter change, and a **draft figure caption** that updates with the settings. |

Controls: index (ONI / RONI / Split), ENSO source (CPC official, ERSST.v6-extended, LMRv2.1), core, proxy
(δ¹⁸O, accumulation, dust, 9 major ions, d-excess), ENSO season (DJF, Aug–Jul, Jan–Dec, custom window incl. wrap across
the year), lag −2…+2 yr (lag > 0 = ENSO leads the ice), detrending, analysis window, event threshold, running window,
miscount rate. Every panel exports to **SVG** and **PNG (3×)**; panels c and e also export **CSV**.

## Data (all downloaded from NOAA; nothing is placeholder)

| Record | Variables | Years | Source |
|---|---|---|---|
| Quelccaya Summit Dome 2003 core | δ¹⁸O, accumulation, dust, F⁻ Cl⁻ SO₄²⁻ NO₃⁻ Na⁺ NH₄⁺ K⁺ Mg²⁺ Ca²⁺ | 226–2009 | Thompson et al. 2013, NCEI `quelccaya2013-noaa.txt` |
| Quelccaya 2018 firn core | δ¹⁸O, 7 ions (sub-annual → Aug–Jul years) | 2004–2018 | NCEI `mayewski2023/` |
| Huascarán Col / Summit composites | δ¹⁸O, d-excess, accumulation, dust | 1960–2019 | Weber et al. 2023; dust Weber et al. 2026 |
| Huascarán 1993 core 2 | δ¹⁸O, particles, NO₃⁻ (thermal Aug–Jul years) | 1894–1993 | Thompson et al. 1995 |
| Illimani 1999 deep core | NH₄⁺ | 1800–1998 | Kellerhals et al. 2010 |
| Illimani 2017 firn core | δ¹⁸O, 7 ions (sub-annual → Aug–Jul years) | 2000–2017 | Lindau et al. 2021 |
| ONI, RONI | official monthly (3-month) series | 1950–2026 | NOAA CPC `oni.ascii.txt`, `RONI.ascii.txt` (ERSST.v6) |
| ERSST.v6 grid | ONI and RONI recomputed with CPC's algorithm | 1850–2026 | NOAA PSL |
| LMRv2.1 | Niño-3.4 ensemble + SST fields → reconstructed ONI/RONI | 1–2000 | Tardif et al. 2019, NCEI |

**Not available at annual resolution in the NCEI ice-core archive (searched Sep 2026):** Sajama (only 5-m / 100-yr
means), Coropuna, Hualcán, Chimborazo. They appear on the map as triangles but are not analysed. To add one if you
obtain it (e.g. from the Byrd Center), supply a tab-separated file `data/raw/<site>.txt` with columns
`age_CE<TAB><variable>…` (one row per year, `NA` for missing), then add a block in `pipeline/build_dataset.py`
modelled on the Illimani entry (id, site, lat, lon, elevation, reference, variables) and re-run the build.

## Verification

* **RONI definition reproduced:** our ERSST.v6 recomputation matches CPC RONI with r = 0.9994, RMSE = 0.03 °C
  (919 months); ONI with CPC's sliding base: r = 0.999, RMSE = 0.038 °C. This validates the pre-1950 extension.
* **Statistics library:** 24 quantities cross-checked against SciPy/NumPy (relative error < 10⁻⁶).
* **End-to-end:** an independent Python implementation from the raw files (`tests/verify_endtoend.py`) reproduces the
  app's n, r, ρ, N_eff, p, composite means/counts/p and the reclassified years **exactly** for four configurations
  (different cores, sources, seasons, lags, detrending, thresholds); LMR ensemble summaries also match
  (`tests/verify_lmr.py`).
* Smoke test: 260 control combinations rendered without errors; layout checked at desktop and phone widths; all
  panels serialize to valid SVG and non-blank PNG.

Run: `tests/run_all.sh`. Rebuild data: `pipeline/fetch_data.sh && python3 pipeline/build_dataset.py`.

## Findings from this dataset (so far)

* **1960–2019, CPC, DJF, lag 0:** Huascarán Col/Summit δ¹⁸O r ≈ 0.58–0.61, Huascarán 1993 core 0.53, Illimani firn
  0.66, Quelccaya firn 0.55 (all p < 0.05 after N_eff correction); Quelccaya 2003 core only 0.07. ΔR (RONI − ONI) is
  **+0.00 to +0.03 for every core and never significant**, because ONI and RONI correlate at ≈0.98 in this period.
  Seven ENSO years change class (e.g. 1992/93: neutral under ONI, El Niño under RONI because of post-Pinatubo cooling).
* **Decomposition (panel g):** for δ¹⁸O the ENSO-relative term dominates (Huascarán Col β = 0.65, p < 0.001;
  tropical-mean β = −0.08, n.s.): the interannual δ¹⁸O signal is dynamical, *not* tropical-mean — **the central
  hypothesis is not supported for δ¹⁸O at interannual scale**.
* **Illimani NH₄⁺** (a published temperature proxy) is the exception: with ERSST.v6 1851–1998 its tropical-mean
  coefficient is significant (annual β = 0.53, p < 10⁻⁶), yet ΔR remains small (−0.03 annual, +0.04 DJF, n.s.).
* **Before 1950:** Quelccaya δ¹⁸O coupling is strong ~1870–1920 (21-yr running r ≈ 0.5–0.7), collapses in the
  1920s–40s and is weak after ~1980. With LMRv2.1, ΔR is negative (−0.03 to −0.12), but (i) reconstructed-RONI members
  are noisier than ONI members, which attenuates per-member correlations, and (ii) Quelccaya is assimilated by LMR
  (PAGES2k), so treat LMR-based ΔR as an upper bound on the RONI penalty, not as independent evidence.

## Draft caption (150 words)

**Figure 1 | Andean ice-core δ¹⁸O–ENSO coupling is largely insensitive to removing tropical-mean warming.**
**a**, NOAA NCEI Andean ice cores with annual resolution (circles; fill, Pearson *r* with DJF ONI, 1960–2019; heavy
outline, *p* < 0.05 after effective-sample-size correction) and drill sites lacking annual archives (triangles); inset,
elevation versus latitude. **b**, Huascarán Col δ¹⁸O (black, z-score) with DJF ONI (purple) and RONI (green); shading,
El Niño (orange) and La Niña (blue) years. **c**, Correlation statistics. **d**, El Niño minus La Niña composites (95%
CI). **e**, ΔR = *r*(RONI) − *r*(ONI) per core (Meng–Rosenthal–Rubin test) and reclassified events. **f**, 21-yr
running correlation for Quelccaya δ¹⁸O using ERSST.v6-derived indices (dashed, 1851–2000) and 100 LMRv2.1
reconstruction members (shading, 5–95%). **g**, Standardized regression on ENSO-relative and tropical-mean SST.
**h**, Correlation spread under 1% yr⁻¹ layer-counting error. RONI subtracts the 20°S–20°N mean SST anomaly from
Niño-3.4, rescaled to Niño-3.4 variance.

## Layout

```
app/                 the interactive figure (index.html, css/, js/app.js, js/stats.js, data/dataset.js, vendor/d3)
pipeline/            fetch_data.sh, build_dataset.py, ersst_indices.py (RONI/ONI from ERSST.v6)
data/raw/            downloaded NOAA files (large/ holds ERSST.v6 and LMRv2.1 NetCDF, ≈2 GB)
tests/               stats + end-to-end + LMR verification (run_all.sh)
docs/                GitHub Pages site (built by site/build_site.py)
site/                website shell (src/) and build script
blender/             3-D scene: data prep, Blender build script, renders, web viewer
```

## 3-D infographic

A Blender-built 3-D version with an interactive web viewer lives in `blender/` — open `blender/web/index.html`; see `blender/README.md`.
