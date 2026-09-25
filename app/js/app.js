/* Andean ice cores × ENSO (ONI vs RONI) — interactive figure.
 * State → compute (Stats) → render (d3/SVG). Every panel is a self-contained SVG so it can be
 * exported to SVG/PNG exactly as displayed. */
(function () {
  'use strict';
  const D = window.DATA, S = window.Stats;
  const MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
  const IDX = ['ONI', 'RONI'];
  // Okabe–Ito colour-blind-safe palette
  const C = { proxy: '#1b1b1b', ONI: '#CC79A7', RONI: '#009E73', EN: '#D55E00', LN: '#0072B2', N: '#9a9a9a', grid: '#ececec', T20: '#E69F00' };
  const rColor = d3.scaleDiverging(t => d3.interpolateRdBu(1 - t)).domain([-0.6, 0, 0.6]).clamp(true);
  const FIG_CSS = `text{font-family:"Helvetica Neue",Helvetica,Arial,sans-serif}text:not([fill]){fill:#1b1b1b}
    .ax text{font-size:11px}.ax path,.ax line{stroke:#555;shape-rendering:crispEdges}
    .gridl line{stroke:#ececec;shape-rendering:crispEdges}.gridl path{display:none}
    .lab{font-size:11.5px}.ttl{font-size:12px;font-weight:700}.sm{font-size:10px}.sm:not([fill]){fill:#555}.it{font-style:italic}
    .zero{stroke:#999;stroke-dasharray:2 2}`;
  const SRC_LABEL = { cpc: 'NOAA CPC', ersst: 'ERSST.v6 (recomputed)', lmr: 'LMRv2.1 reconstruction' };
  const VAR_LABEL = { d18O: 'δ¹⁸O', accum: 'accumulation', dust: 'dust', ions: 'major ions', dexcess: 'd-excess' };

  // ------------------------------------------------------------------ state
  const state = {
    mode: 'ONI', source: 'cpc', core: 'HUACOL', variable: 'd18O', ion: 'NO3', season: 'djf', cm0: 12, cm1: 3,
    lag: 0, detrend: false, y0: 1960, y1: 2019, thr: 0.5, runW: 21, theta: 0.01
  };
  const cores = D.cores;
  const SHORT = { QUE13: 'Q03', QUE18: 'Q18', HUACOL: 'HC', HUASUM: 'HS', HUA93: 'H93', ILL99: 'I99', ILL17: 'I17' };
  cores.forEach(c => c.short = SHORT[c.id] || c.id.slice(-3));
  const coreById = Object.fromEntries(cores.map(c => [c.id, c]));

  // ------------------------------------------------------------------ index access
  const lmrMembers = { ONI: D.lmr.oni.map(a => Float64Array.from(a, v => v / 100)), RONI: D.lmr.roni.map(a => Float64Array.from(a, v => v / 100)) };
  const lmrSigma = { ONI: lmrMembers.ONI.map(a => S.sd(Array.from(a))), RONI: lmrMembers.RONI.map(a => S.sd(Array.from(a))) };
  const LMR_END = D.lmr.start + D.lmr.oni[0].length - 1;

  function monthlySeries(source, name) {        // name: oni | roni | t20 | n34rel
    if (source === 'cpc' && (name === 'oni' || name === 'roni')) return { start: D.cpc.start, a: D.cpc[name] };
    return { start: D.ersst.start, a: D.ersst[name] };
  }
  function monthsFor(Y, st) {
    const out = [];
    if (st.season === 'djf') return [[Y, 1]];
    if (st.season === 'ann') { for (let m = 1; m <= 12; m++) out.push([Y, m]); return out; }
    let m0 = st.season === 'hydro' ? 8 : st.cm0, m1 = st.season === 'hydro' ? 7 : st.cm1;
    if (m0 <= m1) { for (let m = m0; m <= m1; m++) out.push([Y, m]); }
    else { for (let m = m0; m <= 12; m++) out.push([Y - 1, m]); for (let m = 1; m <= m1; m++) out.push([Y, m]); }
    return out;
  }
  function seasonal(ser, Y, st) {
    let s = 0, n = 0;
    for (const [y, m] of monthsFor(Y, st)) {
      const v = ser.a[(y - ser.start) * 12 + m - 1];
      if (v === null || v === undefined) return null; s += v; n++;
    }
    return s / n;
  }
  const seasonName = st => st.season === 'djf' ? 'DJF' : st.season === 'ann' ? 'Jan–Dec' : st.season === 'hydro' ? 'Aug–Jul' :
    `${MON[st.cm0 - 1]}${st.cm0 > st.cm1 ? '(−1)' : ''}–${MON[st.cm1 - 1]}`;

  // annual index value (deterministic) or array of member values (LMR)
  function indexAt(name, Y, st) {
    if (st.source === 'lmr') {
      if (Y < D.lmr.start || Y > LMR_END) return null;
      return lmrMembers[name].map(a => a[Y - D.lmr.start]);
    }
    return seasonal(monthlySeries(st.source, name.toLowerCase()), Y, st);
  }
  function t20At(Y, st) {
    if (st.source === 'lmr') { const v = D.lmr.t20[Y - D.lmr.start]; return v === undefined ? null : v; }
    return seasonal(monthlySeries('ersst', 't20'), Y, st);
  }
  function relAt(Y, st) {        // un-rescaled relative Niño3.4 (for the decomposition)
    if (st.source === 'lmr') { const m = indexAt('RONI', Y, st); return m ? d3.median(m) : null; }
    return seasonal(monthlySeries('ersst', 'n34rel'), Y, st);
  }

  // tropical-mean SST anomaly relative to the same centred 30-yr base as ONI (for explaining ONI−RONI gaps)
  const VOLC = [{ n: 'Tambora', y: 1815 }, { n: 'Krakatau', y: 1883 }, { n: 'Santa María', y: 1902 }, { n: 'Agung', y: 1963 }, { n: 'El Chichón', y: 1982 }, { n: 'Pinatubo', y: 1991 }];
  const t20ClimCache = new Map();
  function t20Sliding(Y, m) {
    const ser = monthlySeries('ersst', 't20'), y0 = ser.start;
    let b = 5 * Math.floor((Y - 1) / 5) + 1; if (Y === 1950) b = 1951;
    let s0 = b - 15, e0 = b + 14; if (e0 > 2020) { s0 = 1991; e0 = 2020; } if (s0 < y0) { s0 = y0; e0 = y0 + 29; }
    const key = s0 + ':' + m;
    if (!t20ClimCache.has(key)) { const v = []; for (let y = s0; y <= e0; y++) { const x = ser.a[(y - y0) * 12 + m - 1]; if (x !== null) v.push(x); } t20ClimCache.set(key, S.mean(v)); }
    const x = ser.a[(Y - y0) * 12 + m - 1]; return x === null || x === undefined ? NaN : x - t20ClimCache.get(key);
  }
  // CPC episode rule: ≥5 consecutive overlapping seasons beyond ±thr
  const epiCache = new Map();
  function episodeFlags(source, name, thr) {
    const key = source + name + thr; if (epiCache.has(key)) return epiCache.get(key);
    const ser = monthlySeries(source, name.toLowerCase()), a = ser.a, n = a.length, f = new Int8Array(n);
    for (const sg of [1, -1]) {
      let i = 0;
      while (i < n) {
        if (a[i] !== null && sg * a[i] >= thr - 1e-9) {
          let j = i; while (j + 1 < n && a[j + 1] !== null && sg * a[j + 1] >= thr - 1e-9) j++;
          if (j - i + 1 >= 5) for (let k = i; k <= j; k++) f[k] = sg;
          i = j + 1;
        } else i++;
      }
    }
    const out = { start: ser.start, f }; epiCache.set(key, out); return out;
  }
  function episodes(source, name, thr, y0, y1) {
    const E = episodeFlags(source, name, thr), ser = monthlySeries(source, name.toLowerCase()), out = [];
    let i = 0; const n = E.f.length;
    while (i < n) {
      if (E.f[i] !== 0) {
        let j = i; while (j + 1 < n && E.f[j + 1] === E.f[i]) j++;
        const ys = E.start + Math.floor(i / 12), ye = E.start + Math.floor(j / 12);
        let pk = 0; for (let k = i; k <= j; k++) if (Math.abs(ser.a[k]) > Math.abs(pk)) pk = ser.a[k];
        if (ye >= y0 && ys <= y1) out.push({ sign: E.f[i], i0: i, i1: j, t0: E.start + i / 12, t1: E.start + (j + 1) / 12, peak: pk, ys, ye, m0: i % 12, m1: j % 12 });
        i = j + 1;
      } else i++;
    }
    return out;
  }
  // ENSO class of "ENSO year" Y: episode status of the DJF season centred on Jan Y
  // (LMR: annual value beyond ±thr·σ, evaluated per ensemble member)
  function yearClass(name, Y, st) {
    if (st.source === 'lmr') {
      const m = indexAt(name, Y, st); if (!m) return null;
      return m.map((v, k) => v >= st.thr * lmrSigma[name][k] ? 1 : v <= -st.thr * lmrSigma[name][k] ? -1 : 0);
    }
    const E = episodeFlags(st.source, name, st.thr), i = (Y - E.start) * 12;
    if (i < 0 || i >= E.f.length) return null;
    const ser = monthlySeries(st.source, name.toLowerCase()); if (ser.a[i] === null) return null;
    return E.f[i];
  }

  // ------------------------------------------------------------------ proxy access & alignment
  function proxyOf(core, st) {
    const v = core.vars[st.variable]; if (!v) return null;
    if (st.variable === 'ions') { const s = v[st.ion]; return s ? { ...s, key: 'ions:' + st.ion } : null; }
    return { ...v, key: st.variable };
  }
  function indexRange(st) {
    if (st.source === 'lmr') return [D.lmr.start + 1, LMR_END];
    if (st.source === 'cpc') return [D.cpc.start + 1, D.cpc.start + Math.floor((D.cpc.oni.length - 1) / 12)];
    return [D.ersst.start + 1, D.ersst.start + Math.floor((D.ersst.oni.length - 2) / 12)];
  }

  // Align proxy(Y) with index(Y − lag), within [y0, y1]. Returns z-scored, optionally detrended arrays.
  function align(core, name, st, opts = {}) {
    const P = proxyOf(core, st); if (!P) return null;
    const y0 = opts.y0 ?? st.y0, y1 = opts.y1 ?? st.y1, lag = opts.lag ?? st.lag;
    const Y = [], p = [], x = [];
    let pv = P.values;
    if (opts.proxyValues) pv = opts.proxyValues;
    for (let i = 0; i < P.years.length; i++) {
      const yr = P.years[i]; if (yr < y0 || yr > y1 || pv[i] === null) continue;
      const iv = indexAt(name, yr - lag, st); if (iv === null) continue;
      Y.push(yr); p.push(pv[i]); x.push(iv);
    }
    if (Y.length < 8) return { n: Y.length, Y, P };
    const ens = st.source === 'lmr';
    let pp = p.slice(), xx;
    if (st.detrend) pp = S.detrend(Y, pp).y;
    pp = S.zscore(pp);
    if (ens) {
      const M = x[0].length; xx = [];
      for (let k = 0; k < M; k++) { let col = x.map(v => v[k]); if (st.detrend) col = S.detrend(Y, col).y; xx.push(col); }
    } else { xx = st.detrend ? S.detrend(Y, x).y : x; }
    return { n: Y.length, Y, pRaw: p, p: pp, x: xx, ens, P };
  }

  function summarizeEns(list, key) { return { med: d3.median(list, d => d[key]), lo: S.quantile(list.map(d => d[key]), 0.025), hi: S.quantile(list.map(d => d[key]), 0.975) }; }

  function correlate(A) {
    if (!A || A.n < 8) return null;
    if (!A.ens) return S.corrStats(A.x, A.p);
    const per = A.x.map(col => S.corrStats(col, A.p));
    // total uncertainty: reconstruction ensemble ⊕ sampling (Fisher z with N_eff)
    const rnd = S.mulberry32(11), pooled = [];
    per.forEach(s => { for (let j = 0; j < 20; j++) { const u1 = rnd(), u2 = rnd(); const g = Math.sqrt(-2 * Math.log(u1 + 1e-12)) * Math.cos(2 * Math.PI * u2); pooled.push(Math.tanh(Math.atanh(s.r) + g / Math.sqrt(Math.max(s.neff - 3, 1)))); } });
    const r = summarizeEns(per, 'r');
    return {
      n: A.n, r: r.med, r_lo: r.lo, r_hi: r.hi, rho: d3.median(per, d => d.rho), neff: d3.median(per, d => d.neff),
      ar1x: d3.median(per, d => d.ar1x), ar1y: per[0].ar1y, p: d3.median(per, d => d.p), p_rho: d3.median(per, d => d.p_rho),
      p_naive: d3.median(per, d => d.p_naive), ci: [S.quantile(pooled, 0.025), S.quantile(pooled, 0.975)],
      fracSig: per.filter(d => d.p < 0.05).length / per.length, members: per.map(d => d.r), ens: true
    };
  }

  function analyze(core, st) {
    const out = {};
    for (const nm of IDX) { const A = align(core, nm, st); out[nm] = { A, s: correlate(A) }; }
    // ΔR significance (dependent correlations share the proxy)
    const a = out.ONI.A, b = out.RONI.A;
    if (out.ONI.s && out.RONI.s && a.n === b.n) {
      const xo = a.ens ? a.x.map((_, i) => 0) : a.x;
      let r12;
      if (a.ens) { const mo = a.Y.map((_, i) => d3.median(a.x, c => c[i])), mr = b.Y.map((_, i) => d3.median(b.x, c => c[i])); r12 = S.pearson(mo, mr); }
      else r12 = S.pearson(a.x, b.x);
      const ne = Math.min(out.ONI.s.neff, out.RONI.s.neff);
      out.diff = { dr: out.RONI.s.r - out.ONI.s.r, r12, ...S.compareDepCorr(out.RONI.s.r, out.ONI.s.r, r12, ne), neff: ne };
      void xo;
    }
    return out;
  }

  function composites(A, name, st) {
    if (!A || A.n < 8) return null;
    const lag = st.lag;
    if (!A.ens) {
      const g = { 1: [], 0: [], [-1]: [] };
      A.Y.forEach((y, i) => { const c = yearClass(name, y - lag, st); if (c !== null) g[c].push(A.p[i]); });
      const w = S.welch(g[1], g[-1]);
      return {
        EN: S.groupSummary(g[1]), N: S.groupSummary(g[0]), LN: S.groupSummary(g[-1]),
        diff: g[1].length && g[-1].length ? S.mean(g[1]) - S.mean(g[-1]) : NaN, p_welch: w.p, p_perm: S.permTest(g[1], g[-1], 4000)
      };
    }
    // ensemble: member-wise classification
    const M = A.x.length, res = { EN: [], N: [], LN: [], diff: [], p: [], nEN: [], nLN: [] };
    const cls = A.Y.map(y => yearClass(name, y - lag, st));
    for (let k = 0; k < M; k++) {
      const g = { 1: [], 0: [], [-1]: [] };
      A.Y.forEach((y, i) => { if (cls[i]) g[cls[i][k]].push(A.p[i]); });
      if (g[1].length < 2 || g[-1].length < 2) continue;
      res.EN.push(S.groupSummary(g[1])); res.N.push(S.groupSummary(g[0])); res.LN.push(S.groupSummary(g[-1]));
      res.diff.push(S.mean(g[1]) - S.mean(g[-1])); res.p.push(S.welch(g[1], g[-1]).p); res.nEN.push(g[1].length); res.nLN.push(g[-1].length);
    }
    if (!res.diff.length) return null;
    const pool = arr => {           // member spread ⊕ sampling SE
      const m = arr.map(d => d.mean), se = d3.median(arr, d => d.se) || 0, mm = d3.median(m);
      const sdm = m.length > 1 ? S.sd(m) : 0, tot = Math.sqrt(sdm * sdm + se * se);
      return { n: Math.round(d3.median(arr, d => d.n)), mean: mm, se: tot, ci: [mm - 1.96 * tot, mm + 1.96 * tot] };
    };
    return { EN: pool(res.EN), N: pool(res.N), LN: pool(res.LN), diff: d3.median(res.diff), p_welch: d3.median(res.p), p_perm: NaN, ens: true, fracSig: res.p.filter(p => p < 0.05).length / res.p.length };
  }

  // ------------------------------------------------------------------ helpers
  const fmt = (v, d = 2) => (v === null || v === undefined || !Number.isFinite(v)) ? '—' : (v < 0 ? '−' : '') + Math.abs(v).toFixed(d);
  const fmtP = p => !Number.isFinite(p) ? '—' : p < 0.001 ? '<0.001' : p.toFixed(3);
  const eqP = p => !Number.isFinite(p) ? '= —' : p < 0.001 ? '< 0.001' : '= ' + p.toFixed(3);
  const sigStar = p => !Number.isFinite(p) ? '' : p < 0.01 ? '**' : p < 0.05 ? '*' : '';
  const clsName = c => c === 1 ? 'El Niño' : c === -1 ? 'La Niña' : 'neutral';
  const clsKey = c => c === 1 ? 'EN' : c === -1 ? 'LN' : 'N';
  const tip = d3.select('#tooltip');
  function showTip(ev, html) { tip.html(html).style('opacity', 1).style('left', Math.min(ev.clientX + 14, innerWidth - 290) + 'px').style('top', (ev.clientY + 12) + 'px'); }
  const hideTip = () => tip.style('opacity', 0);
  function newSvg(sel, W, H) {
    const el = d3.select(sel); el.selectAll('*').remove();
    const svg = el.append('svg').attr('viewBox', `0 0 ${W} ${H}`).attr('width', W).attr('height', H)
      .attr('xmlns', 'http://www.w3.org/2000/svg').attr('font-family', 'Helvetica Neue, Helvetica, Arial, sans-serif');
    svg.append('style').text(FIG_CSS);
    svg.append('rect').attr('width', W).attr('height', H).attr('fill', '#fff');
    return svg;
  }
  const widthOf = sel => Math.max(320, document.querySelector(sel).clientWidth);
  function proxyLabel(core, st) { const P = proxyOf(core, st); return P ? `${P.label} (${P.units})` : '—'; }
  function varTitle(st) { return st.variable === 'ions' ? (ionLabel(st.ion) || 'ion') : VAR_LABEL[st.variable]; }
  function ionLabel(k) { for (const c of cores) if (c.vars.ions && c.vars.ions[k]) return c.vars.ions[k].label; return k; }

  // ------------------------------------------------------------------ controls
  function initControls() {
    const cs = d3.select('#ctl-core');
    const bySite = d3.group(cores, c => c.site);
    for (const [site, list] of bySite) {
      const og = cs.append('optgroup').attr('label', site);
      list.forEach(c => og.append('option').attr('value', c.id).text(c.name));
    }
    ['ctl-cm0', 'ctl-cm1'].forEach(id => MON.forEach((m, i) => d3.select('#' + id).append('option').attr('value', i + 1).text(m)));
    const ions = new Set(); cores.forEach(c => c.vars.ions && Object.keys(c.vars.ions).forEach(k => ions.add(k)));
    [...ions].forEach(k => d3.select('#ctl-ion').append('option').attr('value', k).text(ionLabel(k)));

    d3.selectAll('#ctl-mode button').on('click', function () { state.mode = this.dataset.v; update(); });
    d3.selectAll('#ctl-detrend button').on('click', function () { state.detrend = this.dataset.v === '1'; update(); });
    d3.select('#ctl-source').on('change', function () {
      state.source = this.value;
      if (state.source === 'lmr') { state.season = 'ann'; }
      fitWindow(true); update();
    });
    d3.select('#ctl-core').on('change', function () { selectCore(this.value); });
    d3.select('#ctl-var').on('change', function () { state.variable = this.value; ensureVariable(); fitWindow(true); update(); });
    d3.select('#ctl-ion').on('change', function () { state.ion = this.value; ensureVariable(); fitWindow(true); update(); });
    d3.select('#ctl-season').on('change', function () { state.season = this.value; update(); });
    d3.select('#ctl-cm0').on('change', function () { state.cm0 = +this.value; update(); });
    d3.select('#ctl-cm1').on('change', function () { state.cm1 = +this.value; update(); });
    d3.select('#ctl-lag').on('input', function () { state.lag = +this.value; update(); });
    d3.select('#ctl-thr').on('input', function () { state.thr = +this.value; update(); });
    d3.select('#ctl-runw').on('input', function () { state.runW = +this.value; update(); });
    d3.select('#ctl-theta').on('input', function () { state.theta = +this.value; update(); });
    d3.select('#ctl-y0').on('change', function () { state.y0 = +this.value; clampWindow(); update(); });
    d3.select('#ctl-y1').on('change', function () { state.y1 = +this.value; clampWindow(); update(); });
    d3.select('#btn-fullwin').on('click', () => { fitWindow(true); update(); });
    d3.select('#btn-copycap').on('click', () => navigator.clipboard && navigator.clipboard.writeText(document.getElementById('caption').innerText));
    // export buttons
    d3.selectAll('.ex').each(function () {
      const t = this.dataset.target, box = d3.select(this);
      box.append('button').text('SVG').attr('title', 'Download panel as SVG').on('click', () => exportPanel(t, 'svg'));
      box.append('button').text('PNG').attr('title', 'Download panel as PNG (3×)').on('click', () => exportPanel(t, 'png'));
      if (t === 'fig-c' || t === 'fig-e') box.append('button').text('CSV').attr('title', 'Download underlying numbers').on('click', () => exportCSV(t));
    });
    let rt; window.addEventListener('resize', () => { clearTimeout(rt); rt = setTimeout(update, 150); });
  }
  function selectCore(id) {
    state.core = id; const c = coreById[id];
    if (!proxyOf(c, state)) { state.variable = Object.keys(c.vars)[0]; }
    ensureVariable(); fitWindow(true); update();
  }
  function ensureVariable() {
    const c = coreById[state.core];
    if (state.variable === 'ions' && c.vars.ions && !c.vars.ions[state.ion]) state.ion = Object.keys(c.vars.ions)[0];
  }
  function overlap(core, st) {
    const P = proxyOf(core, st); if (!P) return null;
    const [a, b] = indexRange(st);
    const lo = Math.max(P.years[0], a + 2), hi = Math.min(P.years[P.years.length - 1], b - 2);
    return hi - lo >= 7 ? [lo, hi] : null;
  }
  function fitWindow(force) {
    const ov = overlap(coreById[state.core], state);
    if (!ov) return;
    if (force || state.y0 < ov[0] || state.y1 > ov[1]) { state.y0 = ov[0]; state.y1 = ov[1]; }
    // the LMR ensemble is uninformative about interannual timing far back; default to the last 400 yr for long cores
    if (force && state.source === 'lmr' && state.y1 - state.y0 > 400) state.y0 = state.y1 - 399;
  }
  function clampWindow() { if (state.y1 < state.y0 + 9) state.y1 = state.y0 + 9; }
  function syncControls() {
    d3.selectAll('#ctl-mode button').classed('on', function () { return this.dataset.v === state.mode; });
    d3.selectAll('#ctl-detrend button').classed('on', function () { return (this.dataset.v === '1') === state.detrend; });
    d3.select('#ctl-source').property('value', state.source);
    d3.select('#ctl-core').property('value', state.core);
    const c = coreById[state.core];
    d3.select('#ctl-var').selectAll('option').property('disabled', function () { return !c.vars[this.value]; });
    d3.select('#ctl-var').property('value', state.variable);
    d3.select('#ctl-ion').style('display', state.variable === 'ions' ? null : 'none')
      .selectAll('option').property('disabled', function () { return !(c.vars.ions && c.vars.ions[this.value]); });
    d3.select('#ctl-ion').property('value', state.ion);
    const sea = d3.select('#ctl-season'); sea.property('value', state.season);
    sea.selectAll('option').property('disabled', function () { return state.source === 'lmr' && this.value !== 'ann'; });
    d3.select('#custom-win').style('display', state.season === 'custom' ? null : 'none');
    d3.select('#ctl-cm0').property('value', state.cm0); d3.select('#ctl-cm1').property('value', state.cm1);
    d3.select('#lag-val').text((state.lag > 0 ? '+' : state.lag < 0 ? '−' : '') + Math.abs(state.lag) + ' yr');
    d3.select('#ctl-lag').property('value', state.lag);
    d3.select('#thr-val').text('±' + state.thr.toFixed(2) + (state.source === 'lmr' ? ' σ' : ' °C'));
    d3.select('#runw-val').text(state.runW + ' yr');
    d3.select('#theta-val').text((state.theta * 100).toFixed(1) + '% yr⁻¹');
    d3.select('#ctl-y0').property('value', state.y0); d3.select('#ctl-y1').property('value', state.y1);
  }

  // ------------------------------------------------------------------ (a) map
  function renderMap(results) {
    const W = widthOf('#fig-a'), split = state.mode === 'SPLIT', narrow = W < 1000;
    const blockW = split && !narrow ? W / 2 : narrow ? W : W * 0.58;
    const mapWest = blockW - Math.min(150, Math.max(100, blockW * 0.2)) - 60;       // width available to the map itself
    const H = Math.min(640, Math.max(narrow ? 300 : 480, mapWest * 28 / 24 + 84));  // 28° lat × 24° lon extent
    const listH = 110 + cores.length * 34;
    const totalH = split ? (narrow ? 2 * H : H) : (narrow ? H + listH : H);
    const svg = newSvg('#fig-a', W, totalH);
    const blocks = split ? (narrow ? [['ONI', 0, 0, W], ['RONI', 0, H, W]] : [['ONI', 0, 0, W / 2], ['RONI', W / 2, 0, W / 2]])
      : [[state.mode, 0, 0, narrow ? W : W * 0.58]];
    blocks.forEach(([nm, x0, y0, w]) => mapBlock(svg.append('g').attr('transform', `translate(0,${y0})`), nm, x0, w, H, results));
    if (!split) narrow ? siteList(svg.append('g').attr('transform', `translate(0,${H})`), 8, W - 12, listH, results) : siteList(svg, W * 0.6, W * 0.4 - 4, H, results);
    d3.select('#note-a').html(`Circles: annual-resolution NCEI cores (fill = Pearson <i>r</i> between ${varTitle(state)} and the ${seasonName(state)} index at lag ${state.lag} yr over ${state.y0}–${state.y1}; heavy outline = p&lt;0.05 after N<sub>eff</sub> correction; open = variable not measured or &lt;8 overlapping years). Triangles: Andean drill sites without annual-resolution NCEI data. Click a core to open its record. Grey relief: ETOPO1; line: 3000 m contour.`);
  }
  function mapBlock(svg, nm, x0, w, H, results) {
    const ext = D.map.ext, g = svg.append('g').attr('transform', `translate(${x0},0)`);
    const insetW = Math.min(150, Math.max(100, w * 0.2)), mapW = w - insetW - 26, top = 26, bottom = 58;
    const proj = d3.geoEquirectangular().fitExtent([[34, top], [mapW, H - bottom]],
      { type: 'MultiPoint', coordinates: [[ext.lon0, ext.lat0], [ext.lon1, ext.lat1]] });
    const path = d3.geoPath(proj);
    const [px0, py0] = proj([ext.lon0, ext.lat1]), [px1, py1] = proj([ext.lon1, ext.lat0]);
    const cid = 'clip' + nm + Math.round(x0);
    g.append('clipPath').attr('id', cid).append('rect').attr('x', px0).attr('y', py0).attr('width', px1 - px0).attr('height', py1 - py0);
    const mg = g.append('g').attr('clip-path', `url(#${cid})`);
    mg.append('rect').attr('x', px0).attr('y', py0).attr('width', px1 - px0).attr('height', py1 - py0).attr('fill', '#eef3f6');
    mg.append('image').attr('href', D.map.relief).attr('x', px0).attr('y', py0).attr('width', px1 - px0).attr('height', py1 - py0).attr('preserveAspectRatio', 'none');
    mg.append('path').datum(d3.geoGraticule().step([5, 5])()).attr('d', path).attr('fill', 'none').attr('stroke', '#fff').attr('stroke-width', 0.6).attr('opacity', 0.7);
    mg.selectAll('path.c').data(D.map.countries).join('path').attr('d', path).attr('fill', 'none').attr('stroke', '#777').attr('stroke-width', 0.6);
    D.map.contour3000.forEach(line => mg.append('path').datum({ type: 'LineString', coordinates: line }).attr('d', path).attr('fill', 'none').attr('stroke', '#5a4a35').attr('stroke-width', 0.5).attr('opacity', 0.55));
    g.append('rect').attr('x', px0).attr('y', py0).attr('width', px1 - px0).attr('height', py1 - py0).attr('fill', 'none').attr('stroke', '#333');
    // ticks
    for (let lo = -80; lo <= -60; lo += 5) { const [x] = proj([lo, 0]); g.append('text').attr('class', 'sm').attr('x', x).attr('y', py1 + 12).attr('text-anchor', 'middle').text(`${-lo}°W`); }
    for (let la = -20; la <= 0; la += 5) { const [, y] = proj([0, la]); g.append('text').attr('class', 'sm').attr('x', px0 - 4).attr('y', y + 3).attr('text-anchor', 'end').text(la === 0 ? '0°' : `${-la}°S`); }
    g.append('text').attr('class', 'ttl').attr('x', px0).attr('y', 16).text(`${nm} · ${varTitle(state)} · ${seasonName(state)} · lag ${state.lag} yr`).attr('fill', C[nm]);
    // other (non-analysed) sites
    D.other_sites.forEach(s => {
      const [x, y] = proj([s.lon, s.lat]);
      g.append('path').attr('d', d3.symbol(d3.symbolTriangle, 55)()).attr('transform', `translate(${x},${y})`).attr('fill', '#fff').attr('stroke', '#666').attr('stroke-width', 1.2)
        .on('mousemove', ev => showTip(ev, `<b>${s.site}</b> (${s.country})<br>${Math.abs(s.lat).toFixed(2)}°S, ${s.elev} m<br>${s.status}`)).on('mouseleave', hideTip);
      const L = { 'Hualcán': [7, 27, 'start'], Sajama: [7, 14, 'start'], Coropuna: [-8, 10, 'end'], Chimborazo: [7, 3, 'start'] }[s.site] || [7, 3, 'start'];
      if (px1 - px0 > 260) g.append('text').attr('class', 'sm it').attr('x', x + L[0]).attr('y', y + L[1]).attr('text-anchor', L[2]).text(`${s.site} ${s.elev} m`);
    });
    // analysed cores, clustered by site
    const bySite = d3.group(cores, c => c.site);
    for (const [site, list] of bySite) {
      const [sx, sy] = proj([list[0].lon, list[0].lat]);
      g.append('circle').attr('cx', sx).attr('cy', sy).attr('r', 2).attr('fill', '#111');
      const k = list.length, R = Math.max(7, Math.min(11, (px1 - px0) / 30));
      list.forEach((c, i) => {
        const cx = sx - R - 5 - (k - 1 - i) * (2 * R + 3), cy = sy;
        const res = results[c.id]?.[nm]?.s;
        g.append('line').attr('x1', sx).attr('y1', sy).attr('x2', cx).attr('y2', cy).attr('stroke', '#111').attr('stroke-width', 0.6);
        const circ = g.append('circle').attr('cx', cx).attr('cy', cy).attr('r', R).style('cursor', 'pointer')
          .attr('fill', res ? rColor(res.r) : '#fff').attr('stroke', res ? '#111' : '#888')
          .attr('stroke-width', res && res.p < 0.05 ? 2.6 : 0.8).attr('stroke-dasharray', res ? null : '2 2');
        if (c.id === state.core) g.append('circle').attr('cx', cx).attr('cy', cy).attr('r', R + 4).attr('fill', 'none').attr('stroke', '#0b5d8f').attr('stroke-width', 1.6);
        g.append('text').attr('x', cx).attr('y', cy + 3.5).attr('text-anchor', 'middle').attr('font-size', 8).attr('font-weight', 700)
          .attr('fill', res && Math.abs(res.r) > 0.35 ? '#fff' : '#111').style('pointer-events', 'none').text(c.short);
        circ.on('click', () => selectCore(c.id))
          .on('mousemove', ev => showTip(ev, `<b>${c.name}</b><br>${Math.abs(c.lat).toFixed(2)}°S ${Math.abs(c.lon).toFixed(2)}°W · ${c.elev} m a.s.l.<br>` +
            (res ? `${nm}: r = ${fmt(res.r)} (p<sub>eff</sub> ${eqP(res.p)}, n = ${res.n}, N<sub>eff</sub> = ${res.neff.toFixed(0)})` : 'No overlapping data for this variable/window') + '<br><i>click to open record</i>'))
          .on('mouseleave', hideTip);
      });
      g.append('text').attr('class', 'lab').attr('font-weight', 700).attr('x', sx + 6).attr('y', sy - 7).text(site);
      g.append('text').attr('class', 'sm').attr('x', sx + 6).attr('y', sy + 12).text(`${d3.extent(list, c => c.elev).join('–').replace(/^(\d+)–\1$/, '$1')} m`);
    }
    // latitude–elevation inset (shares the map's latitude axis)
    const ix0 = px1 + 18, ix1 = Math.min(w - 8, ix0 + insetW);
    const xe = d3.scaleLinear().domain([3, 7]).range([ix0, ix1]);
    const yl = la => proj([ext.lon0, la])[1];
    g.append('rect').attr('x', ix0).attr('y', py0).attr('width', ix1 - ix0).attr('height', py1 - py0).attr('fill', '#fafafa').attr('stroke', '#333');
    const crest = d3.area().x0(xe(3)).x1(d => xe(Math.max(3, d[1] / 1000))).y(d => yl(d[0])).curve(d3.curveBasis);
    g.append('path').datum(D.map.crest).attr('d', crest).attr('fill', '#d9d4c7').attr('stroke', '#8b826d').attr('stroke-width', 0.6);
    g.append('g').attr('class', 'ax').attr('transform', `translate(0,${py1})`).call(d3.axisBottom(xe).ticks(4).tickSize(3));
    g.append('text').attr('class', 'sm').attr('x', (ix0 + ix1) / 2).attr('y', py1 + 26).attr('text-anchor', 'middle').text('Elevation (km a.s.l.)');
    g.append('text').attr('class', 'sm').attr('x', ix0 + 3).attr('y', py0 + 10).text('Andean crest');
    D.other_sites.forEach(s => g.append('path').attr('d', d3.symbol(d3.symbolTriangle, 36)()).attr('transform', `translate(${xe(s.elev / 1000)},${yl(s.lat)})`).attr('fill', '#fff').attr('stroke', '#666'));
    cores.forEach(c => {
      const res = results[c.id]?.[nm]?.s;
      g.append('circle').attr('cx', xe(c.elev / 1000)).attr('cy', yl(c.lat)).attr('r', 4.5).attr('fill', res ? rColor(res.r) : '#fff')
        .attr('stroke', res ? '#111' : '#888').attr('stroke-width', res && res.p < 0.05 ? 1.8 : 0.7).style('cursor', 'pointer').on('click', () => selectCore(c.id))
        .on('mousemove', ev => showTip(ev, `<b>${c.name}</b><br>${c.elev} m, ${Math.abs(c.lat).toFixed(2)}°S`)).on('mouseleave', hideTip);
    });
    // colour bar
    const cbW = Math.max(60, Math.min(220, mapW - 60)), cbx = px0 + 6, cby = H - 22;
    const gid = 'grad' + nm + Math.round(x0), gr = g.append('linearGradient').attr('id', gid);
    d3.range(0, 1.01, 0.1).forEach(t => gr.append('stop').attr('offset', t).attr('stop-color', rColor(-0.6 + 1.2 * t)));
    g.append('rect').attr('x', cbx).attr('y', cby).attr('width', cbW).attr('height', 8).attr('fill', `url(#${gid})`).attr('stroke', '#555').attr('stroke-width', 0.5);
    const cbs = d3.scaleLinear().domain([-0.6, 0.6]).range([cbx, cbx + cbW]);
    g.append('g').attr('class', 'ax').attr('transform', `translate(0,${cby + 8})`).call(d3.axisBottom(cbs).ticks(5).tickSize(3).tickFormat(d3.format('.1f'))).select('.domain').remove();
    g.append('text').attr('class', 'sm').attr('x', cbx + cbW + 8).attr('y', cby + 8).text(`Pearson r (vs ${nm})`);
  }
  function siteList(svg, x0, w, H, results) {
    const g = svg.append('g').attr('transform', `translate(${x0},0)`);
    g.append('text').attr('class', 'ttl').attr('x', 0).attr('y', 16).text('Cores analysed (annual resolution, NOAA NCEI)');
    const nar = w < 520, cols = nar ? [0, w * 0.40, w * 0.56, w * 0.70, w * 0.86] : [0, w * 0.44, w * 0.58, w * 0.71, w * 0.85];
    const hy = 38;
    ['Core', 'Lat.', 'Elev.', 'r ONI', 'r RONI'].forEach((t, i) => g.append('text').attr('class', 'sm').attr('font-weight', 700).attr('x', cols[i] + (i ? 40 : 0)).attr('text-anchor', i ? 'end' : 'start').attr('y', hy).text(t));
    g.append('line').attr('x1', 0).attr('x2', w).attr('y1', hy + 5).attr('y2', hy + 5).attr('stroke', '#111');
    const rowH = Math.min(34, (H - 130) / cores.length);
    cores.forEach((c, i) => {
      const y = hy + 8 + (i + 0.7) * rowH, row = g.append('g').style('cursor', 'pointer').on('click', () => selectCore(c.id));
      if (c.id === state.core) row.append('rect').attr('x', -4).attr('y', y - rowH * 0.62).attr('width', w + 4).attr('height', rowH).attr('fill', '#eef4f8');
      const P = proxyOf(c, state);
      row.append('text').attr('class', 'lab').attr('x', 0).attr('y', y).attr('font-weight', c.id === state.core ? 700 : 400).text(nar ? `${c.short} · ${c.site}` : c.name.replace(/ \(.*\)$/, '').replace(' composite', ''));
      row.append('text').attr('class', 'sm').attr('x', 0).attr('y', y + 11).text(P ? `${P.years[0]}–${P.years[P.years.length - 1]} · ${c.id}` : `${varTitle(state)} not measured · ${c.id}`);
      row.append('text').attr('class', 'lab').attr('x', cols[1] + 40).attr('y', y).attr('text-anchor', 'end').text(`${Math.abs(c.lat).toFixed(1)}°S`);
      row.append('text').attr('class', 'lab').attr('x', cols[2] + 40).attr('y', y).attr('text-anchor', 'end').text(c.elev);
      IDX.forEach((nm, j) => {
        const s = results[c.id]?.[nm]?.s, xx = cols[3 + j] + 40;
        if (s) {
          row.append('rect').attr('x', xx - 44).attr('y', y - 11).attr('width', 46).attr('height', 15).attr('rx', 2).attr('fill', rColor(s.r));
          row.append('text').attr('class', 'lab').attr('x', xx - 2).attr('y', y).attr('text-anchor', 'end').attr('font-weight', s.p < 0.05 ? 700 : 400)
            .attr('fill', Math.abs(s.r) > 0.35 ? '#fff' : '#111').text(fmt(s.r) + sigStar(s.p));
        } else row.append('text').attr('class', 'sm').attr('x', xx - 2).attr('y', y).attr('text-anchor', 'end').text('n/a');
      });
    });
    const fy = H - 44;
    g.append('text').attr('class', 'sm').attr('x', 0).attr('y', fy).text(`* p<0.05, ** p<0.01 (two-sided, N_eff-corrected). Source: ${SRC_LABEL[state.source]}.`);
    g.append('text').attr('class', 'sm').attr('x', 0).attr('y', fy + 13).text(`Window ${state.y0}–${state.y1} CE${state.detrend ? ', linearly detrended' : ''}. Lag > 0: ENSO leads the ice.`);
    if (state.source === 'lmr') g.append('text').attr('class', 'sm').attr('x', 0).attr('y', fy + 26).attr('fill', '#8a4b00').text('LMR: r shown is the ensemble median; see panel c for the propagated uncertainty.');
  }

  // ------------------------------------------------------------------ (b) time series
  function renderTimeseries(res) {
    const core = coreById[state.core], W = widthOf('#fig-b');
    const panels = state.mode === 'SPLIT' ? ['ONI', 'RONI'] : [state.mode];
    const PH = 250, H = 64 + panels.length * PH + 30;
    const svg = newSvg('#fig-b', W, H);
    d3.select('#title-b').text(`${core.name}: ${proxyLabel(core, state)} vs. ${seasonName(state)} ONI and RONI`);
    const A0 = res.ONI.A;
    const m = { l: 56, r: 62 };
    if (!A0 || A0.n < 8) { svg.append('text').attr('x', 20).attr('y', 40).attr('class', 'lab').text('Fewer than 8 overlapping years for this variable, source and window — widen the window or change source.'); return; }
    const x = d3.scaleLinear().domain([state.y0 - 0.5, state.y1 + 0.5]).range([m.l, W - m.r]);
    // event strips (monthly episodes; annual probability for LMR)
    IDX.forEach((nm, k) => {
      const y = 22 + k * 15;
      svg.append('text').attr('class', 'sm').attr('x', m.l - 6).attr('y', y + 8).attr('text-anchor', 'end').attr('font-weight', 700).attr('fill', C[nm]).text(nm);
      svg.append('rect').attr('x', m.l).attr('y', y).attr('width', W - m.l - m.r).attr('height', 11).attr('fill', '#f3f3f3');
      if (state.source === 'lmr') {
        for (let Y = state.y0; Y <= state.y1; Y++) {
          const c = yearClass(nm, Y, state); if (!c) continue;
          const pEN = c.filter(v => v === 1).length / c.length, pLN = c.filter(v => v === -1).length / c.length;
          if (pEN > 0.05) svg.append('rect').attr('x', x(Y - 0.5)).attr('y', y).attr('width', x(Y + 0.5) - x(Y - 0.5)).attr('height', 11).attr('fill', C.EN).attr('opacity', pEN);
          if (pLN > 0.05) svg.append('rect').attr('x', x(Y - 0.5)).attr('y', y).attr('width', x(Y + 0.5) - x(Y - 0.5)).attr('height', 11).attr('fill', C.LN).attr('opacity', pLN);
        }
      } else {
        episodes(state.source, nm, state.thr, state.y0 - 1, state.y1 + 1).forEach(e => {
          // the episode month t maps to proxy year t + lag on this axis (DJF of Jan Y belongs to year Y)
          const a = Math.max(x.domain()[0], e.t0 + state.lag - 0.5), b = Math.min(x.domain()[1], e.t1 + state.lag - 0.5);
          if (b <= a) return;
          svg.append('rect').attr('x', x(a)).attr('y', y).attr('width', x(b) - x(a)).attr('height', 11).attr('fill', e.sign > 0 ? C.EN : C.LN)
            .on('mousemove', ev => showTip(ev, `${nm} ${e.sign > 0 ? 'El Niño' : 'La Niña'} episode<br>${MON[e.m0]} ${e.ys} – ${MON[e.m1]} ${e.ye}<br>peak ${fmt(e.peak, 1)} °C`)).on('mouseleave', hideTip);
        });
      }
    });
    svg.append('text').attr('class', 'sm').attr('x', W - m.r + 4).attr('y', 38).text(state.source === 'lmr' ? 'P(event)' : 'episodes');
    panels.forEach((shadeIdx, pi) => {
      const top = 64 + pi * PH, bot = top + PH - 40;
      const g = svg.append('g');
      const A = res[shadeIdx].A;
      const pz = A.p, Y = A.Y;
      const zmax = Math.ceil(Math.max(2.5, d3.max(pz, Math.abs)));
      const yz = d3.scaleLinear().domain([-zmax, zmax]).range([bot, top + 8]);   // symmetric: zero aligned with index axis
      // index values (plotted at proxy year)
      const idxVals = {};
      IDX.forEach(nm => {
        const B = res[nm].A;
        idxVals[nm] = B.ens ? B.Y.map((yy, i) => ({ Y: yy, v: d3.median(B.x, c => c[i]), lo: d3.quantile(B.x.map(c => c[i]).sort(d3.ascending), 0.05), hi: d3.quantile(B.x.map(c => c[i]).sort(d3.ascending), 0.95) }))
          : B.Y.map((yy, i) => ({ Y: yy, v: B.x[i] }));
      });
      const allIdx = IDX.flatMap(nm => idxVals[nm].flatMap(d => [d.v, d.lo ?? d.v, d.hi ?? d.v]));
      const ext = Math.max(1.5, d3.max(allIdx, Math.abs));
      const yi = d3.scaleLinear().domain([-Math.ceil(ext * 2) / 2, Math.ceil(ext * 2) / 2]).range([bot, top + 8]);
      // year shading by class of the shading index
      Y.forEach(yy => {
        const c = yearClass(shadeIdx, yy - state.lag, state); if (c === null) return;
        if (Array.isArray(c)) {
          const pEN = c.filter(v => v === 1).length / c.length, pLN = c.filter(v => v === -1).length / c.length;
          if (pEN > 0.5) g.append('rect').attr('x', x(yy - 0.5)).attr('y', top).attr('width', x(yy + 0.5) - x(yy - 0.5)).attr('height', bot - top).attr('fill', C.EN).attr('opacity', 0.16 * pEN);
          if (pLN > 0.5) g.append('rect').attr('x', x(yy - 0.5)).attr('y', top).attr('width', x(yy + 0.5) - x(yy - 0.5)).attr('height', bot - top).attr('fill', C.LN).attr('opacity', 0.16 * pLN);
        } else if (c !== 0) g.append('rect').attr('x', x(yy - 0.5)).attr('y', top).attr('width', x(yy + 0.5) - x(yy - 0.5)).attr('height', bot - top).attr('fill', c > 0 ? C.EN : C.LN).attr('opacity', 0.15);
      });
      g.append('g').attr('class', 'gridl').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(yz).ticks(5).tickSize(-(W - m.l - m.r)).tickFormat(''));
      g.append('line').attr('class', 'zero').attr('x1', m.l).attr('x2', W - m.r).attr('y1', yz(0)).attr('y2', yz(0));
      // index lines
      IDX.forEach(nm => {
        const dat = idxVals[nm], emph = state.mode === 'SPLIT' ? nm === shadeIdx : nm === state.mode;
        if (state.mode === 'SPLIT' && nm !== shadeIdx) return;
        if (dat[0] && dat[0].lo !== undefined) g.append('path').datum(dat).attr('d', d3.area().x(d => x(d.Y)).y0(d => yi(d.lo)).y1(d => yi(d.hi))).attr('fill', C[nm]).attr('opacity', 0.18);
        g.append('path').datum(dat).attr('d', d3.line().x(d => x(d.Y)).y(d => yi(d.v)).defined(d => d.v !== null))
          .attr('fill', 'none').attr('stroke', C[nm]).attr('stroke-width', emph ? 2 : 1.2).attr('stroke-dasharray', emph ? null : '4 3');
      });
      // proxy
      const pd = Y.map((yy, i) => ({ Y: yy, z: pz[i], raw: A.pRaw[i] }));
      g.append('path').datum(pd).attr('d', d3.line().x(d => x(d.Y)).y(d => yz(d.z))).attr('fill', 'none').attr('stroke', C.proxy).attr('stroke-width', 1.5);
      g.selectAll('circle.pp').data(pd).join('circle').attr('cx', d => x(d.Y)).attr('cy', d => yz(d.z)).attr('r', Y.length > 150 ? 0 : 2.2).attr('fill', C.proxy);
      // axes
      g.append('g').attr('class', 'ax').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(yz).ticks(5));
      g.append('g').attr('class', 'ax').attr('transform', `translate(${W - m.r},0)`).call(d3.axisRight(yi).ticks(5));
      g.append('g').attr('class', 'ax').attr('transform', `translate(0,${bot})`).call(d3.axisBottom(x).ticks(Math.min(12, W / 80)).tickFormat(d3.format('d')));
      g.append('text').attr('class', 'lab').attr('transform', `translate(14,${(top + bot) / 2}) rotate(-90)`).attr('text-anchor', 'middle').text(`${varTitle(state)} anomaly (σ)${state.detrend ? ', detrended' : ''}`);
      g.append('text').attr('class', 'lab').attr('transform', `translate(${W - 14},${(top + bot) / 2}) rotate(90)`).attr('text-anchor', 'middle')
        .text(`${state.mode === 'SPLIT' ? shadeIdx : 'ONI, RONI'} ${seasonName(state)} (${state.source === 'lmr' ? 'annual, °C' : '°C'})`);
      g.append('text').attr('class', 'lab').attr('x', (m.l + W - m.r) / 2).attr('y', bot + 30).attr('text-anchor', 'middle').text(`Ice year (CE)${state.lag ? ` — index shifted: ENSO year Y − (${state.lag}) plotted at Y` : ''}`);
      // legend
      const lg = g.append('g').attr('transform', `translate(${m.l + 8},${top + 4})`);
      const items = [['Proxy (z)', C.proxy, null]].concat((state.mode === 'SPLIT' ? [shadeIdx] : IDX).map(nm => [nm, C[nm], (state.mode !== 'SPLIT' && nm !== state.mode) ? '4 3' : null]))
        .concat([[`${shadeIdx} El Niño yr`, C.EN, 'box'], [`${shadeIdx} La Niña yr`, C.LN, 'box']]);
      let lx = 0;
      items.forEach(([t, col, dash]) => {
        if (dash === 'box') lg.append('rect').attr('x', lx).attr('y', 1).attr('width', 12).attr('height', 9).attr('fill', col).attr('opacity', 0.3);
        else lg.append('line').attr('x1', lx).attr('x2', lx + 14).attr('y1', 6).attr('y2', 6).attr('stroke', col).attr('stroke-width', 2).attr('stroke-dasharray', dash);
        lg.append('text').attr('class', 'sm').attr('x', lx + 18).attr('y', 9.5).text(t); lx += 24 + t.length * 5.6;
      });
      const s = res[shadeIdx].s;
      if (s) g.append('text').attr('class', 'lab').attr('x', W - m.r - 6).attr('y', top + 14).attr('text-anchor', 'end').attr('font-weight', 700)
        .text(`r(${shadeIdx}) = ${fmt(s.r)}${s.ens ? ` [${fmt(s.ci[0])}, ${fmt(s.ci[1])}]` : ''}, p_eff ${eqP(s.p)}, n = ${s.n}`);
      // hover
      g.append('rect').attr('x', m.l).attr('y', top).attr('width', W - m.l - m.r).attr('height', bot - top).attr('fill', 'transparent')
        .on('mousemove', ev => {
          const yy = Math.round(x.invert(d3.pointer(ev)[0])), i = Y.indexOf(yy); if (i < 0) return hideTip();
          const cls = IDX.map(nm => { const c = yearClass(nm, yy - state.lag, state); return `${nm}: ${fmt(idxVals[nm][i]?.v)} °C · ${Array.isArray(c) ? `P(EN)=${(c.filter(v => v === 1).length / c.length).toFixed(2)}` : clsName(c)}`; });
          showTip(ev, `<b>${yy}</b> — proxy ${fmt(A.pRaw[i], 2)} ${proxyOf(core, state).units} (z = ${fmt(pz[i])})<br>${cls.join('<br>')}`);
        }).on('mouseleave', hideTip);
    });
    d3.select('#core-card').html(`<b>${core.name}</b> — ${Math.abs(core.lat).toFixed(2)}°S, ${Math.abs(core.lon).toFixed(2)}°W, ${core.elev} m a.s.l. · ${core.ref}${core.doi ? ` · doi:<a href="https://doi.org/${core.doi}" target="_blank" rel="noopener">${core.doi}</a>` : ''} · <a href="${core.url}" target="_blank" rel="noopener">NCEI data</a><br>Dating: ${core.dating} Year convention: ${core.year_def}.` +
      (state.source === 'lmr' && core.in_lmr ? `<br><span class="warn"><b>Circularity warning:</b> Quelccaya δ¹⁸O is in the PAGES2k (2017) network assimilated by LMRv2.1, so LMR-based correlations for this core are not independent.</span>` : ''));
  }

  // ------------------------------------------------------------------ (c) statistics + scatter
  function renderStats(res) {
    const rows = [
      ['n (years)', s => s.n, 0], ['Pearson r', s => s.r, 2, 'p'], ['95% CI (Fisher, N_eff)', s => `[${fmt(s.ci[0])}, ${fmt(s.ci[1])}]`],
      ['p naive (n)', s => fmtP(s.p_naive)], ['Lag-1 autocorr. proxy / index', s => `${fmt(s.ar1y)} / ${fmt(s.ar1x)}`],
      ['N_eff (Bretherton 1999)', s => s.neff, 1], ['p (N_eff)', s => fmtP(s.p), null, 'p'], ['Spearman ρ', s => s.rho, 2], ['p_ρ (N_eff)', s => fmtP(s.p_rho)]
    ];
    if (state.source === 'lmr') rows.splice(3, 0, ['Ensemble-only 95% range', s => `[${fmt(s.r_lo)}, ${fmt(s.r_hi)}]`], ['Members with p<0.05', s => (100 * s.fracSig).toFixed(0) + '%']);
    const sel = state.mode;
    let h = `<table class="st"><thead><tr><th>${varTitle(state)} vs ${seasonName(state)} index</th>${IDX.map(nm => `<th class="${sel === nm ? 'sel' : ''}" style="color:${C[nm]}">${nm}</th>`).join('')}</tr></thead><tbody>`;
    for (const [lab, f, d] of rows) {
      h += `<tr><td>${lab}</td>${IDX.map(nm => { const s = res[nm].s; if (!s) return `<td>—</td>`; const v = f(s); return `<td class="${sel === nm ? 'sel' : ''} ${lab.startsWith('p (') && s.p < 0.05 ? 'sig' : ''}">${typeof v === 'number' ? fmt(v, d) : v}</td>`; }).join('')}</tr>`;
    }
    if (res.diff) h += `<tr><td>ΔR = r<sub>RONI</sub> − r<sub>ONI</sub></td><td colspan="2" style="text-align:center">${fmt(res.diff.dr, 3)} (r<sub>ONI,RONI</sub> = ${fmt(res.diff.r12)}; MRR Z = ${fmt(res.diff.Z)}, p ${eqP(res.diff.p)})</td></tr>`;
    h += `</tbody></table><div class="small">${state.source === 'lmr' ? 'LMR: values are ensemble medians over 100 reconstruction members; CI combines reconstruction spread and sampling error (Fisher z, N_eff). ' : ''}Significance of ΔR: Meng, Rosenthal & Rubin (1992) test for dependent correlations sharing the proxy.</div>`;
    d3.select('#stats-table').html(h);
    // scatter
    const W = widthOf('#fig-c'), H = 280, m = { l: 46, r: 12, t: 16, b: 38 };
    const svg = newSvg('#fig-c', W, H);
    const pts = [];
    IDX.forEach(nm => { const A = res[nm].A; if (!A || A.n < 8) return; A.Y.forEach((y, i) => pts.push({ nm, y, px: A.p[i], ix: A.ens ? d3.median(A.x, c => c[i]) : A.x[i] })); });
    if (!pts.length) return;
    const x = d3.scaleLinear().domain(d3.extent(pts, d => d.ix)).nice().range([m.l, W - m.r]);
    const y = d3.scaleLinear().domain(d3.extent(pts, d => d.px)).nice().range([H - m.b, m.t]);
    svg.append('g').attr('class', 'gridl').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5).tickSize(-(W - m.l - m.r)).tickFormat(''));
    const shown = state.mode === 'SPLIT' ? IDX : [state.mode, ...IDX.filter(n => n !== state.mode)];
    shown.slice().reverse().forEach(nm => {
      const P = pts.filter(d => d.nm === nm), main = state.mode === 'SPLIT' || nm === state.mode;
      svg.selectAll(null).data(P).join('circle').attr('cx', d => x(d.ix)).attr('cy', d => y(d.px)).attr('r', 3)
        .attr('fill', main ? C[nm] : 'none').attr('stroke', C[nm]).attr('opacity', main ? 0.75 : 0.45)
        .on('mousemove', (ev, d) => showTip(ev, `${d.y}: ${nm} ${fmt(d.ix)} °C, proxy z ${fmt(d.px)}`)).on('mouseleave', hideTip);
      const b = S.pearson(P.map(d => d.ix), P.map(d => d.px)) * S.sd(P.map(d => d.px)) / S.sd(P.map(d => d.ix));
      const mx = S.mean(P.map(d => d.ix)), my = S.mean(P.map(d => d.px)), [a0, a1] = x.domain();
      svg.append('line').attr('x1', x(a0)).attr('x2', x(a1)).attr('y1', y(my + b * (a0 - mx))).attr('y2', y(my + b * (a1 - mx))).attr('stroke', C[nm]).attr('stroke-width', main ? 2 : 1.2).attr('stroke-dasharray', main ? null : '4 3');
    });
    svg.append('g').attr('class', 'ax').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(6));
    svg.append('g').attr('class', 'ax').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5));
    svg.append('text').attr('class', 'lab').attr('x', (m.l + W - m.r) / 2).attr('y', H - 6).attr('text-anchor', 'middle').text(`${seasonName(state)} index (°C), ENSO year Y − ${state.lag}`);
    svg.append('text').attr('class', 'lab').attr('transform', `translate(12,${(H - m.b + m.t) / 2}) rotate(-90)`).attr('text-anchor', 'middle').text(`${varTitle(state)} (σ)`);
    IDX.forEach((nm, i) => { const s = res[nm].s; if (s) svg.append('text').attr('class', 'lab').attr('x', m.l + 8).attr('y', m.t + 12 + i * 14).attr('fill', C[nm]).attr('font-weight', 700).text(`${nm}: r = ${fmt(s.r)}, ρ = ${fmt(s.rho)}`); });
  }

  // ------------------------------------------------------------------ (d) composites
  function renderComposites(res) {
    const comp = { ONI: composites(res.ONI.A, 'ONI', state), RONI: composites(res.RONI.A, 'RONI', state) };
    const W = widthOf('#fig-d'), H = 320, m = { l: 50, r: 12, t: 40, b: 44 };
    const svg = newSvg('#fig-d', W, H);
    if (!comp.ONI && !comp.RONI) { svg.append('text').attr('class', 'lab').attr('x', 20).attr('y', 40).text('Insufficient overlap for composites.'); d3.select('#note-d').text(''); return comp; }
    const groups = ['EN', 'N', 'LN'], gname = { EN: 'El Niño years', N: 'Neutral years', LN: 'La Niña years' };
    const x0 = d3.scaleBand().domain(groups).range([m.l, W - m.r]).paddingInner(0.25).paddingOuter(0.1);
    const x1 = d3.scaleBand().domain(IDX).range([0, x0.bandwidth()]).padding(0.08);
    const vals = IDX.flatMap(nm => comp[nm] ? groups.flatMap(g => comp[nm][g].ci.filter(Number.isFinite)) : []);
    const y = d3.scaleLinear().domain([Math.min(-0.8, d3.min(vals)), Math.max(0.8, d3.max(vals))]).nice().range([H - m.b, m.t]);
    svg.append('g').attr('class', 'gridl').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6).tickSize(-(W - m.l - m.r)).tickFormat(''));
    groups.forEach(g => IDX.forEach(nm => {
      const c = comp[nm]; if (!c) return; const s = c[g]; if (!s.n) return;
      const bx = x0(g) + x1(nm), bw = x1.bandwidth(), fill = g === 'EN' ? C.EN : g === 'LN' ? C.LN : C.N;
      svg.append('rect').attr('x', bx).attr('width', bw).attr('y', y(Math.max(0, s.mean))).attr('height', Math.abs(y(s.mean) - y(0)))
        .attr('fill', fill).attr('opacity', nm === 'ONI' ? 0.45 : 0.9).attr('stroke', C[nm]).attr('stroke-width', 2.2);
      if (Number.isFinite(s.ci[0])) {
        const cx = bx + bw / 2;
        svg.append('line').attr('x1', cx).attr('x2', cx).attr('y1', y(s.ci[0])).attr('y2', y(s.ci[1])).attr('stroke', '#111');
        [0, 1].forEach(k => svg.append('line').attr('x1', cx - 4).attr('x2', cx + 4).attr('y1', y(s.ci[k])).attr('y2', y(s.ci[k])).attr('stroke', '#111'));
      }
      svg.append('text').attr('class', 'sm').attr('x', bx + bw / 2).attr('y', H - m.b + 26).attr('text-anchor', 'middle').text(`${nm} n=${s.n}`);
    }));
    svg.append('line').attr('x1', m.l).attr('x2', W - m.r).attr('y1', y(0)).attr('y2', y(0)).attr('stroke', '#333');
    svg.append('g').attr('class', 'ax').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6));
    groups.forEach(g => svg.append('text').attr('class', 'lab').attr('x', x0(g) + x0.bandwidth() / 2).attr('y', H - m.b + 13).attr('text-anchor', 'middle').attr('font-weight', 700).text(gname[g]));
    svg.append('text').attr('class', 'lab').attr('transform', `translate(13,${(H - m.b + m.t) / 2}) rotate(-90)`).attr('text-anchor', 'middle').text(`Mean ${varTitle(state)} anomaly (σ)`);
    IDX.forEach((nm, i) => { const c = comp[nm]; if (!c) return; svg.append('text').attr('class', 'lab').attr('x', m.l + 4).attr('y', 12 + i * 14).attr('fill', C[nm]).attr('font-weight', 700).text(`${nm}: EN − LN = ${fmt(c.diff)} σ (Welch p ${eqP(c.p_welch)}${Number.isFinite(c.p_perm) ? `, perm. p ${eqP(c.p_perm)}` : ''})`); });
    d3.select('#note-d').html(`Years classified by the ${state.source === 'lmr' ? `annual-mean index beyond ±${state.thr}σ, per ensemble member` : `CPC episode rule (≥5 consecutive overlapping 3-month seasons beyond ±${state.thr} °C) evaluated for DJF of ENSO year Y − ${state.lag}`}. Bars: mean proxy z-anomaly; whiskers: 95% CI (${state.source === 'lmr' ? 'ensemble spread ⊕ sampling error' : 't-based'}). Faded/solid bars outlined in ONI/RONI colour.`);
    return comp;
  }

  // ------------------------------------------------------------------ (e) what changed
  function renderWhatChanged(all, res) {
    const W = widthOf('#fig-e'), list = cores.filter(c => all[c.id] && (all[c.id].ONI.s || all[c.id].RONI.s));
    const rowH = 30, m = { l: 190, r: 118, t: 34, b: 44 }, H = m.t + m.b + Math.max(1, list.length) * rowH;
    const svg = newSvg('#fig-e', W, H);
    const rs = list.flatMap(c => IDX.map(nm => all[c.id][nm].s?.r).filter(Number.isFinite));
    const x = d3.scaleLinear().domain([Math.min(-0.5, d3.min(rs) - 0.05), Math.max(0.5, d3.max(rs) + 0.05)]).nice().range([m.l, W - m.r - 10]).clamp(true);
    svg.append('g').attr('class', 'gridl').attr('transform', `translate(0,${m.t})`).call(d3.axisTop(x).ticks(6).tickSize(-(H - m.t - m.b)).tickFormat(''));
    svg.append('line').attr('x1', x(0)).attr('x2', x(0)).attr('y1', m.t).attr('y2', H - m.b).attr('stroke', '#555');
    svg.append('g').attr('class', 'ax').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(6));
    svg.append('text').attr('class', 'lab').attr('x', (m.l + W - m.r) / 2).attr('y', H - 10).attr('text-anchor', 'middle').text(`Pearson r (${varTitle(state)} vs ${seasonName(state)} index, lag ${state.lag} yr)`);
    svg.append('text').attr('class', 'ttl').attr('x', W - m.r + 8).attr('y', 22).text('ΔR (MRR p)');
    [['ONI', 'none'], ['RONI', C.RONI]].forEach(([nm, f], i) => {
      svg.append('circle').attr('cx', m.l + i * 70 + 6).attr('cy', 16).attr('r', 5).attr('fill', nm === 'ONI' ? '#fff' : f).attr('stroke', C[nm]).attr('stroke-width', 2);
      svg.append('text').attr('class', 'lab').attr('x', m.l + i * 70 + 15).attr('y', 20).text(nm);
    });
    svg.append('text').attr('class', 'sm').attr('x', m.l + 140).attr('y', 20).text('arrow ONI → RONI; bar = 95% CI; ● filled = p<0.05');
    const csv = [['core', 'n', 'r_ONI', 'p_ONI', 'r_RONI', 'p_RONI', 'dR', 'p_dR']];
    list.forEach((c, i) => {
      const y = m.t + (i + 0.5) * rowH, o = all[c.id].ONI.s, r = all[c.id].RONI.s, d = all[c.id].diff;
      if (c.id === state.core) svg.append('rect').attr('x', 4).attr('y', y - rowH / 2).attr('width', W - 8).attr('height', rowH).attr('fill', '#eef4f8');
      const lab = svg.append('g').style('cursor', 'pointer').on('click', () => selectCore(c.id));
      lab.append('text').attr('class', 'lab').attr('x', 10).attr('y', y - 1).attr('font-weight', c.id === state.core ? 700 : 400).text(c.name.replace(/ \(.*\)$/, ''));
      lab.append('text').attr('class', 'sm').attr('x', 10).attr('y', y + 10).text(`${c.id} · n = ${o ? o.n : r ? r.n : 0}`);
      [[o, 'ONI', -5], [r, 'RONI', 5]].forEach(([s, nm, dy]) => {
        if (!s) return;
        svg.append('line').attr('x1', x(s.ci[0])).attr('x2', x(s.ci[1])).attr('y1', y + dy).attr('y2', y + dy).attr('stroke', C[nm]).attr('stroke-width', 1.2).attr('opacity', 0.7);
      });
      if (o && r) svg.append('line').attr('x1', x(o.r)).attr('x2', x(r.r)).attr('y1', y - 5).attr('y2', y + 5).attr('stroke', '#333').attr('stroke-width', 1).attr('marker-end', null);
      [[o, 'ONI', -5], [r, 'RONI', 5]].forEach(([s, nm, dy]) => {
        if (!s) return;
        svg.append('circle').attr('cx', x(s.r)).attr('cy', y + dy).attr('r', 5).attr('fill', s.p < 0.05 ? C[nm] : '#fff').attr('stroke', C[nm]).attr('stroke-width', 2)
          .on('mousemove', ev => showTip(ev, `${c.name}<br>${nm}: r = ${fmt(s.r)} [${fmt(s.ci[0])}, ${fmt(s.ci[1])}], p_eff ${eqP(s.p)}`)).on('mouseleave', hideTip);
      });
      if (d) svg.append('text').attr('class', 'lab').attr('x', W - m.r + 8).attr('y', y + 4).attr('font-weight', d.p < 0.05 ? 700 : 400)
        .text(`${d.dr >= 0 ? '+' : '−'}${Math.abs(d.dr).toFixed(3)} (${fmtP(d.p)})`);
      csv.push([c.id, o?.n, o?.r, o?.p, r?.r, r?.p, d?.dr, d?.p]);
    });
    if (!list.length) svg.append('text').attr('class', 'lab').attr('x', m.l).attr('y', m.t + 18).text('No core has this variable in the selected window.');
    renderWhatChanged.csv = csv;
    // reclassified years
    const rc = reclassified();
    const T = rc.rows;
    let html = '';
    if (state.source === 'lmr') {
      html = `<table class="reclass"><thead><tr><th>ENSO yr</th><th>P(EN) ONI → RONI</th><th>P(LN) ONI → RONI</th><th>Tropical-mean anomaly</th></tr></thead><tbody>` +
        T.slice(0, 40).map(r => `<tr><td>${r.Y}</td><td>${r.pENo.toFixed(2)} → ${r.pENr.toFixed(2)}</td><td>${r.pLNo.toFixed(2)} → ${r.pLNr.toFixed(2)}</td><td>${fmt(r.t20)} °C</td></tr>`).join('') + '</tbody></table>' +
        (T.length > 40 ? `<div class="small">…and ${T.length - 40} more years.</div>` : '') + (!T.length ? '<div class="small">No year changes class probability by ≥0.25 in this window.</div>' : '');
    } else {
      html = `<table class="reclass"><thead><tr><th>ENSO yr (DJF)</th><th>ONI</th><th>RONI</th><th>Trop. mean*</th><th>Why</th></tr></thead><tbody>` +
        T.map(r => `<tr><td>${r.Y - 1}/${String(r.Y).slice(-2)}</td><td><span class="chip ${clsKey(r.co)}">${clsKey(r.co)}</span> ${fmt(r.vo, 1)}</td><td><span class="chip ${clsKey(r.cr)}">${clsKey(r.cr)}</span> ${fmt(r.vr, 1)}</td><td>${fmt(r.t20)} °C</td><td>${r.why}</td></tr>`).join('') + '</tbody></table><div class="small">*20°S–20°N SST anomaly (ERSST.v6, DJF) relative to the same centred 30-yr base as ONI.</div>' +
        (!T.length ? '<div class="small">No ENSO year in this window changes class between ONI and RONI at this threshold.</div>' : '');
    }
    d3.select('#reclass-e').html(`<div class="small" style="margin-top:6px"><b>Reclassified ENSO years, ${state.y0}–${state.y1}</b> (${SRC_LABEL[state.source]}${state.source === 'lmr' ? '' : `, DJF, ±${state.thr} °C, ≥5 overlapping seasons`})</div>` + html);
    d3.select('#why-e').html(whyText(all, rc));
    return rc;
  }
  function reclassified() {
    const rows = [];
    if (state.source === 'lmr') {
      for (let Y = state.y0; Y <= state.y1; Y++) {
        const o = yearClass('ONI', Y, state), r = yearClass('RONI', Y, state); if (!o || !r) continue;
        const f = (c, v) => c.filter(x => x === v).length / c.length;
        const row = { Y, pENo: f(o, 1), pENr: f(r, 1), pLNo: f(o, -1), pLNr: f(r, -1), t20: D.lmr.t20[Y - D.lmr.start] };
        if (Math.abs(row.pENo - row.pENr) >= 0.25 || Math.abs(row.pLNo - row.pLNr) >= 0.25) rows.push(row);
      }
      return { rows, n: rows.length };
    }
    const so = monthlySeries(state.source, 'oni'), sr = monthlySeries(state.source, 'roni');
    for (let Y = Math.max(state.y0, so.start + 1); Y <= state.y1; Y++) {
      const co = yearClass('ONI', Y, state), cr = yearClass('RONI', Y, state);
      if (co === null || cr === null || co === cr) continue;
      const vo = so.a[(Y - so.start) * 12], vr = sr.a[(Y - sr.start) * 12], t20 = t20Sliding(Y, 1);
      const d = vo - vr, volc = VOLC.find(v => Y - v.y >= 0 && Y - v.y <= 2);
      let why;
      if (Math.abs(d) < 0.15) why = `Values nearly equal (${fmt(vo, 1)} vs ${fmt(vr, 1)} °C); the ≥5-season duration rule tips the class.`;
      else if (d > 0 && t20 > 0.08) why = `Tropics warm (${fmt(t20)} °C vs 30-yr climatology): part of the Niño-3.4 anomaly is basin-wide, so RONI is ${fmt(d, 1)} °C lower.`;
      else if (d < 0 && t20 < -0.08) why = `Tropics cool (${fmt(t20)} °C vs 30-yr climatology)${volc ? ` after the ${volc.n} eruption (${volc.y})` : ''}: relative to the tropics Niño-3.4 is ${fmt(-d, 1)} °C warmer, so RONI is higher.`;
      else why = `Tropical mean near normal (${fmt(t20)} °C); the ${fmt(Math.abs(d), 1)} °C gap comes from RONI’s variance rescaling (×${D.ersst.k}) and ONI’s step-wise 30-yr base.`;
      rows.push({ Y, co, cr, vo, vr, t20, why });
    }
    return { rows, n: rows.length };
  }
  function whyText(all, rc) {
    const list = cores.filter(c => all[c.id]?.diff);
    const dr = list.map(c => all[c.id].diff.dr), sig = list.filter(c => all[c.id].diff.p < 0.05);
    const weaker = list.filter(c => Math.abs(all[c.id].RONI.s.r) < Math.abs(all[c.id].ONI.s.r));
    let t = '';
    if (list.length) {
      t += `<p>Across ${list.length} core${list.length > 1 ? 's' : ''} with ${varTitle(state)} in ${state.y0}–${state.y1}, switching from ONI to RONI changes r by ${fmt(d3.min(dr), 2)} to ${fmt(d3.max(dr), 2)} (median ${fmt(d3.median(dr), 2)}); |r| decreases for ${weaker.length} of ${list.length}. ` +
        (sig.length ? `The change is significant (MRR test, p<0.05) for ${sig.map(c => c.id).join(', ')}.` : `No change is statistically significant (MRR test), because ONI and RONI are themselves highly correlated (r ≈ ${fmt(d3.median(list, c => all[c.id].diff.r12))}) over this window.`) + '</p>';
    }
    t += `<p><b>Why the indices differ.</b> ONI measures Niño-3.4 against a sliding 30-year climatology; RONI measures Niño-3.4 against the simultaneous 20°S–20°N mean. The two diverge whenever the whole tropical ocean is anomalously warm or cool — e.g. after strong El Niños (tropical-wide warming lags Niño-3.4 by ~3–5 months), during volcanic cooling, and increasingly with anthropogenic warming. When a proxy responds to tropical-mean tropospheric temperature (a plausible pathway for Andean δ¹⁸O via upstream Amazonian rain-out and the temperature of the free troposphere), that shared variance is <i>removed</i> from RONI, so its correlation with the proxy should drop. A rise in |r| under RONI instead indicates that the proxy tracks the east–west SST gradient (Walker circulation, easterly moisture flux over the Altiplano) rather than basin-wide warmth.</p>`;
    t += `<p><b>${rc.n} ENSO year${rc.n === 1 ? '' : 's'} change class</b> in this window${state.source === 'lmr' ? ' (event probability shifts ≥ 0.25 across the reconstruction ensemble)' : ''}. Each reclassification moves that year between composite groups in panel d, which is why composite differences can change even when r barely moves.</p>`;
    return t;
  }

  // ------------------------------------------------------------------ (f) running correlation
  function runningCorr(core, nm, st, src) {
    const s2 = { ...st, source: src || st.source, y0: -1e9, y1: 1e9 };
    const A = align(core, nm, s2); if (!A || A.n < st.runW) return [];
    const h = (st.runW - 1) / 2, out = [];
    const YY = A.Y;
    for (let c = YY[0] + h; c <= YY[YY.length - 1] - h; c++) {
      const idx = []; for (let i = 0; i < YY.length; i++) if (YY[i] >= c - h && YY[i] <= c + h) idx.push(i);
      if (idx.length < 0.8 * st.runW) continue;
      const p = S.zscore(idx.map(i => A.p[i]));
      if (A.ens) {
        const rr = A.x.map(col => S.pearson(idx.map(i => col[i]), p)).filter(Number.isFinite);
        out.push({ c, r: d3.median(rr), lo: S.quantile(rr, 0.05), hi: S.quantile(rr, 0.95), neff: S.neff(idx.map(i => A.x[0][i]), p).neff });
      } else {
        const xx = idx.map(i => A.x[i]);
        out.push({ c, r: S.pearson(xx, p), neff: S.neff(xx, p).neff });
      }
    }
    return out;
  }
  function renderRunning() {
    const core = coreById[state.core], W = widthOf('#fig-f'), H = 300, m = { l: 52, r: 150, t: 20, b: 40 };
    const svg = newSvg('#fig-f', W, H);
    const R = { ONI: runningCorr(core, 'ONI', state), RONI: runningCorr(core, 'RONI', state) };
    let X = null;
    if (state.source === 'lmr') X = { ONI: runningCorr(core, 'ONI', state, 'ersst'), RONI: runningCorr(core, 'RONI', state, 'ersst') };
    const all = [...R.ONI, ...R.RONI, ...(X ? [...X.ONI, ...X.RONI] : [])];
    if (!all.length) { svg.append('text').attr('class', 'lab').attr('x', 20).attr('y', 40).text(`Record overlap shorter than the ${state.runW}-yr window.`); d3.select('#note-f').text(''); return R; }
    const x = d3.scaleLinear().domain(d3.extent(all, d => d.c)).range([m.l, W - m.r]);
    const y = d3.scaleLinear().domain([-1, 1]).range([H - m.b, m.t]);
    const rc = S.rCrit(state.runW), neMed = d3.median(all, d => d.neff), rcE = S.rCrit(Math.max(4, neMed));
    svg.append('rect').attr('x', Math.max(m.l, x(state.y0))).attr('width', Math.max(0, Math.min(W - m.r, x(state.y1)) - Math.max(m.l, x(state.y0)))).attr('y', m.t).attr('height', H - m.t - m.b).attr('fill', '#f4f1ea');
    svg.append('g').attr('class', 'gridl').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(8).tickSize(-(W - m.l - m.r)).tickFormat(''));
    [rc, -rc].forEach(v => svg.append('line').attr('x1', m.l).attr('x2', W - m.r).attr('y1', y(v)).attr('y2', y(v)).attr('stroke', '#777').attr('stroke-dasharray', '5 3'));
    [rcE, -rcE].forEach(v => svg.append('line').attr('x1', m.l).attr('x2', W - m.r).attr('y1', y(v)).attr('y2', y(v)).attr('stroke', '#777').attr('stroke-dasharray', '1 2'));
    svg.append('line').attr('x1', m.l).attr('x2', W - m.r).attr('y1', y(0)).attr('y2', y(0)).attr('stroke', '#333');
    if (state.source !== 'lmr' && x.domain()[0] < 1950 && state.source === 'ersst') svg.append('line').attr('x1', x(1950)).attr('x2', x(1950)).attr('y1', m.t).attr('y2', H - m.b).attr('stroke', '#999');
    IDX.forEach(nm => {
      const d = R[nm]; if (!d.length) return;
      if (d[0].lo !== undefined) svg.append('path').datum(d).attr('d', d3.area().x(v => x(v.c)).y0(v => y(v.lo)).y1(v => y(v.hi))).attr('fill', C[nm]).attr('opacity', 0.2);
      svg.append('path').datum(d).attr('d', d3.line().x(v => x(v.c)).y(v => y(v.r)).defined(v => Number.isFinite(v.r))).attr('fill', 'none').attr('stroke', C[nm]).attr('stroke-width', state.mode === nm || state.mode === 'SPLIT' ? 2.2 : 1.4);
      if (X && X[nm].length) svg.append('path').datum(X[nm]).attr('d', d3.line().x(v => x(v.c)).y(v => y(v.r))).attr('fill', 'none').attr('stroke', C[nm]).attr('stroke-width', 1).attr('stroke-dasharray', '3 2');
    });
    // Δ running
    const dmap = new Map(R.ONI.map(d => [d.c, d.r])), dd = R.RONI.filter(d => dmap.has(d.c)).map(d => ({ c: d.c, v: d.r - dmap.get(d.c) }));
    svg.append('path').datum(dd).attr('d', d3.line().x(v => x(v.c)).y(v => y(v.v))).attr('fill', 'none').attr('stroke', '#111').attr('stroke-width', 0.9).attr('opacity', 0.7);
    svg.append('g').attr('class', 'ax').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(10).tickFormat(d3.format('d')));
    svg.append('g').attr('class', 'ax').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(8));
    svg.append('text').attr('class', 'lab').attr('x', (m.l + W - m.r) / 2).attr('y', H - 6).attr('text-anchor', 'middle').text(`Centre year of ${state.runW}-yr window (CE)`);
    svg.append('text').attr('class', 'lab').attr('transform', `translate(13,${(H - m.b + m.t) / 2}) rotate(-90)`).attr('text-anchor', 'middle').text('Running r');
    const lg = svg.append('g').attr('transform', `translate(${W - m.r + 12},${m.t + 6})`);
    const items = [['ONI', C.ONI, null], ['RONI', C.RONI, null], ['ΔR (RONI − ONI)', '#111', null], [`p=0.05 (n=${state.runW})`, '#777', '5 3'], [`p=0.05 (N_eff≈${neMed.toFixed(0)})`, '#777', '1 2']];
    if (X) items.push(['ERSST.v6 (instr.)', '#555', '3 2'], ['LMR 5–95%', '#bbb', 'band']);
    items.forEach(([t, col, dash], i) => {
      if (dash === 'band') lg.append('rect').attr('x', 0).attr('y', i * 16).attr('width', 18).attr('height', 9).attr('fill', col);
      else lg.append('line').attr('x1', 0).attr('x2', 18).attr('y1', i * 16 + 5).attr('y2', i * 16 + 5).attr('stroke', col).attr('stroke-width', 2).attr('stroke-dasharray', dash);
      lg.append('text').attr('class', 'sm').attr('x', 23).attr('y', i * 16 + 9).text(t);
    });
    svg.append('text').attr('class', 'sm').attr('x', Math.max(m.l, x(state.y0)) + 4).attr('y', m.t + 11).text('analysis window');
    const note = state.source === 'lmr'
      ? `Shaded bands: 5–95% range of running r across 100 LMRv2.1 ensemble members (reconstruction uncertainty only; add sampling uncertainty via the dotted N_eff threshold). Dashed thin lines: same computation with instrumental ERSST.v6 indices (1850–), an independent check on the reconstruction. LMR skill vs ERSST.v6 (${D.lmr.skill.period}): r = ${D.lmr.skill.r_oni} (ONI), ${D.lmr.skill.r_roni} (RONI).`
      : `Running Pearson r between ${varTitle(state)} and the ${seasonName(state)} index in ${state.runW}-yr centred windows (≥80% coverage), full record overlap; the analysis window is shaded. Dashed: two-sided p=0.05 threshold for n=${state.runW}; dotted: threshold using the median N_eff of the windows. Black: ΔR(t).`;
    d3.select('#note-f').html(note + (state.source === 'cpc' ? ' <b>Switch source to ERSST.v6 or LMRv2.1 to extend before 1950.</b>' : ''));
    return R;
  }

  // ------------------------------------------------------------------ (g) decomposition
  function decomposition(core, st) {
    const P = proxyOf(core, st); if (!P) return null;
    const Y = [], p = [], rel = [], t = [];
    for (let i = 0; i < P.years.length; i++) {
      const yr = P.years[i]; if (yr < st.y0 || yr > st.y1 || P.values[i] === null) continue;
      const a = relAt(yr - st.lag, st), b = t20At(yr - st.lag, st); if (a === null || b === null) continue;
      Y.push(yr); p.push(P.values[i]); rel.push(a); t.push(b);
    }
    if (Y.length < 12) return null;
    let pp = p, rr = rel, tt = t;
    if (st.detrend) { pp = S.detrend(Y, p).y; rr = S.detrend(Y, rel).y; tt = S.detrend(Y, t).y; }
    const ne = Math.min(S.neff(pp, rr).neff, S.neff(pp, tt).neff);
    return { n: Y.length, neff: ne, ...S.regress2(pp, rr, tt, ne), rT: S.pearson(pp, tt), pT: S.pCorr(S.pearson(pp, tt), S.neff(pp, tt).neff) };
  }
  function renderDecomp() {
    const W = widthOf('#fig-g'), list = cores.map(c => ({ c, d: decomposition(c, state) })).filter(o => o.d);
    const rowH = 30, m = { l: 150, r: 16, t: 40, b: 42 }, H = m.t + m.b + Math.max(1, list.length) * rowH;
    const svg = newSvg('#fig-g', W, H);
    if (!list.length) { svg.append('text').attr('class', 'lab').attr('x', 10).attr('y', 30).text('Not enough overlap (≥12 yr) with tropical-mean SST for any core.'); d3.select('#note-g').text(''); return null; }
    const ext = d3.max(list, o => Math.max(Math.abs(o.d.b1) + 2 * o.d.se1, Math.abs(o.d.b2) + 2 * o.d.se2));
    const lim = Math.min(1.5, Math.max(0.6, ext));
    const x = d3.scaleLinear().domain([-lim, lim]).nice().range([m.l, W - m.r]).clamp(true);
    svg.append('g').attr('class', 'gridl').attr('transform', `translate(0,${m.t})`).call(d3.axisTop(x).ticks(6).tickSize(-(H - m.t - m.b)).tickFormat(''));
    svg.append('line').attr('x1', x(0)).attr('x2', x(0)).attr('y1', m.t).attr('y2', H - m.b).attr('stroke', '#555');
    const items = [['β ENSO-relative (Niño3.4 − tropics)', C.RONI, 'b1', 'se1', 'p1', -5], ['β tropical-mean SST (20°S–20°N)', C.T20, 'b2', 'se2', 'p2', 5]];
    items.forEach(([t, col], i) => {
      svg.append('rect').attr('x', m.l + i * (W - m.l) / 2).attr('y', 10).attr('width', 12).attr('height', 10).attr('fill', col);
      svg.append('text').attr('class', 'sm').attr('x', m.l + 16 + i * (W - m.l) / 2).attr('y', 19).text(t);
    });
    list.forEach((o, i) => {
      const y = m.t + (i + 0.5) * rowH, d = o.d, tq = S.tQuantile(0.975, Math.max(d.df, 1));
      if (o.c.id === state.core) svg.append('rect').attr('x', 2).attr('y', y - rowH / 2).attr('width', W - 4).attr('height', rowH).attr('fill', '#eef4f8');
      svg.append('text').attr('class', 'lab').attr('x', 8).attr('y', y - 1).attr('font-weight', o.c.id === state.core ? 700 : 400).text(o.c.id + ' ' + o.c.site);
      svg.append('text').attr('class', 'sm').attr('x', 8).attr('y', y + 10).text(`R² = ${fmt(d.R2)}, n = ${d.n}, N_eff ≈ ${d.neff.toFixed(0)}`);
      items.forEach(([, col, b, se, pk, dy]) => {
        svg.append('line').attr('x1', x(d[b] - tq * d[se])).attr('x2', x(d[b] + tq * d[se])).attr('y1', y + dy).attr('y2', y + dy).attr('stroke', col).attr('stroke-width', 1.5);
        svg.append('circle').attr('cx', x(d[b])).attr('cy', y + dy).attr('r', 4.5).attr('fill', d[pk] < 0.05 ? col : '#fff').attr('stroke', col).attr('stroke-width', 2)
          .on('mousemove', ev => showTip(ev, `${o.c.name}<br>β = ${fmt(d[b])} ± ${fmt(tq * d[se])} (95%), p ${eqP(d[pk])}<br>partial r = ${fmt(b === 'b1' ? d.pr1 : d.pr2)}; predictor r = ${fmt(d.r12)}, VIF = ${fmt(d.vif, 1)}`)).on('mouseleave', hideTip);
      });
    });
    svg.append('g').attr('class', 'ax').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(6));
    svg.append('text').attr('class', 'lab').attr('x', (m.l + W - m.r) / 2).attr('y', H - 8).attr('text-anchor', 'middle').text(`Standardized coefficient (proxy σ per predictor σ), ${varTitle(state)}`);
    const sel = list.find(o => o.c.id === state.core);
    d3.select('#note-g').html(`Multiple regression of each proxy on the two components that ONI mixes: the ENSO-relative signal (Niño-3.4 minus tropical mean, i.e. un-rescaled RONI) and the tropical-mean SST anomaly (both ${seasonName(state)}, ${state.source === 'lmr' ? 'LMRv2.1 ensemble median' : 'ERSST.v6'}; lag ${state.lag} yr). ` +
      `A significant tropical-mean β (filled orange) supports the hypothesis that the ice records basin-wide tropical warmth that RONI removes. Whiskers: 95% CI with N_eff-based degrees of freedom.` +
      (sel ? ` <b>${sel.c.id}:</b> β<sub>rel</sub> = ${fmt(sel.d.b1)} (p ${eqP(sel.d.p1)}), β<sub>trop</sub> = ${fmt(sel.d.b2)} (p ${eqP(sel.d.p2)}).` : ''));
    return { list, sel };
  }

  // ------------------------------------------------------------------ (h) banded age model
  function renderBAM(res) {
    const core = coreById[state.core], P = proxyOf(core, state), W = widthOf('#fig-h'), H = 260, m = { l: 44, r: 14, t: 22, b: 40 };
    const svg = newSvg('#fig-h', W, H);
    if (!P || !res.ONI.s) { svg.append('text').attr('class', 'lab').attr('x', 10).attr('y', 30).text('No overlapping record.'); d3.select('#note-h').text(''); return null; }
    // perturb the full record from the core top downward (youngest → oldest); the index series is
    // fixed (ensemble median for LMR), so only the proxy is re-sampled in each realization
    const order = P.years.map((y, i) => i).sort((a, b) => P.years[b] - P.years[a]);
    const vals = order.map(i => P.values[i]);
    const pos = new Map(order.map((i, j) => [P.years[i], j]));
    const rnd = S.mulberry32(42), NIT = 400, dist = { ONI: [], RONI: [] };
    const base = {};
    IDX.forEach(nm => {
      const A = res[nm].A; if (!A || A.n < 8) return;
      let xi = A.ens ? A.Y.map((_, i) => d3.median(A.x, c => c[i])) : A.x.slice();
      base[nm] = { Y: A.Y, x: xi, j: A.Y.map(y => pos.get(y)) };
    });
    for (let k = 0; k < NIT; k++) {
      const pv = S.bamPerturb(vals, state.theta, rnd);
      IDX.forEach(nm => {
        const b = base[nm]; if (!b) return;
        let p = b.j.map(j => pv[j]); if (p.some(v => v === null)) return;
        if (state.detrend) p = S.detrend(b.Y, p).y;
        dist[nm].push(S.pearson(b.x, p));
      });
    }
    const all = [...dist.ONI, ...dist.RONI, res.ONI.s.r, res.RONI.s?.r ?? 0];
    const x = d3.scaleLinear().domain([Math.min(-0.3, d3.min(all)), Math.max(0.3, d3.max(all))]).nice().range([m.l, W - m.r]);
    const bins = d3.bin().domain(x.domain()).thresholds(x.ticks(30));
    const hs = IDX.map(nm => bins(dist[nm]));
    const y = d3.scaleLinear().domain([0, d3.max(hs.flat(), b => b.length) || 1]).nice().range([H - m.b, m.t]);
    IDX.forEach((nm, i) => {
      svg.selectAll(null).data(hs[i]).join('rect').attr('x', b => x(b.x0) + 0.5).attr('width', b => Math.max(0, x(b.x1) - x(b.x0) - 1)).attr('y', b => y(b.length)).attr('height', b => y(0) - y(b.length))
        .attr('fill', C[nm]).attr('opacity', 0.45);
      const s = res[nm].s; if (!s) return;
      svg.append('line').attr('x1', x(s.r)).attr('x2', x(s.r)).attr('y1', m.t).attr('y2', H - m.b).attr('stroke', C[nm]).attr('stroke-width', 2);
    });
    svg.append('g').attr('class', 'ax').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(8));
    svg.append('g').attr('class', 'ax').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(4));
    svg.append('text').attr('class', 'lab').attr('x', (m.l + W - m.r) / 2).attr('y', H - 6).attr('text-anchor', 'middle').text(`r under perturbed age models (${NIT} realizations)`);
    svg.append('text').attr('class', 'lab').attr('transform', `translate(12,${(H - m.b + m.t) / 2}) rotate(-90)`).attr('text-anchor', 'middle').text('Count');
    const q = nm => dist[nm].length ? [S.quantile(dist[nm], 0.025), d3.median(dist[nm]), S.quantile(dist[nm], 0.975)] : [NaN, NaN, NaN];
    IDX.forEach((nm, i) => { const [a, b, c] = q(nm); svg.append('text').attr('class', 'lab').attr('x', m.l + 6).attr('y', m.t + 4 + i * 14).attr('fill', C[nm]).attr('font-weight', 700).text(`${nm}: obs ${fmt(res[nm].s?.r)}; perturbed median ${fmt(b)} [${fmt(a)}, ${fmt(c)}]`); });
    const depth = P.years[P.years.length - 1] - state.y0;
    d3.select('#note-h').html(`Banded age model (Comboul et al. 2014): each annual layer is independently missed or double-counted with probability ${(state.theta * 100).toFixed(1)}%, and errors accumulate from the core top (${P.years[P.years.length - 1]} CE). At the start of the window (${depth} layers down) the expected |age error| is ≈ ${Math.sqrt(2 * state.theta * depth).toFixed(1)} yr (1σ). Vertical lines: observed r. A distribution collapsing toward 0 means the ENSO signal cannot survive realistic dating error at this depth.`);
    return { ONI: q('ONI'), RONI: q('RONI') };
  }

  // ------------------------------------------------------------------ interpretation & caption
  function lagScan(core) {
    const out = [];
    for (let L = -2; L <= 2; L++) { const s = correlate(align(core, state.mode === 'RONI' ? 'RONI' : 'ONI', state, { lag: L })); out.push({ L, r: s ? s.r : NaN, p: s ? s.p : NaN }); }
    return out;
  }
  function renderInterp(res, comp, dec, bam, rc) {
    const core = coreById[state.core], nm = state.mode === 'SPLIT' ? 'ONI' : state.mode, s = res[nm].s, so = res.ONI.s, sr = res.RONI.s;
    const el = d3.select('#interp');
    if (!s) { el.html('<p>Fewer than 8 overlapping years: no statistics are computed. Widen the window, change the ENSO source (ERSST.v6 extends to 1850; LMRv2.1 to 1 CE) or choose a variable measured in this core.</p>'); return; }
    const P = proxyOf(core, state), strength = a => a < 0.1 ? 'negligible' : a < 0.3 ? 'weak' : a < 0.5 ? 'moderate' : 'strong';
    const ps = [];
    ps.push(`<span class="k">Coupling</span><p>${core.name} ${varTitle(state)} (${P.units}) and the ${seasonName(state)} ${nm} (${SRC_LABEL[state.source]}) share a <b>${strength(Math.abs(s.r))} ${s.r >= 0 ? 'positive' : 'negative'}</b> relationship over ${state.y0}–${state.y1} (r = ${fmt(s.r)}, ρ = ${fmt(s.rho)}, n = ${s.n}). ` +
      (s.p < 0.05 ? `It is significant after accounting for serial correlation (p ${eqP(s.p)}).` : `It is <b>not</b> significant once serial correlation is considered (p ${eqP(s.p)}).`) +
      (s.n - s.neff > 2 ? ` Persistence in the ${Math.abs(s.ar1y) > Math.abs(s.ar1x) ? 'proxy' : 'index'} (lag-1 autocorrelation ${fmt(Math.abs(s.ar1y) > Math.abs(s.ar1x) ? s.ar1y : s.ar1x)}) reduces the effective sample size from ${s.n} to ${s.neff.toFixed(0)}${s.p_naive < 0.05 && s.p >= 0.05 ? ', which is what removes the nominal significance' : ''}.` : ' Both series are close to white noise at lag 1, so N_eff ≈ n.') +
      (varTitle(state) === 'δ¹⁸O' ? ` Physically, ${s.r > 0 ? 'positive r means isotopically enriched (less negative) ice in warm-ENSO years, consistent with reduced upstream rain-out and weaker convection over the Amazon and Altiplano during El Niño' : 'negative r means depleted ice in warm-ENSO years, opposite to the canonical Altiplano response; check the dating and the season'}.` : '') + '</p>');
    if (so && sr && res.diff) {
      const d = res.diff;
      ps.push(`<span class="k">ONI vs RONI</span><p>Under RONI the correlation ${Math.abs(sr.r) < Math.abs(so.r) ? 'weakens' : 'strengthens'} (r ${fmt(so.r)} → ${fmt(sr.r)}, ΔR = ${fmt(d.dr, 3)}). ` +
        (d.p < 0.05 ? `This difference is statistically significant (MRR Z = ${fmt(d.Z)}, p ${eqP(d.p)}).` : `The difference is not significant (p ${eqP(d.p)}): the two indices correlate at ${fmt(d.r12)} here, so a larger record or a stronger proxy would be needed to separate them.`) +
        (Math.abs(sr.r) < Math.abs(so.r) ? ' The drop is what the hypothesis predicts if part of the proxy–ONI covariance is carried by tropical-mean temperature, which RONI subtracts.' : ' A rise is the opposite of the hypothesis: this proxy follows the zonal Pacific SST gradient more than basin-wide warmth, and removing the tropical mean cleans the ENSO signal.') + '</p>');
    }
    if (dec && dec.sel) {
      const q = dec.sel.d;
      ps.push(`<span class="k">Decomposition (panel g)</span><p>With both components in one regression, the ENSO-relative term has β = ${fmt(q.b1)} (p ${eqP(q.p1)}) and the tropical-mean term β = ${fmt(q.b2)} (p ${eqP(q.p2)}); together they explain R² = ${fmt(q.R2)} of the variance. ` +
        (q.p2 < 0.05 && q.p1 >= 0.05 ? `The tropical-mean term carries the signal, consistent with the hypothesis${q.R2 < 0.15 ? ', although the explained variance is small' : ''}.` :
          q.p2 < 0.05 && q.p1 < 0.05 ? 'Both terms matter: the ice records ENSO dynamics <i>and</i> tropical-mean warmth, so ONI and RONI capture different fractions of it.' :
            q.p1 < 0.05 ? 'Only the ENSO-relative term is significant: here the coupling is dynamical rather than a tropical-mean thermodynamic signal.' :
              'Neither term is individually significant at this sample size; the predictors are correlated (r = ' + fmt(q.r12) + ', VIF = ' + fmt(q.vif, 1) + '), which inflates the standard errors.') + '</p>');
    }
    const other = d3.select('#ctl-detrend').empty() ? null : correlate(align(core, nm, { ...state, detrend: !state.detrend }));
    if (other) ps.push(`<span class="k">Detrending</span><p>${state.detrend ? 'With' : 'Without'} linear detrending r = ${fmt(s.r)}; ${state.detrend ? 'without' : 'with'} it r = ${fmt(other.r)}. ` +
      (Math.abs(other.r - s.r) < 0.05 ? 'Shared trends contribute little to the coupling.' : state.detrend ? 'The difference is the part of the covariance carried by common long-term trends (e.g. 20th-century warming and δ¹⁸O enrichment) rather than interannual ENSO variability.' : 'Some of the apparent coupling is carried by common long-term trends; detrending isolates interannual variability. Note that RONI already removes the tropical-mean trend by construction, whereas ONI’s sliding climatology removes it only in 30-yr steps.') + '</p>');
    const ls = lagScan(core), best = ls.filter(d => Number.isFinite(d.r)).sort((a, b) => Math.abs(b.r) - Math.abs(a.r))[0];
    if (best) ps.push(`<span class="k">Lag and season</span><p>Scanning lags −2…+2 yr for ${nm}: ${ls.map(d => `${d.L > 0 ? '+' : d.L < 0 ? '−' : ''}${Math.abs(d.L)}: ${fmt(d.r)}`).join(', ')}. The strongest |r| is at ${best.L > 0 ? '+' : ''}${best.L} yr${best.L === 0 ? ', consistent with same-wet-season deposition' : best.L > 0 ? ' (ENSO leading the ice), which may reflect the Aug–Jul vs calendar-year convention or ±1-yr layer-counting offsets' : ' (ice leading ENSO), physically implausible as a cause, so it most likely signals a dating offset'}. ` +
      `The ${seasonName(state)} aggregation ${state.season === 'djf' ? 'isolates the ENSO peak that coincides with the Andean wet season (≈70–80% of annual snowfall falls in Nov–Mar), maximizing the physical match with the annual layer' : state.season === 'ann' ? 'averages the index over the calendar year, blending the decay of one event with the onset of the next' : state.season === 'hydro' ? 'matches the thermal (Aug–Jul) year used to date the cores' : 'uses a custom window; compare with DJF to test seasonal specificity'}.</p>`);
    if (comp && comp[nm]) ps.push(`<span class="k">Composites</span><p>El Niño years average ${fmt(comp[nm].EN.mean)} σ and La Niña years ${fmt(comp[nm].LN.mean)} σ (n = ${comp[nm].EN.n} and ${comp[nm].LN.n}); the difference of ${fmt(comp[nm].diff)} σ has p ${eqP(comp[nm].p_welch)} (Welch)${Number.isFinite(comp[nm].p_perm) ? ` and ${fmtP(comp[nm].p_perm)} (permutation)` : ''}. ${rc.n} year${rc.n === 1 ? '' : 's'} change class between ONI and RONI in this window.</p>`);
    if (bam) ps.push(`<span class="k">Dating robustness</span><p>With ${(state.theta * 100).toFixed(1)}% layer miscounts per year, ${nm} r ranges ${fmt(bam[nm][0])} to ${fmt(bam[nm][2])} (95%) around a median of ${fmt(bam[nm][1])}. ${(() => { const q = Math.abs(bam[nm][1]) / Math.max(1e-6, Math.abs(s.r)); return q >= 0.8 ? 'The coupling is robust to this level of dating error.' : q >= 0.5 ? `Dating error of this size erodes about ${Math.round(100 * (1 - q))}% of the correlation; the sign survives but the magnitude is uncertain.` : 'The coupling largely disappears under realistic age uncertainty — interpret this section cautiously, and treat the observed r as an upper-bound estimate of a dating-sensitive signal.'; })()}</p>`);
    if (state.source === 'lmr') ps.push(`<span class="k">Reconstruction uncertainty</span><p>Correlations use 100 LMRv2.1 members; reported CIs combine reconstruction spread and sampling error (${(100 * s.fracSig).toFixed(0)}% of members give p<0.05). LMR is an annual-mean product, so DJF-specific analysis is not possible before 1850; its interannual skill is limited (r ≈ ${D.lmr.skill.r_oni} vs instrumental Niño-3.4, ${D.lmr.skill.period}).${core.in_lmr ? ' <b>Quelccaya was assimilated into LMR (PAGES2k), so these correlations are partly circular.</b>' : ''}</p>`);
    el.html(ps.join(''));
  }
  function renderCaption(all, res) {
    const nm = 'ONI', core = coreById[state.core], list = cores.filter(c => all[c.id]?.diff), s = res;
    const med = list.length ? d3.median(list, c => all[c.id].diff.dr) : NaN;
    const nsig = list.filter(c => all[c.id].diff.p < 0.05).length;
    const src = state.source === 'lmr' ? 'LMRv2.1 multi-proxy reconstruction, 100 members' : state.source === 'ersst' ? 'ERSST.v6-derived indices (1850–present)' : 'NOAA CPC indices';
    void nm;
    const cap = `<b>Figure 1 | ENSO coupling of Andean ice-core ${varTitle(state)} under ONI and RONI.</b> ` +
      `<b>a</b>, NOAA NCEI Andean ice-core sites (${cores.length} annually resolved records, ${new Set(cores.map(c => c.site)).size} sites) coloured by Pearson <i>r</i> with the ${seasonName(state)} index (${src}), ${state.y0}–${state.y1} CE, lag ${state.lag} yr${state.detrend ? ', detrended' : ''}; inset, site elevation versus latitude. ` +
      `<b>b</b>, ${core.name} ${varTitle(state)} (z-score) with ONI and RONI; shading marks El Niño/La Niña years (±${state.thr} ${state.source === 'lmr' ? 'σ' : '°C'}). ` +
      `<b>c</b>, Correlations with significance from effective sample size (r = ${fmt(s.ONI.s?.r)} and ${fmt(s.RONI.s?.r)}). ` +
      `<b>d</b>, El Niño–La Niña composites (95% CI). ` +
      `<b>e</b>, Change in <i>r</i> from ONI to RONI (median ΔR = ${fmt(med)}; ${nsig} of ${list.length} significant) and reclassified events. ` +
      `<b>f</b>, ${state.runW}-yr running correlation${state.source === 'lmr' ? ', 5–95% ensemble range' : ''}. ` +
      `<b>g</b>, Regression on ENSO-relative and tropical-mean SST. ` +
      `<b>h</b>, Sensitivity to ${(state.theta * 100).toFixed(1)}% yr⁻¹ layer-counting error. ` +
      `RONI subtracts the 20°S–20°N mean SST anomaly from Niño-3.4; p-values use autocorrelation-adjusted effective sample sizes, and ΔR is tested for dependent correlations${state.source === 'lmr' ? '; uncertainties combine reconstruction ensemble spread and sampling error' : ''}.`;
    d3.select('#caption').html(cap);
    const words = d3.select('#caption').node().innerText.split(/\s+/).filter(Boolean).length;
    d3.select('#cap-count').text(`${words} words · updates with the settings`);
  }

  function renderMethods() {
    const e = D.ersst, sk = D.lmr.skill;
    d3.select('#methods').html(`
      <p><b>ENSO indices.</b> Official NOAA CPC ONI and RONI (<code>oni.ascii.txt</code>, <code>RONI.ascii.txt</code>; ERSST.v6; last value ${D.meta.cpc_last}). Monthly series are assigned to the centre month of each 3-month season (DJF → January). To extend before 1950 we recomputed both indices from the ERSST.v6 grid with CPC’s algorithm: Niño-3.4 (5°N–5°S, 170°–120°W) and the cos(lat)-weighted 20°S–20°N ocean mean; ONI with centred 30-yr base periods updated every 5 yr; RONI = 3-month mean of (Niño-3.4 − tropical mean) anomalies (1991–2020 base) × k, k = ${e.k} matching the Niño-3.4 variance. Verification against CPC (1950–2026, n = ${e.val_oni.n} months): ONI r = ${e.val_oni.r}, RMSE = ${e.val_oni.rmse} °C; RONI r = ${e.val_roni.r}, RMSE = ${e.val_roni.rmse} °C.</p>
      <p><b>Pre-instrumental ENSO.</b> LMRv2.1 (Tardif et al. 2019) offline data assimilation of PAGES2k proxies (corals, tree rings, ice cores, etc.), annual means 1–2000 CE. Reconstructed ONI = member Niño-3.4 with CPC’s centred 30-yr base periods (as for ONI); reconstructed RONI = member Niño-3.4 minus the 20°S–20°N mean SST of the same Monte-Carlo iteration, plus a random draw N(0, σ<sub>T</sub> = ${sk.sigma_T} °C) for the unarchived within-iteration spread of the tropical mean, then variance-rescaled. σ<sub>T</sub> is estimated from the Niño-3.4 box, where box-mean and grid-point spreads are both known (${sk.n_box} independent error regions), scaled by ocean area and doubled for conservatism; the spread between iterations adds ${sk.mc_spread_T} °C. 100 of the 2000 members (20 iterations × 100) are sampled. Skill vs ERSST.v6 annual means (${sk.period}): ONI r = ${sk.r_oni}, RONI r = ${sk.r_roni}, tropical mean r = ${sk.r_t20}; mean ensemble 1σ = ${sk.spread_oni} °C (ONI), ${sk.spread_roni} °C (RONI). Uncertainty is propagated by computing every statistic for each member and pooling with sampling error (Fisher z with N<sub>eff</sub>). Because RONI members are noisier, per-member correlations with RONI are slightly more attenuated than with ONI; compare against the ERSST.v6 overlap (panel f) before interpreting small ΔR. <span class="warn">Quelccaya δ¹⁸O is part of PAGES2k and may be assimilated by LMR (circularity).</span></p>
      <p><b>Ice cores.</b> Only South American Andean cores archived by NOAA NCEI Paleoclimatology (ice-core data type) with annual resolution overlapping the index record: Quelccaya Summit Dome 2003 core (δ¹⁸O, accumulation, dust, 9 ions; 226–2009 CE), Quelccaya 2018 firn core (2004–2018), Huascarán Col and Summit composites (δ¹⁸O, d-excess, accumulation, dust; 1960–2019), Huascarán 1993 core 2 (δ¹⁸O, particles, NO₃⁻; 1894–1993), Illimani 1999 deep core (NH₄⁺; 1800–1998) and Illimani 2017 firn core (2000–2017). Sub-annual firn samples are averaged into Aug–Jul years. Sajama is archived only at 5-m/100-yr resolution; no Coropuna, Hualcán or Chimborazo records exist in the NCEI ice-core archive (searched Sep 2026); these sites are shown but not analysed.</p>
      <p><b>Statistics.</b> Proxy year Y is paired with ENSO year Y − lag (lag &gt; 0: ENSO leads). Optional linear detrending of each series over the analysis window. Pearson r and Spearman ρ; two-sided p from Student-t with N<sub>eff</sub> = n(1 − r₁r₂)/(1 + r₁r₂) (Bretherton et al. 1999), capped at n; Fisher-z 95% CI with N<sub>eff</sub>. ΔR tested with Meng–Rosenthal–Rubin (1992). Composites: ENSO year classified by the CPC episode rule (≥5 consecutive overlapping seasons beyond ±threshold) for DJF; Welch t-test and 4000-shuffle permutation test for El Niño − La Niña. All routines were cross-checked against SciPy/NumPy (24 quantities, max relative error &lt;10⁻⁶; <code>tests/</code>).</p>
      <p><b>Reproduce.</b> <code>python3 pipeline/build_dataset.py</code> regenerates <code>app/data/dataset.js</code> from <code>data/raw/</code> (download script: <code>pipeline/fetch_data.sh</code>). Open <code>app/index.html</code> in any modern browser; no server needed.</p>`);
  }

  // ------------------------------------------------------------------ export
  function exportPanel(id, kind) {
    const svg = document.querySelector('#' + id + ' svg'); if (!svg) return;
    const clone = svg.cloneNode(true);
    const W = +svg.getAttribute('width'), H = +svg.getAttribute('height');
    clone.setAttribute('width', W); clone.setAttribute('height', H);
    const src = '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(clone);
    const name = `fig_${id.replace('fig-', '')}_${state.core}_${state.variable}${state.variable === 'ions' ? state.ion : ''}_${state.mode}_${state.source}`;
    if (kind === 'svg') return download(new Blob([src], { type: 'image/svg+xml' }), name + '.svg');
    const img = new Image(), scale = 3;
    img.onload = () => {
      const cv = document.createElement('canvas'); cv.width = W * scale; cv.height = H * scale;
      const ctx = cv.getContext('2d'); ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, cv.width, cv.height); ctx.drawImage(img, 0, 0, cv.width, cv.height);
      cv.toBlob(b => download(b, name + '.png'), 'image/png');
    };
    img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(src);
  }
  function exportCSV(id) {
    let rows;
    if (id === 'fig-e') rows = renderWhatChanged.csv;
    else {
      const core = coreById[state.core], A = { ONI: align(core, 'ONI', state), RONI: align(core, 'RONI', state) };
      rows = [['year', 'proxy_raw', 'proxy_z', 'ONI', 'RONI', 'class_ONI', 'class_RONI']];
      if (A.ONI && A.ONI.n >= 8) A.ONI.Y.forEach((y, i) => {
        const j = A.RONI.Y.indexOf(y), cO = yearClass('ONI', y - state.lag, state), cR = yearClass('RONI', y - state.lag, state);
        const v = (B, k) => B.ens ? d3.median(B.x, c => c[k]) : B.x[k];
        rows.push([y, A.ONI.pRaw[i], A.ONI.p[i], v(A.ONI, i), j >= 0 ? v(A.RONI, j) : '', Array.isArray(cO) ? '' : cO, Array.isArray(cR) ? '' : cR]);
      });
    }
    const txt = rows.map(r => r.map(v => typeof v === 'number' ? +v.toFixed(5) : (v ?? '')).join(',')).join('\n');
    download(new Blob([txt], { type: 'text/csv' }), `${id}_${state.core}_${state.source}.csv`);
  }
  function download(blob, name) { const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name; document.body.appendChild(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500); }

  // ------------------------------------------------------------------ main update
  function update() {
    syncControls();
    const all = {};
    cores.forEach(c => { if (proxyOf(c, state)) all[c.id] = analyze(c, state); });
    const res = all[state.core] || analyze(coreById[state.core], state);
    renderMap(all);
    renderTimeseries(res);
    renderStats(res);
    const comp = renderComposites(res);
    const rc = renderWhatChanged(all, res);
    renderRunning();
    const dec = renderDecomp();
    const bam = renderBAM(res);
    renderInterp(res, comp, dec, bam, rc);
    renderCaption(all, res);
    window.__APP = { state, all, res, comp, dec, bam, rc };   // for inspection/testing
  }

  d3.select('#built').text(D.meta.built);
  initControls(); renderMethods(); fitWindow(true); update();
  window.__setState = (o) => { Object.assign(state, o); update(); };
})();
