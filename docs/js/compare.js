/* "Your data" tab: load a user time series and compare it with the historical records in window.DATA.
 * Everything runs in the browser; the file never leaves the viewer's machine.
 * Conventions follow the Analysis tab: user year Y is paired with ENSO year Y − lag (DJF = Jan Y season),
 * p-values use N_eff (Bretherton et al. 1999), ΔR uses Meng–Rosenthal–Rubin (1992),
 * the screen controls the false-discovery rate (Benjamini–Hochberg 1995). */
(function () {
  'use strict';
  const D = window.DATA, S = window.Stats;
  const $ = id => document.getElementById(id);
  const MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
  const C = { user: '#1b1b1b', ONI: '#CC79A7', RONI: '#009E73', EN: '#D55E00', LN: '#0072B2', T: '#E69F00', grid: '#ececec' };
  const FIG_CSS = `text{font-family:"IBM Plex Sans","Helvetica Neue",Arial,sans-serif}text:not([fill]){fill:#1b1b1b}
    .ax text{font-size:11px}.ax path,.ax line{stroke:#666;shape-rendering:crispEdges}.gridl line{stroke:#ececec}.gridl path{display:none}
    .lab{font-size:11.5px}.sm{font-size:10px}.sm:not([fill]){fill:#555}`;
  const st = {
    raw: null, rows: [], cols: [], yearCol: null, monthCol: null, valCol: null, bp: false, missing: '-999',
    agg: 'hydro', label: 'Your series', units: '', user: null,
    source: 'cpc', season: 'djf', cm0: 12, cm1: 3, lag: 0, detrend: false, y0: null, y1: null, runW: 21, thr: 0.5, target: null
  };

  // ------------------------------------------------------------------ parsing
  const MISSING = new Set(['', 'na', 'nan', 'null', 'none', '-', '--', '.']);
  function parseText(text) {
    const lines = text.replace(/\r/g, '').split('\n').filter(l => l.trim() && !/^\s*#/.test(l));
    if (!lines.length) throw new Error('No data rows found (lines starting with # are treated as comments).');
    const cand = ['\t', ',', ';', /\s+/];
    let best = null;
    for (const d of cand) {
      const n = lines.slice(0, 20).map(l => l.trim().split(d).length);
      const mode = n.sort((a, b) => a - b)[Math.floor(n.length / 2)];
      if (mode > 1 && (!best || mode > best.n)) best = { d, n: mode };
    }
    if (!best) throw new Error('Could not detect columns. Use comma-, tab- or space-separated values with at least a year and a value column.');
    const split = l => l.trim().split(best.d).map(s => s.trim().replace(/^"|"$/g, ''));
    let first = split(lines[0]);
    const numeric = s => s !== '' && !isNaN(+s);
    const hasHeader = first.some(s => !numeric(s) && !/^\d{4}-\d{2}/.test(s));
    const cols = hasHeader ? first.map((s, i) => s || `col${i + 1}`) : first.map((_, i) => `col${i + 1}`);
    const rows = lines.slice(hasHeader ? 1 : 0).map(split).filter(r => r.length >= 2);
    return { cols, rows };
  }
  function guessColumns() {
    const { cols, rows } = st;
    const colVals = i => rows.map(r => r[i]);
    const isYearLike = i => { const v = colVals(i).filter(s => s !== undefined && s !== ''); return v.length && v.every(s => /^\d{4}-\d{2}/.test(s) || (!isNaN(+s) && +s > -10000 && +s < 60000)); };
    st.yearCol = cols.findIndex((c, i) => /year|age|yr|date|time|ce\b|bp\b/i.test(c) && isYearLike(i));
    if (st.yearCol < 0) st.yearCol = cols.findIndex((c, i) => isYearLike(i));
    if (st.yearCol < 0) st.yearCol = 0;
    st.monthCol = cols.findIndex((c, i) => i !== st.yearCol && /^(month|mon|mth|mm)$/i.test(c));
    const numericCol = i => i !== st.yearCol && i !== st.monthCol && colVals(i).some(s => !isNaN(+s) && s !== '');
    const bookkeeping = c => /^(sample|samp|id|index|idx|no|n|num|count|depth|top_?depth|bottom_?depth|bot_?depth|mid_?depth|.*depth.*|top_year|bottom_year|age_?err.*|.*_err|.*sigma.*)$/i.test(c);
    st.valCol = cols.findIndex((c, i) => numericCol(i) && !bookkeeping(c));
    if (st.valCol < 0) st.valCol = cols.findIndex((c, i) => numericCol(i));
    st.bp = /bp/i.test(cols[st.yearCol] || '');
    st.label = cols[st.valCol] || 'Your series';
  }

  // → Map(year → mean value), plus a description of the time resolution
  function buildSeries() {
    const miss = new Set(st.missing.split(',').map(s => s.trim()).filter(Boolean));
    const pts = [];
    for (const r of st.rows) {
      let ys = r[st.yearCol], vs = r[st.valCol];
      if (ys === undefined || vs === undefined) continue;
      if (MISSING.has(String(vs).toLowerCase()) || miss.has(String(vs).trim())) continue;
      const v = +vs; if (!Number.isFinite(v)) continue;
      let t, month = null;
      const iso = /^(\d{4})-(\d{2})/.exec(ys);
      if (iso) { t = +iso[1]; month = +iso[2]; }
      else { t = +ys; if (!Number.isFinite(t)) continue; if (st.bp) t = 1950 - t; }
      if (st.monthCol >= 0 && r[st.monthCol] !== undefined && r[st.monthCol] !== '') month = +r[st.monthCol];
      pts.push({ t, month, v });
    }
    if (pts.length < 5) throw new Error('Fewer than 5 valid values after removing missing entries.');
    const monthly = pts.some(p => p.month), fractional = !monthly && pts.some(p => Math.abs(p.t - Math.round(p.t)) > 1e-6);
    const buckets = new Map(), add = (y, v) => { if (!buckets.has(y)) buckets.set(y, []); buckets.get(y).push(v); };
    let res = 'annual';
    for (const p of pts) {
      if (monthly) {
        res = 'monthly';
        add(st.agg === 'hydro' && p.month >= 8 ? Math.floor(p.t) + 1 : Math.floor(p.t), p.v);
      } else if (fractional) {
        res = 'sub-annual (decimal years)';
        add(st.agg === 'hydro' ? Math.floor(p.t - 7 / 12 + 1) : Math.floor(p.t), p.v);
      } else add(Math.round(p.t), p.v);
    }
    // sub-annual data: drop incomplete years (fewer than 60% of the typical number of samples per year)
    let dropped = [];
    if (res !== 'annual') {
      const med = d3.median([...buckets.values()], b => b.length);
      for (const [y, b] of [...buckets]) if (b.length < 0.6 * med) { buckets.delete(y); dropped.push(y); }
    }
    const years = [...buckets.keys()].sort((a, b) => a - b);
    const map = new Map(years.map(y => [y, S.mean(buckets.get(y))]));
    return { map, years, res, nRaw: pts.length, perYear: d3.median([...buckets.values()], b => b.length), dropped };
  }

  // ------------------------------------------------------------------ historical series
  const lmrMed = {};
  ['oni', 'roni'].forEach(k => {
    const M = D.lmr[k], n = M[0].length; lmrMed[k] = new Float64Array(n);
    for (let t = 0; t < n; t++) { const col = M.map(a => a[t] / 100).sort((a, b) => a - b); lmrMed[k][t] = col[Math.floor(col.length / 2)]; }
  });
  const lmrSig = { oni: S.sd(Array.from(lmrMed.oni)), roni: S.sd(Array.from(lmrMed.roni)) };
  function monthsFor(Y) {
    const s = st.season; if (s === 'djf') return [[Y, 1]];
    if (s === 'ann') return MON.map((_, i) => [Y, i + 1]);
    const [m0, m1] = s === 'hydro' ? [8, 7] : [st.cm0, st.cm1], out = [];
    if (m0 <= m1) for (let m = m0; m <= m1; m++) out.push([Y, m]);
    else { for (let m = m0; m <= 12; m++) out.push([Y - 1, m]); for (let m = 1; m <= m1; m++) out.push([Y, m]); }
    return out;
  }
  function seasonal(block, key, Y) {
    let s = 0, n = 0;
    for (const [y, m] of monthsFor(Y)) { const v = block[key][(y - block.start) * 12 + m - 1]; if (v === null || v === undefined) return null; s += v; n++; }
    return s / n;
  }
  const seasonName = () => ({ djf: 'DJF', ann: 'Jan–Dec', hydro: 'Aug–Jul' })[st.season] || `${MON[st.cm0 - 1]}–${MON[st.cm1 - 1]}`;
  // ENSO targets for the current source
  function ensoGet(name, Y) {             // name: oni | roni | t20
    if (st.source === 'lmr') {
      const arr = name === 't20' ? D.lmr.t20 : lmrMed[name]; const v = arr[Y - D.lmr.start]; return v === undefined ? null : v;
    }
    const block = (st.source === 'cpc' && name !== 't20') ? D.cpc : D.ersst;
    return seasonal(block, name, Y);
  }
  function lmrMembers(name, Y) { const i = Y - D.lmr.start; if (i < 0 || i >= D.lmr[name][0].length) return null; return D.lmr[name].map(a => a[i] / 100); }
  const SRC_LABEL = { cpc: 'NOAA CPC (1950–)', ersst: 'ERSST.v6, CPC algorithm (1850–)', lmr: 'LMRv2.1 reconstruction (1–2000)' };
  function targets() {
    const out = [];
    const src = SRC_LABEL[st.source], sea = st.source === 'lmr' ? 'annual' : seasonName();
    out.push({ id: 'enso:oni', group: 'ENSO indices', label: `ONI · ${sea} · ${src}`, short: 'ONI', units: '°C', get: Y => ensoGet('oni', Y) });
    out.push({ id: 'enso:roni', group: 'ENSO indices', label: `RONI · ${sea} · ${src}`, short: 'RONI', units: '°C', get: Y => ensoGet('roni', Y) });
    out.push({ id: 'enso:t20', group: 'ENSO indices', label: `Tropical-mean SST 20°S–20°N · ${sea} · ${st.source === 'lmr' ? 'LMRv2.1' : 'ERSST.v6'}`, short: 'Trop. SST', units: '°C', get: Y => ensoGet('t20', Y) });
    D.cores.forEach(c => {
      const add = (key, v) => { const m = new Map(v.years.map((y, i) => [y, v.values[i]])); out.push({ id: `core:${c.id}:${key}`, group: 'Andean ice cores (NOAA NCEI)', label: `${c.name} — ${v.label}`, short: `${c.id} ${v.label}`, units: v.units, get: Y => (m.has(Y) ? m.get(Y) : null) }); };
      Object.entries(c.vars).forEach(([k, v]) => { if (k === 'ions') Object.entries(v).forEach(([ik, iv]) => add('ions:' + ik, iv)); else add(k, v); });
    });
    return out;
  }
  function epiFlag(name, Y) {            // ENSO class of ENSO year Y
    if (st.source === 'lmr') { const v = lmrMed[name][Y - D.lmr.start]; if (v === undefined) return null; return v >= st.thr * lmrSig[name] ? 1 : v <= -st.thr * lmrSig[name] ? -1 : 0; }
    const block = st.source === 'cpc' ? D.cpc : D.ersst, a = block[name], i0 = (Y - block.start) * 12;
    if (i0 < 0 || i0 >= a.length || a[i0] === null) return null;
    for (const sg of [1, -1]) {
      if (sg * a[i0] < st.thr - 1e-9) continue;
      let i = i0, j = i0;
      while (i - 1 >= 0 && a[i - 1] !== null && sg * a[i - 1] >= st.thr - 1e-9) i--;
      while (j + 1 < a.length && a[j + 1] !== null && sg * a[j + 1] >= st.thr - 1e-9) j++;
      if (j - i + 1 >= 5) return sg;
    }
    return 0;
  }

  // ------------------------------------------------------------------ alignment & statistics
  function align(getT, lag = st.lag, y0 = st.y0, y1 = st.y1) {
    const Y = [], u = [], t = [];
    for (const y of st.user.years) {
      if (y < y0 || y > y1) continue;
      const tv = getT(y - lag); if (tv === null || tv === undefined || !Number.isFinite(tv)) continue;
      Y.push(y); u.push(st.user.map.get(y)); t.push(tv);
    }
    if (Y.length < 8) return { n: Y.length, Y };
    let uu = u, tt = t;
    if (st.detrend) { uu = S.detrend(Y, u).y; tt = S.detrend(Y, t).y; }
    return { n: Y.length, Y, u: uu, t: tt, uRaw: u, tRaw: t };
  }
  function ensembleStats(name, A) {      // LMR: correlation for each of the 100 members + sampling error
    const rs = [], ne = [];
    for (let k = 0; k < D.lmr[name].length; k++) {
      let tt = A.Y.map(y => D.lmr[name][k][y - st.lag - D.lmr.start] / 100);
      if (st.detrend) tt = S.detrend(A.Y, tt).y;
      const s = S.corrStats(tt, A.u); rs.push(s.r); ne.push(s.neff);
    }
    const rnd = S.mulberry32(3), pooled = [];
    rs.forEach((r, k) => { for (let j = 0; j < 20; j++) { const g = Math.sqrt(-2 * Math.log(rnd() + 1e-12)) * Math.cos(2 * Math.PI * rnd()); pooled.push(Math.tanh(Math.atanh(r) + g / Math.sqrt(Math.max(ne[k] - 3, 1)))); } });
    return { lo: S.quantile(pooled, 0.025), hi: S.quantile(pooled, 0.975), fracSig: rs.filter((r, k) => S.pCorr(r, ne[k]) < 0.05).length / rs.length };
  }
  function bh(ps) {                       // Benjamini–Hochberg q-values
    const idx = ps.map((p, i) => [p, i]).sort((a, b) => a[0] - b[0]), m = ps.length, q = new Array(m);
    let prev = 1;
    for (let k = m - 1; k >= 0; k--) { const [p, i] = idx[k]; prev = Math.min(prev, p * m / (k + 1)); q[i] = prev; }
    return q;
  }

  // ------------------------------------------------------------------ formatting / svg helpers
  const f = (v, d = 2) => (v == null || !Number.isFinite(v)) ? '—' : (v < 0 ? '−' : '') + Math.abs(v).toFixed(d);
  const fp = p => !Number.isFinite(p) ? '—' : p < 0.001 ? '< 0.001' : p.toFixed(3);
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]);
  function svgIn(id, W, H) {
    const el = d3.select('#' + id); el.selectAll('*').remove();
    const s = el.append('svg').attr('viewBox', `0 0 ${W} ${H}`).attr('width', W).attr('height', H).attr('xmlns', 'http://www.w3.org/2000/svg');
    s.append('style').text(FIG_CSS); s.append('rect').attr('width', W).attr('height', H).attr('fill', '#fff'); return s;
  }
  const wOf = id => Math.max(320, $(id).clientWidth);
  const tip = d3.select('#cmp-tip');
  const showTip = (ev, h) => tip.html(h).style('opacity', 1).style('left', Math.min(ev.clientX + 14, innerWidth - 280) + 'px').style('top', ev.clientY + 12 + 'px');
  const hideTip = () => tip.style('opacity', 0);

  // ------------------------------------------------------------------ rendering
  function renderPreview() {
    const box = $('cmp-preview');
    if (!st.rows.length) { box.innerHTML = ''; return; }
    const opt = (i, sel) => st.cols.map((c, k) => `<option value="${k}" ${k === sel ? 'selected' : ''}>${esc(c)}</option>`).join('');
    $('cmp-ycol').innerHTML = opt(0, st.yearCol);
    $('cmp-mcol').innerHTML = '<option value="-1">— none —</option>' + opt(0, st.monthCol);
    $('cmp-vcol').innerHTML = opt(0, st.valCol);
    $('cmp-bp').checked = st.bp;
    const head = st.cols.map((c, k) => `<th class="${k === st.yearCol ? 'y' : k === st.valCol ? 'v' : ''}">${esc(c)}</th>`).join('');
    const body = st.rows.slice(0, 6).map(r => '<tr>' + st.cols.map((_, k) => `<td>${esc(r[k] ?? '')}</td>`).join('') + '</tr>').join('');
    box.innerHTML = `<div class="tbl-wrap"><table class="prev"><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table></div>
      <div class="small">${st.rows.length} rows · ${st.cols.length} columns · <span class="k-y">year column</span> · <span class="k-v">value column</span></div>`;
  }

  function analyse() {
    const msg = $('cmp-msg');
    if (!st.rows.length) { $('cmp-results').hidden = true; return; }
    try { st.user = buildSeries(); } catch (e) { msg.textContent = e.message; msg.className = 'msg err'; $('cmp-results').hidden = true; return; }
    const U = st.user, yr = [U.years[0], U.years[U.years.length - 1]];
    msg.className = 'msg ok';
    msg.innerHTML = `<b>${esc(st.label)}</b>: ${U.years.length} years (${yr[0]}–${yr[1]} CE) from ${U.nRaw} values · resolution: ${U.res}` +
      (U.res !== 'annual' ? ` → averaged into ${st.agg === 'hydro' ? 'Aug–Jul (hydrological) years, labelled by the year of the wet-season January' : 'calendar years'} (≈${Math.round(U.perYear)} samples per year)` : '') + '.' +
      (U.dropped.length ? ` Dropped ${U.dropped.length} incomplete year${U.dropped.length > 1 ? 's' : ''} (${U.dropped.slice(0, 4).join(', ')}${U.dropped.length > 4 ? '…' : ''}).` : '');
    // default window = overlap with the chosen ENSO source
    const srcRange = st.source === 'cpc' ? [1951, 2026] : st.source === 'ersst' ? [1851, 2026] : [1, 2000];
    if (st.y0 === null || st._auto) { st.y0 = Math.max(yr[0], srcRange[0]); st.y1 = Math.min(yr[1], srcRange[1]); st._auto = true; }
    $('cmp-y0').value = st.y0; $('cmp-y1').value = st.y1;
    $('cmp-results').hidden = false;
    const T = targets();
    if (!st.target || !T.find(t => t.id === st.target)) st.target = 'enso:oni';
    const tsel = $('cmp-target');
    tsel.innerHTML = [...d3.group(T, t => t.group)].map(([g, list]) => `<optgroup label="${esc(g)}">${list.map(t => `<option value="${t.id}" ${t.id === st.target ? 'selected' : ''}>${esc(t.label)}</option>`).join('')}</optgroup>`).join('');
    const enso = renderEnso(T);
    renderLag(T); renderRunning(T); renderComposite();
    const scr = renderScreen(T);
    renderTarget(T);
    renderInterp(enso, scr);
  }

  // A. ENSO coupling: ONI vs RONI -------------------------------------------------------------
  function renderEnso(T) {
    const out = {};
    ['oni', 'roni'].forEach(k => { const A = align(T.find(t => t.id === 'enso:' + k).get); out[k] = { A, s: A.n >= 8 ? S.corrStats(A.t, A.u) : null }; if (out[k].s && st.source === 'lmr') out[k].e = ensembleStats(k, A); });
    const o = out.oni.s, r = out.roni.s;
    let diff = null;
    if (o && r) { const r12 = S.pearson(out.oni.A.t, out.roni.A.t); diff = { dr: r.r - o.r, r12, ...S.compareDepCorr(r.r, o.r, r12, Math.min(o.neff, r.neff)) }; }
    out.diff = diff;
    const row = (lab, fn) => `<tr><td>${lab}</td><td>${o ? fn(o, out.oni.e) : '—'}</td><td>${r ? fn(r, out.roni.e) : '—'}</td></tr>`;
    $('cmp-enso-table').innerHTML = !o ? `<p class="small">Fewer than 8 overlapping years with ${SRC_LABEL[st.source]} in ${st.y0}–${st.y1}. Try another ENSO source or a wider window.</p>` :
      `<table class="st"><thead><tr><th>${esc(st.label)} vs ${st.source === 'lmr' ? 'annual' : seasonName()} index</th><th style="color:${C.ONI}">ONI</th><th style="color:${C.RONI}">RONI</th></tr></thead><tbody>
      ${row('n (years)', s => s.n)}${row('Pearson r', s => `<b>${f(s.r)}</b>`)}${row('95% CI', (s, e) => e ? `[${f(e.lo)}, ${f(e.hi)}]*` : `[${f(s.ci[0])}, ${f(s.ci[1])}]`)}
      ${row('p (N_eff)', s => (s.p < 0.05 ? '<b>' : '') + fp(s.p) + (s.p < 0.05 ? '</b>' : ''))}${row('N_eff / lag-1 autocorr. (yours)', s => `${s.neff.toFixed(1)} / ${f(s.ar1y)}`)}
      ${row('Spearman ρ (p)', s => `${f(s.rho)} (${fp(s.p_rho)})`)}
      <tr><td>ΔR = r<sub>RONI</sub> − r<sub>ONI</sub></td><td colspan="2">${diff ? `${f(diff.dr, 3)} (MRR p ${diff.p < 0.001 ? '< 0.001' : '= ' + diff.p.toFixed(3)}; r<sub>ONI,RONI</sub> = ${f(diff.r12)})` : '—'}</td></tr></tbody></table>
      ${st.source === 'lmr' ? '<div class="small">* LMR: CI combines the spread of 100 reconstruction members with sampling error; r is for the ensemble-median index.</div>' : ''}`;
    // time series
    const W = wOf('cmp-enso-fig'), H = 270, m = { l: 48, r: 56, t: 34, b: 34 }, svg = svgIn('cmp-enso-fig', W, H);
    if (!o) return out;
    const A = out.oni.A, uz = S.zscore(A.u), x = d3.scaleLinear().domain([st.y0 - 0.5, st.y1 + 0.5]).range([m.l, W - m.r]);
    const zmax = Math.ceil(Math.max(2.5, d3.max(uz, Math.abs))), yz = d3.scaleLinear().domain([-zmax, zmax]).range([H - m.b, m.t]);
    const imax = Math.ceil(2 * Math.max(1.5, d3.max([...out.oni.A.t, ...(out.roni.A.t || [])], Math.abs))) / 2, yi = d3.scaleLinear().domain([-imax, imax]).range([H - m.b, m.t]);
    A.Y.forEach(y => { const c = epiFlag('oni', y - st.lag); if (c) svg.append('rect').attr('x', x(y - 0.5)).attr('width', x(y + 0.5) - x(y - 0.5)).attr('y', m.t).attr('height', H - m.t - m.b).attr('fill', c > 0 ? C.EN : C.LN).attr('opacity', 0.13); });
    svg.append('g').attr('class', 'gridl').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(yz).ticks(5).tickSize(-(W - m.l - m.r)).tickFormat(''));
    ['oni', 'roni'].forEach(k => { const B = out[k].A; if (B.n < 8) return; svg.append('path').datum(B.Y.map((y, i) => [y, B.t[i]])).attr('d', d3.line().x(d => x(d[0])).y(d => yi(d[1]))).attr('fill', 'none').attr('stroke', C[k.toUpperCase()]).attr('stroke-width', 1.5).attr('stroke-dasharray', k === 'roni' ? '4 3' : null); });
    svg.append('path').datum(A.Y.map((y, i) => [y, uz[i]])).attr('d', d3.line().x(d => x(d[0])).y(d => yz(d[1]))).attr('fill', 'none').attr('stroke', C.user).attr('stroke-width', 1.8);
    svg.append('g').attr('class', 'ax').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(Math.min(12, W / 80)).tickFormat(d3.format('d')));
    svg.append('g').attr('class', 'ax').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(yz).ticks(5));
    svg.append('g').attr('class', 'ax').attr('transform', `translate(${W - m.r},0)`).call(d3.axisRight(yi).ticks(5));
    svg.append('text').attr('class', 'lab').attr('transform', `translate(13,${(H + m.t - m.b) / 2}) rotate(-90)`).attr('text-anchor', 'middle').text(`${st.label} (σ)`);
    svg.append('text').attr('class', 'lab').attr('transform', `translate(${W - 12},${(H + m.t - m.b) / 2}) rotate(90)`).attr('text-anchor', 'middle').text(`index (°C)${st.lag ? `, ENSO year Y − ${st.lag}` : ''}`);
    const lg = [[st.label, C.user, null], ['ONI', C.ONI, null], ['RONI', C.RONI, '4 3'], ['ONI El Niño yr', C.EN, 'box'], ['ONI La Niña yr', C.LN, 'box']];
    let lx = m.l; lg.forEach(([t, c, d]) => { if (d === 'box') svg.append('rect').attr('x', lx).attr('y', 11).attr('width', 12).attr('height', 9).attr('fill', c).attr('opacity', 0.3); else svg.append('line').attr('x1', lx).attr('x2', lx + 14).attr('y1', 16).attr('y2', 16).attr('stroke', c).attr('stroke-width', 2).attr('stroke-dasharray', d); svg.append('text').attr('class', 'sm').attr('x', lx + 18).attr('y', 20).text(t); lx += 26 + t.length * 5.6; });
    svg.append('rect').attr('x', m.l).attr('y', m.t).attr('width', W - m.l - m.r).attr('height', H - m.t - m.b).attr('fill', 'transparent')
      .on('mousemove', ev => { const y = Math.round(x.invert(d3.pointer(ev)[0])), i = A.Y.indexOf(y); if (i < 0) return hideTip(); showTip(ev, `<b>${y}</b> — ${esc(st.label)}: ${f(A.uRaw[i], 3)} ${esc(st.units)}<br>ONI ${f(out.oni.A.t[i])} · RONI ${f(out.roni.A.t?.[i])} °C`); }).on('mouseleave', hideTip);
    return out;
  }

  // B. lag scan -------------------------------------------------------------------------------
  function renderLag(T) {
    const W = wOf('cmp-lag-fig'), H = 220, m = { l: 44, r: 12, t: 22, b: 36 }, svg = svgIn('cmp-lag-fig', W, H);
    const L = d3.range(-5, 6), rows = [];
    ['oni', 'roni'].forEach(k => L.forEach(l => { const A = align(T.find(t => t.id === 'enso:' + k).get, l); if (A.n >= 8) { const s = S.corrStats(A.t, A.u); rows.push({ k, l, r: s.r, p: s.p, n: A.n }); } }));
    if (!rows.length) return;
    const x0 = d3.scaleBand().domain(L).range([m.l, W - m.r]).padding(0.2), x1 = d3.scaleBand().domain(['oni', 'roni']).range([0, x0.bandwidth()]).padding(0.1);
    const ext = Math.max(0.3, d3.max(rows, d => Math.abs(d.r)) + 0.05), y = d3.scaleLinear().domain([-ext, ext]).nice().range([H - m.b, m.t]);
    svg.append('g').attr('class', 'gridl').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5).tickSize(-(W - m.l - m.r)).tickFormat(''));
    rows.forEach(d => svg.append('rect').attr('x', x0(d.l) + x1(d.k)).attr('width', x1.bandwidth()).attr('y', y(Math.max(0, d.r))).attr('height', Math.abs(y(d.r) - y(0)))
      .attr('fill', C[d.k.toUpperCase()]).attr('opacity', d.p < 0.05 ? 1 : 0.35).attr('stroke', d.l === st.lag ? '#111' : 'none')
      .on('mousemove', ev => showTip(ev, `${d.k.toUpperCase()}, lag ${d.l}: r = ${f(d.r)}, p = ${fp(d.p)}, n = ${d.n}`)).on('mouseleave', hideTip));
    svg.append('line').attr('x1', m.l).attr('x2', W - m.r).attr('y1', y(0)).attr('y2', y(0)).attr('stroke', '#333');
    svg.append('g').attr('class', 'ax').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x0).tickFormat(d => (d > 0 ? '+' : '') + d));
    svg.append('g').attr('class', 'ax').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5));
    svg.append('text').attr('class', 'lab').attr('x', (m.l + W - m.r) / 2).attr('y', H - 4).attr('text-anchor', 'middle').text('Lag (yr; > 0 = ENSO leads your series)');
    svg.append('text').attr('class', 'sm').attr('x', m.l).attr('y', 12).text('Solid bars: p < 0.05 (N_eff); outlined: current lag');
    const best = rows.slice().sort((a, b) => Math.abs(b.r) - Math.abs(a.r))[0];
    $('cmp-lag-note').textContent = `Strongest |r| at lag ${best.l > 0 ? '+' : ''}${best.l} yr with ${best.k.toUpperCase()} (r = ${f(best.r)}, p = ${fp(best.p)}). Scanning 11 lags × 2 indices inflates the chance of a spurious peak; treat isolated peaks at |lag| > 1 with caution.`;
  }

  // C. running correlation --------------------------------------------------------------------
  function renderRunning(T) {
    const W = wOf('cmp-run-fig'), H = 230, m = { l: 44, r: 12, t: 22, b: 36 }, svg = svgIn('cmp-run-fig', W, H), h = (st.runW - 1) / 2;
    const series = {};
    ['oni', 'roni'].forEach(k => {
      const A = align(T.find(t => t.id === 'enso:' + k).get, st.lag, -1e9, 1e9), out = [];
      if (A.n >= st.runW) for (let c = A.Y[0] + h; c <= A.Y[A.Y.length - 1] - h; c++) {
        const idx = A.Y.map((y, i) => (y >= c - h && y <= c + h) ? i : -1).filter(i => i >= 0);
        if (idx.length >= 0.8 * st.runW) out.push({ c, r: S.pearson(idx.map(i => A.t[i]), idx.map(i => A.u[i])) });
      }
      series[k] = out;
    });
    const all = [...series.oni, ...series.roni];
    if (!all.length) { svg.append('text').attr('class', 'lab').attr('x', 12).attr('y', 30).text(`Overlap shorter than the ${st.runW}-yr window.`); return; }
    const x = d3.scaleLinear().domain(d3.extent(all, d => d.c)).range([m.l, W - m.r]), y = d3.scaleLinear().domain([-1, 1]).range([H - m.b, m.t]);
    svg.append('g').attr('class', 'gridl').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5).tickSize(-(W - m.l - m.r)).tickFormat(''));
    const rc = S.rCrit(st.runW); [rc, -rc].forEach(v => svg.append('line').attr('x1', m.l).attr('x2', W - m.r).attr('y1', y(v)).attr('y2', y(v)).attr('stroke', '#888').attr('stroke-dasharray', '5 3'));
    svg.append('line').attr('x1', m.l).attr('x2', W - m.r).attr('y1', y(0)).attr('y2', y(0)).attr('stroke', '#333');
    ['oni', 'roni'].forEach(k => svg.append('path').datum(series[k]).attr('d', d3.line().x(d => x(d.c)).y(d => y(d.r))).attr('fill', 'none').attr('stroke', C[k.toUpperCase()]).attr('stroke-width', 2).attr('stroke-dasharray', k === 'roni' ? '4 3' : null));
    svg.append('g').attr('class', 'ax').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(8).tickFormat(d3.format('d')));
    svg.append('g').attr('class', 'ax').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5));
    svg.append('text').attr('class', 'lab').attr('x', (m.l + W - m.r) / 2).attr('y', H - 4).attr('text-anchor', 'middle').text(`Centre year of ${st.runW}-yr window (full overlap)`);
    svg.append('text').attr('class', 'sm').attr('x', m.l).attr('y', 12).text(`Dashed grey: p = 0.05 for n = ${st.runW}. Pink: ONI · green dashed: RONI`);
  }

  // D. composites -----------------------------------------------------------------------------
  function renderComposite() {
    const W = wOf('cmp-comp-fig'), H = 230, m = { l: 44, r: 12, t: 26, b: 40 }, svg = svgIn('cmp-comp-fig', W, H);
    const res = {};
    ['oni', 'roni'].forEach(k => {
      const g = { 1: [], 0: [], [-1]: [] }, A = align(Y => ensoGet(k, Y)); if (A.n < 8) return;
      const z = S.zscore(A.u); A.Y.forEach((y, i) => { const c = epiFlag(k, y - st.lag); if (c !== null) g[c].push(z[i]); });
      res[k] = { EN: S.groupSummary(g[1]), N: S.groupSummary(g[0]), LN: S.groupSummary(g[-1]), p: S.welch(g[1], g[-1]).p, perm: S.permTest(g[1], g[-1], 3000) };
    });
    if (!res.oni) return;
    const G = ['EN', 'N', 'LN'], name = { EN: 'El Niño', N: 'Neutral', LN: 'La Niña' };
    const x0 = d3.scaleBand().domain(G).range([m.l, W - m.r]).padding(0.25), x1 = d3.scaleBand().domain(['oni', 'roni']).range([0, x0.bandwidth()]).padding(0.08);
    const vals = Object.values(res).flatMap(r => G.flatMap(g => r[g].ci.filter(Number.isFinite)));
    const y = d3.scaleLinear().domain([Math.min(-0.8, d3.min(vals)), Math.max(0.8, d3.max(vals))]).nice().range([H - m.b, m.t]);
    svg.append('g').attr('class', 'gridl').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5).tickSize(-(W - m.l - m.r)).tickFormat(''));
    G.forEach(g => ['oni', 'roni'].forEach(k => {
      const s = res[k]?.[g]; if (!s || !s.n) return; const bx = x0(g) + x1(k), bw = x1.bandwidth();
      svg.append('rect').attr('x', bx).attr('width', bw).attr('y', y(Math.max(0, s.mean))).attr('height', Math.abs(y(s.mean) - y(0))).attr('fill', g === 'EN' ? C.EN : g === 'LN' ? C.LN : '#9a9a9a').attr('opacity', k === 'oni' ? 0.45 : 0.9).attr('stroke', C[k.toUpperCase()]).attr('stroke-width', 2);
      if (Number.isFinite(s.ci[0])) svg.append('line').attr('x1', bx + bw / 2).attr('x2', bx + bw / 2).attr('y1', y(s.ci[0])).attr('y2', y(s.ci[1])).attr('stroke', '#111');
      svg.append('text').attr('class', 'sm').attr('x', bx + bw / 2).attr('y', H - m.b + 26).attr('text-anchor', 'middle').text(`n=${s.n}`);
    }));
    svg.append('line').attr('x1', m.l).attr('x2', W - m.r).attr('y1', y(0)).attr('y2', y(0)).attr('stroke', '#333');
    svg.append('g').attr('class', 'ax').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5));
    G.forEach(g => svg.append('text').attr('class', 'lab').attr('x', x0(g) + x0.bandwidth() / 2).attr('y', H - m.b + 13).attr('text-anchor', 'middle').attr('font-weight', 700).text(name[g]));
    svg.append('text').attr('class', 'sm').attr('x', m.l).attr('y', 12).text(`EN − LN: ONI p = ${fp(res.oni.p)} (perm. ${fp(res.oni.perm)})${res.roni ? `; RONI p = ${fp(res.roni.p)} (perm. ${fp(res.roni.perm)})` : ''} · faded = ONI, solid = RONI · mean z ± 95% CI`);
  }

  // E. screen against every historical series --------------------------------------------------
  function renderScreen(T) {
    const rows = [];
    T.forEach(t => { const A = align(t.get); if (A.n >= 10) { const s = S.corrStats(A.t, A.u); rows.push({ t, n: A.n, r: s.r, p: s.p, neff: s.neff, y0: A.Y[0], y1: A.Y[A.Y.length - 1], same: Math.abs(s.r) > 0.9995 }); } });
    // identical series (e.g. your file is a copy of an archived record) are shown but excluded from the FDR family
    const fam = rows.filter(r => !r.same), q = bh(fam.map(r => r.p)); fam.forEach((r, i) => r.q = q[i]); rows.filter(r => r.same).forEach(r => r.q = NaN);
    rows.sort((a, b) => Math.abs(b.r) - Math.abs(a.r));
    const nsig = rows.filter(r => r.q < 0.05).length, nsame = rows.filter(r => r.same).length;
    $('cmp-screen').innerHTML = !rows.length ? '<p class="small">No historical series overlaps your data by ≥ 10 years in this window.</p>' :
      `<div class="small">${rows.length} series overlap by ≥ 10 years · ${nsig} significant after false-discovery-rate control (q < 0.05)${nsame ? ` · ${nsame} identical to your series (excluded from the FDR test)` : ''} · click a row to compare in detail</div>
      <div class="tbl-wrap"><table class="st screen"><thead><tr><th>Historical series</th><th>n</th><th>overlap</th><th>r</th><th>p (N_eff)</th><th>q (FDR)</th></tr></thead><tbody>` +
      rows.map(r => `<tr data-id="${r.t.id}" class="${r.t.id === st.target ? 'on' : ''} ${r.same ? 'same' : ''}"><td>${esc(r.t.label)}${r.same ? ' — identical to your series' : ''}</td><td>${r.n}</td><td>${r.y0}–${r.y1}</td>
        <td><span class="rbar" style="--w:${Math.min(100, Math.abs(r.r) * 100)}%;--c:${r.r >= 0 ? '#b2182b' : '#2166ac'}"></span>${f(r.r)}</td><td>${fp(r.p)}</td><td class="${r.q < 0.05 ? 'sig' : ''}">${r.same ? '—' : fp(r.q)}</td></tr>`).join('') + '</tbody></table></div>';
    document.querySelectorAll('#cmp-screen tr[data-id]').forEach(tr => tr.onclick = () => { st.target = tr.dataset.id; $('cmp-target').value = st.target; analyse(); $('cmp-target-card').scrollIntoView({ behavior: 'smooth', block: 'start' }); });
    return { rows, nsig };
  }

  // F. any target: time series + scatter --------------------------------------------------------
  function renderTarget(T) {
    const t = T.find(x => x.id === st.target), A = align(t.get);
    const W = wOf('cmp-tgt-fig'), H = 300, svg = svgIn('cmp-tgt-fig', W, H);
    if (A.n < 8) { svg.append('text').attr('class', 'lab').attr('x', 12).attr('y', 30).text(`Fewer than 8 overlapping years with ${t.label}.`); $('cmp-tgt-stats').innerHTML = ''; return; }
    const s = S.corrStats(A.t, A.u), uz = S.zscore(A.u), tz = S.zscore(A.t);
    const split = W > 700, wTs = split ? W * 0.64 : W, hTs = split ? H : H * 0.5;
    const m = { l: 44, r: 14, t: 26, b: 32 }, x = d3.scaleLinear().domain(d3.extent(A.Y)).range([m.l, wTs - m.r]);
    const zm = Math.ceil(Math.max(2.5, d3.max([...uz, ...tz], Math.abs))), y = d3.scaleLinear().domain([-zm, zm]).range([hTs - m.b, m.t]);
    svg.append('g').attr('class', 'gridl').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5).tickSize(-(wTs - m.l - m.r)).tickFormat(''));
    svg.append('path').datum(A.Y.map((yy, i) => [yy, tz[i]])).attr('d', d3.line().x(d => x(d[0])).y(d => y(d[1]))).attr('fill', 'none').attr('stroke', C.T).attr('stroke-width', 1.6);
    svg.append('path').datum(A.Y.map((yy, i) => [yy, uz[i]])).attr('d', d3.line().x(d => x(d[0])).y(d => y(d[1]))).attr('fill', 'none').attr('stroke', C.user).attr('stroke-width', 1.8);
    svg.append('g').attr('class', 'ax').attr('transform', `translate(0,${hTs - m.b})`).call(d3.axisBottom(x).ticks(8).tickFormat(d3.format('d')));
    svg.append('g').attr('class', 'ax').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5));
    svg.append('text').attr('class', 'sm').attr('x', m.l).attr('y', 14).text(`Black: ${st.label} · orange: ${t.short} (both z-scores${st.detrend ? ', detrended' : ''})`);
    // scatter
    const sx0 = split ? wTs + 24 : m.l, sy0 = split ? 0 : hTs, sw = split ? W - wTs - 34 : W - m.l - m.r, sh = split ? H : H * 0.5;
    const g = svg.append('g').attr('transform', `translate(${sx0},${sy0})`);
    const px = d3.scaleLinear().domain(d3.extent(A.t)).nice().range([30, sw - 6]), py = d3.scaleLinear().domain(d3.extent(A.u)).nice().range([sh - m.b, m.t]);
    g.selectAll('circle').data(A.Y.map((yy, i) => [A.t[i], A.u[i], yy])).join('circle').attr('cx', d => px(d[0])).attr('cy', d => py(d[1])).attr('r', 2.8).attr('fill', C.T).attr('opacity', 0.75)
      .on('mousemove', (ev, d) => showTip(ev, `${d[2]}: ${f(d[0], 3)} vs ${f(d[1], 3)}`)).on('mouseleave', hideTip);
    const b = S.pearson(A.t, A.u) * S.sd(A.u) / S.sd(A.t), mx = S.mean(A.t), my = S.mean(A.u), [a0, a1] = px.domain();
    g.append('line').attr('x1', px(a0)).attr('x2', px(a1)).attr('y1', py(my + b * (a0 - mx))).attr('y2', py(my + b * (a1 - mx))).attr('stroke', '#111').attr('stroke-width', 1.5);
    g.append('g').attr('class', 'ax').attr('transform', `translate(0,${sh - m.b})`).call(d3.axisBottom(px).ticks(4));
    g.append('g').attr('class', 'ax').attr('transform', 'translate(30,0)').call(d3.axisLeft(py).ticks(4));
    g.append('text').attr('class', 'sm').attr('x', 30).attr('y', 14).text(`r = ${f(s.r)}, ρ = ${f(s.rho)}`);
    $('cmp-tgt-stats').innerHTML = `n = ${s.n} (${A.Y[0]}–${A.Y[A.Y.length - 1]}) · r = <b>${f(s.r)}</b> [${f(s.ci[0])}, ${f(s.ci[1])}] · p<sub>N_eff</sub> = ${fp(s.p)} (N_eff = ${s.neff.toFixed(1)}; naive p = ${fp(s.p_naive)}) · Spearman ρ = ${f(s.rho)} · slope = ${f(b, 3)} ${esc(st.units || 'units')} per ${esc(t.units)} · lag ${st.lag} yr${st.detrend ? ' · detrended' : ''}`;
    renderTarget.last = { t, A };
  }

  // G. interpretation --------------------------------------------------------------------------
  function renderInterp(enso, scr) {
    const o = enso.oni.s, r = enso.roni.s, d = enso.diff, p = [];
    if (o) {
      const strength = a => a < 0.1 ? 'negligible' : a < 0.3 ? 'weak' : a < 0.5 ? 'moderate' : 'strong';
      p.push(`<p><b>ENSO coupling.</b> Over ${st.y0}–${st.y1}, ${esc(st.label)} shows a ${strength(Math.abs(o.r))} ${o.r >= 0 ? 'positive' : 'negative'} correlation with the ${st.source === 'lmr' ? 'annual' : seasonName()} ONI (r = ${f(o.r)}, n = ${o.n}), ${o.p < 0.05 ? `significant after accounting for autocorrelation (p ${o.p < 0.001 ? '< 0.001' : '= ' + o.p.toFixed(3)})` : `not significant once autocorrelation is considered (p = ${fp(o.p)})`}.` +
        (o.n - o.neff > 3 ? ` Persistence in your series (lag-1 autocorrelation ${f(o.ar1y)}) reduces the effective sample size from ${o.n} to ${o.neff.toFixed(0)}.` : '') + '</p>');
      if (r && d) p.push(`<p><b>ONI vs RONI.</b> Under RONI r = ${f(r.r)} (ΔR = ${f(d.dr, 3)}). ${d.p < 0.05 ? `The difference is significant (p = ${fp(d.p)}): ${Math.abs(r.r) < Math.abs(o.r) ? 'part of the ONI coupling is carried by tropical-mean SST, which RONI removes' : 'removing tropical-mean SST sharpens the ENSO signal in your record'}.` : `The difference is not significant (p = ${fp(d.p)}); ONI and RONI correlate at ${f(d.r12)} over this window.`}</p>`);
    }
    if (scr.rows.length) {
      const same = scr.rows.filter(r => r.same), top = scr.rows.find(r => !r.same);
      const qs = q => q < 0.001 ? 'q < 0.001' : `q = ${q.toFixed(3)}`;
      p.push(`<p><b>Screen.</b> ${same.length ? `Your series is identical to ${esc(same[0].t.label)} over the overlap, so it is set aside. ` : ''}` +
        (top ? `The strongest other match among ${scr.rows.length - same.length} historical series is ${esc(top.t.label)} (r = ${f(top.r)}, ${qs(top.q)}). ${scr.nsig ? `${scr.nsig} ${scr.nsig === 1 ? 'series remains' : 'series remain'} significant after false-discovery-rate control.` : 'None survives false-discovery-rate control; with this many comparisons, some nominal p < 0.05 hits are expected by chance.'}` : '') + '</p>');
    }
    p.push(`<p class="small">Your data stay in this browser tab; nothing is uploaded.</p>`);
    $('cmp-interp').innerHTML = p.join('');
  }

  // ------------------------------------------------------------------ export
  function exportSvg(id, kind) {
    const svg = document.querySelector(`#${id} svg`); if (!svg) return;
    const src = new XMLSerializer().serializeToString(svg), W = +svg.getAttribute('width'), H = +svg.getAttribute('height');
    const url = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(src), name = `yourdata_${id.replace('cmp-', '').replace('-fig', '')}`;
    if (kind === 'svg') return Site.exportPanel({ kind: 'image', src: url, name: name + '.svg' });
    const img = new Image(); img.onload = () => { const c = document.createElement('canvas'); c.width = W * 3; c.height = H * 3; const x = c.getContext('2d'); x.fillStyle = '#fff'; x.fillRect(0, 0, c.width, c.height); x.drawImage(img, 0, 0, c.width, c.height); Site.exportPanel({ kind: 'image', src: c.toDataURL('image/png'), name: name + '.png' }); }; img.src = url;
  }
  function exportCsv() {
    const T = targets(), get = id => T.find(t => t.id === id).get;
    const lines = [['year', st.label.replace(/,/g, ' '), 'ONI', 'RONI', 'tropical_mean_SST', 'ONI_class', 'RONI_class'].join(',')];
    for (const y of st.user.years) {
      if (y < st.y0 || y > st.y1) continue;
      const v = k => { const x = get('enso:' + k)(y - st.lag); return x === null || x === undefined ? '' : +x.toFixed(4); };
      lines.push([y, st.user.map.get(y), v('oni'), v('roni'), v('t20'), epiFlag('oni', y - st.lag) ?? '', epiFlag('roni', y - st.lag) ?? ''].join(','));
    }
    Site.exportPanel({ kind: 'text', text: `# ${st.label} aligned with ${SRC_LABEL[st.source]} (${st.source === 'lmr' ? 'annual' : seasonName()}), lag ${st.lag} yr (ENSO year = year − lag)\n` + lines.join('\n'), name: 'yourdata_aligned.csv' });
  }

  // ------------------------------------------------------------------ loading
  function load(text, name, keepLabel) {
    try {
      const p = parseText(text); st.raw = { name, keepLabel }; st.cols = p.cols; st.rows = p.rows; st.y0 = null; st._auto = true;
      guessColumns(); if (keepLabel) st.label = keepLabel; $('cmp-label').value = st.label;
      $('cmp-file-name').textContent = name; renderPreview(); analyse();
    } catch (e) { const m = $('cmp-msg'); m.textContent = e.message; m.className = 'msg err'; $('cmp-results').hidden = true; }
  }
  function exampleCsv(coreId, key, ion) {
    const c = D.cores.find(c => c.id === coreId), v = ion ? c.vars.ions[ion] : c.vars[key];
    return `# Example exported from ${c.name} (${c.ref}), NOAA NCEI\nyear,${v.label.replace(/[^\w]+/g, '_')}\n` + v.years.map((y, i) => `${y},${v.values[i] ?? 'NA'}`).join('\n');
  }
  const EXAMPLES = {
    quel_acc: ['Quelccaya accumulation (m w.e.)', () => exampleCsv('QUE13', 'accum'), 'Quelccaya accumulation'],
    ill_nh4: ['Illimani NH₄⁺ (1800–1998)', () => exampleCsv('ILL99', 'ions', 'NH4'), 'Illimani NH4'],
    hua_dex: ['Huascarán Summit d-excess', () => exampleCsv('HUASUM', 'dexcess'), 'Huascarán d-excess']
  };

  function init() {
    const drop = $('cmp-drop');
    $('cmp-file').onchange = e => { const fl = e.target.files[0]; if (fl) fl.text().then(t => load(t, fl.name)); };
    ['dragenter', 'dragover'].forEach(ev => drop.addEventListener(ev, e => { e.preventDefault(); drop.classList.add('hover'); }));
    ['dragleave', 'drop'].forEach(ev => drop.addEventListener(ev, e => { e.preventDefault(); drop.classList.remove('hover'); }));
    drop.addEventListener('drop', e => { const fl = e.dataTransfer.files[0]; if (fl) fl.text().then(t => load(t, fl.name)); });
    $('cmp-paste-go').onclick = () => load($('cmp-paste').value, 'pasted text');
    Object.entries(EXAMPLES).forEach(([k, [lab, fn, name]]) => { const b = document.createElement('button'); b.className = 'chipbtn'; b.textContent = lab; b.onclick = () => load(fn(), 'example: ' + lab, name); $('cmp-examples').appendChild(b); });
    const re = () => { if (st.rows.length) { renderPreview(); analyse(); } };
    $('cmp-ycol').onchange = e => { st.yearCol = +e.target.value; re(); };
    $('cmp-mcol').onchange = e => { st.monthCol = +e.target.value; re(); };
    $('cmp-vcol').onchange = e => { st.valCol = +e.target.value; st.label = st.cols[st.valCol]; $('cmp-label').value = st.label; re(); };
    $('cmp-bp').onchange = e => { st.bp = e.target.checked; st._auto = true; re(); };
    $('cmp-missing').onchange = e => { st.missing = e.target.value; re(); };
    $('cmp-agg').onchange = e => { st.agg = e.target.value; re(); };
    $('cmp-label').onchange = e => { st.label = e.target.value || 'Your series'; re(); };
    $('cmp-units').onchange = e => { st.units = e.target.value; re(); };
    $('cmp-source').onchange = e => { st.source = e.target.value; st._auto = true; $('cmp-season').disabled = st.source === 'lmr'; re(); };
    $('cmp-season').onchange = e => { st.season = e.target.value; $('cmp-custom').hidden = st.season !== 'custom'; re(); };
    MON.forEach((mn, i) => ['cmp-cm0', 'cmp-cm1'].forEach(id => $(id).insertAdjacentHTML('beforeend', `<option value="${i + 1}">${mn}</option>`)));
    $('cmp-cm0').value = st.cm0; $('cmp-cm1').value = st.cm1;
    $('cmp-cm0').onchange = e => { st.cm0 = +e.target.value; re(); }; $('cmp-cm1').onchange = e => { st.cm1 = +e.target.value; re(); };
    $('cmp-lag').oninput = e => { st.lag = +e.target.value; $('cmp-lag-v').textContent = (st.lag > 0 ? '+' : '') + st.lag + ' yr'; re(); };
    $('cmp-detrend').onchange = e => { st.detrend = e.target.checked; re(); };
    $('cmp-y0').onchange = e => { st.y0 = +e.target.value; st._auto = false; re(); };
    $('cmp-y1').onchange = e => { st.y1 = +e.target.value; st._auto = false; re(); };
    $('cmp-runw').oninput = e => { st.runW = +e.target.value; $('cmp-runw-v').textContent = st.runW + ' yr'; re(); };
    $('cmp-target').onchange = e => { st.target = e.target.value; re(); };
    document.querySelectorAll('[data-cmp-export]').forEach(b => b.onclick = () => exportSvg(b.dataset.cmpExport, b.dataset.kind));
    $('cmp-csv').onclick = exportCsv;
    let rt; addEventListener('resize', () => { clearTimeout(rt); rt = setTimeout(() => { if (Site.current === 'compare') re(); }, 200); });
    if (window.Site) Site.onTab(t => { if (t === 'compare') requestAnimationFrame(re); });
    // open in a working state with an example
    load(EXAMPLES.quel_acc[1](), 'example: ' + EXAMPLES.quel_acc[0], EXAMPLES.quel_acc[2]);
  }
  init();
})();
