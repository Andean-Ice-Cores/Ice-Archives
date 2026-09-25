/* Statistics for the ENSO–ice-core explorer. Pure functions, no DOM.
 * Works in the browser (window.Stats) and in Node (module.exports) for testing.
 *
 * Key references
 *  - Effective sample size for correlation of two AR(1) series:
 *      Bretherton et al. (1999) J. Climate 12:1990, eq. 31:  N* = N (1 − ρ1ρ2)/(1 + ρ1ρ2)
 *  - Comparing two dependent correlations sharing one variable:
 *      Meng, Rosenthal & Rubin (1992) Psychol. Bull. 111:172
 *  - Banded age model (layer-counting errors): Comboul et al. (2014) Clim. Past 10:825
 */
(function (root) {
  'use strict';

  const isNum = v => v !== null && v !== undefined && Number.isFinite(v);
  const mean = a => a.reduce((s, v) => s + v, 0) / a.length;
  function sd(a) { const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
  const zscore = a => { const m = mean(a), s = sd(a); return a.map(v => (v - m) / s); };

  function pearson(x, y) {
    const n = x.length, mx = mean(x), my = mean(y);
    let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < n; i++) { const a = x[i] - mx, b = y[i] - my; sxy += a * b; sxx += a * a; syy += b * b; }
    return sxy / Math.sqrt(sxx * syy);
  }
  function ranks(a) {           // average ranks for ties (as scipy.stats.rankdata)
    const idx = a.map((v, i) => [v, i]).sort((p, q) => p[0] - q[0]);
    const r = new Array(a.length);
    for (let i = 0; i < idx.length;) {
      let j = i; while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++;
      const rk = (i + j) / 2 + 1; for (let k = i; k <= j; k++) r[idx[k][1]] = rk; i = j + 1;
    }
    return r;
  }
  const spearman = (x, y) => pearson(ranks(x), ranks(y));
  function ar1(x) {             // lag-1 autocorrelation (biased estimator, as commonly used for N*)
    const m = mean(x); let num = 0, den = 0;
    for (let i = 0; i < x.length; i++) { den += (x[i] - m) ** 2; if (i) num += (x[i] - m) * (x[i - 1] - m); }
    return num / den;
  }
  function neff(x, y) {
    const n = x.length, r1 = ar1(x), r2 = ar1(y), p = r1 * r2;
    let ne = n * (1 - p) / (1 + p);
    ne = Math.min(n, Math.max(3, ne));   // never inflate beyond N (conservative); keep ≥3
    return { neff: ne, r1, r2 };
  }

  // ---- special functions ------------------------------------------------------------
  function lgamma(z) {          // Lanczos
    const g = 7, c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
      -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
    if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lgamma(1 - z);
    z -= 1; let x = c[0]; for (let i = 1; i < g + 2; i++) x += c[i] / (z + i);
    const t = z + g + 0.5; return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
  }
  function betacf(a, b, x) {
    const MAXIT = 300, EPS = 3e-14, FPMIN = 1e-300;
    let qab = a + b, qap = a + 1, qam = a - 1, c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN; d = 1 / d; let h = d;
    for (let m = 1; m <= MAXIT; m++) {
      const m2 = 2 * m; let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d; const del = d * c; h *= del; if (Math.abs(del - 1) < EPS) break;
    }
    return h;
  }
  function ibeta(x, a, b) {     // regularized incomplete beta I_x(a,b)
    if (x <= 0) return 0; if (x >= 1) return 1;
    const bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(1 - x));
    return x < (a + 1) / (a + b + 2) ? bt * betacf(a, b, x) / a : 1 - bt * betacf(b, a, 1 - x) / b;
  }
  const pT2 = (t, df) => ibeta(df / (df + t * t), df / 2, 0.5);              // two-sided Student-t p
  function erf(x) {             // Abramowitz–Stegun 7.1.26 refined via series for small |x|
    const s = Math.sign(x); x = Math.abs(x);
    if (x < 2.5) {              // Taylor series is very accurate here
      let sum = x, term = x, n = 0;
      while (Math.abs(term) > 1e-16 * Math.abs(sum) && n < 200) { n++; term *= -x * x / n; sum += term / (2 * n + 1); }
      return s * 2 / Math.sqrt(Math.PI) * sum;
    }
    // continued fraction for erfc
    let f = 0; for (let k = 60; k >= 1; k--) f = k / 2 / (x + f);
    return s * (1 - Math.exp(-x * x) / Math.sqrt(Math.PI) / (x + f));
  }
  const pNorm2 = z => 1 - erf(Math.abs(z) / Math.SQRT2);                      // two-sided normal p
  function qNorm(p) {           // inverse normal CDF (Acklam)
    const a = [-3.969683028665376e1, 2.209460984245205e2, -2.759285104469687e2, 1.383577518672690e2, -3.066479806614716e1, 2.506628277459239];
    const b = [-5.447609879822406e1, 1.615858368580409e2, -1.556989798598866e2, 6.680131188771972e1, -1.328068155288572e1];
    const c = [-7.784894002430293e-3, -3.223964580411365e-1, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783];
    const d = [7.784695709041462e-3, 3.224671290700398e-1, 2.445134137142996, 3.754408661907416];
    const pl = 0.02425; let q, r;
    if (p < pl) { q = Math.sqrt(-2 * Math.log(p)); return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
    if (p > 1 - pl) { q = Math.sqrt(-2 * Math.log(1 - p)); return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
    q = p - 0.5; r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  }

  // ---- correlation with autocorrelation-corrected significance ------------------------
  function pCorr(r, n) {        // two-sided p for correlation r with (effective) sample size n
    if (n <= 2 || !isNum(r)) return NaN;
    const df = n - 2; if (Math.abs(r) >= 1) return 0;
    return pT2(r * Math.sqrt(df / (1 - r * r)), df);
  }
  function fisherCI(r, n, level = 0.95) {
    const z = Math.atanh(r), se = 1 / Math.sqrt(Math.max(n - 3, 1e-9)), q = qNorm(0.5 + level / 2);
    return [Math.tanh(z - q * se), Math.tanh(z + q * se)];
  }
  function rCrit(n, alpha = 0.05) {   // |r| needed for p<alpha given (effective) n — bisection
    let lo = 0, hi = 0.9999; for (let i = 0; i < 60; i++) { const m = (lo + hi) / 2; pCorr(m, n) > alpha ? lo = m : hi = m; }
    return hi;
  }
  function corrStats(x, y) {
    const n = x.length; if (n < 5) return null;
    const r = pearson(x, y), rho = spearman(x, y), ne = neff(x, y);
    return {
      n, r, rho, neff: ne.neff, ar1x: ne.r1, ar1y: ne.r2,
      p_naive: pCorr(r, n), p: pCorr(r, ne.neff), p_rho: pCorr(rho, ne.neff),
      ci: fisherCI(r, ne.neff)
    };
  }

  // Meng–Rosenthal–Rubin test: H0 r_xy1 == r_xy2 where y1,y2 correlated at r12, sample size n
  function compareDepCorr(r1, r2, r12, n) {
    const z1 = Math.atanh(r1), z2 = Math.atanh(r2), rm2 = (r1 * r1 + r2 * r2) / 2;
    const f = Math.min(1, (1 - r12) / (2 * (1 - rm2)));
    const h = (1 - f * rm2) / (1 - rm2);
    const Z = (z1 - z2) * Math.sqrt((n - 3) / (2 * (1 - r12) * h));
    return { Z, p: pNorm2(Z) };
  }

  function detrend(t, y) {
    const mt = mean(t), my = mean(y); let num = 0, den = 0;
    for (let i = 0; i < t.length; i++) { num += (t[i] - mt) * (y[i] - my); den += (t[i] - mt) ** 2; }
    const b = num / den; return { y: y.map((v, i) => v - my - b * (t[i] - mt)), slope: b };
  }

  // ---- composites ------------------------------------------------------------------------
  function welch(a, b) {
    if (a.length < 2 || b.length < 2) return { t: NaN, df: NaN, p: NaN };
    const va = sd(a) ** 2 / a.length, vb = sd(b) ** 2 / b.length;
    const t = (mean(a) - mean(b)) / Math.sqrt(va + vb);
    const df = (va + vb) ** 2 / (va * va / (a.length - 1) + vb * vb / (b.length - 1));
    return { t, df, p: pT2(t, df) };
  }
  function mulberry32(seed) { return function () { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
  function permTest(a, b, nperm = 5000, seed = 7) {   // two-sided test of mean(a) − mean(b)
    if (a.length < 2 || b.length < 2) return NaN;
    const all = a.concat(b), na = a.length, obs = Math.abs(mean(a) - mean(b)), rnd = mulberry32(seed);
    let cnt = 0; const w = all.slice();
    for (let k = 0; k < nperm; k++) {
      for (let i = w.length - 1; i > 0; i--) { const j = Math.floor(rnd() * (i + 1)); [w[i], w[j]] = [w[j], w[i]]; }
      let sa = 0; for (let i = 0; i < na; i++) sa += w[i];
      const sb = all.reduce((s, v) => s + v, 0) - sa;
      if (Math.abs(sa / na - sb / (w.length - na)) >= obs - 1e-12) cnt++;
    }
    return (cnt + 1) / (nperm + 1);
  }
  function groupSummary(a) {
    if (!a.length) return { n: 0, mean: NaN, se: NaN, ci: [NaN, NaN] };
    const m = mean(a), se = a.length > 1 ? sd(a) / Math.sqrt(a.length) : NaN;
    const tq = a.length > 1 ? tQuantile(0.975, a.length - 1) : NaN;
    return { n: a.length, mean: m, se, ci: [m - tq * se, m + tq * se] };
  }
  function tQuantile(p, df) {   // inverse Student t by bisection on the CDF
    let lo = 0, hi = 50; const target = 2 * (1 - p);
    for (let i = 0; i < 80; i++) { const m = (lo + hi) / 2; pT2(m, df) > target ? lo = m : hi = m; }
    return (lo + hi) / 2;
  }

  // ---- two-predictor standardized regression (decomposition) ------------------------------
  // y ~ b1 x1 + b2 x2 (all z-scored). Returns betas, SEs (using effective n), partial correlations.
  function regress2(y, x1, x2, nEff) {
    const zy = zscore(y), z1 = zscore(x1), z2 = zscore(x2), n = y.length;
    const r1 = pearson(zy, z1), r2 = pearson(zy, z2), r12 = pearson(z1, z2);
    const det = 1 - r12 * r12;
    const b1 = (r1 - r12 * r2) / det, b2 = (r2 - r12 * r1) / det;
    const R2 = b1 * r1 + b2 * r2;
    const ne = Math.min(n, nEff || n), df = Math.max(ne - 3, 1);
    const se = Math.sqrt((1 - R2) / (df * det));
    const pr1 = (r1 - r2 * r12) / Math.sqrt((1 - r2 * r2) * (1 - r12 * r12));
    const pr2 = (r2 - r1 * r12) / Math.sqrt((1 - r1 * r1) * (1 - r12 * r12));
    return {
      b1, b2, se1: se, se2: se, p1: pT2(b1 / se, df), p2: pT2(b2 / se, df), R2, r1, r2, r12,
      pr1, pr2, vif: 1 / det, df
    };
  }

  // ---- banded age model (Comboul et al. 2014) -----------------------------------------------
  // series ordered youngest→oldest; each year independently has P=theta of a missed layer and
  // P=theta of a double-counted layer; errors accumulate down-core.
  function bamPerturb(vals, theta, rnd) {
    const n = vals.length, out = new Array(n); let off = 0;
    for (let i = 0; i < n; i++) {
      const u = rnd(); if (u < theta) off -= 1; else if (u < 2 * theta) off += 1;
      const j = Math.min(n - 1, Math.max(0, i + off)); out[i] = vals[j];
    }
    return out;
  }

  function quantile(a, q) {
    const s = a.filter(isNum).sort((x, y) => x - y); if (!s.length) return NaN;
    const pos = (s.length - 1) * q, lo = Math.floor(pos), hi = Math.ceil(pos);
    return s[lo] + (s[hi] - s[lo]) * (pos - lo);
  }

  const Stats = {
    isNum, mean, sd, zscore, pearson, spearman, ranks, ar1, neff, pCorr, fisherCI, rCrit, corrStats,
    compareDepCorr, detrend, welch, permTest, groupSummary, tQuantile, regress2, bamPerturb,
    quantile, mulberry32, qNorm, pNorm2, pT2, ibeta
  };
  if (typeof module !== 'undefined' && module.exports) module.exports = Stats; else root.Stats = Stats;
})(typeof window !== 'undefined' ? window : globalThis);
