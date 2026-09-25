// Emits JS statistics on a fixed set of test vectors for comparison with scipy.
const S = require('../app/js/stats.js');
const fs = require('fs');
const src = fs.readFileSync(__dirname + '/../app/data/dataset.js', 'utf8');
const DATA = JSON.parse(src.slice(src.indexOf('=') + 1).trim().replace(/;$/, ''));
// Quelccaya d18O vs CPC ONI DJF (Jan value), 1950-2009
const c = DATA.cores.find(c => c.id === 'QUE13').vars.d18O;
const x = [], y = [], t = [];
for (let i = 0; i < c.years.length; i++) {
  const Y = c.years[i]; if (Y < 1951 || Y > 2009) continue;
  const v = DATA.cpc.oni[(Y - 1950) * 12]; if (v == null) continue;
  x.push(c.values[i]); y.push(v); t.push(Y);
}
const cs = S.corrStats(x, y);
const a = [0.3, -1.2, 0.8, 1.9, -0.4, 0.2, 1.1], b = [-0.5, -1.4, 0.1, -0.9, -2.2, 0.4];
const w = S.welch(a, b);
const reg = S.regress2(x, y, t, x.length);
const out = {
  n: cs.n, r: cs.r, rho: cs.rho, ar1x: cs.ar1x, ar1y: cs.ar1y, neff: cs.neff, p_naive: cs.p_naive, p: cs.p, p_rho: cs.p_rho, ci: cs.ci,
  welch: w, tq: S.tQuantile(0.975, 7.3), pnorm: S.pNorm2(1.96), qn: S.qNorm(0.975),
  mrr: S.compareDepCorr(0.45, 0.30, 0.95, 60), detr_slope: S.detrend(t, x).slope,
  reg: { b1: reg.b1, b2: reg.b2, R2: reg.R2, pr1: reg.pr1, pr2: reg.pr2 },
  x, y, t, a, b
};
console.log(JSON.stringify(out));
