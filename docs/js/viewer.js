/* Interactive viewer for the Blender-built ENSO × Andean ice-core scene.
 * Geometry comes from blender/out/andes_enso.glb (embedded in scene_data.js); column heights, colours and
 * significance discs are driven here from the statistics in window.ENSO3D.
 * Blender (x, y, z) → glTF/three (x, z, −y). Column height in Blender units = |r| × H. */
(function () {
  'use strict';
  const D = window.ENSO3D, H = 3.2;
  const RDBU = [[-0.7, '#053061'], [-0.35, '#4393c3'], [0, '#f7f7f7'], [0.35, '#d6604d'], [0.7, '#67001f']];
  const state = { period: 'full', show: 'both', year: 1930, sel: 'HUACOL', playing: false, rotate: false };
  const coreById = Object.fromEntries(D.cores.map(c => [c.id, c]));

  // ---------------------------------------------------------------- colour scale
  const lerp = (a, b, t) => a + (b - a) * t;
  function rgb(hex) { const n = parseInt(hex.slice(1), 16); return [(n >> 16) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255]; }
  function colorFor(v, dom = 0.7) {
    const x = Math.max(-0.7, Math.min(0.7, v * 0.7 / dom));
    for (let i = 0; i < RDBU.length - 1; i++) {
      const [a, ca] = RDBU[i], [b, cb] = RDBU[i + 1];
      if (x <= b) { const t = (x - a) / (b - a), A = rgb(ca), B = rgb(cb); return new THREE.Color().setRGB(lerp(A[0], B[0], t), lerp(A[1], B[1], t), lerp(A[2], B[2], t)).convertSRGBToLinear(); }
    }
    return new THREE.Color(RDBU[4][1]).convertSRGBToLinear();
  }
  document.getElementById('cbar').style.background = `linear-gradient(90deg, ${RDBU.map(([, c], i) => `${c} ${i * 25}%`).join(',')})`;

  // ---------------------------------------------------------------- three.js setup
  const stage = document.getElementById('stage');
  const SW = () => stage.clientWidth || innerWidth, SH = () => stage.clientHeight || innerHeight;
  const renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
  renderer.setPixelRatio(Math.min(2, devicePixelRatio)); renderer.setSize(SW(), SH());
  renderer.outputEncoding = THREE.sRGBEncoding; renderer.toneMapping = THREE.ACESFilmicToneMapping; renderer.toneMappingExposure = 1.05;
  stage.appendChild(renderer.domElement);
  const labels = new THREE.CSS2DRenderer(); labels.setSize(SW(), SH());
  labels.domElement.style.position = 'absolute'; labels.domElement.style.inset = '0'; labels.domElement.style.pointerEvents = 'none';
  document.getElementById('labels').appendChild(labels.domElement);
  const scene = new THREE.Scene(); scene.background = new THREE.Color('#0e1116'); scene.fog = new THREE.Fog('#0e1116', 55, 110);
  const camera = new THREE.PerspectiveCamera(38, SW() / SH(), 0.1, 400);
  camera.position.set(-1.5, 27, 38);
  const controls = new THREE.OrbitControls(camera, renderer.domElement);
  controls.target.set(0, -0.8, -0.8); controls.enableDamping = true; controls.maxPolarAngle = Math.PI * 0.49; controls.minDistance = 4; controls.maxDistance = 80;
  scene.add(new THREE.HemisphereLight('#cfd9e6', '#1a1f28', 0.55));
  const sun = new THREE.DirectionalLight('#ffffff', 1.6); sun.position.set(-20, 30, -20); scene.add(sun);
  const fill = new THREE.DirectionalLight('#9fb8d8', 0.35); fill.position.set(20, 15, 25); scene.add(fill);

  // ---------------------------------------------------------------- load the Blender scene
  const bin = Uint8Array.from(atob(window.SCENE_GLB_B64), ch => ch.charCodeAt(0));
  const P = {};              // id → {ONI, RONI, halo, base}
  new THREE.GLTFLoader().parse(bin.buffer, '', gltf => {
    scene.add(gltf.scene);
    gltf.scene.traverse(o => {
      if (!o.isMesh && !o.isObject3D) return;
      const m = o.name.match(/^P_(\w+?)_(ONI|RONI)$/);
      if (m) {
        (P[m[1]] ||= {})[m[2]] = o;
        o.material = new THREE.MeshStandardMaterial({ roughness: 0.3, metalness: 0.05 });
        o.userData = { id: m[1], idx: m[2], h: o.scale.y }; P[m[1]].base = o.position.y;
      }
      const h = o.name.match(/^Halo_(\w+)$/); if (h) { (P[h[1]] ||= {}).halo = o; o.material = new THREE.MeshBasicMaterial({ color: '#cfe3ff', transparent: true, opacity: 0.2, depthWrite: false, side: THREE.DoubleSide }); }
      if (o.name === 'Terrain' && o.material) { o.material.roughness = 0.95; o.material.metalness = 0; }
      if (/^(Borders|Stalk_)/.test(o.name) && o.material) o.material.emissiveIntensity = 0.6;
    });
    addLabels(gltf.scene); addYearMarker();
    document.getElementById('loading').remove();
    setShow(state.show); setPeriod(state.period); select(state.sel); flyTo('overview', true);
  }, err => { document.getElementById('loading').textContent = 'Could not load scene: ' + err; });

  // ---------------------------------------------------------------- labels
  const coreLabels = {};
  function label(text, cls, obj, dy) {
    const div = document.createElement('div'); div.className = cls; div.textContent = text;
    const l = new THREE.CSS2DObject(div); l.position.set(0, dy, 0); obj.add(l); return div;
  }
  function addLabels(root) {
    const bySite = {};
    D.cores.forEach(c => (bySite[c.site] ||= []).push(c));
    for (const [site, list] of Object.entries(bySite)) {
      const s = root.getObjectByName('Site_' + site); if (!s) continue;
      const a = new THREE.Object3D(); a.position.copy(s.position).add(new THREE.Vector3(-0.3 * list.length, 3.2, 0)); root.add(a);
      label(site, 'lbl-site', a, 0);
    }
    D.cores.forEach(c => {
      const p = P[c.id]; if (!p) return;
      const a = new THREE.Object3D(); a.position.set((p.ONI.position.x + p.RONI.position.x) / 2, p.base - 0.15, p.ONI.position.z + 0.45); root.add(a);
      const div = label(c.short, 'lbl-core', a, 0); div.style.pointerEvents = 'auto'; div.style.cursor = 'pointer';
      div.onclick = () => select(c.id); coreLabels[c.id] = div;
    });
    D.other_sites.forEach(s => { const o = root.getObjectByName('Other_' + s.site); if (o) label(`${s.site} ▵`, 'lbl-other', o, 0.45); });
  }

  // year marker on the time wall (Blender: x = −11 + 22 (y−1960)/(2019−1960), plane at y = 15.6 → three z = −15.6)
  let marker;
  function addYearMarker() {
    marker = new THREE.Mesh(new THREE.BoxGeometry(0.08, 4.4, 0.05), new THREE.MeshBasicMaterial({ color: '#7cc4ff' }));
    marker.position.set(0, 3.45, -15.5); scene.add(marker);
    // wall annotations (Blender text objects are not part of the glTF)
    for (let y = 1960; y <= 2020; y += 10) {
      const a = new THREE.Object3D(); a.position.set(-11 + 22 * (y - 1960) / 59, 0.95, -15.5); scene.add(a); label(String(y), 'lbl-other', a, 0);
    }
    const t = new THREE.Object3D(); t.position.set(0, 5.95, -15.5); scene.add(t);
    const d = label('', 'lbl-core', t, 0);
    d.innerHTML = 'Huascarán Col δ¹⁸O <b style="color:#fff">(white)</b> vs DJF <b style="color:#CC79A7">ONI</b> and <b style="color:#009E73">RONI</b>, 1960–2019 · bars: ONI <span style="color:#D55E00">El Niño</span> / <span style="color:#0072B2">La Niña</span> years';
  }

  // ---------------------------------------------------------------- state → scene
  function valuesFor(c) {
    if (state.period === 'full') return { ONI: c.full.ONI.r, RONI: c.full.RONI.r, rc: c.full.rcrit, ok: true };
    const v = c.running[String(state.year)];
    return v ? { ONI: v[0], RONI: v[1], rc: D.meta.rcrit21, ok: true } : { ONI: 0, RONI: 0, rc: D.meta.rcrit21, ok: false };
  }
  const target = {};           // pillar → {h, color}
  function apply() {
    D.cores.forEach(c => {
      const p = P[c.id]; if (!p || !p.ONI) return; const v = valuesFor(c);
      ['ONI', 'RONI'].forEach(k => {
        let h, col, vis = v.ok;
        if (state.show === 'delta') {
          const d = v.RONI - v.ONI; vis = v.ok && k === 'RONI'; h = Math.abs(d) * H * 10; col = colorFor(d, 0.1);
        } else { vis = v.ok && (state.show === 'both' || state.show === k); h = Math.abs(v[k]) * H; col = colorFor(v[k]); }
        target[p[k].uuid] = { obj: p[k], h: vis ? Math.max(0.01, h) : 0.001, col };
        if (!p[k].userData.init) { p[k].material.color.copy(col); p[k].userData.init = true; }
      });
      if (p.halo) { p.halo.visible = v.ok && state.show !== 'delta'; p.halo.userData.y = p.base + v.rc * H; }
    });
    marker.visible = state.period === 'run' && state.year >= 1960 && state.year <= 2019;
    marker.position.x = -11 + 22 * (state.year - 1960) / 59;
    updateUI();
  }
  function tick() {
    for (const t of Object.values(target)) {
      const o = t.obj; o.scale.y += (t.h - o.scale.y) * 0.18;
      o.material.color.lerp(t.col, 0.2); o.material.emissive.copy(o.material.color).multiplyScalar(state.sel === o.userData.id ? 0.55 : 0.22);
    }
    Object.values(P).forEach(p => { if (p.halo && p.halo.userData.y !== undefined) p.halo.position.y += (p.halo.userData.y - p.halo.position.y) * 0.18; });
    if (flight) {
      flight.t = Math.min(1, flight.t + 0.035); const e = 1 - Math.pow(1 - flight.t, 3);
      camera.position.lerpVectors(flight.p0, flight.p1, e); controls.target.lerpVectors(flight.t0, flight.t1, e); if (flight.t >= 1) flight = null;
    }
    controls.autoRotate = state.rotate; controls.autoRotateSpeed = 0.6; controls.update();
    if (!document.getElementById('explore').hidden) { renderer.render(scene, camera); labels.render(scene, camera); }
    requestAnimationFrame(tick);
  }
  requestAnimationFrame(tick);

  // ---------------------------------------------------------------- camera presets
  let flight = null;
  const CAMS = {
    overview: [[-1.5, 27, 38], [0, -0.8, -0.8]],
    huascaran: [[-11.5, 9, 10.5], [-6, 1.2, -0.9]],
    south: [[-1.0, 13, 20], [0.2, 1.0, 4.6]],
    wall: [[0, 5.5, 4], [0, 3.6, -15.6]]
  };
  function flyTo(name, instant) {
    const [p, t] = CAMS[name];
    flight = { t: instant ? 1 : 0, p0: camera.position.clone(), p1: new THREE.Vector3(...p), t0: controls.target.clone(), t1: new THREE.Vector3(...t) };
    if (instant) { camera.position.set(...p); controls.target.set(...t); flight = null; }
  }
  document.querySelectorAll('[data-cam]').forEach(b => b.onclick = () => flyTo(b.dataset.cam));

  // ---------------------------------------------------------------- picking
  const ray = new THREE.Raycaster(), mouse = new THREE.Vector2(), tip = document.getElementById('tip');
  function pick(ev) {
    const rc = renderer.domElement.getBoundingClientRect(); mouse.set((ev.clientX - rc.left) / rc.width * 2 - 1, -(ev.clientY - rc.top) / rc.height * 2 + 1); ray.setFromCamera(mouse, camera);
    const objs = Object.values(P).flatMap(p => [p.ONI, p.RONI]).filter(o => o && o.scale.y > 0.01);
    const hit = ray.intersectObjects(objs, false)[0]; return hit && hit.object;
  }
  renderer.domElement.addEventListener('pointermove', ev => {
    const o = pick(ev); renderer.domElement.style.cursor = o ? 'pointer' : '';
    if (!o) { tip.style.opacity = 0; return; }
    const c = coreById[o.userData.id], v = valuesFor(c);
    tip.innerHTML = state.show === 'delta' ? `<b>${c.short}</b> ΔR = ${f(v.RONI - v.ONI, 3)}` : `<b>${c.short}</b> ${o.userData.idx}: r = ${f(v[o.userData.idx])}`;
    tip.style.opacity = 1; tip.style.left = ev.clientX + 14 + 'px'; tip.style.top = ev.clientY + 10 + 'px';
  });
  let down; renderer.domElement.addEventListener('pointerdown', ev => down = [ev.clientX, ev.clientY]);
  renderer.domElement.addEventListener('pointerup', ev => {
    if (!down || Math.hypot(ev.clientX - down[0], ev.clientY - down[1]) > 4) return;
    const o = pick(ev); if (o) select(o.userData.id);
  });

  // ---------------------------------------------------------------- UI
  const f = (v, d = 2) => (v == null || !isFinite(v)) ? '—' : (v < 0 ? '−' : '') + Math.abs(v).toFixed(d);
  const fp = p => p < 0.001 ? '< 0.001' : '= ' + p.toFixed(3);
  function seg(id, fn) { document.querySelectorAll(`#${id} button`).forEach(b => b.onclick = () => fn(b.dataset.v)); }
  function mark(id, v) { document.querySelectorAll(`#${id} button`).forEach(b => b.classList.toggle('on', b.dataset.v === v)); }
  function setPeriod(v) { state.period = v; mark('period', v); document.getElementById('runBox').style.display = v === 'run' ? '' : 'none'; if (v !== 'run') stop(); apply(); }
  function setShow(v) { state.show = v; mark('show', v); apply(); }
  seg('period', setPeriod); seg('show', setShow);
  const yr = document.getElementById('yr');
  yr.oninput = () => { state.year = +yr.value; apply(); };
  let timer = null;
  function stop() { clearInterval(timer); timer = null; document.getElementById('play').textContent = '▶ Play'; }
  document.getElementById('play').onclick = () => {
    if (timer) return stop();
    document.getElementById('play').textContent = '❚❚ Pause';
    timer = setInterval(() => { state.year = state.year >= 2015 ? 1861 : state.year + 1; yr.value = state.year; apply(); }, 180);
  };
  document.getElementById('rotate').onclick = e => { state.rotate = !state.rotate; e.target.classList.toggle('on', state.rotate); e.target.textContent = state.rotate ? 'Stop rotation' : 'Auto-rotate'; };
  document.getElementById('shot').onclick = () => { renderer.render(scene, camera); Site.exportPanel({ kind: 'image', src: renderer.domElement.toDataURL('image/png'), name: `andes_enso_${state.period}_${state.show}${state.period === 'run' ? '_' + state.year : ''}.png` }); };
  function resize() { camera.aspect = SW() / SH(); camera.updateProjectionMatrix(); renderer.setSize(SW(), SH()); labels.setSize(SW(), SH()); }
  addEventListener('resize', resize);
  if (window.Site) Site.onTab(t => { if (t === 'explore') requestAnimationFrame(resize); else stop(); });
  addEventListener('keydown', e => {
    if ((window.Site && Site.current !== 'explore') || /INPUT|SELECT|TEXTAREA/.test(document.activeElement.tagName)) return;
    if (e.key === ' ') { e.preventDefault(); if (state.period !== 'run') setPeriod('run'); document.getElementById('play').click(); }
    if (e.key === 'ArrowRight' || e.key === 'ArrowLeft') { if (state.period !== 'run') setPeriod('run'); state.year = Math.max(1861, Math.min(2015, state.year + (e.key === 'ArrowRight' ? 1 : -1))); yr.value = state.year; apply(); }
  });

  function select(id) {
    state.sel = id; Object.entries(coreLabels).forEach(([k, d]) => d.classList.toggle('sel', k === id)); renderCard(); apply();
  }

  function updateUI() {
    document.getElementById('year').textContent = state.year;
    document.getElementById('yearNote').textContent = `window ${state.year - 10}–${state.year + 10} · ERSST.v6`;
    const delta = state.show === 'delta';
    document.getElementById('cbarLab').textContent = delta ? 'ΔR = r(RONI) − r(ONI)  (column height ×10)' : 'Pearson r (proxy vs DJF index)';
    document.getElementById('cticks').innerHTML = delta ? '<span>−0.10</span><span>0</span><span>+0.10</span>' : '<span>−0.7</span><span>0</span><span>+0.7</span>';
    document.getElementById('key').innerHTML = delta ? 'One column per core: height = |ΔR| × 10, colour = sign (red: RONI couples more strongly; blue: less).'
      : `Column height = |r|, colour = r. <span class="sw" style="background:var(--oni)"></span><b>ONI</b> cap (west) · <span class="sw" style="background:var(--roni)"></span><b>RONI</b> cap (east). Translucent disc = |r| needed for p < 0.05 (${state.period === 'run' ? 'n = 21' : 'N_eff of each record'}). Proxy: δ¹⁸O (Illimani 1999: NH₄⁺).`;
    renderStory(); renderCard();
  }

  function renderStory() {
    const el = document.getElementById('story'), rows = D.cores.map(c => ({ c, v: valuesFor(c) })).filter(o => o.v.ok);
    if (!rows.length) { el.innerHTML = '<p>No record overlaps this 21-yr window.</p>'; return; }
    const sig = rows.filter(o => Math.max(Math.abs(o.v.ONI), Math.abs(o.v.RONI)) >= o.v.rc);
    const dr = rows.map(o => o.v.RONI - o.v.ONI), md = dr.slice().sort((a, b) => a - b)[Math.floor(dr.length / 2)];
    const best = rows.slice().sort((a, b) => Math.abs(b.v.ONI) - Math.abs(a.v.ONI))[0];
    const when = state.period === 'full' ? 'Over each record’s overlap with 1950–present' : `In the 21-yr window centred on ${state.year}`;
    let t = `<p><b>${when}</b>, ${sig.length} of ${rows.length} records ${sig.length === 1 ? 'is' : 'are'} significantly coupled to ENSO (column above its disc); the strongest is ${best.c.name} (r = ${f(best.v.ONI)} with ONI).</p>`;
    t += `<p>Switching ONI → RONI changes r by ${f(Math.min(...dr))} to ${f(Math.max(...dr))} (median ${f(md)}). ` +
      (Math.abs(md) < 0.03 ? 'The two indices rank the cores identically here: removing the tropical-mean SST barely alters the ice-core ENSO signal.' :
        md < 0 ? 'Coupling weakens under RONI, as expected if the ice also records basin-wide tropical warmth.' :
          'Coupling strengthens under RONI: the cores follow the zonal Pacific SST gradient more than basin-wide warmth.') + '</p>';
    if (state.period === 'run') {
      const q = coreById.QUE13.running[String(state.year)];
      if (q) t += `<p class="muted">Quelccaya’s coupling ${Math.abs(q[0]) >= D.meta.rcrit21 ? 'is significant in this window' : 'is not significant in this window'} (r = ${f(q[0])}); it was strongest around 1870–1920 and weakened after the 1920s.</p>`;
    }
    el.innerHTML = t;
  }

  function spark(c) {
    const W = 300, Hh = 90, ys = Object.keys(c.running).map(Number); if (!ys.length) return '<div class="muted">Record too short for 21-yr running correlations.</div>';
    const x = y => 8 + (W - 16) * (y - 1861) / (2015 - 1861), yv = v => Hh / 2 - v * (Hh / 2 - 6);
    const line = k => ys.map((y, i) => `${i ? 'L' : 'M'}${x(y).toFixed(1)},${yv(c.running[y][k]).toFixed(1)}`).join('');
    const rc = D.meta.rcrit21;
    return `<svg viewBox="0 0 ${W} ${Hh}" width="100%" style="display:block;margin-top:4px">
      <rect x="0" y="0" width="${W}" height="${Hh}" fill="rgba(255,255,255,.03)"/>
      <line x1="8" x2="${W - 8}" y1="${yv(0)}" y2="${yv(0)}" stroke="#555"/>
      <line x1="8" x2="${W - 8}" y1="${yv(rc)}" y2="${yv(rc)}" stroke="#777" stroke-dasharray="3 3"/>
      <line x1="8" x2="${W - 8}" y1="${yv(-rc)}" y2="${yv(-rc)}" stroke="#777" stroke-dasharray="3 3"/>
      <path d="${line(0)}" fill="none" stroke="#CC79A7" stroke-width="1.6"/><path d="${line(1)}" fill="none" stroke="#009E73" stroke-width="1.6"/>
      ${state.period === 'run' ? `<line x1="${x(state.year)}" x2="${x(state.year)}" y1="0" y2="${Hh}" stroke="#7cc4ff"/>` : ''}
      ${[1880, 1920, 1960, 2000].map(y => `<text x="${x(y)}" y="${Hh - 3}" fill="#8a94a2" font-size="9" text-anchor="middle">${y}</text>`).join('')}
      <text x="10" y="11" fill="#8a94a2" font-size="9">running r (21 yr) · dashed: p = 0.05</text></svg>`;
  }
  function renderCard() {
    const c = coreById[state.sel]; if (!c) return; const F = c.full, v = valuesFor(c);
    document.getElementById('card').innerHTML = `
      <div class="lab" style="margin-top:0">Selected record</div>
      <div style="font-weight:700;font-size:14.5px">${c.name}</div>
      <div class="muted">${Math.abs(c.lat).toFixed(2)}°S, ${Math.abs(c.lon).toFixed(2)}°W · ${c.elev} m a.s.l. · ${c.proxy} · ${c.years[0]}–${c.years[1]}</div>
      <div class="muted" style="margin-top:6px">Full period vs CPC DJF index, ${F.y0}–${F.y1} (n = ${F.n})</div>
      <table class="st"><tr><th></th><th style="color:var(--oni)">ONI</th><th style="color:var(--roni)">RONI</th></tr>
        <tr><td>Pearson r</td><td>${f(F.ONI.r)}</td><td>${f(F.RONI.r)}</td></tr>
        <tr><td>p (N_eff)</td><td>${fp(F.ONI.p)}</td><td>${fp(F.RONI.p)}</td></tr>
        <tr><td>N_eff</td><td>${F.ONI.neff}</td><td>${F.RONI.neff}</td></tr>
        <tr><td>ΔR = RONI − ONI</td><td colspan="2">${f(F.dR, 3)} (MRR p ${fp(F.p_dR)})</td></tr>
        ${state.period === 'run' ? `<tr><td>21-yr window ${state.year}</td><td>${v.ok ? f(v.ONI) : '—'}</td><td>${v.ok ? f(v.RONI) : '—'}</td></tr>` : ''}
      </table>
      ${spark(c)}
      <div class="muted" style="margin-top:6px">${c.ref}. ${c.id === 'QUE13' ? 'Near-zero r since 1950 reflects a post-1920s breakdown in coupling (see running r).' : ''}</div>
      <div class="muted" style="margin-top:8px">Full statistics, composites and pre-1950 reconstructions: <a href="#analysis" data-tab="analysis">Analysis tab</a> · <a href="#about" data-tab="about">Blender animation</a></div>`;
  }
})();
