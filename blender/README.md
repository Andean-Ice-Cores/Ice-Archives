# 3-D infographic (Blender)

A Blender-built 3-D version of the ONI-vs-RONI ice-core figure, with an interactive web viewer.

## Open it

* **Interactive:** double-click `blender/web/index.html` (offline; three.js is vendored and the Blender scene is embedded).
  Orbit/zoom with the mouse; click a column or core label for its statistics and running-correlation sparkline;
  switch **Full period ↔ Running 21-yr** (year slider, ▶ Play, or ←/→ and space); show **ONI + RONI / ONI / RONI / ΔR**;
  camera presets (Overview, Huascarán, Quelccaya·Illimani, Time wall); PNG snapshot.
* **Blender:** open `blender/out/andes_enso.blend` in Blender ≥ 3.6. Scrub the timeline: frame *f* = running-window
  centre year 1860 + *f* (1861–2015). Allow auto-run scripts to get the live year counter.
* **Renders:** `blender/out/hero_full.png` (1920×1080, full-period correlations) and `blender/out/running.mp4`
  (1861–2015 running correlations, 13 s).

## What is encoded

* Terrain: ETOPO1, 0.1°, vertical exaggeration ≈ 12× (0.32 units km⁻¹, 1 unit = 1°).
* Each core = two columns: **ONI (pink cap, west)** and **RONI (green cap, east)**; height = |r|, colour = r (RdBu, ±0.7).
* Translucent disc = |r| needed for p < 0.05 (full period: from each record's N_eff; running: n = 21 → |r| ≥ 0.433).
* Proxy = δ¹⁸O except Illimani 1999 (NH₄⁺). Full period = official CPC DJF ONI/RONI over each core's overlap with
  1950–present; running = ERSST.v6 DJF ONI/RONI recomputed with CPC's algorithm (validated r = 0.999/0.9994 vs CPC).
* Time wall (north): Huascarán Col δ¹⁸O (z) vs DJF ONI and RONI 1960–2019, with ONI El Niño/La Niña year bars.

Numbers are identical to the 2-D app (`app/index.html`, e.g. Huascarán Col r = 0.582 ONI / 0.607 RONI).

## Rebuild

```bash
python3 blender/prepare_data.py
/Applications/Blender.app/Contents/MacOS/Blender -b -P blender/build_scene.py -- --still --glb --video
python3 blender/build_web.py
```
