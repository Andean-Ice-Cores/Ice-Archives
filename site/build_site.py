"""Assemble the website (site/dist/) from the 2-D app, the Blender viewer and the shell in site/src/.

The sandboxed host blocks file downloads, so export buttons open an in-page export panel instead,
and the two apps are scoped (CSS under #analysis / #explore) to live on one page as tabs.
Run: python3 site/build_site.py
"""
import os, re, shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, 'site', 'src')
DIST = os.path.join(ROOT, 'site', 'dist')
shutil.rmtree(DIST, ignore_errors=True)
for d in ('vendor', 'data', 'js', 'css', 'media'):
    os.makedirs(os.path.join(DIST, d))
cp = lambda a, b: shutil.copy(os.path.join(ROOT, a), os.path.join(DIST, b))

# ---- libraries, data, media
cp('app/vendor/d3.min.js', 'vendor/d3.min.js')
for f in ('three.min.js', 'GLTFLoader.js', 'OrbitControls.js', 'CSS2DRenderer.js'):
    cp(f'blender/web/vendor/{f}', f'vendor/{f}')
cp('app/data/dataset.js', 'data/dataset.js')
cp('blender/web/scene_data.js', 'data/scene_data.js')
cp('app/js/stats.js', 'js/stats.js')
cp('blender/out/running.mp4', 'media/running.mp4')
cp('blender/out/hero_full.png', 'media/hero_full.png')
shutil.copy(os.path.join(SRC, 'shell.js'), os.path.join(DIST, 'js/shell.js'))
shutil.copy(os.path.join(SRC, 'site.css'), os.path.join(DIST, 'css/site.css'))
shutil.copy(os.path.join(SRC, 'explore.css'), os.path.join(DIST, 'css/explore.css'))
shutil.copy(os.path.join(SRC, 'compare.css'), os.path.join(DIST, 'css/compare.css'))
shutil.copy(os.path.join(SRC, 'compare.js'), os.path.join(DIST, 'js/compare.js'))


def patch(s, old, new, count=1):
    assert old in s, 'patch target not found: ' + old[:70]
    return s.replace(old, new, count)


# ---- 2-D app
a = open(os.path.join(ROOT, 'app/js/app.js')).read()
a = re.sub(r"  function download\(blob, name\) \{.*?\n", """  function download(blob, name) {       // hosted build: show in the export panel (downloads are sandboxed)
    if (/^image\\//.test(blob.type)) { const r = new FileReader(); r.onload = () => Site.exportPanel({ kind: 'image', src: r.result, name }); r.readAsDataURL(blob); }
    else blob.text().then(text => Site.exportPanel({ kind: 'text', text, name }));
  }
""", a, count=1, flags=re.S)
a = patch(a, "  window.__setState = (o) => { Object.assign(state, o); update(); };",
          "  window.__setState = (o) => { Object.assign(state, o); update(); };\n  if (window.Site) Site.onTab(t => { if (t === 'analysis') requestAnimationFrame(update); });")
open(os.path.join(DIST, 'js/app.js'), 'w').write(a)

# ---- 3-D viewer: size to its container, only render/keys while visible, export panel, tab links
v = open(os.path.join(ROOT, 'blender/web/viewer.js')).read()
v = patch(v, "  const stage = document.getElementById('stage');",
          "  const stage = document.getElementById('stage');\n  const SW = () => stage.clientWidth || innerWidth, SH = () => stage.clientHeight || innerHeight;")
v = patch(v, "renderer.setSize(innerWidth, innerHeight);\n", "renderer.setSize(SW(), SH());\n")
v = patch(v, "labels.setSize(innerWidth, innerHeight);\n  labels.domElement.style.position = 'fixed';",
          "labels.setSize(SW(), SH());\n  labels.domElement.style.position = 'absolute';")
v = patch(v, "new THREE.PerspectiveCamera(38, innerWidth / innerHeight, 0.1, 400)", "new THREE.PerspectiveCamera(38, SW() / SH(), 0.1, 400)")
v = patch(v, "mouse.set(ev.clientX / innerWidth * 2 - 1, -ev.clientY / innerHeight * 2 + 1);",
          "const rc = renderer.domElement.getBoundingClientRect(); mouse.set((ev.clientX - rc.left) / rc.width * 2 - 1, -(ev.clientY - rc.top) / rc.height * 2 + 1);")
v = re.sub(r"  addEventListener\('resize', \(\) => \{.*?\}\);\n",
           "  function resize() { camera.aspect = SW() / SH(); camera.updateProjectionMatrix(); renderer.setSize(SW(), SH()); labels.setSize(SW(), SH()); }\n"
           "  addEventListener('resize', resize);\n  if (window.Site) Site.onTab(t => { if (t === 'explore') requestAnimationFrame(resize); else stop(); });\n", v, count=1, flags=re.S)
v = patch(v, "  addEventListener('keydown', e => {\n", "  addEventListener('keydown', e => {\n    if ((window.Site && Site.current !== 'explore') || /INPUT|SELECT|TEXTAREA/.test(document.activeElement.tagName)) return;\n")
v = patch(v, "    renderer.render(scene, camera); labels.render(scene, camera); requestAnimationFrame(tick);",
          "    if (!document.getElementById('explore').hidden) { renderer.render(scene, camera); labels.render(scene, camera); }\n    requestAnimationFrame(tick);")
v = re.sub(r"  document.getElementById\('shot'\).onclick = \(\) => \{.*?\};\n",
           "  document.getElementById('shot').onclick = () => { renderer.render(scene, camera); Site.exportPanel({ kind: 'image', src: renderer.domElement.toDataURL('image/png'), name: `andes_enso_${state.period}_${state.show}${state.period === 'run' ? '_' + state.year : ''}.png` }); };\n", v, count=1, flags=re.S)
v = patch(v, '<a href="../../app/index.html">2-D analysis app</a> · <a href="../out/running.mp4">Blender animation</a>',
          '<a href="#analysis" data-tab="analysis">Analysis tab</a> · <a href="#about" data-tab="about">Blender animation</a>')
open(os.path.join(DIST, 'js/viewer.js'), 'w').write(v)

# ---- 2-D CSS scoped under #analysis
css = open(os.path.join(ROOT, 'app/css/style.css')).read()


def scope_sel(sel):
    sel = sel.strip()
    if sel in (':root', 'body', 'html'):
        return '#analysis'
    if sel == '*':
        return '#analysis, #analysis *'
    return '#analysis ' + sel


def scope_block(text):
    out, i = [], 0
    while i < len(text):
        j = text.find('{', i)
        if j < 0: break
        head = text[i:j].strip()
        if head.startswith('@media'):
            depth, k = 1, j + 1
            while depth:
                depth += {'{': 1, '}': -1}.get(text[k], 0); k += 1
            out.append(head + ' {\n' + scope_block(text[j + 1:k - 1]) + '}\n'); i = k
        else:
            k = text.find('}', j)
            out.append(', '.join(scope_sel(s) for s in head.split(',')) + ' ' + text[j:k + 1] + '\n'); i = k + 1
    return ''.join(out)


scoped = scope_block(re.sub(r'/\*.*?\*/', '', css, flags=re.S))
scoped += """
#analysis { color-scheme: light; background: var(--paper); color: var(--ink); min-height: calc(100vh - var(--barH)); }
#analysis .controls { top: calc(env(safe-area-inset-top, 0px) + var(--barH)); }
#analysis a { color: var(--accent); }
#analysis a.load-btn { display: inline-block; font: 600 12.5px var(--sans); color: #fff; background: var(--accent); border-radius: 4px; padding: 5px 10px; text-decoration: none; }
"""
open(os.path.join(DIST, 'css/analysis.css'), 'w').write(scoped)

# ---- page
app_html = open(os.path.join(ROOT, 'app/index.html')).read()
body = app_html[app_html.index('<header class="masthead">'):app_html.index('</main>') + len('</main>')]
body = body.replace('<button class="ghost" id="btn-fullwin" title="Use full overlap">Full</button>\n  </div>',
                    '<button class="ghost" id="btn-fullwin" title="Use full overlap">Full</button>\n  </div>\n  <div class="ctl"><label>Your own record</label><a href="#compare" data-tab="compare" class="load-btn">⤒ Load your data</a></div>', 1)
body += '\n<div id="tooltip" class="tooltip" role="status"></div>\n'
shell = open(os.path.join(SRC, 'shell.html')).read()
page = shell.replace('<!--ANALYSIS-->', body).replace('<!--COMPARE-->', open(os.path.join(SRC, 'compare.html')).read())
open(os.path.join(DIST, 'index.html'), 'w').write(page)          # fragment for claude.ai (host adds the skeleton)
# standalone copy for GitHub Pages: python3 site/build_site.py --pages docs
import sys
if '--pages' in sys.argv:
    dest = os.path.join(ROOT, sys.argv[sys.argv.index('--pages') + 1])
    shutil.rmtree(dest, ignore_errors=True); shutil.copytree(DIST, dest)
    open(os.path.join(dest, 'index.html'), 'w').write(
        '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
        + page.split('<header class="bar">')[0] + '</head>\n<body>\n<header class="bar">'
        + page.split('<header class="bar">', 1)[1] + '\n</body>\n</html>\n')
    open(os.path.join(dest, '.nojekyll'), 'w').write('')
    print('GitHub Pages site written to', dest)
total = sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(DIST) for f in fs)
print('site/dist built:', round(total / 1e6, 1), 'MB')
