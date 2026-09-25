"""Build the 3-D ENSO × Andean ice-core infographic in Blender (headless).

  /Applications/Blender.app/Contents/MacOS/Blender -b -P blender/build_scene.py -- [--still] [--video] [--glb]

Outputs (blender/out/):
  andes_enso.blend   scene; scrub the timeline (frame = running-window centre year 1861–2015)
  hero_full.png      still: full-period correlations (CPC 1950–present)
  running.mp4        animation: 21-yr running correlations (ERSST.v6 1850–present)
  andes_enso.glb     geometry for the web viewer (blender/web/index.html)

Coordinates: x = lon + 72, y = lat + 10 (1 unit = 1°), z = elevation × 0.4 per km (≈15× exaggeration).
Each core = two columns (ONI north, RONI south); height = |r| × H, colour = r on a diverging RdBu scale.
The translucent disc marks the |r| needed for p < 0.05.
"""
import bpy, bmesh, json, math, os, sys
from mathutils import Vector

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'blender', 'out'); os.makedirs(OUT, exist_ok=True)
ARGS = sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else ['--still', '--video', '--glb']
D = json.load(open(os.path.join(ROOT, 'blender', 'enso3d.json')))
_src = open(os.path.join(ROOT, 'app', 'data', 'dataset.js')).read()
MAP = json.loads(_src[_src.index('=') + 1:].strip().rstrip(';'))['map']

ZS, H, Y0 = 0.32, 3.2, 1860          # km→units, pillar height per unit r, frame 1 = 1861
X = lambda lon: lon + 72.0
Yc = lambda lat: lat + 10.0
C = dict(ONI=(0.800, 0.475, 0.655), RONI=(0.0, 0.620, 0.451), EN=(0.835, 0.369, 0.0), LN=(0.0, 0.447, 0.698))
RDBU = [(-0.7, (0.020, 0.188, 0.380)), (-0.35, (0.263, 0.576, 0.765)), (0.0, (0.969, 0.969, 0.969)),
        (0.35, (0.839, 0.376, 0.302)), (0.7, (0.404, 0.0, 0.122))]

bpy.ops.wm.read_factory_settings(use_empty=True)
scene = bpy.context.scene
srgb2lin = lambda c: tuple(((x + 0.055) / 1.055) ** 2.4 if x > 0.04045 else x / 12.92 for x in c)


def mat(name, color, rough=0.6, emit=0.0, alpha=1.0, metal=0.0):
    m = bpy.data.materials.new(name); m.use_nodes = True
    b = m.node_tree.nodes['Principled BSDF']
    b.inputs['Base Color'].default_value = (*srgb2lin(color), 1)
    b.inputs['Roughness'].default_value = rough; b.inputs['Metallic'].default_value = metal
    if emit:
        b.inputs['Emission'].default_value = (*srgb2lin(color), 1); b.inputs['Emission Strength'].default_value = emit
    if alpha < 1:
        b.inputs['Alpha'].default_value = alpha; m.blend_method = 'BLEND'; m.shadow_method = 'NONE'; m.use_backface_culling = False
    return m


def link(obj, col):
    col.objects.link(obj); return obj


col_geo = bpy.data.collections.new('Geometry'); scene.collection.children.link(col_geo)
col_txt = bpy.data.collections.new('Labels'); scene.collection.children.link(col_txt)

# ------------------------------------------------------------------ terrain block
T = D['terrain']; lats, lons, Z = T['lats'], T['lons'], T['z']
nr, nc = len(lats), len(lons)
def zval(e): return e / 1000 * ZS if e > 0 else max(e, -5000) / 1000 * 0.04
def terrain_z(lon, lat):             # bilinear sample
    fx = (lon - lons[0]) / (lons[1] - lons[0]); fy = (lat - lats[0]) / (lats[1] - lats[0])
    i0, j0 = max(0, min(nr - 2, int(fy))), max(0, min(nc - 2, int(fx))); ty, tx = fy - i0, fx - j0
    z = lambda i, j: zval(Z[i][j])
    return (z(i0, j0) * (1 - tx) * (1 - ty) + z(i0, j0 + 1) * tx * (1 - ty) + z(i0 + 1, j0) * (1 - tx) * ty + z(i0 + 1, j0 + 1) * tx * ty)

verts, faces, cols = [], [], []
def hyps(e):
    if e <= 0:
        t = min(1, -e / 4000); return (0.09 - 0.04 * t, 0.15 - 0.06 * t, 0.22 - 0.07 * t)
    stops = [(0, (0.25, 0.27, 0.27)), (1500, (0.40, 0.40, 0.37)), (3500, (0.58, 0.55, 0.50)), (4800, (0.78, 0.78, 0.78)), (5600, (0.97, 0.98, 1.0))]
    for (a, ca), (b, cb) in zip(stops, stops[1:]):
        if e <= b:
            t = (e - a) / (b - a); return tuple(ca[k] + (cb[k] - ca[k]) * t for k in range(3))
    return stops[-1][1]
for i in range(nr):
    for j in range(nc):
        verts.append((X(lons[j]), Yc(lats[i]), zval(Z[i][j]))); cols.append(hyps(Z[i][j]))
for i in range(nr - 1):
    for j in range(nc - 1):
        a = i * nc + j; faces.append((a, a + 1, a + nc + 1, a + nc))
# skirt (diorama sides)
BASE = -0.7
edge = [(0, j) for j in range(nc)] + [(i, nc - 1) for i in range(1, nr)] + [(nr - 1, j) for j in range(nc - 2, -1, -1)] + [(i, 0) for i in range(nr - 2, 0, -1)]
bidx = []
for (i, j) in edge:
    v = verts[i * nc + j]; verts.append((v[0], v[1], BASE)); cols.append((0.05, 0.06, 0.08)); bidx.append((i * nc + j, len(verts) - 1))
for k in range(len(bidx)):
    (t0, b0), (t1, b1) = bidx[k], bidx[(k + 1) % len(bidx)]
    faces.append((t1, t0, b0, b1))
faces.append(tuple(b for _, b in bidx)[::-1])
me = bpy.data.meshes.new('Terrain'); me.from_pydata(verts, [], faces); me.update()
ca = me.color_attributes.new('Col', 'FLOAT_COLOR', 'POINT')
for k, c in enumerate(cols): ca.data[k].color = (*srgb2lin(c), 1)
for p in me.polygons: p.use_smooth = True
terrain = link(bpy.data.objects.new('Terrain', me), col_geo)
tm = bpy.data.materials.new('TerrainMat'); tm.use_nodes = True
nt = tm.node_tree; ab = nt.nodes.new('ShaderNodeVertexColor'); ab.layer_name = 'Col'
bs = nt.nodes['Principled BSDF']; bs.inputs['Roughness'].default_value = 0.92; bs.inputs['Specular'].default_value = 0.2
nt.links.new(ab.outputs['Color'], bs.inputs['Base Color']); me.materials.append(tm)

# ------------------------------------------------------------------ borders & 3000 m contour (thin tubes draped on terrain)
def tube(name, lines, radius, material, lift=0.012, coll=col_geo):
    cu = bpy.data.curves.new(name, 'CURVE'); cu.dimensions = '3D'; cu.bevel_depth = radius; cu.bevel_resolution = 1
    for line in lines:
        pts = [(X(lo), Yc(la), max(terrain_z(lo, la), 0) + lift) for lo, la in line]
        if len(pts) < 2: continue
        sp = cu.splines.new('POLY'); sp.points.add(len(pts) - 1)
        for p, q in zip(sp.points, pts): p.co = (*q, 1)
    ob = link(bpy.data.objects.new(name, cu), coll); cu.materials.append(material)
    return ob
ext = MAP['ext']
def clip(ring):
    out, cur = [], []
    for lo, la in ring:
        if ext['lon0'] <= lo <= ext['lon1'] and ext['lat0'] <= la <= ext['lat1']: cur.append((lo, la))
        elif cur: out.append(cur); cur = []
    if cur: out.append(cur)
    return out
blines = []
for f in MAP['countries']:
    for poly in f['geometry']['coordinates']:
        for ring in poly: blines += clip(ring)
m_border = mat('BorderMat', (0.85, 0.87, 0.9), rough=0.5, emit=0.25)
tube('Borders', blines, 0.012, m_border)

# ------------------------------------------------------------------ pillars
def pillar_mesh(name, radius, seg=40):
    bm = bmesh.new(); bmesh.ops.create_cone(bm, cap_ends=True, segments=seg, radius1=radius, radius2=radius, depth=1.0)
    for v in bm.verts: v.co.z += 0.5
    m = bpy.data.meshes.new(name); bm.to_mesh(m); bm.free()
    for p in m.polygons: p.use_smooth = True
    return m
def ring_mesh(name, r_in, r_out, seg=48, h=0.03):
    bm = bmesh.new(); bmesh.ops.create_cone(bm, cap_ends=True, segments=seg, radius1=r_out, radius2=r_out, depth=h)
    m = bpy.data.meshes.new(name); bm.to_mesh(m); bm.free(); return m

# pillar material: colour from the object's custom property "r" (Attribute node, Object type) → RdBu ramp
pm = bpy.data.materials.new('PillarMat'); pm.use_nodes = True; nt = pm.node_tree
att = nt.nodes.new('ShaderNodeAttribute'); att.attribute_type = 'OBJECT'; att.attribute_name = 'r'
mr = nt.nodes.new('ShaderNodeMapRange'); mr.inputs['From Min'].default_value = -0.7; mr.inputs['From Max'].default_value = 0.7
rp = nt.nodes.new('ShaderNodeValToRGB'); el = rp.color_ramp.elements
el[0].position, el[0].color = 0, (*srgb2lin(RDBU[0][1]), 1); el[1].position, el[1].color = 1, (*srgb2lin(RDBU[-1][1]), 1)
for v, c in RDBU[1:-1]:
    e = el.new((v + 0.7) / 1.4); e.color = (*srgb2lin(c), 1)
bs = nt.nodes['Principled BSDF']; bs.inputs['Roughness'].default_value = 0.25; bs.inputs['Emission Strength'].default_value = 0.35
nt.links.new(att.outputs['Fac'], mr.inputs['Value']); nt.links.new(mr.outputs['Result'], rp.inputs['Fac'])
nt.links.new(rp.outputs['Color'], bs.inputs['Base Color']); nt.links.new(rp.outputs['Color'], bs.inputs['Emission'])

m_halo = mat('HaloMat', (0.85, 0.92, 1.0), rough=0.1, emit=0.8, alpha=0.16)
m_ring = {k: mat(f'Ring{k}', C[k], rough=0.3, emit=1.2) for k in ('ONI', 'RONI')}
m_site = mat('SiteMat', (1, 1, 1), rough=0.3, emit=1.0)
m_other = mat('OtherSiteMat', (0.7, 0.72, 0.75), rough=0.4, emit=0.2)
m_stalk = mat('StalkMat', (0.9, 0.9, 0.95), rough=0.4, emit=0.4)

PILLARS = {}          # id -> {'ONI': obj, 'RONI': obj, 'halo': obj, 'base': z}
sites = {}
for c in D['cores']: sites.setdefault(c['site'], []).append(c)
for site, lst in sites.items():
    sx, sy = X(lst[0]['lon']), Yc(lst[0]['lat']); sz = terrain_z(lst[0]['lon'], lst[0]['lat'])
    dot = link(bpy.data.objects.new(f'Site_{site}', pillar_mesh(f'SiteM_{site}', 0.07, 16)), col_geo)
    dot.location = (sx, sy, sz); dot.scale.z = 0.15; dot.data.materials.append(m_site)
    for i, c in enumerate(lst):
        px, py = sx - 0.75 * (i + 1), sy
        base = max(terrain_z(c['lon'] - 0.75 * (i + 1), c['lat']), sz - 0.6, 0.05)
        # ground stalk from site to pair
        st = tube(f'Stalk_{c["id"]}', [[(c['lon'], c['lat']), (c['lon'] - 0.75 * (i + 1), c['lat'])]], 0.012, m_stalk, lift=0.03)
        entry = {'base': base}
        for k, dx in (('ONI', -0.15), ('RONI', 0.15)):
            ob = link(bpy.data.objects.new(f'P_{c["id"]}_{k}', pillar_mesh(f'PM_{c["id"]}_{k}', 0.12)), col_geo)
            ob.location = (px + dx, py, base); ob.data.materials.append(pm)
            r = c['full'][k]['r']; ob['r'] = r; ob.scale.z = max(0.01, abs(r)) * H
            cap = link(bpy.data.objects.new(f'Cap_{c["id"]}_{k}', ring_mesh(f'CM_{c["id"]}_{k}', 0, 0.125, h=0.02)), col_geo)
            cap.parent = ob; cap.location = (0, 0, 1.0); cap.data.materials.append(m_ring[k])
            rg = link(bpy.data.objects.new(f'R_{c["id"]}_{k}', ring_mesh(f'RM_{c["id"]}_{k}', 0.12, 0.17)), col_geo)
            rg.location = (px + dx, py, base + 0.015); rg.data.materials.append(m_ring[k])
            entry[k] = ob
        # significance halo at |r| for p = 0.05 given N_eff of the full period
        rc = c['full']['rcrit']
        hl = link(bpy.data.objects.new(f'Halo_{c["id"]}', ring_mesh(f'HM_{c["id"]}', 0, 0.34, h=0.01)), col_geo)
        hl.location = (px, py, base + rc * H); hl.data.materials.append(m_halo); hl['rcrit_full'] = rc
        entry['halo'] = hl; PILLARS[c['id']] = entry
for s in D['other_sites']:
    ob = bpy.data.objects.new(f'Other_{s["site"]}', None)
    bm = bmesh.new(); bmesh.ops.create_cone(bm, cap_ends=True, segments=3, radius1=0.16, radius2=0, depth=0.3)
    m = bpy.data.meshes.new(f'OtherM_{s["site"]}'); bm.to_mesh(m); bm.free()
    ob = link(bpy.data.objects.new(f'Other_{s["site"]}', m), col_geo)
    ob.location = (X(s['lon']), Yc(s['lat']), terrain_z(s['lon'], s['lat']) + 0.15); m.materials.append(m_other)

# ------------------------------------------------------------------ time wall: Huascarán Col δ18O vs ONI/RONI 1960–2019
W = D['wall']; yrs = W['years']; WY = 15.6; WZ0 = 3.4; WS = 0.62; WL = 1.0
wx = lambda y: -11 + 22 * (y - yrs[0]) / (yrs[-1] - yrs[0])
panel = bpy.data.meshes.new('WallPanel'); panel.from_pydata([(-11.6, WY + 0.05, 0.2 + WL), (11.6, WY + 0.05, 0.2 + WL), (11.6, WY + 0.05, 4.7 + WL), (-11.6, WY + 0.05, 4.7 + WL)], [], [(0, 1, 2, 3)])
link(bpy.data.objects.new('WallPanel', panel), col_geo).data.materials.append(mat('WallMat', (0.07, 0.08, 0.1), rough=0.9, alpha=0.85))
m_en, m_ln = mat('ENMat', C['EN'], emit=0.5, alpha=0.45), mat('LNMat', C['LN'], emit=0.5, alpha=0.45)
bw = 22 / (len(yrs) - 1)
for y, cl in zip(yrs, W['clsONI']):
    if cl == 0: continue
    bm = bmesh.new(); bmesh.ops.create_cube(bm, size=1); m = bpy.data.meshes.new(f'EvM_{y}'); bm.to_mesh(m); bm.free()
    ob = link(bpy.data.objects.new(f'Ev_{y}', m), col_geo)
    ob.location = (wx(y), WY + 0.02, (0.35 + 4.55) / 2 + WL); ob.scale = (bw * 0.92, 0.02, 4.55 - 0.35)
    m.materials.append(m_en if cl > 0 else m_ln)
def _series(name, vals, color, radius, emit):
    cu = bpy.data.curves.new(name, 'CURVE'); cu.dimensions = '3D'; cu.bevel_depth = radius; cu.bevel_resolution = 2
    sp = cu.splines.new('POLY'); sp.points.add(len(yrs) - 1)
    for p, y, v in zip(sp.points, yrs, vals): p.co = (wx(y), WY - 0.05, WZ0 + WS * v, 1)
    ob = link(bpy.data.objects.new(name, cu), col_geo); cu.materials.append(mat(name + 'Mat', color, rough=0.3, emit=emit)); return ob
_series('WallONI', W['ONI'], C['ONI'], 0.035, 1.5)
_series('WallRONI', W['RONI'], C['RONI'], 0.035, 1.5)
_series('WallProxy', W['proxy'], (1, 1, 1), 0.05, 2.0)
zl = bpy.data.meshes.new('WallZero'); zl.from_pydata([(-11.2, WY - 0.03, WZ0 - 0.006), (11.2, WY - 0.03, WZ0 - 0.006), (11.2, WY - 0.03, WZ0 + 0.006), (-11.2, WY - 0.03, WZ0 + 0.006)], [], [(0, 1, 2, 3)])
link(bpy.data.objects.new('WallZero', zl), col_geo).data.materials.append(mat('ZeroMat', (0.6, 0.6, 0.6), emit=0.3))

# ------------------------------------------------------------------ labels (Blender renders only; the web viewer uses HTML labels)
FONT = None
for fp in ('/System/Library/Fonts/Supplemental/Arial Unicode.ttf', '/System/Library/Fonts/Supplemental/Arial.ttf'):
    if os.path.exists(fp): FONT = bpy.data.fonts.load(fp); break
m_txt = mat('TextMat', (0.95, 0.96, 0.98), rough=0.5, emit=1.2)
m_txt2 = mat('TextMat2', (0.7, 0.73, 0.78), rough=0.5, emit=0.6)
def text(name, body, loc, size, material=m_txt, align='CENTER', rot=(math.radians(62), 0, 0)):
    cu = bpy.data.curves.new(name, 'FONT'); cu.body = body; cu.size = size; cu.align_x = align
    if FONT: cu.font = FONT
    ob = link(bpy.data.objects.new(name, cu), col_txt); ob.location = loc; ob.rotation_euler = rot; cu.materials.append(material)
    return ob
for site, lst in sites.items():
    text(f'TXT_{site}', site, (X(lst[0]['lon']) - 0.4 * len(lst), Yc(lst[0]['lat']) + 0.6, 3.9), 0.6)
    for i, c in enumerate(lst):
        text(f'TXT_{c["id"]}', c['short'], (X(c['lon']) - 0.75 * (i + 1), Yc(c['lat']) - 0.6, PILLARS[c['id']]['base'] + 0.05), 0.34, m_txt2)
for s in D['other_sites']:
    text(f'TXT_o_{s["site"]}', s['site'], (X(s['lon']) + 0.4, Yc(s['lat']), terrain_z(s['lon'], s['lat']) + 0.4), 0.32, m_txt2, 'LEFT')
up = (math.radians(90), 0, 0)
text('TXT_wall', 'Huascarán Col δ18O (white) vs DJF ONI (pink) and RONI (green), 1960–2019 · bars: ONI El Niño / La Niña years', (0, WY - 0.1, WL - 0.55), 0.36, m_txt2, rot=up)
for y in range(1960, 2021, 10):
    text(f'TXT_w{y}', str(y), (wx(y), WY - 0.1, WL - 0.05), 0.34, m_txt2, rot=up)

# ------------------------------------------------------------------ camera, lights, world, render settings
cam = bpy.data.objects.new('Camera', bpy.data.cameras.new('Camera')); scene.collection.objects.link(cam); scene.camera = cam
cam.data.lens = 36; cam.location = (-1.5, -38.0, 27.0)
tgt = bpy.data.objects.new('CamTarget', None); scene.collection.objects.link(tgt); tgt.location = (0.0, 0.8, -0.8)
tc = cam.constraints.new('TRACK_TO'); tc.target = tgt; tc.track_axis = 'TRACK_NEGATIVE_Z'; tc.up_axis = 'UP_Y'
for ob in col_txt.objects:
    if ob.rotation_euler[0] == math.radians(62):
        c2 = ob.constraints.new('DAMPED_TRACK'); c2.target = cam; c2.track_axis = 'TRACK_Z'
sun = bpy.data.objects.new('Sun', bpy.data.lights.new('Sun', 'SUN')); scene.collection.objects.link(sun)
sun.data.energy = 3.2; sun.rotation_euler = (math.radians(50), 0, math.radians(135)); sun.data.angle = math.radians(3)
fill = bpy.data.objects.new('Fill', bpy.data.lights.new('Fill', 'AREA')); scene.collection.objects.link(fill)
fill.data.energy = 4000; fill.data.size = 30; fill.location = (10, -25, 25); fill.rotation_euler = (math.radians(45), 0, math.radians(20))
world = bpy.data.worlds.new('World'); scene.world = world; world.use_nodes = True
world.node_tree.nodes['Background'].inputs['Color'].default_value = (*srgb2lin((0.055, 0.066, 0.086)), 1)
world.node_tree.nodes['Background'].inputs['Strength'].default_value = 0.6
scene.render.engine = 'BLENDER_EEVEE'
ee = scene.eevee; ee.use_bloom = True; ee.bloom_intensity = 0.035; ee.use_gtao = True; ee.gtao_distance = 1.2
ee.use_soft_shadows = True; ee.shadow_cascade_size = '4096'; ee.taa_render_samples = 64
scene.view_settings.view_transform = 'Filmic'; scene.view_settings.look = 'Medium High Contrast'
scene.render.resolution_x, scene.render.resolution_y = 1920, 1080

# HUD title & year counter (children of the camera)
def hud(name, body, xy, size, material=m_txt, align='LEFT'):
    ob = text(name, body, (0, 0, 0), size, material, align, rot=(0, 0, 0)); ob.parent = cam
    ob.location = (xy[0], xy[1], -10); return ob
hud('HUD_title', 'ENSO coupling of Andean ice cores: ONI vs RONI', (-4.75, 2.5), 0.2)
hud_sub = hud('HUD_sub', 'Full period · NOAA CPC DJF indices, 1950–present', (-4.75, 2.28), 0.1, m_txt2)
hud('HUD_key', 'column height = |r|, colour = r  ·  pink cap = ONI, green cap = RONI  ·  translucent disc = |r| needed for p < 0.05', (-4.75, -2.1), 0.085, m_txt2)
cb = bpy.data.meshes.new('HUD_cbar'); cb.from_pydata([(0, 0, 0), (2.4, 0, 0), (2.4, 0.09, 0), (0, 0.09, 0)], [], [(0, 1, 2, 3)])
cbo = bpy.data.objects.new('HUD_cbar', cb); col_txt.objects.link(cbo); cbo.parent = cam; cbo.location = (-4.75, -2.42, -10)
cm = bpy.data.materials.new('CbarMat'); cm.use_nodes = True; n2 = cm.node_tree
tcn = n2.nodes.new('ShaderNodeTexCoord'); sep = n2.nodes.new('ShaderNodeSeparateXYZ'); rp2 = n2.nodes.new('ShaderNodeValToRGB')
rp2.color_ramp.elements[0].color = (*srgb2lin(RDBU[0][1]), 1); rp2.color_ramp.elements[1].color = (*srgb2lin(RDBU[-1][1]), 1)
for v, c in RDBU[1:-1]:
    e = rp2.color_ramp.elements.new((v + 0.7) / 1.4); e.color = (*srgb2lin(c), 1)
em = n2.nodes.new('ShaderNodeEmission'); em.inputs['Strength'].default_value = 1.0
n2.links.new(tcn.outputs['Generated'], sep.inputs[0]); n2.links.new(sep.outputs['X'], rp2.inputs['Fac'])
n2.links.new(rp2.outputs['Color'], em.inputs['Color']); n2.links.new(em.outputs[0], n2.nodes['Material Output'].inputs['Surface'])
cb.materials.append(cm)
for v, xx in ((-0.7, -4.75), (0.0, -3.55), (0.7, -2.35)):
    hud(f'HUD_cb{v}', f'{v:+.1f}'.replace('+0.0', '0'), (xx, -2.58), 0.09, m_txt2, 'CENTER')
hud('HUD_cbl', 'Pearson r (proxy vs index)', (-2.2, -2.4), 0.1, m_txt2)
hud_year = hud('HUD_year', '', (4.75, -2.6), 0.3, m_txt, 'RIGHT')

bpy.context.view_layer.update()
for o in bpy.data.objects: o.select_set(False)
curves = [o for o in col_geo.objects if o.type == 'CURVE']
for o in curves: o.select_set(True)
bpy.context.view_layer.objects.active = curves[0]
bpy.ops.object.convert(target='MESH')

# ------------------------------------------------------------------ outputs
def set_state(year=None):
    """year=None → full-period values; otherwise 21-yr running values centred on `year`."""
    for c in D['cores']:
        e = PILLARS[c['id']]; run = c['running'].get(str(year)) if year else None
        for j, k in enumerate(('ONI', 'RONI')):
            r = c['full'][k]['r'] if year is None else (run[j] if run else 0.0)
            e[k]['r'] = r; e[k].scale.z = max(0.004, abs(r)) * H if (year is None or run) else 0.004
        e['halo'].location.z = e['base'] + (e['halo']['rcrit_full'] if year is None else D['meta']['rcrit21']) * H
        e['halo'].hide_render = bool(year and not run)

if '--glb' in ARGS:
    set_state(None)
    for o in bpy.data.objects: o.select_set(o.name in col_geo.objects and o.type == 'MESH')
    bpy.ops.export_scene.gltf(filepath=os.path.join(OUT, 'andes_enso.glb'), use_selection=True, export_colors=True,
                              export_extras=True, export_apply=True, export_yup=True)
    print('GLB OK')

if '--still' in ARGS:
    set_state(None); hud_year.data.body = 'Full period'
    scene.render.filepath = os.path.join(OUT, 'hero_full.png'); bpy.ops.render.render(write_still=True); print('STILL OK')

# animation keyframes: frame f ↔ centre year Y0 + f
YR = D['meta']['years']; scene.frame_start, scene.frame_end = YR[0] - Y0, YR[1] - Y0
for f in range(scene.frame_start, scene.frame_end + 1):
    set_state(Y0 + f)
    for c in D['cores']:
        e = PILLARS[c['id']]
        for k in ('ONI', 'RONI'):
            e[k].keyframe_insert('["r"]', frame=f); e[k].keyframe_insert('scale', index=2, frame=f)
        e['halo'].keyframe_insert('location', index=2, frame=f); e['halo'].keyframe_insert('hide_render', frame=f)
# year counter via a registered text-block handler (also runs when the .blend is opened with auto-run)
handler_src = '''import bpy
def _enso_year(scene, *a):
    ob = bpy.data.objects.get("HUD_year")
    if ob: ob.data.body = f"{1860 + scene.frame_current}  ·  21-yr window"
    s = bpy.data.objects.get("HUD_sub")
    if s: s.data.body = "21-yr running correlation · ERSST.v6 DJF ONI/RONI (CPC algorithm), 1850–present"
bpy.app.handlers.frame_change_pre[:] = [h for h in bpy.app.handlers.frame_change_pre if getattr(h, "__name__", "") != "_enso_year"]
bpy.app.handlers.frame_change_pre.append(_enso_year)
'''
tx = bpy.data.texts.new('enso_year_handler.py'); tx.write(handler_src); tx.use_module = True
exec(compile(handler_src, 'enso_year_handler.py', 'exec'))
scene.frame_set(1930 - Y0)
bpy.ops.wm.save_as_mainfile(filepath=os.path.join(OUT, 'andes_enso.blend'))
print('BLEND OK')

if '--video' in ARGS:
    scene.render.resolution_x, scene.render.resolution_y = 1280, 720
    ee.taa_render_samples = 24
    scene.render.fps = 12
    scene.render.image_settings.file_format = 'FFMPEG'; scene.render.ffmpeg.format = 'MPEG4'
    scene.render.ffmpeg.codec = 'H264'; scene.render.ffmpeg.constant_rate_factor = 'HIGH'
    scene.render.filepath = os.path.join(OUT, 'running.mp4')
    bpy.ops.render.render(animation=True); print('VIDEO OK')
