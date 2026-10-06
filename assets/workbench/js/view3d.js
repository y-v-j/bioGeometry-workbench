// Two linked three.js views: the structure (cartoon / trace / spheres / sticks + hetero groups)
// and the computational geometry (hull, Delaunay, α-shape, Voronoi + pocket / contact overlays).
import * as THREE from 'three';
import { TrackballControls } from 'three/addons/controls/TrackballControls.js';
import { LineSegments2 } from 'three/addons/lines/LineSegments2.js';
import { LineSegmentsGeometry } from 'three/addons/lines/LineSegmentsGeometry.js';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';
import { store, KIND, PALETTE, ELEMENT_COLOR, GROUP_COLORS, BIND_COLOR, RAMP, CAT10, LIGAND_CARBON, LIGAND_BOND, LIGAND_METAL, ION_COLOR, ION_DEFAULT, WATER_COLOR, hetColor, renderLegend, toast, busy, fmt } from './store.js';

const d3 = window.d3;
const VDW = { C: 1.7, N: 1.55, O: 1.52, S: 1.8, P: 1.8, Se: 1.9, Mg: 1.73, Zn: 1.39, Fe: 1.94, Ca: 2.31, Mn: 1.97,
  Cl: 1.75, Na: 2.27, K: 2.75, Cu: 1.4, Co: 1.92, Ni: 1.63, Br: 1.85, I: 1.98, F: 1.47, X: 1.7 };
const tmpC = new THREE.Color();
const WHITE = new THREE.Color('#ffffff');
const DIM_TO = new THREE.Color('#34344a');      // unselected residues fade towards the background
const MUTE = new THREE.Color('#c8c8d8');
// how strongly the selection colour is applied: kept moderate so a selection reads clearly without glare
const SEL_MIX = { cartoon: 0.35, point: 0.45, star: 0.15, layer: 0.3 };     // a light tint; the outline carries the selection
// unselected things lose their hue (so no scheme colour can be mistaken for the selection) but keep some light
function fade(c) {
  const l = 0.3 * c.r + 0.59 * c.g + 0.11 * c.b;
  return c.setRGB(c.r + (l - c.r) * 0.78, c.g + (l - c.g) * 0.78, c.b + (l - c.b) * 0.78).lerp(DIM_TO, 0.35).multiplyScalar(0.92);
}
// ligand / ion / water atoms: the same radii in both views
const hetRadius = (k, el) => (k === KIND.ion ? (VDW[el] || 1.6) * 0.62 : k === KIND.water ? 0.34 : 0.4);
// the selection accent follows store.setSelColor; cached per value
let selCache = ['', new THREE.Color()];
const selColor = () => { if (selCache[0] !== PALETTE.sel) selCache = [PALETTE.sel, new THREE.Color(PALETTE.sel)]; return selCache[1]; };
const MESH_EDGE = '#42cbee';                     // Delaunay edges, one colour (nodes: chain colour, A = #5eabd9)
const STAR_EDGE = '#aae8f7';                     // edges of the tetrahedra a selection takes part in: a lighter shade
const STAR_WIDEN = 0.75;                         // … and this many px wider than the mesh (same opacity)
const MAX_STAR_EDGES = 250000;                   // at most this many selection edges beyond the mesh cutoff

// ------------------------------------------------------------------ shared materials
// Shaded sphere impostors. 'world' spheres scale with the zoom (atoms in the structure view); 'dot' spheres
// keep a fixed size in pixels (geometry points, like Plotly markers: they never crowd the view when zoomed).
// Edges are anti-aliased (alpha-to-coverage with MSAA) for a crisp outline.
const SPHERE_FRAG = /* glsl */`
  varying vec3 vColor;
  #include <clipping_planes_pars_fragment>
  void main() {
    #include <clipping_planes_fragment>
    vec2 c = gl_PointCoord * 2.0 - 1.0; float r = length(c);
    float aa = max(fwidth(r), 1e-3); float alpha = 1.0 - smoothstep(1.0 - aa, 1.0, r);
    if (alpha <= 0.0) discard;
    float r2 = min(r * r, 1.0);
    vec3 n = vec3(c.x, -c.y, sqrt(1.0 - r2)); vec3 L = normalize(vec3(0.35, 0.55, 1.0));
    float diff = max(dot(n, L), 0.0); float spec = pow(max(dot(reflect(-L, n), vec3(0.0, 0.0, 1.0)), 0.0), 28.0);
    float rim = pow(1.0 - n.z, 3.0) * 0.2;
    gl_FragColor = vec4(vColor * (0.34 + 0.70 * diff) + 0.18 * spec + rim * vColor, alpha);
    #include <colorspace_fragment>
  }`;                                  // colours are linear (THREE.Color); without the sRGB step they came out too dark
function impostorMaterial(uniforms, dot = false) {
  return new THREE.ShaderMaterial({
    uniforms,
    vertexShader: /* glsl */`
      attribute float aRadius; attribute vec3 aColor; varying vec3 vColor; uniform float uScale; uniform float uPx;
      #include <clipping_planes_pars_vertex>
      void main() {
        vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
        vColor = aColor;
        ${dot ? 'gl_PointSize = aRadius * uPx;            // aRadius = diameter in CSS pixels'
              : 'gl_PointSize = clamp(2.0 * aRadius * uScale / -mvPosition.z, 1.5, 512.0);'}
        gl_Position = projectionMatrix * mvPosition;
        #include <clipping_planes_vertex>
      }`,
    fragmentShader: SPHERE_FRAG,
    clipping: true,
    alphaToCoverage: true,
  });
}
const dotMaterial = (uniforms) => impostorMaterial(uniforms, true);
// thick, anti-aliased line segments (WebGL lines are always 1 px): pos / col hold two xyz / rgb per segment
function fatLines(pos, col, width, opacity = 1) {
  const g = new LineSegmentsGeometry();
  g.setPositions(pos); g.setColors(col);
  const m = new LineMaterial({ linewidth: width, vertexColors: true, transparent: opacity < 1, opacity,
    depthWrite: opacity >= 1, alphaToCoverage: opacity >= 1 });
  m.userData.fat = true;
  return new LineSegments2(g, m);
}
function recolorFat(ls, col) {
  const buf = ls.geometry.attributes.instanceColorStart.data;
  buf.array.set(col); buf.needsUpdate = true;
}
// bonds as instanced half-cylinders, each half coloured like its atom
function bondMesh(A, pairs, colorOf, radius) {
  const n = pairs.length / 2;
  const geo = new THREE.CylinderGeometry(radius, radius, 1, 7, 1, true);
  geo.translate(0, 0.5, 0); geo.rotateX(Math.PI / 2);          // unit length along +z from the origin
  const mesh = new THREE.InstancedMesh(geo, new THREE.MeshPhongMaterial({ shininess: 50, specular: 0x333344 }), n * 2);
  const m = new THREE.Matrix4(), pa = new THREE.Vector3(), pb = new THREE.Vector3(), mid = new THREE.Vector3();
  const q = new THREE.Quaternion(), sc = new THREE.Vector3(), z = new THREE.Vector3(0, 0, 1), dir = new THREE.Vector3(), c = new THREE.Color();
  for (let k = 0; k < n; k++) {
    const a = pairs[2 * k], b = pairs[2 * k + 1];
    pa.fromArray(A, 3 * a); pb.fromArray(A, 3 * b); mid.addVectors(pa, pb).multiplyScalar(0.5);
    dir.subVectors(pb, pa); const len = dir.length(); dir.divideScalar(len || 1);
    q.setFromUnitVectors(z, dir); sc.set(1, 1, len / 2);
    m.compose(pa, q, sc); mesh.setMatrixAt(2 * k, m); mesh.setColorAt(2 * k, c.set(colorOf(a)));
    m.compose(mid, q, sc); mesh.setMatrixAt(2 * k + 1, m); mesh.setColorAt(2 * k + 1, c.set(colorOf(b)));
  }
  mesh.instanceMatrix.needsUpdate = true;
  if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true;
  return mesh;
}
// ligands (ball-and-stick), ions and coordinated waters: one builder and one painter serve both 3D views, so
// they are drawn with the same colours, sizes, bonds and lights in the structure and the geometry view
function buildHetBallStick(uniforms) {
  const b = store.S.buf, A = b.atom_pos, els = store.S.stats.elements, idx = [];
  for (let i = 0; i < b.atom_res.length; i++) if (b.res_kind[b.atom_res[i]] >= KIND.ligand) idx.push(i);
  if (!idx.length) return null;
  const pos = new Float32Array(idx.length * 3), rad = new Float32Array(idx.length), col = new Float32Array(idx.length * 3);
  idx.forEach((a, i) => {
    const k = b.res_kind[b.atom_res[a]], el = els[b.atom_elem[a]];
    pos.set(A.subarray(3 * a, 3 * a + 3), 3 * i);
    rad[i] = hetRadius(k, el);
    tmpC.set(hetColor(k, el));
    col[3 * i] = tmpC.r; col[3 * i + 1] = tmpC.g; col[3 * i + 2] = tmpC.b;
  });
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  g.setAttribute('aColor', new THREE.BufferAttribute(col.slice(), 3));
  g.setAttribute('aRadius', new THREE.BufferAttribute(rad.slice(), 1));
  const pts = new THREE.Points(g, impostorMaterial(uniforms));
  const B = b.bonds, pairs = [];
  for (let k = 0; k < B.length / 2; k++) if (b.res_kind[b.atom_res[B[2 * k]]] === KIND.ligand) pairs.push(B[2 * k], B[2 * k + 1]);
  const bonds = pairs.length ? bondMesh(A, pairs, () => LIGAND_BOND, 0.13) : null;
  return { pts, bonds, atoms: Uint32Array.from(idx), col, rad };
}
// selected ligands / ions / waters grow (their colour stays), the rest fade with the dim toggle, hover lightens
function paintHetBallStick(H) {
  const b = store.S.buf, any = store.sel.size > 0, { pts, atoms, col, rad } = H;
  const ca = pts.geometry.attributes.aColor, ra = pts.geometry.attributes.aRadius;
  for (let i = 0; i < atoms.length; i++) {
    const r = b.atom_res[atoms[i]], sel = store.sel.has(r);
    tmpC.setRGB(col[3 * i], col[3 * i + 1], col[3 * i + 2]);
    if (any && !sel && store.dim) fade(tmpC);
    if (store.hover === r) tmpC.lerp(WHITE, 0.4);
    ca.array[3 * i] = tmpC.r; ca.array[3 * i + 1] = tmpC.g; ca.array[3 * i + 2] = tmpC.b;
    ra.array[i] = sel ? rad[i] * 1.2 : rad[i];
  }
  ca.needsUpdate = true; ra.needsUpdate = true;
}
// the same lights for the main scene and the overlay scene (ligands drawn on top in the geometry view)
function addLights(scene) {
  scene.add(new THREE.AmbientLight(0xffffff, 0.5));
  const key = new THREE.DirectionalLight(0xffffff, 1.5), rim = new THREE.DirectionalLight(0xc4b5fd, 0.35);
  scene.add(key, rim);
  return [key, rim];
}
// tetrahedra shrunken towards their centroid (so neighbours stay distinguishable), flat shaded
function tetMesh(P, tets, list, shrink, colorOf, opacity) {
  const n = list.length, pos = new Float32Array(n * 36), col = new Float32Array(n * 36);
  const v = [new THREE.Vector3(), new THREE.Vector3(), new THREE.Vector3(), new THREE.Vector3()], cen = new THREE.Vector3(), c = new THREE.Color();
  const FACES = [[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]];
  for (let i = 0; i < n; i++) {
    const t = list[i];
    cen.set(0, 0, 0);
    for (let k = 0; k < 4; k++) { v[k].fromArray(P, 3 * tets[4 * t + k]); cen.add(v[k]); }
    cen.multiplyScalar(0.25);
    for (let k = 0; k < 4; k++) v[k].lerp(cen, 1 - shrink);
    c.set(colorOf(t, i));
    let o = 36 * i;
    for (const f of FACES) for (const k of f) { pos[o] = v[k].x; pos[o + 1] = v[k].y; pos[o + 2] = v[k].z; col[o] = c.r; col[o + 1] = c.g; col[o + 2] = c.b; o += 3; }
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  g.setAttribute('color', new THREE.BufferAttribute(col, 3));
  g.computeVertexNormals();
  return new THREE.Mesh(g, meshMat(opacity, true));
}
const meshMat = (opacity, flat = false) => new THREE.MeshPhongMaterial({
  vertexColors: true, transparent: opacity < 1, opacity, side: THREE.DoubleSide, shininess: 45, specular: 0x222233,
  depthWrite: opacity >= 0.95, flatShading: flat, polygonOffset: true, polygonOffsetFactor: 1,
});

function disposeGroup(g) {
  g.traverse((o) => {
    if (o.geometry) o.geometry.dispose();
    if (o.material) (Array.isArray(o.material) ? o.material : [o.material]).forEach((m) => m.dispose());
  });
  g.clear();
}

// ------------------------------------------------------------------ base view
class View3D {
  constructor(hostId, name) {
    this.name = name;
    this.host = document.getElementById(hostId);
    this.renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true, alpha: false });
    // supersample on low-density screens: thin lines and sphere edges stay crisp (cf. Plotly's rendering)
    this.renderer.setPixelRatio(Math.min(2, Math.max(1.5, window.devicePixelRatio || 1)));
    this.renderer.setClearColor(0x161622, 1);
    this.renderer.localClippingEnabled = true;
    this.renderer.domElement.className = 'gl';
    this.host.prepend(this.renderer.domElement);
    this.scene = new THREE.Scene();
    this.scene.fog = null;                // no depth fog: it softened the picture
    this.camera = new THREE.PerspectiveCamera(38, 1, 0.5, 20000);
    this.camera.position.set(0, 0, 200);
    this.controls = new TrackballControls(this.camera, this.renderer.domElement);
    Object.assign(this.controls, { rotateSpeed: 3.2, zoomSpeed: 1.5, panSpeed: 0.9, dynamicDampingFactor: 0.2 });
    [this.key, this.rim] = addLights(this.scene);
    // drawn after the main scene with a fresh depth buffer: always visible, still depth-correct among themselves
    this.top = new THREE.Scene();
    [this.topKey, this.topRim] = addLights(this.top);
    this.uniforms = { uScale: { value: 800 }, uPx: { value: this.renderer.getPixelRatio() } };
    this.cssSize = [1, 1];
    this.dirty = 2; this.tweening = 0; this.spin = false; this.slab = false;
    this.slabPlanes = [new THREE.Plane(), new THREE.Plane()];
    this.clip = [];
    this.peer = null; this.linked = false; this.syncing = false;
    this.frames = 0; this.lastFps = performance.now(); this.fps = 0;
    this.controls.addEventListener('change', () => { this.dirty = 2; this.syncPeer(); });
    new ResizeObserver(() => this.resize()).observe(this.host);
    this.tip = this.host.querySelector('.vtip');
    this.renderer.domElement.addEventListener('pointermove', (e) => this.onMove(e));
    this.renderer.domElement.addEventListener('pointerleave', () => { if (this.tip) this.tip.hidden = true; store.setHover(null, this.name); });
    this.renderer.domElement.addEventListener('pointerdown', (e) => { this.downAt = [e.clientX, e.clientY]; });
    this.renderer.domElement.addEventListener('click', (e) => this.onClick(e));
    this.renderer.domElement.addEventListener('dblclick', (e) => this.onDblClick(e));
  }
  request() { this.dirty = Math.max(this.dirty, 1); }
  resize() {
    const w = this.host.clientWidth, h = this.host.clientHeight;
    if (!w || !h) return;
    this.renderer.setSize(w, h, false);
    this.cssSize = [w, h];
    this.camera.aspect = w / h; this.camera.updateProjectionMatrix(); this.controls.handleResize();
    this.uniforms.uScale.value = this.renderer.domElement.height / (2 * Math.tan(THREE.MathUtils.degToRad(this.camera.fov) / 2));
    this.request();
  }
  frame(extent, animate = true) {
    const r = extent || 50, dist = r / Math.sin(THREE.MathUtils.degToRad(this.camera.fov / 2)) * 1.02;
    this.camera.near = Math.max(0.5, dist / 400); this.camera.far = dist * 8; this.camera.updateProjectionMatrix();
    if (this.scene.fog) { this.scene.fog.near = dist * 0.65; this.scene.fog.far = dist * 2.3; }
    this.flyTo(new THREE.Vector3(), new THREE.Vector3(0, 0, dist), animate ? 700 : 0);
  }
  flyTo(target, position, ms = 900) {
    if (!ms) { this.controls.target.copy(target); this.camera.position.copy(position); this.camera.up.set(0, 1, 0); this.request(); return; }
    const t0 = this.controls.target.clone(), p0 = this.camera.position.clone();
    this.tweening++;
    d3.select({}).transition().duration(ms).ease(d3.easeCubicInOut).tween('fly', () => (t) => {
      this.controls.target.lerpVectors(t0, target, t); this.camera.position.lerpVectors(p0, position, t); this.request(); this.syncPeer();
    }).on('end interrupt', () => { this.tweening = Math.max(0, this.tweening - 1); this.request(); });
  }
  focusResidues(rs) {
    if (!rs.length || !store.S) return;
    const P = store.S.buf.res_pos, c = new THREE.Vector3();
    rs.forEach((r) => { c.x += P[3 * r]; c.y += P[3 * r + 1]; c.z += P[3 * r + 2]; });
    c.divideScalar(rs.length);
    let rad = 0;
    rs.forEach((r) => { rad = Math.max(rad, Math.hypot(P[3 * r] - c.x, P[3 * r + 1] - c.y, P[3 * r + 2] - c.z)); });
    const dir = this.camera.position.clone().sub(this.controls.target).normalize();
    this.flyTo(c, c.clone().add(dir.multiplyScalar(Math.max(22, rad * 3.0))));
  }
  syncPeer() {
    if (!this.linked || !this.peer || this.syncing) return;
    this.peer.syncing = true;
    this.peer.camera.position.copy(this.camera.position); this.peer.camera.up.copy(this.camera.up);
    this.peer.controls.target.copy(this.controls.target); this.peer.camera.lookAt(this.controls.target);
    this.peer.request();
    this.peer.syncing = false;
  }
  setSlab(on) {
    this.slab = on;
    this.clip = on ? this.slabPlanes : [];
    const upd = (o) => { if (o.material) (Array.isArray(o.material) ? o.material : [o.material]).forEach((m) => { m.clippingPlanes = this.clip; m.needsUpdate = true; }); };
    this.scene.traverse(upd); this.top.traverse(upd);
    this.request();
  }
  updateSlab() {
    if (!this.slab) return;
    const d = this.controls.target.clone().sub(this.camera.position).normalize();
    const mid = d.dot(this.controls.target), half = Math.max(8, (store.S ? store.S.stats.extent : 40) * 0.25);
    this.slabPlanes[0].set(d, -(mid - half)); this.slabPlanes[1].set(d.clone().negate(), mid + half);
  }
  applyClip(root) { root.traverse((o) => { if (o.material) (Array.isArray(o.material) ? o.material : [o.material]).forEach((m) => { m.clippingPlanes = this.clip; }); }); }
  tick(now) {
    if (this.spin) {
      const off = this.camera.position.clone().sub(this.controls.target);
      off.applyAxisAngle(this.camera.up, 0.004);
      this.camera.position.copy(this.controls.target).add(off); this.dirty = 1; this.syncPeer();
    }
    this.controls.update();
    if (!this.dirty && !this.tweening) return false;
    this.dirty = Math.max(0, this.dirty - 1);
    this.key.position.copy(this.camera.position);
    this.rim.position.copy(this.camera.position).multiplyScalar(-1);
    this.updateSlab();
    this.renderAll();
    this.frames++;
    if (now - this.lastFps > 1000) { this.fps = Math.round(this.frames * 1000 / (now - this.lastFps)); this.frames = 0; this.lastFps = now; }
    return true;
  }
  renderAll() {
    const [w, h] = this.cssSize;                 // fat lines measure their width in CSS pixels
    const fat = (o) => { if (o.material && o.material.userData && o.material.userData.fat) o.material.resolution.set(w, h); };
    this.scene.traverse(fat); this.top.traverse(fat);
    this.renderer.render(this.scene, this.camera);
    if (this.top.children.some((o) => !o.isLight)) {
      this.topKey.position.copy(this.key.position); this.topRim.position.copy(this.rim.position);
      this.renderer.autoClear = false;
      this.renderer.clearDepth();
      this.renderer.render(this.top, this.camera);
      this.renderer.autoClear = true;
    }
  }
  screenshot() {
    this.renderAll();
    const a = document.createElement('a');
    a.href = this.renderer.domElement.toDataURL('image/png');
    a.download = `${store.S ? store.S.stats.name : 'view'}_${this.name}.png`;
    a.click();
  }
  // picking: screen-space nearest residue anchor
  pick(clientX, clientY) {
    if (!store.S) return null;
    const rect = this.renderer.domElement.getBoundingClientRect();
    const P = store.S.buf.res_pos, n = P.length / 3, v = new THREE.Vector3();
    let best = null, bestD = 14 * 14, bestZ = Infinity;
    const vis = this.pickable;
    for (let r = 0; r < n; r++) {
      if (vis && !vis[r]) continue;
      v.set(P[3 * r], P[3 * r + 1], P[3 * r + 2]);
      if (this.slab && (this.slabPlanes[0].distanceToPoint(v) < 0 || this.slabPlanes[1].distanceToPoint(v) < 0)) continue;
      v.project(this.camera);
      if (v.z > 1) continue;
      const dx = (v.x + 1) / 2 * rect.width - (clientX - rect.left), dy = (1 - v.y) / 2 * rect.height - (clientY - rect.top);
      const d = dx * dx + dy * dy;
      if (d < bestD - 4 || (d < bestD + 4 && v.z < bestZ)) { best = r; bestD = d; bestZ = v.z; }
    }
    return best;
  }
  onMove(e) {
    if (e.buttons) { if (this.tip) this.tip.hidden = true; return; }
    clearTimeout(this.hoverT);
    this.hoverT = setTimeout(() => {
      const r = this.pick(e.clientX, e.clientY);
      store.setHover(r, this.name);
      if (!this.tip) return;
      if (r === null) { this.tip.hidden = true; return; }
      this.tip.innerHTML = residueTip(r);
      const rc = this.host.getBoundingClientRect();
      this.tip.style.left = `${Math.min(e.clientX - rc.left + 12, rc.width - 230)}px`;
      this.tip.style.top = `${Math.min(e.clientY - rc.top + 12, rc.height - 80)}px`;
      this.tip.hidden = false;
    }, 35);
  }
  onClick(e) {
    if (this.downAt && Math.hypot(e.clientX - this.downAt[0], e.clientY - this.downAt[1]) > 4) return;  // a drag
    const r = this.pick(e.clientX, e.clientY);
    if (r === null) { if (!store.additive(e)) store.clear(this.name); return; }
    store.click(r, e, this.name);
    store.setFocusChain(store.S.buf.res_chain[r], this.name);
  }
  onDblClick(e) {
    const r = this.pick(e.clientX, e.clientY);
    if (r !== null) this.focusResidues([r]); else if (store.sel.size) this.focusResidues([...store.sel]);
  }
}

export function residueTip(r) {
  const S = store.S, b = S.buf;
  const rows = [`<b>${store.label(r)}</b> <span class="k">${S.stats.kinds[b.res_kind[r]]} · ${S.stats.classes[b.res_class[r]]}</span>`];
  if (b.res_kind[r] === KIND.protein) {
    rows.push(`<span class="k">SS</span> ${({ 1: 'helix', 2: 'strand', 0: 'coil' })[b.res_ss[r]] || '–'} &nbsp;<span class="k">φ/ψ</span> ${f0(b.res_phi[r])} / ${f0(b.res_psi[r])}`);
  }
  rows.push(`<span class="k">${S.stats.is_alphafold ? 'pLDDT' : 'B-factor'}</span> ${b.res_b[r].toFixed(1)}`);
  if (store.F) {
    const add = (key, lab, d = 2) => { const v = store.featureValue(key, r); if (typeof v === 'number' ? Number.isFinite(v) : v) rows.push(`<span class="k">${lab}</span> ${typeof v === 'number' ? v.toFixed(d) : v}`); };
    add('rsa', 'RSA'); add('vor_volume', 'Voronoi vol', 1); add('hull_depth', 'depth', 1); add('site', 'site'); add('interface', 'interface'); add('domain', 'domain');
  }
  return rows.join('<br>');
}
const f0 = (v) => (Number.isFinite(v) ? `${v.toFixed(0)}°` : '–');

// ------------------------------------------------------------------ structure view
export class StructureView extends View3D {
  constructor() {
    super('mol-host', 'mol');
    this.groups = { poly: new THREE.Group(), het: new THREE.Group(), sel: new THREE.Group(), bind: new THREE.Group() };
    Object.values(this.groups).forEach((g) => this.scene.add(g));
    this.style = 'cartoon'; this.colorBy = 'chain'; this.dim = false; this.showBinding = true;   // dim: see store.dim
    this.legendEl = document.getElementById('mol-legend');
    store.on('data', () => this.load());
    store.on('selection', () => this.updateSelection());
    store.on('selstyle', () => this.updateSelection());
    store.on('features', () => { if (store.F && store.F.colIndex[this.colorBy]) this.recolor(); });
    store.on('hover', ({ source }) => { if (source !== 'mol') this.updateHover(); });
  }
  load() {
    disposeGroup(this.groups.poly); disposeGroup(this.groups.het); disposeGroup(this.groups.sel); disposeGroup(this.groups.bind);
    this.cartoon = null; this.spheres = null; this.trace = null; this.het = null; this.hetBonds = null; this.bind = null;
    this.build();
    this.buildHetero();
    this.buildBinding();
    this.frame(store.S.stats.extent, true);
    this.updateSelection();
  }
  setStyle(s) { this.style = s; disposeGroup(this.groups.poly); this.cartoon = this.spheres = this.trace = null; this.build(); this.updateSelection(); }
  setColor(c) { this.colorBy = c; this.recolor(); }
  build() {
    const S = store.S, b = S.buf;
    this.resColors = store.residueColors(this.colorBy).map((c) => new THREE.Color(c));
    const poly = new Uint8Array(b.res_kind.length);
    for (let r = 0; r < poly.length; r++) poly[r] = b.res_kind[r] <= KIND.nucleic ? 1 : 0;
    this.pickable = poly;
    if (this.style === 'cartoon' && b.cart_pos.length) this.buildCartoon();
    else if (this.style === 'spheres' || this.style === 'sticks') this.buildAtoms(this.style);
    else this.buildTrace();
    if (this.style === 'cartoon' && b.ladder.length) this.buildLadder();
    this.applyClip(this.groups.poly);
    this.recolor();
  }
  buildCartoon() {
    const S = store.S, b = S.buf, P = b.cart_pos, rings = b.cart_rings, side = S.stats.cartoon.sides, tside = S.stats.cartoon.tube_sides;
    const V = P.length / 3, idx = [];
    let base = 0;
    for (let s = 0; s < rings.length; s++) {
      const n = Math.abs(rings[s]), cols = (rings[s] < 0 ? tside : side) + 1;
      for (let i = 0; i < n - 1; i++) for (let j = 0; j < cols - 1; j++) {
        const a = base + i * cols + j, c = a + cols;
        idx.push(a, c, a + 1, a + 1, c, c + 1);
      }
      base += n * cols;
    }
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(P, 3));
    g.setAttribute('color', new THREE.BufferAttribute(new Float32Array(V * 3), 3));
    g.setIndex(V > 65535 ? new THREE.Uint32BufferAttribute(idx, 1) : new THREE.Uint16BufferAttribute(idx, 1));
    g.computeVertexNormals();
    const m = new THREE.Mesh(g, new THREE.MeshPhongMaterial({ vertexColors: true, shininess: 55, specular: 0x30304a, side: THREE.DoubleSide }));
    this.groups.poly.add(m);
    this.cartoon = { mesh: m, res: b.cart_res };
  }
  buildLadder() {
    const b = store.S.buf, L = b.ladder, A = b.atom_pos, n = L.length / 2;
    const pos = new Float32Array(n * 6), col = new Float32Array(n * 6);
    for (let k = 0; k < n; k++) { pos.set(A.subarray(3 * L[2 * k], 3 * L[2 * k] + 3), 6 * k); pos.set(A.subarray(3 * L[2 * k + 1], 3 * L[2 * k + 1] + 3), 6 * k + 3); }
    const ls = fatLines(pos, col, 2.2);
    this.groups.poly.add(ls);
    this.ladder = { obj: ls, res: Uint32Array.from({ length: n * 2 }, (_, i) => b.atom_res[L[i]]) };
  }
  buildTrace() {
    const b = store.S.buf, T = b.trace, R = b.res_pos, n = T.length / 2;
    const pos = new Float32Array(n * 6);
    for (let k = 0; k < n; k++) { pos.set(R.subarray(3 * T[2 * k], 3 * T[2 * k] + 3), 6 * k); pos.set(R.subarray(3 * T[2 * k + 1], 3 * T[2 * k + 1] + 3), 6 * k + 3); }
    const ls = fatLines(pos, new Float32Array(n * 6), 3);
    this.groups.poly.add(ls);
    // residue beads
    const npoly = [];
    for (let r = 0; r < b.res_kind.length; r++) if (b.res_kind[r] <= KIND.nucleic) npoly.push(r);
    const bp = new Float32Array(npoly.length * 3);
    npoly.forEach((r, i) => bp.set(R.subarray(3 * r, 3 * r + 3), 3 * i));
    const bg = new THREE.BufferGeometry();
    bg.setAttribute('position', new THREE.BufferAttribute(bp, 3));
    bg.setAttribute('aColor', new THREE.BufferAttribute(new Float32Array(npoly.length * 3), 3));
    bg.setAttribute('aRadius', new THREE.BufferAttribute(new Float32Array(npoly.length).fill(0.75), 1));
    const beads = new THREE.Points(bg, impostorMaterial(this.uniforms));
    this.groups.poly.add(beads);
    this.trace = { lines: ls, res: Uint32Array.from({ length: n * 2 }, (_, i) => T[i]), beads, beadRes: Uint32Array.from(npoly) };
  }
  buildAtoms(style) {
    const b = store.S.buf, A = b.atom_pos, els = store.S.stats.elements, idx = [];
    for (let i = 0; i < b.atom_res.length; i++) if (b.res_kind[b.atom_res[i]] <= KIND.nucleic) idx.push(i);
    const pos = new Float32Array(idx.length * 3), rad = new Float32Array(idx.length);
    const scale = style === 'spheres' ? 1.0 : 0.26;
    idx.forEach((a, i) => { pos.set(A.subarray(3 * a, 3 * a + 3), 3 * i); rad[i] = (VDW[els[b.atom_elem[a]]] || 1.7) * scale; });
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    g.setAttribute('aColor', new THREE.BufferAttribute(new Float32Array(idx.length * 3), 3));
    g.setAttribute('aRadius', new THREE.BufferAttribute(rad, 1));
    const pts = new THREE.Points(g, impostorMaterial(this.uniforms));
    this.groups.poly.add(pts);
    this.spheres = { pts, atoms: Uint32Array.from(idx) };
    if (style === 'sticks') {
      const B = b.bonds, keep = [];
      for (let k = 0; k < B.length / 2; k++) if (b.res_kind[b.atom_res[B[2 * k]]] <= KIND.nucleic) keep.push(k);
      const bp = new Float32Array(keep.length * 6);
      keep.forEach((k, i) => { bp.set(A.subarray(3 * B[2 * k], 3 * B[2 * k] + 3), 6 * i); bp.set(A.subarray(3 * B[2 * k + 1], 3 * B[2 * k + 1] + 3), 6 * i + 3); });
      const ls = fatLines(bp, new Float32Array(keep.length * 6), 2.5);
      this.groups.poly.add(ls);
      this.spheres.bonds = { obj: ls, atoms: Uint32Array.from(keep.flatMap((k) => [B[2 * k], B[2 * k + 1]])) };
    }
  }
  buildHetero() {
    const H = buildHetBallStick(this.uniforms);
    if (!H) return;
    this.het = H;
    this.groups.het.add(H.pts);
    if (H.bonds) { this.hetBonds = H.bonds; this.groups.het.add(H.bonds); }
    this.applyClip(this.groups.het);
  }
  // binding contacts as dotted lines: small beads every 0.3 Å between the two atoms, coloured by type
  // (metal coordination, hydrogen bond, water bridge, covalent link)
  buildBinding() {
    disposeGroup(this.groups.bind);
    this.bind = null;
    const b = store.S.buf, Pp = b.bind_pairs, T = b.bind_type, A = b.atom_pos;
    if (!Pp || !Pp.length) return;
    const pos = [], pair = [], step = 0.3, gap = 0.45, va = new THREE.Vector3(), vb = new THREE.Vector3(), d = new THREE.Vector3();
    for (let k = 0; k < T.length; k++) {
      va.fromArray(A, 3 * Pp[2 * k]); vb.fromArray(A, 3 * Pp[2 * k + 1]);
      d.subVectors(vb, va); const len = d.length(); d.divideScalar(len || 1);
      for (let t = gap; t <= len - gap + 1e-6; t += step) { pos.push(va.x + d.x * t, va.y + d.y * t, va.z + d.z * t); pair.push(k); }
    }
    const n = pair.length, g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(new Float32Array(pos), 3));
    g.setAttribute('aColor', new THREE.BufferAttribute(new Float32Array(n * 3), 3));
    g.setAttribute('aRadius', new THREE.BufferAttribute(new Float32Array(n).fill(0.085), 1));
    const pts = new THREE.Points(g, impostorMaterial(this.uniforms));
    this.groups.bind.add(pts);
    this.groups.bind.visible = this.showBinding;
    this.bind = { pts, pair: Uint32Array.from(pair) };
    this.applyClip(this.groups.bind);
    this.paintBinding();
  }
  paintBinding() {
    if (!this.bind) return;
    const b = store.S.buf, Pp = b.bind_pairs, T = b.bind_type, any = store.sel.size > 0, attr = this.bind.pts.geometry.attributes.aColor;
    const cols = BIND_COLOR.map((c) => new THREE.Color(c)), c = new THREE.Color();
    for (let i = 0; i < this.bind.pair.length; i++) {
      const k = this.bind.pair[i];
      c.copy(cols[T[k]] || cols[1]);
      if (any && this.dim && !store.sel.has(b.atom_res[Pp[2 * k]]) && !store.sel.has(b.atom_res[Pp[2 * k + 1]])) c.lerp(DIM_TO, 0.7);
      attr.array[3 * i] = c.r; attr.array[3 * i + 1] = c.g; attr.array[3 * i + 2] = c.b;
    }
    attr.needsUpdate = true;
  }
  setBinding(on) { this.showBinding = on; this.groups.bind.visible = on; this.renderLegend(); this.request(); }
  // colour of residue r after selection / dim rules
  shade(r, out) {
    out.copy(this.resColors[r] || DIM_TO);
    if (store.sel.size) {
      if (store.sel.has(r)) out.lerp(selColor(), SEL_MIX.cartoon).lerp(WHITE, 0.08);
      else if (this.dim) fade(out);
    }
    if (store.hover === r) out.lerp(WHITE, 0.45);
    return out;
  }
  recolor() {
    if (!store.S) return;
    this.resColors = store.residueColors(this.colorBy).map((c) => new THREE.Color(c));
    this.paint();
    this.renderLegend();
  }
  // colour-scheme legend plus a key to the binding lines that are present
  renderLegend() {
    renderLegend(this.legendEl, store.legendFor(this.colorBy));
    const S = store.S, counts = S && S.stats.bind_counts;
    if (!this.showBinding || !counts || !counts.some((n) => n)) return;
    const row = d3.select(this.legendEl).append('div').attr('class', 'bind-key');
    S.stats.bind_types.forEach((name, k) => {
      if (!counts[k]) return;
      const sp = row.append('span').attr('class', 'sw').attr('title', `${counts[k]} ${name}${counts[k] > 1 ? 's' : ''}`);
      sp.append('i').attr('class', 'dots').style('color', BIND_COLOR[k]);
      sp.append('span').text(name.replace('coordination', 'coord.'));
    });
  }
  paint() {
    const c = new THREE.Color(), b = store.S.buf;
    const fillVerts = (attr, resOf) => {
      const a = attr.array;
      let last = -1;
      for (let i = 0; i < resOf.length; i++) {
        const r = resOf[i];
        if (r !== last) { this.shade(r, c); last = r; }
        a[3 * i] = c.r; a[3 * i + 1] = c.g; a[3 * i + 2] = c.b;
      }
      attr.needsUpdate = true;
    };
    if (this.cartoon) fillVerts(this.cartoon.mesh.geometry.attributes.color, this.cartoon.res);
    const fillFat = (obj, resOf) => {
      const col = new Float32Array(resOf.length * 3);
      fillVerts({ array: col }, resOf);
      recolorFat(obj, col);
    };
    if (this.ladder) fillFat(this.ladder.obj, this.ladder.res);
    if (this.trace) { fillFat(this.trace.lines, this.trace.res); fillVerts(this.trace.beads.geometry.attributes.aColor, this.trace.beadRes); }
    if (this.spheres) {
      const atoms = this.spheres.atoms, attr = this.spheres.pts.geometry.attributes.aColor, a = attr.array, els = store.S.stats.elements;
      for (let i = 0; i < atoms.length; i++) {
        const at = atoms[i];
        if (this.colorBy === 'element') { c.set(ELEMENT_COLOR[els[b.atom_elem[at]]] || '#ccc'); if (store.sel.size && !store.sel.has(b.atom_res[at]) && this.dim) c.multiplyScalar(0.28); }
        else this.shade(b.atom_res[at], c);
        a[3 * i] = c.r; a[3 * i + 1] = c.g; a[3 * i + 2] = c.b;
      }
      attr.needsUpdate = true;
      if (this.spheres.bonds) {
        const ba = this.spheres.bonds.atoms, bb = new Float32Array(ba.length * 3);
        for (let i = 0; i < ba.length; i++) { this.shade(b.atom_res[ba[i]], c); bb[3 * i] = c.r; bb[3 * i + 1] = c.g; bb[3 * i + 2] = c.b; }
        recolorFat(this.spheres.bonds.obj, bb);
      }
    }
    this.paintHetero();
    this.paintBinding();
    this.request();
  }
  // ligands / ions / waters: selected ones grow and glow, the rest fade with the dim toggle
  paintHetero() { if (this.het) paintHetBallStick(this.het); }
  updateHover() { if (store.S) this.paint(); }
  updateSelection() {
    if (!store.S) return;
    this.paint();
    // selected residues: a light tint and a crisp outline around their own cartoon (an enlarged back-face
    // shell in the selection colour); the trace style colours its line; side-chain sticks for small selections
    disposeGroup(this.groups.sel);
    const b = store.S.buf, els = store.S.stats.elements, A = b.atom_pos, R = b.res_pos;
    if (!store.sel.size) { this.request(); return; }
    const sc = selColor();
    if (this.cartoon) {
      const g0 = this.cartoon.mesh.geometry, res = this.cartoon.res, idx = g0.index.array;
      const keep = [];
      for (let t = 0; t < idx.length; t += 3) {
        if (store.sel.has(res[idx[t]]) && store.sel.has(res[idx[t + 1]]) && store.sel.has(res[idx[t + 2]])) keep.push(idx[t], idx[t + 1], idx[t + 2]);
      }
      if (keep.length) {
        const P = g0.attributes.position.array, N = g0.attributes.normal.array, out = new Float32Array(P.length), w = 0.28;
        for (let i = 0; i < P.length; i++) out[i] = P[i] + N[i] * w;
        const og = new THREE.BufferGeometry();
        og.setAttribute('position', new THREE.BufferAttribute(out, 3));
        og.setIndex(P.length / 3 > 65535 ? new THREE.Uint32BufferAttribute(keep, 1) : new THREE.Uint16BufferAttribute(keep, 1));
        this.groups.sel.add(new THREE.Mesh(og, new THREE.MeshBasicMaterial({ color: sc, side: THREE.BackSide })));
      }
    } else if (this.trace) {
      const T = b.trace, seg = [], col = [];
      for (let k = 0; k < T.length / 2; k++) {
        const a = T[2 * k], z = T[2 * k + 1];
        if (!store.sel.has(a) || !store.sel.has(z)) continue;
        seg.push(R[3 * a], R[3 * a + 1], R[3 * a + 2], R[3 * z], R[3 * z + 1], R[3 * z + 2]);
        col.push(sc.r, sc.g, sc.b, sc.r, sc.g, sc.b);
      }
      if (seg.length) { const ls = fatLines(new Float32Array(seg), new Float32Array(col), 5); ls.renderOrder = 3; this.groups.sel.add(ls); }
    }
    const carbon = selColor().clone().lerp(MUTE, 0.5).getStyle();
    const atomColor = (a) => (els[b.atom_elem[a]] === 'C' ? carbon : ELEMENT_COLOR[els[b.atom_elem[a]]] || '#ccc');
    let nRes = 0; for (const r of store.sel) if (b.res_kind[r] <= KIND.nucleic) nRes++;
    const atoms = [];
    if (this.style !== 'spheres' && nRes <= 16) {         // binding sites and small picks: side-chain sticks
      for (let i = 0; i < b.atom_res.length && atoms.length < 6000; i++) {
        const r = b.atom_res[i];
        if (store.sel.has(r) && b.res_kind[r] <= KIND.nucleic) atoms.push(i);
      }
    }
    if (atoms.length) {
      const pos = new Float32Array(atoms.length * 3), ac = new Float32Array(atoms.length * 3), rad = new Float32Array(atoms.length).fill(0.2);
      atoms.forEach((a, i) => {
        pos.set(A.subarray(3 * a, 3 * a + 3), 3 * i);
        tmpC.set(atomColor(a)); ac[3 * i] = tmpC.r; ac[3 * i + 1] = tmpC.g; ac[3 * i + 2] = tmpC.b;
      });
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.BufferAttribute(pos, 3)); g.setAttribute('aColor', new THREE.BufferAttribute(ac, 3));
      g.setAttribute('aRadius', new THREE.BufferAttribute(rad, 1));
      this.groups.sel.add(new THREE.Points(g, impostorMaterial(this.uniforms)));
      const aset = new Set(atoms), B = b.bonds, pairs = [];
      for (let k = 0; k < B.length / 2; k++) if (aset.has(B[2 * k]) && aset.has(B[2 * k + 1])) pairs.push(B[2 * k], B[2 * k + 1]);
      if (pairs.length) this.groups.sel.add(bondMesh(A, pairs, atomColor, 0.09));
    }
    this.applyClip(this.groups.sel);
    this.request();
  }
}

// ------------------------------------------------------------------ geometry view
export class GeometryView extends View3D {
  constructor() {
    super('geo-host', 'geo');
    this.groups = { points: new THREE.Group(), cartoon: new THREE.Group(), overlay: new THREE.Group(),
      star: new THREE.Group(), het: new THREE.Group(), selgeo: new THREE.Group() };
    this.layerGroups = { hull: new THREE.Group(), delaunay: new THREE.Group(), alpha: new THREE.Group(), voronoi: new THREE.Group() };
    [...Object.values(this.groups), ...Object.values(this.layerGroups)].forEach((g) => this.scene.add(g));
    this.active = new Set(['delaunay']);
    this.layers = {};            // kind -> layer data
    this.showPoints = true; this.showCartoon = false; this.opacity = 0.5;
    this.showSimplices = true;   // highlight the nodes and edges of the Delaunay tetrahedra touching the selection
    this.pointScale = 1;         // geometry points are fixed-size dots (CSS px): show ▾ → small / medium / large
    this.showOverlay = true;     // pocket / interface overlays from the analyses
    this.star = null;            // {node: Uint8Array per point (2 selected, 1 neighbour), edges, inner, n_tets}
    this.showSelGeo = true;      // draw the selection-geometry overlay when one is computed
    this.legendEl = document.getElementById('geo-legend');
    this.statsEl = document.getElementById('geo-stats');
    this.overlay = null; this.pocketSel = null; this.contactPair = null; this.selPair = null;
    this.topGroup = new THREE.Group(); this.top.add(this.topGroup);
    store.on('data', () => this.load());
    store.on('tab', () => this.buildSelGeo());
    store.on('selection', () => this.updateSelection());
    store.on('selstyle', () => { this.updateSelection(); this.buildSelGeo(); });
    store.on('params', () => this.refreshLayers());
    store.on('selgeo', () => this.buildSelGeo());
  }
  load() {
    Object.values(this.groups).forEach(disposeGroup); Object.values(this.layerGroups).forEach(disposeGroup);
    this.layers = {}; this.overlay = null; this.selPair = null;
    const b = store.S.buf, poly = new Uint8Array(b.res_kind.length);
    for (let i = 0; i < b.pt_res.length; i++) poly[b.pt_res[i]] = 1;
    this.pickable = poly;
    this.star = null;
    this.buildPoints();
    this.buildHeteroGeo();
    this.frame(store.S.stats.extent, true);
    this.refreshLayers();
  }
  setLayer(kind, on) { if (on) this.active.add(kind); else this.active.delete(kind); this.refreshLayers(); }
  async refreshLayers() {
    if (!store.S) return;
    const job = store.job;
    for (const kind of Object.keys(this.layerGroups)) {
      if (store.job !== job) return;          // a newer structure took over; its own refresh runs
      const g = this.layerGroups[kind];
      if (!this.active.has(kind)) { g.visible = false; continue; }
      const params = store.params[kind] || {};
      const chip = document.querySelector(`#geo-toggles [data-layer="${kind}"]`);
      try {
        chip && chip.classList.add('busy');
        const L = await store.requestLayer(kind, params);
        if (!store.S || store.job !== job) return;
        if (!this.layers[kind] || this.layers[kind].key !== L.key) { this.layers[kind] = L; this.buildLayer(kind); }
        g.visible = this.active.has(kind);
      } catch (e) { toast(`${kind}: ${e.message}`); } finally { chip && chip.classList.remove('busy'); }
    }
    this.updateStats(); this.updateLegend(); this.updateSelection();
    this.request();
  }
  buildPoints() {
    disposeGroup(this.groups.points);
    const S = store.S, b = S.buf, n = b.pt_res.length;
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(b.pt_pos, 3));
    g.setAttribute('aColor', new THREE.BufferAttribute(new Float32Array(n * 3), 3));
    g.setAttribute('aRadius', new THREE.BufferAttribute(new Float32Array(n).fill(this.basePx()), 1));
    this.pts = new THREE.Points(g, dotMaterial(this.uniforms));
    this.groups.points.add(this.pts);
    this.groups.points.visible = this.showPoints;
    this.applyClip(this.groups.points);
    this.colorPoints();
  }
  colorPoints() {
    if (!this.pts) return;
    const b = store.S.buf, cols = store.residueColors(store.colorBy).map((c) => new THREE.Color(c));
    const a = this.pts.geometry.attributes.aColor.array, rad = this.pts.geometry.attributes.aRadius.array;
    const base = this.basePx(), any = store.sel.size > 0;
    const sc = selColor(), node = this.star && this.showSimplices ? this.star.node : null;
    for (let i = 0; i < b.pt_res.length; i++) {
      const r = b.pt_res[i], k = b.res_kind[r], sel = store.sel.has(r);
      let rr = base;
      if (k >= KIND.ligand) {                       // drawn as ball-and-stick on top (buildHeteroGeo): no dot
        a[3 * i] = a[3 * i + 1] = a[3 * i + 2] = 0; rad[i] = 0;
        continue;
      }
      tmpC.copy(cols[r]);
      if (any) {        // every node of the tetrahedra the selection takes part in is highlighted
        if (sel) { tmpC.lerp(sc, 0.6); rr = rr * 1.35 + 1; }                  // +1 px: the nodes of those tetrahedra
        else if (node && node[i]) { tmpC.lerp(sc, 0.45); rr = rr * 1.2 + 1; }
        else if (store.dim) { fade(tmpC); rr *= 0.85; }
      }
      a[3 * i] = tmpC.r; a[3 * i + 1] = tmpC.g; a[3 * i + 2] = tmpC.b;
      rad[i] = rr;
    }
    this.pts.geometry.attributes.aColor.needsUpdate = true;
    this.pts.geometry.attributes.aRadius.needsUpdate = true;
  }
  basePx() { return (store.S.stats.level === 'ca' || store.S.stats.level === 'centroid' ? 8.5 : 6.2) * this.pointScale; }
  setPointScale(k) { this.pointScale = k; this.colorPoints(); this.request(); }
  // ligands (ball-and-stick), ions and coordinated waters as in the structure view, drawn on top of the
  // tessellation so edges and nodes never hide them
  buildHeteroGeo() {
    disposeGroup(this.topGroup);
    this.hetExtra = buildHetBallStick(this.uniforms);
    this.hetMetal = false; this.ionElems = [];
    const H = this.hetExtra;
    if (!H) return;
    const sb = store.S.buf, els = store.S.stats.elements, ions = new Set();      // for the legend
    for (const a of H.atoms) {
      const k = sb.res_kind[sb.atom_res[a]], el = els[sb.atom_elem[a]];
      if (k === KIND.ion) ions.add(el); else if (k === KIND.ligand && hetColor(k, el) === LIGAND_METAL) this.hetMetal = true;
    }
    this.ionElems = [...ions].sort().slice(0, 6);
    this.topGroup.add(H.pts);
    if (H.bonds) this.topGroup.add(H.bonds);
    this.applyClip(this.topGroup);
  }
  setCartoon(on) {
    this.showCartoon = on;
    disposeGroup(this.groups.cartoon);
    if (on && store.S && store.S.buf.cart_pos.length) {
      const S = store.S, b = S.buf, P = b.cart_pos, rings = b.cart_rings, side = S.stats.cartoon.sides, tside = S.stats.cartoon.tube_sides;
      const idx = []; let base = 0;
      for (let s = 0; s < rings.length; s++) {
        const n = Math.abs(rings[s]), cols = (rings[s] < 0 ? tside : side) + 1;
        for (let i = 0; i < n - 1; i++) for (let j = 0; j < cols - 1; j++) { const a = base + i * cols + j, c = a + cols; idx.push(a, c, a + 1, a + 1, c, c + 1); }
        base += n * cols;
      }
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.BufferAttribute(P, 3));
      g.setIndex(P.length / 3 > 65535 ? new THREE.Uint32BufferAttribute(idx, 1) : new THREE.Uint16BufferAttribute(idx, 1));
      g.computeVertexNormals();
      this.groups.cartoon.add(new THREE.Mesh(g, new THREE.MeshPhongMaterial({ color: 0x8a8ab0, transparent: true, opacity: 0.28, depthWrite: false, side: THREE.DoubleSide })));
      this.applyClip(this.groups.cartoon);
    }
    this.request();
  }
  // Delaunay edge width (px) and opacity: dense all-atom meshes get finer, lighter lines
  meshLine() { const coarse = store.S.stats.level === 'ca' || store.S.stats.level === 'centroid'; return coarse ? [1.6, 0.85] : [1.1, 0.55]; }
  buildLayer(kind) {
    const g = this.layerGroups[kind]; disposeGroup(g);
    const L = this.layers[kind], P = store.S.buf.pt_pos;
    if (kind === 'delaunay') {
      const e = L.buf.edges, len = L.buf.edge_len, E = len.length, pos = new Float32Array(E * 6);
      for (let k = 0; k < E; k++) { pos.set(P.subarray(3 * e[2 * k], 3 * e[2 * k] + 3), 6 * k); pos.set(P.subarray(3 * e[2 * k + 1], 3 * e[2 * k + 1] + 3), 6 * k + 3); }
      const ls = fatLines(pos, new Float32Array(E * 6), ...this.meshLine());
      L.edgeKeys = null;                                                     // built when a selection needs it
      ls.userData.edges = true;
      g.add(ls);
      L.scale = () => MESH_EDGE;
      if (L.meta.params.show_tets && L.buf.tets) {
        const T = L.buf.tets, nT = T.length / 4, cap = Math.min(nT, +L.meta.params.max_tets || 20000);
        const step = nT / cap, list = Array.from({ length: cap }, (_, i) => Math.floor(i * step));
        const r = L.buf.tet_r, rs = d3.scaleSequential([0, L.meta.params.max_edge * 0.75], RAMP.seq);
        const m = tetMesh(P, T, list, 0.7, (t) => rs(Math.min(r[t], 99)), Math.min(0.55, this.opacity + 0.1));
        m.userData.allTets = true;
        g.add(m);
      }
    } else if (kind === 'voronoi' && L.buf.cell_pos) {
      const site = L.buf.cell_site, n = site.length;
      const geo = new THREE.BufferGeometry();
      geo.setAttribute('position', new THREE.BufferAttribute(L.buf.cell_pos, 3));
      geo.setAttribute('color', new THREE.BufferAttribute(new Float32Array(n * 3), 3));
      geo.computeVertexNormals();
      g.add(new THREE.Mesh(geo, meshMat(this.opacity, true)));
      const v = Array.from(L.buf.site_volume).filter(Number.isFinite).sort(d3.ascending);
      L.scale = d3.scaleSequential(v.length ? [d3.quantileSorted(v, 0.02), d3.quantileSorted(v, 0.98)] : [0, 1], RAMP.seq);
    } else if (kind === 'voronoi' && L.buf.skeleton) {
      const e = L.buf.skeleton, V = L.buf.vertices, E = e.length / 2, pos = new Float32Array(E * 6), col = new Float32Array(E * 6);
      tmpC.set(PALETTE.lavender);
      for (let k = 0; k < E; k++) {
        pos.set(V.subarray(3 * e[2 * k], 3 * e[2 * k] + 3), 6 * k); pos.set(V.subarray(3 * e[2 * k + 1], 3 * e[2 * k + 1] + 3), 6 * k + 3);
        for (let q = 0; q < 6; q += 3) { col[6 * k + q] = tmpC.r; col[6 * k + q + 1] = tmpC.g; col[6 * k + q + 2] = tmpC.b; }
      }
      g.add(fatLines(pos, col, 1.2, 0.7));
    } else if (kind === 'alpha') {
      g.add(surface(L.buf.faces, P, this.opacity));
    } else if (kind === 'hull') {
      const m = surface(L.buf.faces, P, Math.min(this.opacity, 0.22), true);
      g.add(m);
      const idxGeo = new THREE.BufferGeometry();
      idxGeo.setAttribute('position', new THREE.BufferAttribute(P, 3)); idxGeo.setIndex(new THREE.BufferAttribute(L.buf.faces, 1));
      g.add(new THREE.LineSegments(new THREE.EdgesGeometry(idxGeo, 1), new THREE.LineBasicMaterial({ color: PALETTE.teal, transparent: true, opacity: 0.65 })));
    }
    this.applyClip(g);
  }
  // re-tint layers so geometry touching the selection stands out; the tetrahedra touching it are shown by
  // their nodes (corners) and edges
  updateSelection() {
    if (!store.S) return;
    this.computeStar();
    this.colorPoints();
    const b = store.S.buf, ptRes = b.pt_res, any = store.sel.size > 0, sel = store.sel;
    const L = this.layers, c = new THREE.Color(), sc = selColor();
    if (L.delaunay && this.layerGroups.delaunay.children[0]) {
      // the edges of the tetrahedra the selection takes part in leave the mesh while it is shown (buildStar
      // draws them wider, in STAR_EDGE): the mesh's segments are rewritten with the others first and only
      // those are drawn, so no edge is drawn twice; the rest keep their colour
      const D = L.delaunay, e = D.buf.edges, len = D.buf.edge_len, E = len.length, P = b.pt_pos;
      const ls = this.layerGroups.delaunay.children[0], pos = ls.geometry.attributes.instanceStart.data.array, col = new Float32Array(E * 6);
      const star = this.star && this.showSimplices ? this.star : null, keys = this.edgeKeys();
      let n = 0;
      for (let k = 0; k < E; k++) {
        if (star && star.keys.has(keys[k])) continue;
        pos.set(P.subarray(3 * e[2 * k], 3 * e[2 * k] + 3), 6 * n); pos.set(P.subarray(3 * e[2 * k + 1], 3 * e[2 * k + 1] + 3), 6 * n + 3);
        c.set(D.scale(len[k]));
        if (star && store.dim) fade(c).multiplyScalar(0.7);           // otherwise the colours stay as they are
        for (let q = 0; q < 6; q += 3) { col[6 * n + q] = c.r; col[6 * n + q + 1] = c.g; col[6 * n + q + 2] = c.b; }
        n++;
      }
      ls.geometry.attributes.instanceStart.data.needsUpdate = true;
      ls.geometry.instanceCount = n;
      recolorFat(ls, col);
      if (star) star.moved = E - n;
    }
    this.buildStar();
    this.paintHeteroGeo();
    if (L.voronoi && L.voronoi.buf.cell_pos && this.layerGroups.voronoi.children[0]) {
      const site = L.voronoi.buf.cell_site, vol = L.voronoi.buf.site_volume, a = this.layerGroups.voronoi.children[0].geometry.attributes.color;
      for (let k = 0; k < site.length; k++) {
        const s = site[k], isSel = sel.has(ptRes[s]);
        c.set(L.voronoi.scale(vol[s]));
        if (any) { if (isSel) c.lerp(sc, SEL_MIX.layer); else if (store.dim) fade(c); }
        a.array[3 * k] = c.r; a.array[3 * k + 1] = c.g; a.array[3 * k + 2] = c.b;
      }
      a.needsUpdate = true;
    }
    for (const kind of ['alpha', 'hull']) {
      const mesh = this.layerGroups[kind].children[0];
      if (!L[kind] || !mesh) continue;
      const faces = L[kind].buf.faces, a = mesh.geometry.attributes.color, base = new THREE.Color(kind === 'alpha' ? PALETTE.mint : PALETTE.teal);
      for (let f = 0; f < faces.length / 3; f++) {
        let n = 0;
        for (let q = 0; q < 3; q++) if (sel.has(ptRes[faces[3 * f + q]])) n++;
        c.copy(base);
        if (any) { if (n) c.lerp(sc, SEL_MIX.layer * (0.6 + 0.2 * n)); else if (store.dim) c.multiplyScalar(0.6); }
        for (let q = 0; q < 3; q++) { a.array[9 * f + 3 * q] = c.r; a.array[9 * f + 3 * q + 1] = c.g; a.array[9 * f + 3 * q + 2] = c.b; }
      }
      a.needsUpdate = true;
    }
    this.updateStats(); this.updateLegend();
    this.request();
  }
  // the star of the selection in the Delaunay complex: every tetrahedron with a selected corner,
  // its four nodes and six edges (an edge is 'inner' when both ends are selected)
  // key (lower point index × n points + higher) of every shown Delaunay edge, cached per layer
  edgeKeys() {
    const L = this.layers.delaunay, NP = store.S.buf.pt_res.length;
    if (!L.edgeKeys) {
      const e = L.buf.edges, n = e.length / 2, keys = new Float64Array(n);
      for (let k = 0; k < n; k++) { const a = e[2 * k], z = e[2 * k + 1]; keys[k] = Math.min(a, z) * NP + Math.max(a, z); }
      L.edgeKeys = keys; L.edgeKeySet = new Set(keys);
    }
    return L.edgeKeys;
  }
  computeStar() {
    this.star = null;
    const L = this.layers.delaunay;
    if (!store.sel.size || !L || !L.buf.tets || !this.active.has('delaunay')) return;
    const T = L.buf.tets, nT = T.length / 4, ptRes = store.S.buf.pt_res, sel = store.sel, P = ptRes.length;
    const isSel = new Uint8Array(P);
    for (let i = 0; i < P; i++) if (sel.has(ptRes[i])) isSel[i] = 1;
    const node = new Uint8Array(P), seen = new Set(), edges = [];
    let nTets = 0;
    for (let t = 0; t < nT; t++) {
      const v0 = T[4 * t], v1 = T[4 * t + 1], v2 = T[4 * t + 2], v3 = T[4 * t + 3];
      if (!(isSel[v0] || isSel[v1] || isSel[v2] || isSel[v3])) continue;
      nTets++;
      const v = [v0, v1, v2, v3];
      for (const x of v) node[x] = isSel[x] ? 2 : Math.max(node[x], 1);
      for (let i = 0; i < 4; i++) for (let j = i + 1; j < 4; j++) {
        const a = Math.min(v[i], v[j]), z = Math.max(v[i], v[j]), key = a * P + z;
        if (!seen.has(key)) { seen.add(key); edges.push(a, z); }
      }
    }
    let inner = 0;
    for (let k = 0; k < edges.length; k += 2) if (isSel[edges[k]] && isSel[edges[k + 1]]) inner++;
    this.star = { node, edges, keys: seen, inner, nTets, nNodes: node.reduce((s, x) => s + (x ? 1 : 0), 0), moved: 0, extra: 0 };
  }
  buildStar() {
    disposeGroup(this.groups.star);
    if (!this.star || !this.showSimplices) { this.request(); return; }
    // every edge of the tetrahedra the selection takes part in: the mesh's own (taken out of the mesh in
    // updateSelection) and any longer than its cutoff, in STAR_EDGE, STAR_WIDEN px wider, same opacity
    const P = store.S.buf.pt_pos, { edges } = this.star, NP = store.S.buf.pt_res.length, col = new THREE.Color(STAR_EDGE);
    this.edgeKeys();
    const shown = this.layers.delaunay.edgeKeySet, list = [];
    let extra = 0;
    for (let k = 0; k < edges.length; k += 2) {
      if (!shown.has(edges[k] * NP + edges[k + 1])) { if (extra >= MAX_STAR_EDGES) continue; extra++; }
      list.push(edges[k], edges[k + 1]);
    }
    this.star.extra = extra;
    if (list.length) {
      const pos = new Float32Array(list.length * 3), cc = new Float32Array(list.length * 3);
      for (let i = 0; i < list.length; i++) { pos.set(P.subarray(3 * list[i], 3 * list[i] + 3), 3 * i); cc[3 * i] = col.r; cc[3 * i + 1] = col.g; cc[3 * i + 2] = col.b; }
      const [w, opacity] = this.meshLine();
      const ls = fatLines(pos, cc, w + STAR_WIDEN, opacity); ls.renderOrder = 2;
      this.groups.star.add(ls);
    }
    this.applyClip(this.groups.star);
    this.request();
  }
  // ligands / ions / waters: painted exactly as in the structure view (selected ones grow, colours stay)
  paintHeteroGeo() {
    disposeGroup(this.groups.het);
    if (this.hetExtra) paintHetBallStick(this.hetExtra);
  }
  // selection geometry (dock 'Selection' tab): each group's hull and the contacts between groups
  // (the selected atoms themselves are already highlighted)
  buildSelGeo() {
    disposeGroup(this.groups.selgeo);
    const L = store.selGeo;
    if (!L) this.selPair = null;
    if (!store.S) return;
    if (!L || !this.showSelGeo || store.activeTab !== 'selection') { this.updateLegend(); this.request(); return; }
    const P = L.buf.pts, gid = L.buf.pt_group;
    const gcol = (k) => GROUP_COLORS[k % GROUP_COLORS.length];
    if (L.buf.group_hull_faces) {
      const F = L.buf.group_hull_faces, fg = L.buf.group_hull_group;
      const m = surface(F, P, 0.07, true, (f) => gcol(fg[f]));          // a faint envelope per group, no wireframe
      m.material.depthWrite = false;
      this.groups.selgeo.add(m);
    }
    if (L.buf.contacts) {
      const C = L.buf.contacts, K = C.length / 2, pos = [], cc = [], sc = selColor().clone().lerp(MUTE, 0.3);
      for (let k = 0; k < K; k++) {
        const a = C[2 * k], z = C[2 * k + 1];
        if (this.selPair && !((gid[a] === this.selPair[0] && gid[z] === this.selPair[1]) || (gid[a] === this.selPair[1] && gid[z] === this.selPair[0]))) continue;
        pos.push(P[3 * a], P[3 * a + 1], P[3 * a + 2], P[3 * z], P[3 * z + 1], P[3 * z + 2]);
        cc.push(sc.r, sc.g, sc.b, sc.r, sc.g, sc.b);
      }
      if (pos.length) this.groups.selgeo.add(fatLines(new Float32Array(pos), new Float32Array(cc), 2));
    }
    this.applyClip(this.groups.selgeo);
    this.updateLegend();
    this.request();
  }
  setOpacity(op) {
    this.opacity = op;
    for (const kind of ['alpha', 'voronoi']) this.layerGroups[kind].traverse((o) => { if (o.isMesh) { o.material.opacity = op; o.material.transparent = op < 1; o.material.depthWrite = op >= 0.95; } });
    this.request();
  }
  setOverlayVisible(on) { this.showOverlay = on; this.groups.overlay.visible = on; this.request(); }
  setOverlay(L) {
    disposeGroup(this.groups.overlay);
    this.overlay = L;
    this.groups.overlay.visible = this.showOverlay;
    if (!L) { this.request(); return; }
    const P = store.S.buf.pt_pos;
    if (L.kind === 'pockets' && L.buf.faces) {
      const pal = CAT10;
      this.groups.overlay.add(surface(L.buf.faces, P, 0.85, true, (f) => {
        const id = L.buf.face_pocket[f];
        return this.pocketSel === null || this.pocketSel === id ? pal[id % 10] : '#3a3a52';
      }));
    } else if (L.kind === 'contacts') {
      const e = L.buf.res_interchain, R = store.S.buf.res_pos, E = e.length / 2, pos = [], col = [];
      const ch = store.S.buf.res_chain;
      tmpC.set(PALETTE.coral);
      for (let k = 0; k < E; k++) {
        if (this.contactPair) { const x = ch[e[2 * k]], y = ch[e[2 * k + 1]], [p, q] = this.contactPair; if (!((x === p && y === q) || (x === q && y === p))) continue; }
        pos.push(...R.subarray(3 * e[2 * k], 3 * e[2 * k] + 3), ...R.subarray(3 * e[2 * k + 1], 3 * e[2 * k + 1] + 3));
        col.push(tmpC.r, tmpC.g, tmpC.b, tmpC.r, tmpC.g, tmpC.b);
      }
      if (pos.length) this.groups.overlay.add(fatLines(new Float32Array(pos), new Float32Array(col), 2));
    }
    this.applyClip(this.groups.overlay);
    this.request();
  }
  updateStats() {
    const rows = [];
    const L = this.layers, f = d3.format(',.1f'), i = d3.format(',');
    if (this.active.has('delaunay') && L.delaunay) rows.push(`Delaunay <b>${i(L.delaunay.stats.n_tets)}</b> tets · ${i(L.delaunay.stats.n_shown)} edges ≤ ${L.delaunay.meta.params.max_edge} Å`);
    if (this.active.has('voronoi') && L.voronoi) rows.push(`Voronoi <b>${i(L.voronoi.stats.n_bounded)}</b> closed cells · mean ${f(L.voronoi.stats.mean_volume ?? 0)} Å³`);
    if (this.active.has('alpha') && L.alpha) rows.push(`α ${L.alpha.stats.alpha} Å · V <b>${i(Math.round(L.alpha.stats.volume))}</b> Å³ · ${L.alpha.stats.n_components} comp · ${L.alpha.stats.n_cavities} cav`);
    if (this.active.has('hull') && L.hull) rows.push(`hull V <b>${i(Math.round(L.hull.stats.volume))}</b> Å³ · ψ ${L.hull.stats.sphericity.toFixed(3)}`);
    if (this.star && this.showSimplices) rows.push(`<span class="selc">selection</span> <b>${i(this.star.nTets)}</b> tetrahedra · ${i(this.star.nNodes)} nodes · ${i(this.star.edges.length / 2)} edges (${i(this.star.inner)} inside)`);
    this.statsEl.innerHTML = rows.join('<br>');
    const st = document.getElementById('st-layer');
    const t = Object.values(L).filter((x) => this.active.has(x.kind)).map((x) => `${x.kind} ${x.stats.compute_s}s · ${fmt.bytes(x.stats.bytes)}`);
    st.textContent = t.join(' | ');
  }
  updateLegend() {
    const items = [];
    if (this.active.has('hull')) items.push(['hull', PALETTE.teal]);
    if (this.active.has('delaunay')) items.push(['Delaunay edges', MESH_EDGE]);
    if (this.active.has('alpha')) items.push(['α-shape', PALETTE.mint]);
    if (this.active.has('voronoi')) items.push(['Voronoi (by volume)', '#4fb39f']);
    if (this.star && this.showSimplices) items.push(['tetrahedra touching the selection', STAR_EDGE]);
    const st = store.S.stats;
    if (st.n_ligand) items.push(['ligand C', LIGAND_CARBON], ['bonds', LIGAND_BOND]);
    if (this.hetMetal) items.push(['metal centre', LIGAND_METAL]);
    for (const el of this.ionElems || []) items.push([el, ION_COLOR[el] || ION_DEFAULT]);
    if (st.n_water) items.push(['water', WATER_COLOR]);
    if (store.selGeo && this.showSelGeo) store.selGeo.stats.groups.slice(0, 6).forEach((gr, k) => items.push([gr.label, GROUP_COLORS[k % GROUP_COLORS.length]]));
    renderLegend(this.legendEl, items.length ? { type: 'cat', title: `geometry · ${store.S.stats.level} points (${fmt.int(store.S.stats.n_points)})`, items } : null);
  }
}

function surface(faces, P, opacity, flat = false, colorOf = null) {
  const F = faces.length / 3, pos = new Float32Array(F * 9), col = new Float32Array(F * 9);
  for (let f = 0; f < F; f++) for (let k = 0; k < 3; k++) {
    pos.set(P.subarray(3 * faces[3 * f + k], 3 * faces[3 * f + k] + 3), 9 * f + 3 * k);
    if (colorOf) { tmpC.set(colorOf(f)); col[9 * f + 3 * k] = tmpC.r; col[9 * f + 3 * k + 1] = tmpC.g; col[9 * f + 3 * k + 2] = tmpC.b; }
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  g.setAttribute('color', new THREE.BufferAttribute(col, 3));
  g.computeVertexNormals();
  return new THREE.Mesh(g, meshMat(opacity, flat));
}
