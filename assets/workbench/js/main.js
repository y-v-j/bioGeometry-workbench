// bioGeometry workbench — orchestration: messages from the Reflex panel, data loading, toolbars, status bar.
//
// Messages (postMessage {target:'bgworkbench', ...} or window.bgWorkbench.handle):
//   {type:'structure', api, job, key}     load a structure layer
//   {type:'params', params:{delaunay, voronoi, alpha, hull, contacts}}   geometry parameters
//   {type:'analysis', kind, key?, params?}  show an analysis in the dock (pockets, contacts, filtration, ...)
//   {type:'display', opts:{molStyle, molColor, geoOpacity}}
//   {type:'reset'}
// Standalone: index.html?api=http://127.0.0.1:8770&pdb=4HHB | &af=P69905 [&level=ca]
import { store, busy, toast, fmt, SEL_PRESETS } from './store.js';
import { SequenceViewer } from './sequence.js';
import { StructureView, GeometryView } from './view3d.js';
import { RamaPanel } from './rama.js';
import { Dock } from './dock.js';

const d3 = window.d3;
const $ = (id) => document.getElementById(id);

const seq = new SequenceViewer();
const mol = new StructureView();
const geo = new GeometryView();
mol.peer = geo; geo.peer = mol;
mol.linked = geo.linked = true;          // the two 3D views rotate / zoom / pan together (⇄ link)
const rama = new RamaPanel();
const dock = new Dock(geo);

// ------------------------------------------------------------------ loading
async function loadStructure(msg) {
  if (msg.api !== undefined) store.api = msg.api.replace(/\/$/, '');
  busy(true, 'Loading structure…');
  try {
    store.job = msg.job;
    store.layers.clear();
    store.F = null; store.T = null;
    const S = await store.fetchLayer(msg.key);
    store.S = S;
    store.sel = new Set(); store.hover = null; store.anchor = null;
    S.colIndex = {};
    $('empty').hidden = true;
    store.emit('data');
    updateStatus();
    notifyParent({ type: 'loaded', job: msg.job });
    loadFeatures();
  } catch (e) { toast(`Could not load structure: ${e.message}`); } finally { busy(false); }
}

async function loadFeatures() {
  const job = store.job;
  try {
    const F = await store.requestLayer('features');
    if (store.job !== job) return;
    F.colIndex = Object.fromEntries(F.stats.columns.map((c) => [c.key, c]));
    store.F = F;
    store.emit('features');
    updateStatus();
  } catch (e) { toast(`features: ${e.message}`); }
}

async function showAnalysis(msg) {
  if (!store.S) return;
  busy(true, `Running ${msg.kind}…`);
  try {
    const job = store.job;
    const L = msg.key ? await store.fetchLayer(msg.key) : await store.requestLayer(msg.kind, msg.params || {});
    if (store.job === job) dock.showAnalysis(L);
  } catch (e) { toast(`${msg.kind}: ${e.message}`); } finally { busy(false); }
}

function setParams(p) {
  store.params = { ...store.params, ...p };
  store.emit('params');
}

function setDisplay(o) {
  if (o.molStyle) { $('mol-style').value = o.molStyle; mol.setStyle(o.molStyle); }
  if (o.molColor) { $('mol-color').value = o.molColor; store.colorBy = o.molColor; mol.setColor(o.molColor); geo.colorPoints(); geo.request(); }
  if (o.geoOpacity !== undefined) geo.setOpacity(+o.geoOpacity);
  if (o.selColor && SEL_PRESETS[o.selColor]) { $('sel-color').value = o.selColor; store.setSelColor(SEL_PRESETS[o.selColor]); }
}

// messages are applied in order, so a 'params' arriving during a structure download waits for it
let queue = Promise.resolve();
function handle(msg) {
  queue = queue.then(() => handleNow(msg)).catch((e) => toast(String(e.message || e)));
  return queue;
}
function handleNow(msg) {
  if (!msg || typeof msg !== 'object') return undefined;
  switch (msg.type) {
    case 'structure': return loadStructure(msg);
    case 'params': return setParams(msg.params || {});
    case 'analysis': return msg.kind === 'selection' ? dock.measureSelection(true) : showAnalysis(msg);
    case 'display': return setDisplay(msg.opts || {});
    case 'reset': mol.frame(store.S?.stats.extent); geo.frame(store.S?.stats.extent); return undefined;
    case 'resetall': return resetWorkbench(msg);
    case 'select': store.setSelection(msg.residues || [], 'parent'); return undefined;
    default: return undefined;
  }
}
window.addEventListener('message', (ev) => { if (ev.data && ev.data.target === 'bgworkbench') handle(ev.data); });
function notifyParent(msg) { try { if (window.parent !== window) window.parent.postMessage({ source: 'bgworkbench', ...msg }, '*'); } catch (e) { /* ignore */ } }

// ------------------------------------------------------------------ toolbars
$('mol-style').addEventListener('change', (e) => mol.setStyle(e.target.value));
$('mol-color').addEventListener('change', (e) => { store.colorBy = e.target.value; mol.setColor(e.target.value); geo.colorPoints(); geo.request(); seqRecolor(); });
function seqRecolor() { /* sequence keeps its own scheme selector */ }
$('sel-color').addEventListener('change', (e) => store.setSelColor(SEL_PRESETS[e.target.value]));
try {
  const saved = localStorage.getItem('bg.selColor.v2'), name = Object.keys(SEL_PRESETS).find((k) => SEL_PRESETS[k] === saved);
  if (name) { $('sel-color').value = name; store.setSelColor(saved); }
} catch (e) { /* storage unavailable */ }
$('seq-measure').addEventListener('click', () => dock.measureSelection(true));
$('seq-multi').addEventListener('click', (e) => { store.multi = !store.multi; e.currentTarget.classList.toggle('on', store.multi); });
store.on('measure', () => dock.measureSelection(true));
store.on('focus-selection', () => { if (store.sel.size) { mol.focusResidues([...store.sel]); if (!geo.linked) geo.focusResidues([...store.sel]); } });
document.querySelectorAll('#p-mol [data-tog]').forEach((b) => b.addEventListener('click', () => {
  b.classList.toggle('on');
  const on = b.classList.contains('on'), t = b.dataset.tog;
  if (t === 'dim') { store.dim = mol.dim = on; mol.paint(); geo.updateSelection(); seq.draw(); rama.draw(); }
  if (t === 'link') { mol.linked = geo.linked = on; if (on) mol.syncPeer(); }
  if (t === 'spin') mol.spin = on;
  if (t === 'slab') { mol.setSlab(on); if (geo.linked) geo.setSlab(on); }
  if (t === 'bind') mol.setBinding(on);
}));
document.querySelectorAll('#p-mol [data-act]').forEach((b) => b.addEventListener('click', () => {
  if (b.dataset.act === 'reset') mol.frame(store.S?.stats.extent);
  if (b.dataset.act === 'shot') mol.screenshot();
}));
document.querySelectorAll('#geo-toggles [data-layer]').forEach((b) => b.addEventListener('click', () => {
  b.classList.toggle('on');
  geo.setLayer(b.dataset.layer, b.classList.contains('on'));
  notifyParent({ type: 'layers', active: [...geo.active] });
}));
// geometry view: what to draw besides the layers
const GEO_SHOW = [
  ['geometry points', () => geo.showPoints, (on) => { geo.showPoints = on; geo.groups.points.visible = on; geo.request(); }],
  ['point size: small', () => geo.pointScale === 0.75, () => { geo.setPointScale(0.75); return true; }],
  ['point size: medium', () => geo.pointScale === 1, () => { geo.setPointScale(1); return true; }],
  ['point size: large', () => geo.pointScale === 1.35, () => { geo.setPointScale(1.35); return true; }],
  ['analysis overlay (pockets / interfaces)', () => geo.showOverlay, (on) => geo.setOverlayVisible(on)],
  ['cartoon underneath', () => geo.showCartoon, (on) => geo.setCartoon(on)],
  ['tetrahedra touching the selection (nodes + edges)', () => geo.showSimplices, (on) => { geo.showSimplices = on; geo.updateSelection(); }],
  ['selection geometry (⌬ measure: group hulls, contacts)', () => geo.showSelGeo, (on) => { geo.showSelGeo = on; geo.buildSelGeo(); }],
];
$('geo-show').addEventListener('click', (e) => popupMenu(e.currentTarget, GEO_SHOW));
function popupMenu(anchor, items) {
  const menu = d3.select('#menu').html('').attr('hidden', null);
  items.forEach(([text, get, set]) => {
    const l = menu.append('label');
    l.append('input').attr('type', 'checkbox').property('checked', get())
      .on('change', (ev) => { if (set(ev.target.checked)) popupMenu(anchor, items); });   // radio-like items re-render
    l.append('span').text(text);
  });
  const r = anchor.getBoundingClientRect();
  menu.style('left', `${Math.min(r.left, window.innerWidth - 300)}px`).style('top', `${r.bottom + 4}px`);
  const close = (ev) => { if (!document.getElementById('menu').contains(ev.target) && ev.target !== anchor) { menu.attr('hidden', true); window.removeEventListener('pointerdown', close, true); } };
  setTimeout(() => window.addEventListener('pointerdown', close, true), 0);
}
document.querySelectorAll('#geo-toggles [data-act]').forEach((b) => b.addEventListener('click', () => {
  if (b.dataset.act === 'reset') geo.frame(store.S?.stats.extent);
  if (b.dataset.act === 'shot') geo.screenshot();
}));
window.addEventListener('keydown', (e) => {
  if (['INPUT', 'SELECT', 'TEXTAREA'].includes(e.target.tagName)) return;
  const k = e.key.toLowerCase();
  if (k === 'escape') store.clearAll('key');
  if (k === 'm') $('seq-multi').click();
  if (k === 'r') { mol.frame(store.S?.stats.extent); geo.frame(store.S?.stats.extent); }
  if (k === 'f' && store.sel.size) { mol.focusResidues([...store.sel]); geo.focusResidues([...store.sel]); }
  if (k === 's') document.querySelector('#p-mol [data-tog="spin"]').click();
  if (k === 'c') document.querySelector('#p-mol [data-tog="slab"]').click();
  if (k === 'l') document.querySelector('#p-mol [data-tog="link"]').click();
  if (k === 'g') dock.measureSelection(true);
});

// dock resize
d3.select('#dock-resize').call(d3.drag().on('drag', (ev) => {
  const h = Math.max(110, Math.min(window.innerHeight - 300, window.innerHeight - ev.sourceEvent.clientY - 26));
  document.documentElement.style.setProperty('--table-h', `${h}px`);
}));

// ------------------------------------------------------------------ status bar + render loop
store.on('selection', () => updateStatus());
function updateStatus() {
  const S = store.S;
  if (!S) return;
  const s = S.stats;
  const src = s.is_alphafold ? `AlphaFold ${s.meta.af_entry || ''} · pLDDT ${s.meta.mean_plddt ?? '–'}` : `${s.source.split(':')[0].toUpperCase()} ${s.name}${s.meta.resolution ? ` · ${s.meta.resolution} Å` : ''}`;
  $('st-struct').textContent = `${src} · ${fmt.int(s.n_atoms)} atoms · ${fmt.int(s.n_protein + s.n_nucleic)} residues · ${s.n_chains} ch · ${s.n_ligand + s.n_ion} hetero · ${s.n_water} coord. H₂O · ${s.level} ${fmt.int(s.n_points)} pts` +
    (store.F ? '' : ' · features…');
  $('st-sel').textContent = store.selectionSummary();
  $('mol-info').textContent = s.title ? s.title.slice(0, 90) : '';
}
const spark = [];
function drawSpark() {
  const svg = d3.select('#st-spark'), w = 90, h = 14;
  const x = d3.scaleLinear([0, 59], [0, w]), y = d3.scaleLinear([0, 70], [h - 1, 1]);
  svg.selectAll('path').data([spark]).join('path').attr('d', d3.line().x((d, i) => x(i)).y((d) => y(Math.min(70, d))));
}
let lastStat = 0;
function loop(now) {
  requestAnimationFrame(loop);
  mol.tick(now); geo.tick(now);
  if (now - lastStat > 1000) {
    lastStat = now;
    $('st-fps').textContent = `structure ${mol.fps || 'idle'}${mol.fps ? ' fps' : ''} · geometry ${geo.fps || 'idle'}${geo.fps ? ' fps' : ''}`;
    spark.push(Math.max(mol.fps, geo.fps)); if (spark.length > 60) spark.shift();
    drawSpark();
    mol.fps = mol.frames ? mol.fps : 0; geo.fps = geo.frames ? geo.fps : 0;
  }
}
requestAnimationFrame(loop);
busy(false);

// ------------------------------------------------------------------ reset: a clean slate with the loaded structure
async function resetWorkbench(msg = {}) {
  if (!store.S) return;
  store.sel = new Set(); store.hover = null; store.anchor = null; store.selGeo = null;
  Object.assign(geo, { pocketSel: null, contactPair: null, selPair: null, overlay: null, showPoints: true, showSimplices: true,
    showSelGeo: true, showOverlay: true, pointScale: 1 });
  geo.active = new Set(['delaunay']);
  document.querySelectorAll('#geo-toggles [data-layer]').forEach((b) => b.classList.toggle('on', b.dataset.layer === 'delaunay'));
  if (geo.showCartoon) geo.setCartoon(false);
  const tog = (name, on) => { const b = document.querySelector(`#p-mol [data-tog="${name}"]`); if (b) b.classList.toggle('on', on); };
  tog('dim', false); tog('link', true); tog('bind', true); tog('spin', false); tog('slab', false);
  store.dim = mol.dim = false; mol.linked = geo.linked = true; mol.spin = false; mol.showBinding = true;
  store.multi = false; $('seq-multi').classList.remove('on');
  if (mol.slab) mol.setSlab(false);
  if (geo.slab) geo.setSlab(false);
  $('seq-color').value = 'class';
  if (rama.tab !== 'rama') rama.setTab('rama');
  document.querySelector('[data-tog="selonly"]').classList.remove('on');
  if (msg.params) store.params = { ...store.params, ...msg.params };
  if (msg.display) setDisplay(msg.display);
  await loadStructure({ job: store.job, key: store.S.key });        // every view, table and tab rebuilds from scratch
  toast('Workbench reset', 'info');
}
document.addEventListener('keydown', (e) => { if (e.key === 'X' && e.shiftKey && !['INPUT', 'SELECT', 'TEXTAREA'].includes(e.target.tagName)) resetWorkbench(); });

window.bgWorkbench = { handle, store, ready: true, views: { seq, mol, geo, rama, dock }, reset: resetWorkbench };
notifyParent({ type: 'ready' });

// ------------------------------------------------------------------ standalone mode
(async () => {
  const q = new URLSearchParams(location.search);
  if (!q.get('api') || !(q.get('pdb') || q.get('af') || q.get('path'))) return;
  store.api = q.get('api');
  busy(true, 'Opening structure…');
  try {
    const spec = q.get('af') ? { source: 'alphafold', uniprot: q.get('af') }
      : q.get('path') ? { source: 'file', path: q.get('path') } : { source: 'rcsb', pdb_id: q.get('pdb') };
    // mode=allatom | ca: the two analysis modes (polymer level + ligands, ions and coordinated waters as atoms)
    const mode = q.get('mode') || 'allatom';
    Object.assign(spec, { chains: q.get('chains') || '', level: q.get('level') || (mode === 'ca' ? 'ca' : 'heavy'),
      assembly: q.get('assembly') || '', ligands: q.get('ligands') !== '0' });
    const res = await store.apiCall('/bg/jobs', spec);
    setParams(res.summary.defaults);
    await loadStructure({ job: res.summary.job, key: res.structure.key });
    for (const layer of (q.get('layers') || '').split(',').filter(Boolean)) document.querySelector(`#geo-toggles [data-layer="${layer}"]`)?.click();
    if (q.get('select')) {
      const want = q.get('select').split(',');
      const rs = []; for (let r = 0; r < store.R; r++) if (want.includes(String(store.S.buf.res_seq[r])) && store.S.buf.res_chain[r] === store.focusChain) rs.push(r);
      store.setSelection(rs, 'url');
    }
    for (const a of (q.get('analysis') || '').split(',').filter(Boolean)) await showAnalysis({ kind: a });
    if (q.get('tab')) await dock.activate(q.get('tab'));
    window.bgStandaloneDone = true;
  } catch (e) { toast(e.message); window.bgStandaloneDone = `error: ${e.message}`; } finally { busy(false); }
})();
