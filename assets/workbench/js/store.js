// bioGeometry workbench — shared store: data access, the selection bus and colour scales.
//
// Every component reads residues by their index in the structure's residue table (0..R-1).
// The selection is one Set of residue indices; components publish changes through the bus and
// subscribe to 'selection', 'hover', 'focus' (current chain) and 'data' events.
const d3 = window.d3;

export const KIND = { protein: 0, nucleic: 1, ligand: 2, ion: 3, water: 4 };
export const PALETTE = {
  bg: '#191926', card: '#20202f', edge: '#2b2b40', fg: '#e8e8f2', label: '#9a9ab8', dim: '#62627e',
  sky: '#7dd3fc', lavender: '#c4b5fd', rose: '#f9a8d4', mint: '#86efac', teal: '#5eead4', coral: '#fca5a5',
  sel: '#00e5ff',     // selection accent, shared by every view (store.setSelColor); cyan as in the Plotly dashboard
};
// saturated on purpose: every colour scheme is pastel, so the selection never blends into one
export const SEL_PRESETS = { cyan: '#00e5ff', magenta: '#ff5ac8', violet: '#a78bfa', lime: '#a3e635', white: '#f4f4ff' };
// categorical colours for selection-geometry groups (domain, ligand, ion ...)
export const GROUP_COLORS = ['#7dd3fc', '#34d399', '#c4b5fd', '#86efac', '#fca5a5', '#67e8f9', '#f9a8d4', '#f0abfc', '#a5b4fc', '#5eead4'];
// no gold, amber or yellow anywhere: these replace viridis / plasma / cividis / turbo / YlOrRd / RdYlBu / Tableau10
export const RAMP = {
  seq: d3.interpolateRgbBasis(['#2e1065', '#6d28d9', '#a78bfa', '#7dd3fc', '#e0f2fe']),     // low → high (deep violet → pale sky)
  heat: d3.interpolateRgbBasis(['#c7d2fe', '#a78bfa', '#db2777', '#991b1b']),               // B-factor: low → high
  access: d3.interpolateRgbBasis(['#1e3a8a', '#7dd3fc', '#f1f5f9']),                        // burial: buried → exposed
  div: d3.interpolateRgbBasis(['#ef4444', '#f1f5f9', '#3b82f6']),                           // hydropathy: hydrophobic → hydrophilic
  rainbow: d3.interpolateRgbBasis(['#3b82f6', '#22d3ee', '#34d399', '#f472b6', '#ef4444']), // N → C
};
export const CAT10 = ['#4e79a7', '#e15759', '#76b7b2', '#59a14f', '#b07aa1', '#ff9da7', '#a0a0c8', '#86bcb6', '#8cd17d', '#d4a6c8'];
export const SS_COLOR = { 1: '#f9a8d4', 2: '#5eead4', 0: '#7dd3fc', [-1]: '#62627e' };   // H rose, E teal, C sky
export const SS_NAME = { 1: 'helix', 2: 'strand', 0: 'coil', [-1]: '–' };
export const CLASS_COLOR = {
  aliphatic: '#86efac', aromatic: '#c4b5fd', polar: '#7dd3fc', acidic: '#fca5a5', basic: '#a5b4fc',
  special: '#5eead4', purine: '#f9a8d4', pyrimidine: '#f0abfc', ligand: '#c9b26e', ion: '#d0bbff', water: '#d65f5f',
  other: '#9a9ab8',
};
export const KD = { ILE: 4.5, VAL: 4.2, LEU: 3.8, PHE: 2.8, CYS: 2.5, MET: 1.9, ALA: 1.8, GLY: -0.4, THR: -0.7,
  SER: -0.8, TRP: -0.9, TYR: -1.3, PRO: -1.6, HIS: -3.2, GLU: -3.5, GLN: -3.5, ASP: -3.5, ASN: -3.5, LYS: -3.9, ARG: -4.5 };
export const ELEMENT_COLOR = { C: '#b8bccc', N: '#7aa2f7', O: '#f7768e', S: '#c3e88d', P: '#f78c6c', Se: '#f0a0a0',
  Mg: '#9ece6a', Zn: '#7dcfff', Fe: '#e0915a', Ca: '#73daca', Mn: '#bb9af7', Cl: '#9ece6a', Na: '#a9b1d6', K: '#c0caf5',
  Cu: '#d08770', Co: '#d78bce', Ni: '#8bd7a5', Br: '#c47755', I: '#9b59b6', F: '#90e050', X: '#ff79c6' };
// ligands, ions and coordinated waters look identical in both 3D views and differ from every colour scheme.
// HEM and other complex ligands, by role: carbons, bonds (sticks) and the metal centre in warm earth tones
// (the user's exception to 'no gold'); N, O, S … keep element colours so the chemistry stays readable
export const LIGAND_CARBON = '#c9b26e', LIGAND_BOND = '#bd8b64', LIGAND_METAL = '#c15b2a';
const METALS = new Set(['Fe', 'Mg', 'Zn', 'Co', 'Cu', 'Mn', 'Ni', 'Mo', 'W', 'V', 'Cr', 'Cd', 'Hg', 'Pt', 'Pd', 'Ru', 'Rh', 'Ir', 'Os', 'Au', 'Ag', 'Ga', 'Al', 'Ca', 'Na', 'K', 'Li', 'Sr', 'Ba']);
const LIGAND_ELEMENT = { N: '#5b8cff', O: '#ff5252', S: '#c3e88d', P: '#f78c6c', F: '#90e050', Cl: '#a6d854', Br: '#8c8c8c', I: '#e78ac3' };
// ions and waters from seaborn's categorical palettes (deep / muted / pastel / colorblind / Set2), grouped by
// chemistry and chosen to stay clear of the blue tessellation and the warm ligands (CIELAB ΔE > 20 to each)
export const ION_COLOR = {
  Zn: '#d0bbff', Cd: '#cfcfcf', Hg: '#cfcfcf',                                   // group 12: pastel violet / grey
  Mg: '#6acc64', Ca: '#8de5a1', Sr: '#55a868', Ba: '#55a868',                    // alkaline earth: greens
  Na: '#956cb4', K: '#8172b3', Li: '#d0bbff', Rb: '#8172b3', Cs: '#956cb4',      // alkali: purples
  Fe: LIGAND_METAL, Cu: '#ffb482', Mn: '#dc7ec0', Co: '#fab0e4', Ni: '#029e73',  // transition metals (Fe as in haem)
  Cl: '#a6d854', Br: '#8c8c8c', I: '#e78ac3', F: '#66c2a5',                      // halides
};
export const ION_DEFAULT = '#cfcfcf';
export const WATER_COLOR = '#d65f5f';                                            // seaborn muted red: water oxygen
export function hetColor(kind, el) {
  if (kind === 4) return WATER_COLOR;
  if (kind === 3) return ION_COLOR[el] || ION_DEFAULT;
  if (el === 'C') return LIGAND_CARBON;
  return METALS.has(el) ? LIGAND_METAL : LIGAND_ELEMENT[el] || ELEMENT_COLOR[el] || '#ff79c6';
}
// binding contacts (dotted lines in the structure view), indexed by the server's bind_type
export const BIND_COLOR = ['#c4b5fd', '#f9a8d4', '#7dd3fc', '#e5e7eb'];   // metal · H-bond · water bridge · covalent
export const PLDDT_COLOR = (v) => (v > 90 ? '#0053d6' : v > 70 ? '#65cbf3' : v > 50 ? '#c084fc' : '#ff7d45');     // AlphaFold bands, with violet instead of yellow
export const chainColor = (i, n) => {
  const base = ['#5eabd9', '#f9a8d4', '#86efac', '#f0abfc', '#c4b5fd', '#fca5a5', '#67e8f9', '#93c5fd', '#a5b4fc', '#5eead4'];
  // many chains: golden-angle hues, skipping the orange-yellow band (20°-80°)
  return n <= base.length ? base[i % base.length] : d3.hsl(80 + ((i * 137.508) % 300), 0.65, 0.72).formatHex();
};

const DTYPES = { float32: Float32Array, float64: Float64Array, uint32: Uint32Array, int32: Int32Array,
  uint16: Uint16Array, int16: Int16Array, uint8: Uint8Array, int8: Int8Array };

export const fmt = {
  int: d3.format(','), f1: d3.format(',.1f'), f2: d3.format(',.2f'),
  bytes: (b) => (b > 1e6 ? `${(b / 1e6).toFixed(1)} MB` : `${(b / 1e3).toFixed(0)} kB`),
};

// ------------------------------------------------------------------ the store
class Store {
  constructor() {
    this.api = '';
    this.job = null;
    this.S = null;              // structure layer {meta, stats, buf}
    this.F = null;              // residue feature table layer
    this.T = null;              // tetrahedra layer
    this.layers = new Map();    // cache of fetched layers
    this.params = {};           // geometry parameters from the control panel
    this.sel = new Set();
    this.hover = null;          // residue index or null
    this.focusChain = 0;        // chain index shown in the sequence viewer
    this.anchor = null;         // last clicked residue (for shift ranges)
    this.handlers = {};
    this.colorBy = 'chain';
    this.selGeo = null;         // selection-geometry layer for the current selection (dock 'Selection' tab)
    this.multi = false;         // multi-select mode (sequence header ＋ multi, key M): clicks add / toggle
    this.dim = false;           // fade what is not selected (structure toolbar 'dim'); off: colours stay as they are
  }

  setSelColor(c) {
    PALETTE.sel = c;
    document.documentElement.style.setProperty('--sel', c);
    try { localStorage.setItem('bg.selColor.v2', c); } catch (e) { /* private window */ }
    this.emit('selstyle');
  }

  on(evt, fn) { (this.handlers[evt] ||= []).push(fn); }
  emit(evt, payload) { (this.handlers[evt] || []).forEach((fn) => { try { fn(payload); } catch (e) { console.error(evt, e); } }); }

  // ---------------- data
  async apiCall(path, body) {
    const opt = body === undefined ? {} : { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) };
    const r = await fetch(this.api + path, opt);
    const j = await r.json().catch(() => ({ error: `HTTP ${r.status}` }));
    if (!r.ok || j.error) throw new Error(j.error || `HTTP ${r.status}`);
    return j;
  }

  // `job` is pinned by the caller so a request started before a structure switch never mixes jobs
  async fetchLayer(key, job = this.job) {
    const ck = `${job}/${key}`;
    if (this.layers.has(ck)) return this.layers.get(ck);
    const meta = await this.apiCall(`/bg/jobs/${job}/layers/${key}/meta`);
    const names = Object.keys(meta.buffers);
    const arrays = await Promise.all(names.map(async (n) => {
      const r = await fetch(`${this.api}/bg/jobs/${job}/layers/${key}/buf/${n}`);
      if (!r.ok) throw new Error(`buffer ${n}: HTTP ${r.status}`);
      return new DTYPES[meta.buffers[n].dtype](await r.arrayBuffer());
    }));
    const layer = { meta, stats: meta.stats, kind: meta.kind, key, buf: Object.fromEntries(names.map((n, i) => [n, arrays[i]])) };
    this.layers.set(ck, layer);
    if (this.layers.size > 40) this.layers.delete(this.layers.keys().next().value);
    return layer;
  }

  async requestLayer(kind, params = {}) {
    const job = this.job;
    const meta = await this.apiCall(`/bg/jobs/${job}/layers`, { kind, params });
    return this.fetchLayer(meta.key, job);
  }

  // ---------------- residue helpers
  get R() { return this.S ? this.S.buf.res_kind.length : 0; }
  resName(r) { return this.S.stats.res_names[this.S.buf.res_name[r]]; }
  chainId(r) { return this.S.stats.chains[this.S.buf.res_chain[r]]; }
  label(r) {
    const b = this.S.buf, ic = this.S.stats.icodes[b.res_icode[r]] || '';
    return `${this.chainId(r)}:${this.resName(r)}${b.res_seq[r]}${ic}`;
  }
  chainResidues(ci) {
    const out = [], ch = this.S.buf.res_chain;
    for (let r = 0; r < ch.length; r++) if (ch[r] === ci) out.push(r);
    return out;
  }
  // residues of a domain annotation (polymer residues of its chain between r0 and r1)
  domainResidues(d) {
    const b = this.S.buf, out = [];
    for (let r = d.r0; r <= d.r1; r++) if (b.res_kind[r] <= KIND.nucleic) out.push(r);
    return out;
  }
  featureValue(key, r) {
    if (!this.F) return NaN;
    const col = this.F.colIndex[key];
    if (!col) return NaN;
    const v = this.F.buf[`c_${key}`][r];
    return col.type === 'cat' ? col.categories[v] : v;
  }

  // ---------------- selection bus
  setSelection(rs, source = '') {
    this.sel = new Set(rs);
    this.emit('selection', { source });
  }
  toggle(r, source = '') {
    if (this.sel.has(r)) this.sel.delete(r); else this.sel.add(r);
    this.anchor = r;
    this.emit('selection', { source });
  }
  add(rs, source = '') { rs.forEach((r) => this.sel.add(r)); this.emit('selection', { source }); }
  remove(rs, source = '') { rs.forEach((r) => this.sel.delete(r)); this.emit('selection', { source }); }
  // selection rules shared by every component: plain = replace; ctrl/cmd = toggle; shift = add (range);
  // in multi-select mode a plain click toggles and drags / brushes add
  toggling(ev) { return !!(ev && (ev.ctrlKey || ev.metaKey)) || (this.multi && !(ev && ev.shiftKey)); }
  additive(ev) { return this.multi || !!(ev && (ev.shiftKey || ev.ctrlKey || ev.metaKey)); }
  addOrSet(rs, ev, source = '') { if (this.additive(ev)) this.add(rs, source); else this.setSelection(rs, source); }
  // a group (domain, chain, ligand, pocket …) picked: toggled, added or made the selection
  pick(rs, ev, source = '') {
    if (this.toggling(ev)) { if (rs.every((r) => this.sel.has(r))) this.remove(rs, source); else this.add(rs, source); }
    else if (ev && ev.shiftKey) this.add(rs, source);
    else this.setSelection(rs, source);
  }
  // click semantics shared by every component: plain = only this, ctrl/cmd = toggle, shift = range/add
  click(r, ev, source = '') {
    if (r === null || r === undefined) return;
    if (this.toggling(ev)) { this.toggle(r, source); return; }
    if (ev && ev.shiftKey && this.anchor !== null && this.S.buf.res_chain[this.anchor] === this.S.buf.res_chain[r]) {
      const [a, b] = [Math.min(this.anchor, r), Math.max(this.anchor, r)];
      const rs = []; for (let i = a; i <= b; i++) rs.push(i);
      this.add(rs, source); return;
    }
    if (ev && ev.shiftKey) { this.add([r], source); this.anchor = r; return; }
    this.anchor = r;
    this.setSelection([r], source);
  }
  clear(source = '') { this.sel.clear(); this.anchor = null; this.emit('selection', { source }); }
  // the clear button / Esc: the selection, its measurement and the analysis highlights
  clearAll(source = '') { this.selGeo = null; this.clear(source); this.emit('selgeo'); }
  setHover(r, source = '') { if (r === this.hover) return; this.hover = r; this.emit('hover', { source }); }
  setFocusChain(ci, source = '') { if (ci === this.focusChain) return; this.focusChain = ci; this.emit('focus', { source }); }

  selectionSummary() {
    const n = this.sel.size;
    if (!n || !this.S) return 'selection: none';
    const chains = new Set(), het = [];
    for (const r of this.sel) {
      chains.add(this.chainId(r));
      if (this.S.buf.res_kind[r] >= KIND.ligand && het.length < 4) het.push(this.resName(r));
    }
    const first = [...this.sel].sort((a, b) => a - b).slice(0, 6).map((r) => this.label(r)).join(' ');
    return `selection: ${n} residue${n > 1 ? 's' : ''} · ${chains.size} chain${chains.size > 1 ? 's' : ''}` +
      (het.length ? ` · ${het.join(', ')}` : '') + ` · ${first}${n > 6 ? ' …' : ''}`;
  }

  // ---------------- colours (one CSS colour per residue for the active scheme)
  residueColors(scheme) {
    const S = this.S, b = S.buf, R = this.R, nCh = S.stats.chains.length, out = new Array(R);
    const isAF = S.stats.is_alphafold;
    const robust = (arr) => {
      const v = Array.from(arr).filter(Number.isFinite).sort(d3.ascending);
      return v.length ? [d3.quantileSorted(v, 0.02), d3.quantileSorted(v, 0.98)] : [0, 1];
    };
    if (scheme === 'chain') for (let r = 0; r < R; r++) out[r] = chainColor(b.res_chain[r], nCh);
    else if (scheme === 'ss') for (let r = 0; r < R; r++) out[r] = b.res_kind[r] <= 1 ? SS_COLOR[b.res_ss[r]] : CLASS_COLOR[S.stats.kinds[b.res_kind[r]]] || PALETTE.label;
    else if (scheme === 'class') for (let r = 0; r < R; r++) out[r] = CLASS_COLOR[S.stats.classes[b.res_class[r]]] || PALETTE.label;
    else if (scheme === 'hydropathy') {
      const sc = d3.scaleDiverging([4.5, 0, -4.5], RAMP.div);
      for (let r = 0; r < R; r++) { const v = KD[this.resName(r)]; out[r] = v === undefined ? '#6b6b85' : sc(v); }
    } else if (scheme === 'conf') {
      if (isAF) for (let r = 0; r < R; r++) out[r] = PLDDT_COLOR(b.res_b[r]);
      else { const sc = d3.scaleSequential(robust(b.res_b), RAMP.heat); for (let r = 0; r < R; r++) out[r] = sc(b.res_b[r]); }
    } else if (scheme === 'rainbow') {
      const first = new Map(), last = new Map();
      b.res_chain.forEach((c, r) => { if (b.res_kind[r] > 1) return; if (!first.has(c)) first.set(c, r); last.set(c, r); });
      for (let r = 0; r < R; r++) {
        const c = b.res_chain[r], a = first.get(c), z = last.get(c);
        out[r] = a === undefined || b.res_kind[r] > 1 ? CLASS_COLOR[S.stats.kinds[b.res_kind[r]]] || PALETTE.label
          : RAMP.rainbow(z > a ? (r - a) / (z - a) : 0.5);
      }
    } else if (this.F && this.F.colIndex[scheme]) {
      const vals = this.F.buf[`c_${scheme}`];
      const interp = scheme === 'rsa' ? RAMP.access : scheme === 'hull_depth' ? RAMP.seq : RAMP.heat;
      const sc = d3.scaleSequential(robust(vals), interp);
      for (let r = 0; r < R; r++) out[r] = Number.isFinite(vals[r]) ? sc(vals[r]) : '#3a3a52';
    } else for (let r = 0; r < R; r++) out[r] = chainColor(b.res_chain[r], nCh);
    return out;
  }

  legendFor(scheme) {
    const S = this.S;
    if (!S) return null;
    if (scheme === 'chain') {
      const n = S.stats.chains.length;
      return { type: 'cat', title: `chains (${n})`, items: S.stats.chains.slice(0, 16).map((c, i) => [c, chainColor(i, n)]), more: Math.max(0, n - 16) };
    }
    if (scheme === 'ss') return { type: 'cat', title: 'secondary structure', items: [['helix', SS_COLOR[1]], ['strand', SS_COLOR[2]], ['coil', SS_COLOR[0]]] };
    if (scheme === 'class') return { type: 'cat', title: 'residue class', items: Object.entries(CLASS_COLOR).filter(([k]) => S.stats.classes.includes(k)) };
    if (scheme === 'conf' && S.stats.is_alphafold) return { type: 'cat', title: 'pLDDT', items: [['> 90', '#0053d6'], ['70–90', '#65cbf3'], ['50–70', '#c084fc'], ['< 50', '#ff7d45']] };
    if (scheme === 'conf') return { type: 'seq', title: 'B-factor (Å²)', scale: d3.scaleSequential(d3.extent(S.buf.res_b), RAMP.heat) };
    if (scheme === 'hydropathy') return { type: 'seq', title: 'Kyte–Doolittle', scale: d3.scaleDiverging([4.5, 0, -4.5], RAMP.div) };
    if (scheme === 'rainbow') return { type: 'seq', title: 'N → C', scale: d3.scaleSequential([0, 1], RAMP.rainbow), noTicks: true };
    if (scheme === 'element') return { type: 'cat', title: 'element', items: ['C', 'N', 'O', 'S', 'P'].map((e) => [e, ELEMENT_COLOR[e]]) };
    if (this.F && this.F.colIndex[scheme]) {
      const col = this.F.colIndex[scheme], v = Array.from(this.F.buf[`c_${scheme}`]).filter(Number.isFinite).sort(d3.ascending);
      const ext = v.length ? [d3.quantileSorted(v, 0.02), d3.quantileSorted(v, 0.98)] : [0, 1];
      const interp = scheme === 'rsa' ? RAMP.access : scheme === 'hull_depth' ? RAMP.seq : RAMP.heat;
      return { type: 'seq', title: `${col.label}${col.unit ? ` (${col.unit})` : ''}`, scale: d3.scaleSequential(ext, interp) };
    }
    return null;
  }
}

export const store = new Store();
window.bgStore = store;

// ------------------------------------------------------------------ UI helpers
let busyCount = 0;
export function busy(on, msg) {
  busyCount = Math.max(0, busyCount + (on ? 1 : -1));
  document.getElementById('busy').hidden = busyCount === 0;
  if (msg) document.getElementById('busy-msg').textContent = msg;
  const dot = document.getElementById('st-dot'), st = document.getElementById('st-state');
  dot.classList.toggle('busy', busyCount > 0);
  st.textContent = busyCount > 0 ? 'WORKING' : 'READY';
}
// fetch a file and save it through a temporary link: no new tab, no flicker
export async function downloadFile(url, fallbackName = 'download') {
  try {
    const r = await fetch(url);
    if (!r.ok) throw new Error(`HTTP ${r.status}`);
    const m = (r.headers.get('content-disposition') || '').match(/filename="?([^";]+)"?/);
    saveBlob(await r.blob(), m ? m[1] : fallbackName);
  } catch (e) { toast(`download failed: ${e.message}`); }
}
export function saveBlob(blob, name) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob); a.download = name;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 4000);
}

let toastTimer;
export function toast(msg, kind = 'error') {
  const t = document.getElementById('toast');
  t.textContent = msg; t.className = kind; t.hidden = false;
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => (t.hidden = true), kind === 'error' ? 8000 : 3500);
  if (kind === 'error') document.getElementById('st-dot').classList.add('err');
}

export function renderLegend(el, spec) {
  const box = d3.select(el).html('');
  if (!spec) return;
  box.append('div').attr('class', 'lt').text(spec.title);
  if (spec.type === 'cat') {
    const row = box.append('div');
    spec.items.forEach(([n, c]) => { const s = row.append('span').attr('class', 'sw'); s.append('i').style('background', c); s.append('span').text(n); });
    if (spec.more) row.append('span').attr('class', 'sw').text(`+${spec.more}`);
  } else {
    const w = 150, h = 8, dom = spec.scale.domain(), lo = dom[0], hi = dom[dom.length - 1];
    const svg = box.append('svg').attr('width', w).attr('height', spec.noTicks ? h : h + 12);
    const id = `g${Math.random().toString(36).slice(2, 8)}`;
    const grad = svg.append('defs').append('linearGradient').attr('id', id);
    d3.range(0, 1.0001, 0.1).forEach((t) => grad.append('stop').attr('offset', t).attr('stop-color', spec.scale(lo + t * (hi - lo))));
    svg.append('rect').attr('width', w).attr('height', h).attr('rx', 2).attr('fill', `url(#${id})`);
    if (!spec.noTicks) svg.append('g').attr('transform', `translate(0,${h})`)
      .call(d3.axisBottom(d3.scaleLinear([lo, hi], [0, w])).ticks(4).tickSize(2).tickFormat(d3.format('~g')))
      .call((g) => g.select('.domain').remove());
  }
}
