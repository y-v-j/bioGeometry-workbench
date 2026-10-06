// Bottom dock: the Residues / Tetrahedra tables and the analysis tabs, all coupled to the selection.
import { store, KIND, PALETTE, GROUP_COLORS, CAT10, chainColor, toast, busy, fmt, downloadFile, saveBlob } from './store.js';
import { DataTable, GROUP_COLOR } from './table.js';
import { scrollBox } from './scrollbars.js';

const d3 = window.d3;
const TABS = [
  ['residues', 'Residues'], ['tetrahedra', 'Tetrahedra'], ['selection', '⌬ Selection'], ['pockets', 'Pockets'], ['contacts', 'Interfaces'],
  ['profile', 'Profile'], ['filtration', 'α filtration'], ['fourbody', 'Four-body'], ['scaling', 'Scaling'], ['layer', 'Layer stats'],
];
const TABLE_TABS = new Set(['residues', 'tetrahedra']);
// bioGeometry_<id>-<chains>_<atom set>_<content>.csv, like the server's CSV downloads
export function exportName(content, ext = 'csv') {
  const st = store.S.stats, ch = st.chains;
  const tag = ch.length > 26 ? `${ch.length}ch` : (ch.every((c) => c.length === 1) ? ch.join('') : ch.join('+'));
  const set = st.level === 'ca' || st.level === 'centroid' ? 'c-alpha' : 'all-atom';
  return `bioGeometry_${st.name}-${tag}_${set}_${content}.${ext}`;
}

export class Dock {
  constructor(geoView) {
    this.geo = geoView;
    this.body = document.getElementById('dock-body');
    this.tabsEl = document.getElementById('dock-tabs');
    this.info = document.getElementById('table-info');
    this.active = 'residues';
    this.data = {};          // analysis kind -> layer
    this.tables = {};        // tab -> DataTable
    this.tableState = {};    // tab -> saved filter/sort state across rebuilds
    d3.select(this.tabsEl).selectAll('button').data(TABS).join('button').attr('class', 'tab')
      .on('click', (ev, [id]) => this.activate(id));
    this.markTabs();
    document.querySelector('[data-tog="selonly"]').addEventListener('click', (e) => {
      e.currentTarget.classList.toggle('on');
      const t = this.tables[this.active]; if (t) t.setSelOnly(e.currentTarget.classList.contains('on'));
    });
    document.querySelector('[data-act="selfiltered"]').addEventListener('click', () => this.selectFiltered());
    document.querySelector('[data-act="clearfilters"]').addEventListener('click', () => { const t = this.tables[this.active]; if (t) t.clearFilters(); });
    document.querySelector('[data-act="exportview"]').addEventListener('click', () => this.exportView());
    document.querySelector('[data-act="columns"]').addEventListener('click', (e) => this.columnMenu(e));
    document.getElementById('csv-download').addEventListener('change', (e) => {
      const kind = e.target.value; e.target.value = '';
      if (!kind || !store.job) return;
      const path = kind.startsWith('structure:') ? `structure/${kind.slice(10)}` : `csv/${kind}`;
      busy(true, 'Preparing the download…');
      downloadFile(`${store.api}/bg/jobs/${store.job}/${path}`, `${kind}.csv`).finally(() => busy(false));
    });
    new ResizeObserver(() => { if (!TABLE_TABS.has(this.active)) this.redraw(); }).observe(this.body);
    store.on('data', () => this.reset());
    store.on('features', () => { this.buildResidueTable(); if (this.active === 'profile') this.redraw(); });
    store.on('selection', ({ source }) => this.onSelection(source));
    store.on('selstyle', () => { if (this.active === 'selection') this.redraw(); });
    store.on('params', () => {        // the tetrahedra table follows the Delaunay max edge
      const me = store.params.delaunay && +store.params.delaunay.max_edge;
      if (store.T && me && Math.abs(store.T.stats.max_edge - me) > 1e-6) { store.T = null; if (this.active === 'tetrahedra') this.activate('tetrahedra'); }
    });
    store.on('hover', () => { const t = this.tables[this.active]; if (t) t.markRows(); });
  }

  reset() {
    this.data = {}; this.tables = {}; this.tableState = {};
    store.selGeo = null; this.geo.selPair = null; store.emit('selgeo');
    this.buildResidueTable();
    this.activate('residues');
    this.fillCsvMenu();
  }

  async fillCsvMenu() {
    try {
      const kinds = await store.apiCall('/bg/csv-kinds');
      const sel = d3.select('#csv-download');
      sel.selectAll('option:not(:first-child), optgroup').remove();
      const st = sel.append('optgroup').attr('label', 'structure file (loaded model)');
      st.append('option').attr('value', 'structure:cif').text('structure — mmCIF (.cif)');
      st.append('option').attr('value', 'structure:pdb').text('structure — PDB (.pdb)');
      const cg = sel.append('optgroup').attr('label', 'feature tables (CSV)');
      Object.entries(kinds).forEach(([k, desc]) => cg.append('option').attr('value', k).text(`${k.replace(/_/g, ' ')} — ${desc}`));
    } catch (e) { /* offline API */ }
  }

  markTabs() {
    if (store.activeTab !== this.active) { store.activeTab = this.active; store.emit('tab'); }   // overlays follow the tab
    d3.select(this.tabsEl).selectAll('button').classed('on', ([id]) => id === this.active)
      .html(([id, label]) => `${label}${this.data[id] ? '<span class="dot"></span>' : ''}`);
    document.getElementById('table-tools').style.visibility = TABLE_TABS.has(this.active) ? 'visible' : 'hidden';
  }

  async activate(id) {
    this.active = id;
    this.markTabs();
    if (id === 'tetrahedra' && !store.T && store.job) {
      busy(true, 'Computing residue tetrahedra…');
      const job = store.job;
      const maxEdge = store.params.delaunay && store.params.delaunay.max_edge;
      try { const T = await store.requestLayer('tetrahedra', maxEdge ? { max_edge: +maxEdge } : {}); if (store.job === job) store.T = T; } catch (e) { toast(`tetrahedra: ${e.message}`); } finally { busy(false); }
    }
    if (id === 'profile' && !store.F) { this.message('Waiting for the residue feature table…'); return; }
    if (id === 'selection') { this.measureSelection(false); return; }
    this.redraw();
  }

  showAnalysis(L) {
    const tab = L.kind === 'contacts' ? 'contacts' : L.kind;
    this.data[tab] = L;
    if (L.kind === 'pockets' || L.kind === 'contacts') this.geo.setOverlay(L);
    this.activate(tab);
  }

  message(t) { d3.select(this.body).html('').append('div').attr('class', 'vt-empty').text(t); }

  redraw() {
    if (!store.S) return;
    if (this.active === 'residues') return this.mountTable('residues');
    if (this.active === 'tetrahedra') return this.mountTable('tetrahedra');
    if (this.active === 'selection') return this.selection();
    const L = this.data[this.active];
    if (this.active === 'profile') return this.profile();
    if (this.active === 'layer') return this.layerStats();
    if (!L) return this.message('Run this analysis from the panel on the left (Analyses).');
    ({ pockets: () => this.pockets(L), contacts: () => this.contacts(L), filtration: () => this.filtration(L),
      fourbody: () => this.fourbody(L), scaling: () => this.scaling(L) })[this.active]();
  }

  // ---------------------------------------------------------------- tables
  buildResidueTable() {
    if (!store.S) return;
    const S = store.S, b = S.buf, R = store.R;
    let columns, cols;
    if (store.F) {
      columns = store.F.stats.columns;
      cols = Object.fromEntries(columns.map((c) => [c.key, store.F.buf[`c_${c.key}`]]));
    } else {   // provisional table from the structure layer while features compute
      const cat = (key, label, group, cats, codes, desc) => ({ c: { key, label, group, type: 'cat', categories: cats, desc }, v: codes });
      const num = (key, label, group, arr, unit = '', f = '.1f', desc = '') => ({ c: { key, label, group, type: 'num', unit, fmt: f, desc }, v: arr });
      const ssCats = ['C', 'H', 'E', ''], ssCode = Uint16Array.from(b.res_ss, (s) => (s < 0 ? 3 : s));
      const items = [
        cat('chain', 'chain', 'identity', S.stats.chains, Uint16Array.from(b.res_chain)),
        num('resi', 'res #', 'identity', Float32Array.from(b.res_seq), '', 'd'),
        cat('resn', 'residue', 'identity', S.stats.res_names, Uint16Array.from(b.res_name)),
        cat('kind', 'kind', 'identity', S.stats.kinds, Uint16Array.from(b.res_kind)),
        cat('class', 'class', 'biochemistry', S.stats.classes, Uint16Array.from(b.res_class)),
        cat('ss', 'SS', 'backbone', ssCats, ssCode),
        num('phi', 'φ', 'backbone', b.res_phi, '°', '.0f'), num('psi', 'ψ', 'backbone', b.res_psi, '°', '.0f'),
        num('bfactor', S.stats.is_alphafold ? 'pLDDT' : 'B-factor', 'confidence', b.res_b),
      ];
      columns = items.map((x) => x.c);
      cols = Object.fromEntries(items.map((x) => [x.c.key, x.v]));
    }
    this.tableSpec('residues', {
      columns, cols, n: R,
      rowResidues: (i) => [i],
      isSelected: (i) => store.sel.has(i),
      isHovered: (i) => store.hover === i,
      involvesSelection: (i) => store.sel.has(i),
      onRowClick: (i, ev) => { store.click(i, ev, 'table'); store.setFocusChain(b.res_chain[i], 'table'); },
    });
  }

  tableSpec(tab, spec) {
    this.specs = this.specs || {};
    this.specs[tab] = spec;
    if (this.tables[tab]) this.tableState[tab] = this.tables[tab].state();
    delete this.tables[tab];
    if (this.active === tab) this.mountTable(tab);
  }

  mountTable(tab) {
    if (tab === 'tetrahedra' && store.T && !(this.specs && this.specs.tetrahedra && this.specs.tetrahedra._T === store.T)) {
      // the reference-format table: server columns + 'simplex-indices' and 'residues' built here from the vertices
      const T = store.T, r4 = T.buf.res4, tets = T.buf.tets, n = T.stats.n, b = store.S.buf;
      const cols = Object.fromEntries(T.stats.columns.map((c) => [c.key, T.buf[`c_${c.key}`]]));
      const textColumn = (labelOf) => {
        const labels = [], lookup = new Map(), codes = new Uint32Array(n);
        for (let i = 0; i < n; i++) {
          const s = labelOf(i);
          let c = lookup.get(s); if (c === undefined) { c = labels.length; labels.push(s); lookup.set(s, c); } codes[i] = c;
        }
        return [labels, codes];
      };
      const resLabel = (r) => `${b.res_kind[r] === KIND.water ? 'HOH' : store.resName(r)}${b.res_seq[r]}${store.S.stats.icodes[b.res_icode[r]] || ''}`;
      const [idxCats, idxCodes] = textColumn((i) => `${tets[4 * i]}-${tets[4 * i + 1]}-${tets[4 * i + 2]}-${tets[4 * i + 3]}`);
      const [resCats, resCodes] = textColumn((i) => [0, 1, 2, 3].map((k) => resLabel(r4[4 * i + k])).join('-'));
      cols['simplex-indices'] = idxCodes; cols.residues = resCodes;
      const extra = [
        { key: 'simplex-indices', label: 'simplex-indices', group: 'identity', type: 'cat', categories: idxCats, hidden: true, desc: 'indices of the four vertices in the point set' },
        { key: 'residues', label: 'residues', group: 'identity', type: 'cat', categories: resCats, desc: 'the residues of the four vertices, with their numbers' },
      ];
      const at = T.stats.columns.findIndex((c) => c.key === 'residue-class');
      const columns = [...T.stats.columns.slice(0, at), ...extra, ...T.stats.columns.slice(at)];
      const involves = (i) => store.sel.has(r4[4 * i]) || store.sel.has(r4[4 * i + 1]) || store.sel.has(r4[4 * i + 2]) || store.sel.has(r4[4 * i + 3]);
      this.tableSpec('tetrahedra', {
        _T: T, columns, cols, n, selectionFirst: true,
        rowResidues: (i) => [r4[4 * i], r4[4 * i + 1], r4[4 * i + 2], r4[4 * i + 3]],
        isSelected: involves, involvesSelection: involves,
        isHovered: (i) => [0, 1, 2, 3].some((k) => r4[4 * i + k] === store.hover),
        onRowClick: (i, ev) => {
          const rs = [r4[4 * i], r4[4 * i + 1], r4[4 * i + 2], r4[4 * i + 3]];
          store.addOrSet(rs, ev, 'table');
        },
      });
      return;
    }
    const spec = this.specs && this.specs[tab];
    if (!spec) { this.message(tab === 'tetrahedra' ? 'Computing tetrahedra…' : 'No data'); return; }
    if (!this.tables[tab] || !this.body.contains(this.tables[tab].root)) {
      const st = this.tableState[tab] || {};
      st.selOnly = document.querySelector('[data-tog="selonly"]').classList.contains('on');
      this.tables[tab] = new DataTable(this.body, spec, st);
      this.tables[tab].onChange = (t) => this.updateInfo(t);
      if (store.sel.size && !spec.selectionFirst) this.tables[tab].scrollToSelection();
    }
    this.updateInfo(this.tables[tab]);
  }

  updateInfo(t) {
    if (!t || this.tables[this.active] !== t) return;
    const nf = Object.keys(t.filters).length;
    const T = this.active === 'tetrahedra' && store.T ? store.T.stats : null;
    const shown = t.order.length === t.spec.n ? fmt.int(t.spec.n) : `${fmt.int(t.order.length)}/${fmt.int(t.spec.n)}`;
    this.info.textContent = `${shown} ${T ? 'tets' : 'rows'}${nf ? ` · ${nf} filter${nf > 1 ? 's' : ''}` : ''}${t.nInvolved ? ` · ${fmt.int(t.nInvolved)} sel ↑` : ''}` +
      (this.active === 'residues' && !store.F ? ' · features…' : '');
    this.info.title = `${fmt.int(t.order.length)} of ${fmt.int(t.spec.n)} ${T ? 'tetrahedra' : 'rows'} shown` +
      (nf ? `, ${nf} column filter${nf > 1 ? 's' : ''}` : '') +
      (T ? `\nDelaunay tetrahedra with every edge ≤ ${T.max_edge} Å${T.sampled ? ` (sampled from ${fmt.int(T.n_compact)})` : ''}` : '') +
      (t.nInvolved ? `\n${fmt.int(t.nInvolved)} involve the selection and are listed first` : '');
  }

  onSelection(source) {
    if (!store.sel.size && (this.geo.pocketSel !== null || this.geo.contactPair)) {
      this.geo.pocketSel = null; this.geo.contactPair = null;
      if (this.geo.overlay) this.geo.setOverlay(this.geo.overlay);
      if (this.active === 'pockets' && this.data.pockets) this.pockets(this.data.pockets);
    }
    if (this.active === 'selection') {             // the measurement follows the selection while its tab is open
      clearTimeout(this.selTimer);
      this.selTimer = setTimeout(() => this.measureSelection(false), 450);
      return;
    }
    if (store.selGeo && store.selGeo._key !== [...store.sel].sort((a, b) => a - b).join(',')) { store.selGeo = null; store.emit('selgeo'); }
    const t = this.tables[this.active];
    if (!t) { if (this.active === 'profile' || this.active === 'pockets') this.markProfile(); return; }
    if (t.selOnly || t.spec.selectionFirst) { t.apply(); if (t.spec.selectionFirst) t.body.node().scrollTop = 0; } else t.markRows();
    if (source !== 'table' && store.sel.size && !t.spec.selectionFirst) t.scrollToSelection();
  }

  selectFiltered() {
    const t = this.tables[this.active];
    if (!t) return;
    const rs = new Set();
    for (const i of t.order) { for (const r of t.spec.rowResidues(i)) rs.add(r); if (rs.size > 50000) break; }
    store.setSelection([...rs], 'table');
  }

  exportView() {
    const t = this.tables[this.active];
    if (!t) return;
    saveBlob(new Blob([t.toCSV()], { type: 'text/csv' }), exportName(`${this.active}-view`));
  }

  columnMenu(ev) {
    const t = this.tables[this.active];
    if (!t) return;
    const menu = d3.select('#menu').html('').attr('hidden', null);
    const groups = d3.groups(t.spec.columns, (c) => c.group);
    for (const [g, cols] of groups) {
      menu.append('div').attr('class', 'mg').style('color', GROUP_COLOR[g] || PALETTE.label).text(g);
      cols.forEach((c) => {
        const l = menu.append('label').attr('title', c.desc || '');
        l.append('input').attr('type', 'checkbox').property('checked', !t.hidden.has(c.key)).on('change', (e) => t.setHidden(c.key, !e.target.checked));
        l.append('span').text(`${c.label}${c.unit ? ` (${c.unit})` : ''}`);
      });
    }
    const r = ev.currentTarget.getBoundingClientRect();
    menu.style('left', `${Math.min(r.left, window.innerWidth - 240)}px`).style('top', `${Math.max(8, r.top - Math.min(window.innerHeight * 0.6, 420) - 6)}px`);
    const close = (e) => { if (!document.getElementById('menu').contains(e.target)) { menu.attr('hidden', true); window.removeEventListener('pointerdown', close, true); } };
    setTimeout(() => window.addEventListener('pointerdown', close, true), 0);
  }

  // ---------------------------------------------------------------- selection geometry
  // measure the current selection (domains, chains, ligands, ions) and show it in the 'Selection' tab
  async measureSelection(open = true) {
    if (open && this.active !== 'selection') { this.active = 'selection'; this.markTabs(); }
    if (!store.S) return;
    const residues = [...store.sel].sort((a, b) => a - b), key = residues.join(','), job = store.job;
    if (!residues.length) { store.selGeo = null; store.emit('selgeo'); this.selection(); return; }
    if (store.selGeo && store.selGeo._key === key) { this.selection(); return; }
    const token = (this.selToken = (this.selToken || 0) + 1);
    if (this.active === 'selection') this.message(`Measuring ${fmt.int(residues.length)} selected residues…`);
    try {
      const maxEdge = store.params.delaunay && +store.params.delaunay.max_edge;
      const L = await store.requestLayer('selection', maxEdge ? { residues, max_edge: maxEdge } : { residues });
      if (store.job !== job || token !== this.selToken) return;
      L._key = key;
      store.selGeo = L; this.geo.selPair = null;
      store.emit('selgeo');
    } catch (e) { if (token === this.selToken) { toast(`selection geometry: ${e.message}`); store.selGeo = null; store.emit('selgeo'); } }
    if (token === this.selToken && this.active === 'selection') this.selection();
  }

  selection() {
    const row = this.charts(), L = store.selGeo;
    if (!store.sel.size || !L) {
      const c = this.panel(row, 'GEOMETRY OF A SELECTION', '', 'wide3');
      c.append('div').attr('class', 'sel-help').html(store.sel.size ? 'Measuring…' :
        `Pick what to measure: a <b>domain</b> (click its bar under the sequence, or tick it in the list on the right),
         a <b>ligand</b> or <b>ion</b> (tick it in the HETATM list; ◎ adds its binding site), whole chains or residues.
         <b>Shift</b> adds to the selection, <b>ctrl</b> toggles.<br><br>The selection is split into groups (a domain or chain
         segment, every ligand and ion instance) and measured as a whole, per group and per pair of groups: convex hull,
         α-shape, Delaunay complex, gyration shape, solvent accessibility alone and in the structure, atom contacts,
         hydrogen-bond candidates, metal coordination and the surface buried between groups.
         The measurement follows the selection while this tab is open.`);
      return;
    }
    const s = L.stats, S = s.summary, f0 = d3.format(',.0f'), f1 = d3.format(',.1f'), f2 = d3.format('.2f'), pct = d3.format('.0%');
    const v = (x, f = f1) => (x === undefined || x === null || !Number.isFinite(x) ? '–' : f(x));
    const gcol = (k) => GROUP_COLORS[k % GROUP_COLORS.length];
    // whole selection
    const k = this.panel(row, `SELECTION · ${fmt.int(S.n_residues)} RESIDUES · ${fmt.int(S.n_atoms)} ${S.level === 'heavy' ? 'ATOMS' : 'Cα/P + HETERO ATOMS'}`, '', 'scroll');
    const T = S.participating || {};
    this.kv(k, {
      groups: S.n_groups, chains: S.chains.join(' '),
      [`tetrahedra it takes part in (edges ≤ ${T.max_edge} Å)`]: T.n_tets ? `${fmt.int(T.n_tets)} (${fmt.int(T.n_inside)} wholly inside)` : '0',
      'their nodes · edges': T.n_tets ? `${fmt.int(T.n_nodes)} · ${fmt.int(T.n_edges)}` : '–',
      'their volume (Å³) · mean edge (Å)': T.n_tets ? `${v(T.volume, f0)} · ${v(T.mean_edge, f2)}` : '–',
      'neighbouring residues in them': T.n_tets ? `${T.n_neighbour_residues}${T.neighbours.length ? ` · ${T.neighbours.slice(0, 6).join(' ')}${T.n_neighbour_residues > 6 ? ' …' : ''}` : ''}` : '–',
      'hull volume (Å³)': v(S.hull_volume, f0), 'hull area (Å²)': v(S.hull_area, f0), 'hull sphericity ψ': v(S.sphericity, f2),
      [`α-shape volume (α ${S.alpha} Å)`]: v(S.alpha_volume, f0), 'α-shape area (Å²)': v(S.alpha_area, f0),
      'α components · cavities': S.alpha_components !== undefined ? `${S.alpha_components} · ${S.alpha_cavities}` : '–',
      'α-shape / hull volume': v(S.alpha_over_hull, f2),
      'own Delaunay (selection alone) tets · edges': S.n_tets !== undefined ? `${fmt.int(S.n_tets)} · ${fmt.int(S.n_edges)}` : '–',
      'own mean Delaunay edge (Å)': v(S.mean_edge, f2),
      'radius of gyration (Å)': v(S.shape && S.shape.rg), 'shape κ² (0 sphere → 1 rod)': v(S.shape && S.shape.kappa2, f2),
      'principal half-lengths (Å)': S.shape ? S.shape.principal_lengths.map((x) => x.toFixed(1)).join(' · ') : '–',
      'SASA alone (Å²)': v(S.sasa_selection, f0), 'SASA in the structure (Å²)': v(S.sasa_in_structure, f0),
      'buried by the rest (Å²)': `${v(S.sasa_buried_by_rest, f0)}${S.sasa_selection ? ` (${pct(S.sasa_buried_by_rest / S.sasa_selection)})` : ''}`,
      ...(S.note ? { note: S.note } : {}),
    });
    // groups
    const gt = this.panel(row, 'GROUPS · row: fly to it · ⊕ select it (shift adds)', '', 'wide scroll');
    const tb = gt.append('table').attr('class', 'mini nowrap');
    tb.append('thead').append('tr').selectAll('th')
      .data(['', 'group', 'kind', 'res', 'atoms', 'Rg Å', 'hull Å³', 'ψ', 'SASA alone', 'in selection', 'in structure', 'buried · sel', 'buried · struct'])
      .join('th').text((d) => d);
    const gr = tb.append('tbody').selectAll('tr').data(s.groups).join('tr').attr('title', (d) => `${d.label} · chains ${d.chains.join(', ')}`)
      .on('click', (ev, d) => this.geo.focusResidues(d.residues));
    gr.append('td').append('span').attr('class', 'pick').attr('title', 'select this group (shift: add, ctrl: toggle)').text('⊕')
      .on('click', (ev, d) => { ev.stopPropagation(); store.pick(d.residues, ev, 'dock'); });
    gr.append('td').html((d, i) => `<i class="sw" style="background:${gcol(i)}"></i>${d.label}`);
    gr.append('td').text((d) => d.kind);
    gr.selectAll(null).data((d) => [fmt.int(d.n_residues), fmt.int(d.n_atoms), v(d.rg), v(d.hull_volume, f0), v(d.sphericity, f2),
      v(d.sasa_alone, f0), v(d.sasa_in_selection, f0), v(d.sasa_in_structure, f0),
      d.buried_in_selection === undefined ? '–' : pct(d.buried_in_selection), d.buried_in_structure === undefined ? '–' : pct(d.buried_in_structure)])
      .join('td').attr('class', 'num').text((d) => d);
    // pairs
    const note = s.groups.length < 2 ? 'one group: add a ligand, an ion or another domain to measure what lies between them'
      : 'row: show only that pair\'s contacts in the geometry view · ⊕ select the contacting residues';
    const pt = this.panel(row, `PAIRS · contacts ≤ 4 Å · polar ≤ 3.5 Å · metal ≤ 2.8 Å${s.pairs_truncated ? ' · first 16 groups' : ''}`, note, 'wide scroll');
    if (s.pairs.length) {
      const tp = pt.append('table').attr('class', 'mini nowrap');
      tp.append('thead').append('tr').selectAll('th')
        .data(['', 'A', 'B', 'min Å', 'closest atoms', 'centroids Å', 'contacts', 'polar', 'Delaunay', 'BSA Å²', 'interface Å²', 'A buried', 'B buried', 'metal coordination'])
        .join('th').text((d) => d);
      const pr = tp.append('tbody').selectAll('tr').data(s.pairs).join('tr')
        .classed('sel', (d) => this.geo.selPair && this.geo.selPair[0] === d.a && this.geo.selPair[1] === d.b)
        .on('click', (ev, d) => {
          const same = this.geo.selPair && this.geo.selPair[0] === d.a && this.geo.selPair[1] === d.b;
          this.geo.selPair = same ? null : [d.a, d.b];
          this.geo.buildSelGeo();
          if (!same && d.residues_a) this.geo.focusResidues([...d.residues_a, ...d.residues_b]);
          pr.classed('sel', (q) => this.geo.selPair && this.geo.selPair[0] === q.a && this.geo.selPair[1] === q.b);
        });
      pr.append('td').append('span').attr('class', 'pick').attr('title', 'select the contacting residues (shift: add)').text('⊕')
        .on('click', (ev, d) => { ev.stopPropagation(); if (d.residues_a) store.pick([...d.residues_a, ...d.residues_b], ev, 'dock'); });
      pr.selectAll(null).data((d) => [
        `<i class="sw" style="background:${gcol(d.a)}"></i>${s.groups[d.a].label}`, `<i class="sw" style="background:${gcol(d.b)}"></i>${s.groups[d.b].label}`,
        v(d.min_dist, f2), d.min_pair.join(' – '), v(d.centroid_dist), d.n_contacts ?? 0, d.n_polar ?? 0, d.delaunay_contacts,
        v(d.bsa, f0), v(d.interface_area, f0), d.buried_a === undefined || d.buried_a === null ? '–' : pct(d.buried_a),
        d.buried_b === undefined || d.buried_b === null ? '–' : pct(d.buried_b),
        d.metal ? d.metal.map((m) => `${m.atom} ${m.dist.toFixed(2)}`).join(', ') : ''])
        .join('td').attr('class', (d, i) => (i >= 2 && i <= 11 && i !== 3 ? 'num' : null)).html((d) => d)
        .attr('title', (d, i) => (i === 3 || i === 12 ? String(d) : null));
    }
    // buried surface per group
    if (s.groups.some((d) => d.buried_in_selection !== undefined)) {
      const c = this.panel(row, 'SURFACE BURIED PER GROUP');
      c.append('div').attr('class', 'note').html(`<i class="key" style="background:${PALETTE.sel}"></i>by the rest of the selection &nbsp;<i class="key" style="background:${PALETTE.lavender}"></i>in the whole structure`);
      const [w, h] = this.size(c);
      const rows = s.groups.map((d, i) => ({ d, i })).filter((r) => r.d.buried_in_selection !== undefined).slice(0, Math.max(2, Math.floor((h - 26) / 22)));
      const { g, iw, ih } = this.frame(c, w, h, { t: 2, r: 12, b: 22, l: 96 });
      const y = d3.scaleBand(rows.map((r) => r.i), [0, ih]).padding(0.25), x = d3.scaleLinear([0, 1], [0, iw]);
      g.append('g').attr('class', 'grid').attr('transform', `translate(0,${ih})`).call(d3.axisBottom(x).ticks(4).tickSize(-ih).tickFormat(''));
      g.append('g').attr('class', 'axis').attr('transform', `translate(0,${ih})`).call(d3.axisBottom(x).ticks(4).tickFormat(d3.format('.0%')));
      g.append('g').attr('class', 'axis').call(d3.axisLeft(y).tickSize(0).tickFormat((i) => s.groups[i].label.slice(0, 15))).call((q) => q.select('.domain').remove());
      const bh = y.bandwidth() / 2;
      for (const [key, col, off] of [['buried_in_selection', PALETTE.sel, 0], ['buried_in_structure', PALETTE.lavender, bh]]) {
        g.selectAll(null).data(rows).join('rect').attr('y', (r) => y(r.i) + off).attr('height', bh - 1).attr('x', 0)
          .attr('width', (r) => x(Math.max(0, r.d[key]))).attr('fill', key === 'buried_in_selection' ? col : col).attr('rx', 2)
          .append('title').text((r) => `${r.d.label}: ${pct(r.d[key])} of its surface buried ${key === 'buried_in_selection' ? 'by the rest of the selection' : 'in the whole structure'}`);
      }
    }
  }

  // ---------------------------------------------------------------- chart helpers
  charts() { return d3.select(this.body).html('').append('div').attr('class', 'charts'); }
  // a chart panel; 'scroll' panels return a gutter-scrolled box to fill instead of the panel itself
  panel(row, title, note, cls = '') {
    const scroll = /\bscroll\b/.test(cls);
    const c = row.append('div').attr('class', `chart ${cls.replace(/\bscroll\b/, '')}`);
    c.append('h4').text(title);
    if (note) c.append('div').attr('class', 'note').text(note);
    return scroll ? scrollBox(c, { y: true, x: /\bwide\b/.test(cls) }) : c;
  }
  size(c) { const n = c.node(), used = Array.from(n.children).reduce((a, ch) => a + ch.offsetHeight + 2, 0); return [Math.max(200, n.clientWidth), Math.max(80, n.clientHeight - used - 4)]; }
  frame(c, w, h, m = { t: 8, r: 10, b: 24, l: 40 }) {
    const svg = c.append('svg').attr('width', w).attr('height', h);
    return { svg, g: svg.append('g').attr('transform', `translate(${m.l},${m.t})`), iw: w - m.l - m.r, ih: h - m.t - m.b };
  }
  axes(g, x, y, iw, ih, xl, yl) {
    g.append('g').attr('class', 'grid').call(d3.axisLeft(y).ticks(4).tickSize(-iw).tickFormat(''));
    g.append('g').attr('class', 'axis').attr('transform', `translate(0,${ih})`).call(d3.axisBottom(x).ticks(Math.max(2, iw / 70)).tickFormat(d3.format('~s')));
    g.append('g').attr('class', 'axis').call(d3.axisLeft(y).ticks(4).tickFormat(d3.format('~s')));
    if (xl) g.append('text').attr('x', iw).attr('y', ih + 20).attr('text-anchor', 'end').text(xl);
    if (yl) g.append('text').attr('x', 2).attr('y', 6).text(yl);
  }
  hbars(c, rows, { label, value, color = PALETTE.sky, w, h, xl, onClick, diverging, fmtv = d3.format(',.3~f') }) {
    rows = rows.slice(0, Math.max(3, Math.floor((h - 26) / 13)));
    const { g, iw, ih } = this.frame(c, w, h, { t: 2, r: 30, b: 22, l: 104 });
    const y = d3.scaleBand(rows.map(label), [0, ih]).padding(0.2);
    const ext = d3.extent(rows, value);
    const x = d3.scaleLinear(diverging ? [Math.min(0, ext[0]), Math.max(0, ext[1])] : [0, ext[1] || 1], [0, iw]).nice();
    g.append('g').attr('class', 'grid').attr('transform', `translate(0,${ih})`).call(d3.axisBottom(x).ticks(4).tickSize(-ih).tickFormat(''));
    g.append('g').attr('class', 'axis').attr('transform', `translate(0,${ih})`).call(d3.axisBottom(x).ticks(4).tickFormat(d3.format('~s')));
    g.append('g').attr('class', 'axis').call(d3.axisLeft(y).tickSize(0)).call((s) => s.select('.domain').remove());
    if (xl) g.append('text').attr('x', iw).attr('y', ih + 19).attr('text-anchor', 'end').text(xl);
    g.selectAll('rect.bar').data(rows).join('rect').attr('class', 'bar').attr('y', (d) => y(label(d))).attr('height', y.bandwidth())
      .attr('x', (d) => x(Math.min(0, value(d)))).attr('width', (d) => Math.abs(x(value(d)) - x(0)))
      .attr('fill', (d, i) => (typeof color === 'function' ? color(d, i) : color)).attr('rx', 2).style('cursor', onClick ? 'pointer' : null)
      .on('click', (ev, d) => onClick && onClick(d, ev)).append('title').text((d) => `${label(d)}: ${fmtv(value(d))}`);
  }
  kv(c, obj) {
    const el = c.append('div').attr('class', 'kv');
    for (const [k, v] of Object.entries(obj)) { const p = el.append('div').attr('class', 'kvp'); p.append('span').text(k); p.append('span').text(v); }
  }

  // ---------------------------------------------------------------- analysis tabs
  pockets(L) {
    const s = L.stats, row = this.charts();
    const c = this.panel(row, `POCKETS · small α ${s.alpha_small} Å · large α ${s.alpha_large} Å`, 'click a pocket: select its lining residues and fly to it · click it again (or Esc) to deselect');
    const tools = c.append('div').attr('class', 'chart-tools');
    tools.append('button').attr('class', `chip${this.geo.showOverlay ? ' on' : ''}`).text(this.geo.showOverlay ? 'pockets shown' : 'pockets hidden')
      .attr('title', 'show / hide the pocket surfaces in the geometry view')
      .on('click', () => { this.geo.setOverlayVisible(!this.geo.showOverlay); this.pockets(L); });
    c.append('div').attr('class', 'note').text(`Not every structure has pockets, cavities or channels: small, flat or fibrous chains often have none.` +
      (s.polymer_only ? ' Pockets are found in the polymer alone (ligands, ions and waters set aside), so an occupied binding site still shows up.' : ''));
    if (!s.pockets.length) {
      c.append('div').attr('class', 'note').style('color', PALETTE.teal)
        .text('None above the minimum volume — that is a result, not an error. A larger "pocket α large" or a smaller "min pocket Å³" finds shallower ones.');
      return;
    }
    const [w, h] = this.size(c), pal = CAT10;
    const pick = (d, ev) => {
      const again = this.geo.pocketSel === d.id;
      this.geo.pocketSel = again ? null : d.id;
      if (!this.geo.showOverlay) this.geo.setOverlayVisible(true);
      this.geo.setOverlay(L);
      if (again) {                      // clicking the selected pocket again deselects it, residues included
        if (store.additive(ev)) store.remove(d.res_index, 'dock'); else store.clear('dock');
      } else {
        store.addOrSet(d.res_index, ev, 'dock');
        this.geo.focusResidues(d.res_index);
      }
      this.pockets(L);
    };
    this.hbars(c, s.pockets, { label: (d) => `#${d.id + 1}${d.enclosed ? ' (cavity)' : ''}`, value: (d) => d.volume, xl: 'Å³', w, h, onClick: pick,
      color: (d) => (this.geo.pocketSel === null || this.geo.pocketSel === d.id ? pal[d.id % 10] : '#3a3a52') });
    const t = this.panel(row, 'POCKET TABLE', '', 'wide scroll');
    const tb = t.append('table').attr('class', 'mini');
    tb.append('thead').append('tr').selectAll('th').data(['#', 'volume Å³', 'mouth Å²', 'residues', 'chains', 'lining']).join('th').text((d) => d);
    tb.append('tbody').selectAll('tr').data(s.pockets).join('tr').classed('sel', (d) => d.id === this.geo.pocketSel)
      .on('click', (ev, d) => pick(d, ev))
      .selectAll('td').data((d) => [d.id + 1, d3.format(',.0f')(d.volume), d.enclosed ? 'enclosed' : d3.format(',.0f')(d.mouth_area), d.n_residues, d.chains.join(','), d.residues.slice(0, 8).join(' ')])
      .join('td').text((d) => d);
  }

  contacts(L) {
    const s = L.stats, n = s.chains.length, row = this.charts();
    const c = this.panel(row, `CHAIN × CHAIN DELAUNAY CONTACTS (≤ ${s.cutoff} Å)`, 'click a cell: draw that interface and select its residues');
    if (n < 2) c.append('div').attr('class', 'note').text('single chain - no interfaces');
    else {
      const [w, h] = this.size(c), side = Math.max(110, Math.min(w - 50, h - 30)), m = { t: 2, r: 2, b: 26, l: 30 };
      const { g } = this.frame(c, side + m.l + m.r, side + m.t + m.b, m);
      const x = d3.scaleBand(d3.range(n), [0, side]), mat = s.chain_matrix, vmax = d3.max(mat.flat()) || 1;
      const color = d3.scaleSequentialLog([1, vmax], (t) => d3.interpolateRgbBasis(['#20203a', '#7dd3fc', '#c4b5fd', '#f9a8d4'])(t));
      const cells = [];
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) if (mat[i][j] > 0) cells.push({ i, j, v: mat[i][j] });
      g.selectAll('rect').data(cells).join('rect').attr('x', (d) => x(d.j)).attr('y', (d) => x(d.i))
        .attr('width', Math.max(1, x.bandwidth() - (n > 40 ? 0 : 1))).attr('height', Math.max(1, x.bandwidth() - (n > 40 ? 0 : 1)))
        .attr('fill', (d) => color(d.v)).style('cursor', 'pointer')
        .on('click', (ev, d) => this.selectInterface(L, d.i, d.j, ev)).append('title').text((d) => `${s.chains[d.i]} ↔ ${s.chains[d.j]}: ${d.v} contacts`);
      const every = Math.ceil(n / 24);
      g.append('g').selectAll('text').data(d3.range(0, n, every)).join('text').attr('x', (i) => x(i) + x.bandwidth() / 2).attr('y', side + 11).attr('text-anchor', 'middle').text((i) => s.chains[i]);
      g.append('g').selectAll('text').data(d3.range(0, n, every)).join('text').attr('x', -4).attr('y', (i) => x(i) + x.bandwidth() / 2 + 3).attr('text-anchor', 'end').text((i) => s.chains[i]);
    }
    const t = this.panel(row, `INTERFACES · ${fmt.int(s.n_interchain)} inter-chain contacts · ${fmt.int(s.n_interface_residues)} interface residues`, '', 'scroll');
    const tb = t.append('table').attr('class', 'mini');
    tb.append('thead').append('tr').selectAll('th').data(['chain A', 'chain B', 'contacts', 'residues']).join('th').text((d) => d);
    tb.append('tbody').selectAll('tr').data(s.pairs.slice(0, 80)).join('tr')
      .on('click', (ev, d) => this.selectInterface(L, s.chains.indexOf(d.a), s.chains.indexOf(d.b), ev))
      .selectAll('td').data((d) => [d.a, d.b, d.contacts, d.residues]).join('td').text((d) => d);
  }

  selectInterface(L, ci, cj, ev) {
    this.geo.contactPair = [ci, cj];
    this.geo.setOverlay(L);
    const e = L.buf.res_interchain, ch = store.S.buf.res_chain, rs = new Set();
    for (let k = 0; k < e.length / 2; k++) {
      const a = e[2 * k], b = e[2 * k + 1];
      if ((ch[a] === ci && ch[b] === cj) || (ch[a] === cj && ch[b] === ci)) { rs.add(a); rs.add(b); }
    }
    store.addOrSet([...rs], ev, 'dock');
  }

  profile() {
    const F = store.F, S = store.S, row = this.charts();
    const c = this.panel(row, 'PER-RESIDUE PROFILE', 'brush along the sequence to select residues · Esc clears', 'wide');
    const num = F.stats.columns.filter((x) => x.type === 'num' && !['resi', 'charge'].includes(x.key));
    this.metric = num.some((x) => x.key === this.metric) ? this.metric : (num.find((x) => x.key === 'rsa') || num[0]).key;
    const ctrl = c.append('div').style('position', 'absolute').style('right', '0').style('top', '0');
    const sel = ctrl.append('select').on('change', (ev) => { this.metric = ev.target.value; this.profile(); });
    sel.selectAll('option').data(num).join('option').attr('value', (d) => d.key).property('selected', (d) => d.key === this.metric).text((d) => d.label);
    const col = num.find((x) => x.key === this.metric), vals = F.buf[`c_${this.metric}`], n = vals.length;
    const [w, h] = this.size(c);
    const { g, iw, ih } = this.frame(c, w, h);
    const x = d3.scaleLinear([0, n - 1], [0, iw]);
    const fin = Array.from(vals).filter(Number.isFinite);
    const y = d3.scaleLinear([Math.min(0, d3.min(fin) ?? 0), d3.max(fin) ?? 1], [ih, 0]).nice();
    const ch = S.buf.res_chain, bands = [];
    let st = 0;
    for (let r = 1; r <= n; r++) if (r === n || ch[r] !== ch[r - 1]) { bands.push({ c: ch[st], a: st, b: r - 1 }); st = r; }
    g.selectAll('rect.band').data(bands).join('rect').attr('x', (d) => x(d.a)).attr('width', (d) => Math.max(0.5, x(d.b) - x(d.a))).attr('y', ih + 1).attr('height', 3)
      .attr('fill', (d) => chainColor(d.c, S.stats.chains.length));
    this.axes(g, x, y, iw, ih, 'residue index', `${col.label}${col.unit ? ` (${col.unit})` : ''}`);
    const color = GROUP_COLOR[col.group] || PALETTE.mint;
    g.append('path').datum(Array.from(vals)).attr('fill', 'none').attr('stroke', color).attr('stroke-width', n > 3000 ? 0.6 : 1.1)
      .attr('d', d3.line().defined((d) => Number.isFinite(d)).x((d, i) => x(i)).y((d) => y(d)));
    this.profSel = g.append('g');
    this.profScale = { x, y, ih, vals };
    const brush = d3.brushX().extent([[0, 0], [iw, ih]]).on('end', (ev) => {
      if (!ev.selection || !ev.sourceEvent) return;
      const [a, b] = ev.selection.map((p) => Math.round(x.invert(p)));
      const rs = d3.range(Math.max(0, a), Math.min(n - 1, b) + 1);
      g.select('.brush').call(brush.move, null);
      store.addOrSet(rs, ev.sourceEvent, 'dock');
    });
    g.append('g').attr('class', 'brush').call(brush);
    this.markProfile();
    const cc = this.panel(row, 'SPEARMAN ρ BETWEEN FEATURES', 'rank correlation, ties averaged');
    const rows = (F.stats.correlations || []).slice().sort((a, b) => Math.abs(b.spearman) - Math.abs(a.spearman));
    const [w2, h2] = this.size(cc);
    const lab = (k) => (F.stats.columns.find((x) => x.key === k) || { label: k }).label;
    this.hbars(cc, rows, { label: (d) => `${lab(d.a)} · ${lab(d.b)}`, value: (d) => d.spearman, diverging: true, w: Math.min(w2, 440), h: h2,
      color: (d) => (d.spearman > 0 ? PALETTE.mint : PALETTE.coral), fmtv: d3.format('+.3f'), xl: 'ρ' });
  }

  markProfile() {
    if (this.active !== 'profile' || !this.profSel) return;
    const { x, y, ih, vals } = this.profScale;
    const rs = [...store.sel].filter((r) => r < vals.length).sort((a, b) => a - b);
    const runs = [];
    for (const r of rs) { const l = runs[runs.length - 1]; if (l && r === l[1] + 1) l[1] = r; else runs.push([r, r]); }
    this.profSel.selectAll('rect').data(runs.slice(0, 2000)).join('rect').attr('x', (d) => x(d[0]) - 1).attr('width', (d) => Math.max(2, x(d[1]) - x(d[0]) + 2))
      .attr('y', 0).attr('height', ih).attr('fill', PALETTE.sky).attr('fill-opacity', 0.16).attr('pointer-events', 'none');
  }

  filtration(L) {
    const s = L.stats, row = this.charts();
    const alphaNow = (store.params.alpha || {}).alpha;
    for (const [key, title, tf, color] of [
      ['volume', 'VOLUME / HULL VOLUME', (v) => v / s.hull_volume, PALETTE.mint], ['area', 'SURFACE AREA (Å²)', (v) => v, PALETTE.teal],
      ['components', 'CONNECTED COMPONENTS (β₀)', (v) => v, PALETTE.sky], ['cavities', 'ENCLOSED CAVITIES', (v) => v, PALETTE.rose]]) {
      const c = this.panel(row, title);
      const [w, h] = this.size(c);
      const { g, iw, ih } = this.frame(c, w, h);
      const pts = s.alpha.map((a, i) => [a, tf(s[key][i])]);
      const x = d3.scaleLinear(d3.extent(s.alpha), [0, iw]), y = d3.scaleLinear([0, d3.max(pts, (p) => p[1]) || 1], [ih, 0]).nice();
      this.axes(g, x, y, iw, ih, 'α (Å)');
      g.append('path').datum(pts).attr('fill', 'none').attr('stroke', color).attr('stroke-width', 2)
        .attr('d', d3.line().x((p) => x(p[0])).y((p) => y(p[1])).curve(key === 'components' || key === 'cavities' ? d3.curveStepAfter : d3.curveMonotoneX));
      if (alphaNow) g.append('line').attr('x1', x(alphaNow)).attr('x2', x(alphaNow)).attr('y1', 0).attr('y2', ih).attr('stroke', PALETTE.coral).attr('stroke-dasharray', '4 3');
    }
  }

  fourbody(L) {
    const s = L.stats, row = this.charts();
    if (!s.n_tets) { row.append('div').attr('class', 'vt-empty').text('no compact four-residue tetrahedra'); return; }
    const head = this.panel(row, 'DELAUNAY FOUR-BODY STATISTICS');
    this.kv(head, { 'tetrahedra counted': fmt.int(s.n_tets), 'distinct quadruplets': fmt.int(s.n_distinct), 'max edge (Å)': s.max_edge, score: 'log(f_obs / p_multinomial)' });
    for (const [key, title, color] of [['most_enriched', 'MOST ENRICHED (count ≥ 5)', PALETTE.mint], ['most_depleted', 'MOST DEPLETED', PALETTE.coral], ['most_common', 'MOST FREQUENT', PALETTE.sky]]) {
      const c = this.panel(row, title), [w, h] = this.size(c);
      const val = key === 'most_common' ? (d) => d.count : (d) => d.log_likelihood;
      this.hbars(c, s[key] || [], { label: (d) => d.quad, value: val, diverging: key !== 'most_common', color, w, h, xl: key === 'most_common' ? 'count' : 'log-likelihood' });
    }
  }

  scaling(L) {
    const s = L.stats, rows = s.rows, row = this.charts();
    const algs = { delaunay: PALETTE.sky, alpha: PALETTE.mint, voronoi: PALETTE.lavender, hull: PALETTE.teal };
    const c = this.panel(row, `WALL TIME VS POINTS (${s.level}, log–log)`, Object.keys(s.exponents).length ? '' : 'size range too narrow for a meaningful exponent', 'wide');
    const [w, h] = this.size(c), { g, iw, ih } = this.frame(c, w, h);
    const x = d3.scaleLog(d3.extent(rows, (r) => r.n), [0, iw]);
    const y = d3.scaleLog([d3.min(rows, (r) => d3.min(Object.keys(algs), (a) => Math.max(r[a], 1e-4))), d3.max(rows, (r) => d3.max(Object.keys(algs), (a) => r[a]))], [ih, 0]).nice();
    g.append('g').attr('class', 'axis').attr('transform', `translate(0,${ih})`).call(d3.axisBottom(x).ticks(5, '~s'));
    g.append('g').attr('class', 'axis').call(d3.axisLeft(y).ticks(4, '~g'));
    Object.entries(algs).forEach(([a, col], k) => {
      const pts = rows.map((r) => [r.n, Math.max(r[a], 1e-4)]);
      g.append('path').datum(pts).attr('fill', 'none').attr('stroke', col).attr('stroke-width', 2).attr('d', d3.line().x((p) => x(p[0])).y((p) => y(p[1])));
      g.selectAll(null).data(pts).enter().append('circle').attr('cx', (p) => x(p[0])).attr('cy', (p) => y(p[1])).attr('r', 2.5).attr('fill', col);
      g.append('text').attr('x', 6).attr('y', 10 + k * 12).style('fill', col).text(`${a}${s.exponents[a] !== undefined ? `  ~ n^${s.exponents[a].toFixed(2)}` : ''}`);
    });
  }

  layerStats() {
    const row = this.charts(), L = this.geo.layers;
    const kinds = Object.keys(L).filter((k) => this.geo.active.has(k));
    if (!kinds.length) { row.append('div').attr('class', 'vt-empty').text('switch on a geometry layer (hull / Delaunay / α-shape / Voronoi)'); return; }
    for (const k of kinds) {
      const s = L[k].stats, c = this.panel(row, `${k.toUpperCase()} — PARAMETERS & STATISTICS`, '', 'scroll');
      const flat = Object.fromEntries(Object.entries(s).filter(([, v]) => typeof v === 'number' || typeof v === 'string').map(([a, v]) => [a.replace(/_/g, ' '), typeof v === 'number' ? d3.format(',.4~f')(v) : v]));
      this.kv(c, { ...Object.fromEntries(Object.entries(L[k].meta.params).map(([a, v]) => [`param ${a}`, String(v)])), ...flat });
      if (k === 'delaunay') {
        const h = this.panel(row, 'TETRAHEDRON CIRCUMRADII (the α at which each enters the α-shape)'), [w, hh] = this.size(h);
        const vals = Array.from(L[k].buf.hist_circumradius).filter((v) => v < 40);
        const { g, iw, ih } = this.frame(h, w, hh), x = d3.scaleLinear(d3.extent(vals), [0, iw]).nice();
        const bins = d3.bin().domain(x.domain()).thresholds(x.ticks(50))(vals), y = d3.scaleLinear([0, d3.max(bins, (b) => b.length)], [ih, 0]).nice();
        this.axes(g, x, y, iw, ih, 'Å', 'count');
        g.selectAll('rect').data(bins).join('rect').attr('x', (b) => x(b.x0) + 0.5).attr('width', (b) => Math.max(0, x(b.x1) - x(b.x0) - 1)).attr('y', (b) => y(b.length)).attr('height', (b) => ih - y(b.length)).attr('fill', PALETTE.teal);
      }
    }
  }
}
