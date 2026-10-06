// Virtual data table with per-column filters, sorting, mini histograms and selection coupling.
import { store, PALETTE } from './store.js';
import { attachScrollbars } from './scrollbars.js';

const d3 = window.d3;
const ROW_H = 20;
export const GROUP_COLOR = {
  identity: PALETTE.label, biochemistry: PALETTE.rose, backbone: PALETTE.teal, confidence: PALETTE.lavender,
  surface: PALETTE.sky, geometry: PALETTE.mint, interactions: PALETTE.coral, annotation: '#a5b4fc',
  edges: PALETTE.sky, volume: PALETTE.mint, faces: PALETTE.lavender, angles: PALETTE.teal,
};

export class DataTable {
  // spec: {id, columns, cols (key -> typed array), n, rowResidues(i) -> int[], isSelected(i), residueRow (bool)}
  constructor(container, spec, state = {}) {
    this.el = container;
    this.spec = spec;
    this.filters = state.filters || {};
    this.sort = state.sort || null;          // {key, dir}
    this.hidden = state.hidden || new Set(spec.columns.filter((c) => c.hidden).map((c) => c.key));
    this.selOnly = state.selOnly || false;
    this.build();
    this.apply();
  }

  state() { return { filters: this.filters, sort: this.sort, hidden: this.hidden, selOnly: this.selOnly }; }

  visibleColumns() { return this.spec.columns.filter((c) => !this.hidden.has(c.key)); }

  width(c) {
    if (c._w) return c._w;
    let base = c.label.length + (c.unit ? c.unit.length + 1 : 0);
    if (c.type === 'cat') base = Math.max(base, Math.min(22, d3.max(c.categories, (s) => s.length) || 4));
    else base = Math.max(base, 6);
    c._w = Math.max(44, Math.min(230, base * 7 + 18));
    return c._w;
  }

  build() {
    const el = d3.select(this.el).html('');
    const vt = el.append('div').attr('class', 'vt');
    this.root = vt.node();          // the dock re-mounts the table when another tab has replaced this element
    this.head = vt.append('div').attr('class', 'vt-head');
    this.body = vt.append('div').attr('class', 'vt-bodyhost').append('div').attr('class', 'vt-body');
    this.rowsEl = this.body.append('div').attr('class', 'vt-rows');
    this.bars = attachScrollbars(this.body.node(), { x: true, y: true });
    this.body.on('scroll', () => { this.head.node().scrollLeft = this.body.node().scrollLeft; this.renderRows(); });
    this.body.on('mouseleave', () => store.setHover(null, 'table'));
    new ResizeObserver(() => this.renderRows()).observe(this.body.node());
    this.renderHead();
  }

  renderHead() {
    const cols = this.visibleColumns(), total = d3.sum(cols, (c) => this.width(c));
    const head = this.head.html('');
    // group row
    const groups = [];
    for (const c of cols) { const last = groups[groups.length - 1]; if (last && last.g === c.group) last.w += this.width(c); else groups.push({ g: c.group, w: this.width(c) }); }
    head.append('div').attr('class', 'vt-hrow vt-grp').style('width', `${total}px`).selectAll('div').data(groups).join('div')
      .attr('class', 'vt-cell').style('width', (g) => `${g.w}px`).style('color', (g) => GROUP_COLOR[g.g] || PALETTE.label)
      .style('border-bottom', (g) => `1px solid ${GROUP_COLOR[g.g] || PALETTE.edge}`).text((g) => g.g);
    // header row
    const hr = head.append('div').attr('class', 'vt-hrow').style('width', `${total}px`);
    hr.selectAll('div').data(cols).join('div').attr('class', 'vt-cell vt-hcell').style('width', (c) => `${this.width(c)}px`)
      .attr('title', (c) => `${c.label}${c.unit ? ` (${c.unit})` : ''}\n${c.desc || ''}\n\nclick to sort`)
      .html((c) => `${c.label}${c.unit ? `<span class="unit">${c.unit}</span>` : ''}${this.sort && this.sort.key === c.key ? `<span class="arrow">${this.sort.dir > 0 ? '▲' : '▼'}</span>` : ''}`)
      .on('click', (ev, c) => {
        if (!this.sort || this.sort.key !== c.key) this.sort = { key: c.key, dir: 1 };
        else if (this.sort.dir === 1) this.sort.dir = -1; else this.sort = null;
        this.renderHead(); this.apply();
      });
    // histogram row
    const hh = head.append('div').attr('class', 'vt-hrow').style('width', `${total}px`);
    const self = this;
    hh.selectAll('div').data(cols).join('div').attr('class', 'vt-cell').style('width', (c) => `${this.width(c)}px`).style('padding', '0 3px')
      .each(function (c) { if (c.type === 'num') self.histogram(d3.select(this), c); });
    // filter row
    const fr = head.append('div').attr('class', 'vt-hrow').style('width', `${total}px`);
    fr.selectAll('div').data(cols).join('div').attr('class', 'vt-cell vt-fcell').style('width', (c) => `${this.width(c)}px`)
      .each(function (c) {
        const cell = d3.select(this), f = self.filters[c.key] || '';
        if (c.type === 'cat' && c.categories.length <= 14) {
          const s = cell.append('select').classed('active', !!f).on('change', (ev) => { self.setFilter(c.key, ev.target.value); });
          s.append('option').attr('value', '').text('all');
          c.categories.forEach((v) => s.append('option').attr('value', `=${v}`).text(v === '' ? '(blank)' : v));
          s.property('value', f);
        } else {
          cell.append('input').attr('placeholder', c.type === 'num' ? '<, >, a..b' : 'contains').property('value', f).classed('active', !!f)
            .attr('title', c.type === 'num' ? 'numeric filter: >x  <x  >=x  <=x  =x  a..b  blank  !blank' : 'text filter: contains (case-insensitive); prefix ! to exclude; =exact')
            .on('input', debounce((ev) => self.setFilter(c.key, ev.target.value), 220));
        }
      });
  }

  histogram(cell, c) {
    const arr = this.spec.cols[c.key], vals = [];
    for (let i = 0; i < arr.length; i++) if (Number.isFinite(arr[i])) vals.push(arr[i]);
    if (vals.length < 2) return;
    const w = this.width(c) - 6, h = 14;
    const ext = d3.extent(vals);
    if (ext[0] === ext[1]) return;
    const bins = d3.bin().domain(ext).thresholds(16)(vals), ymax = d3.max(bins, (b) => b.length);
    const x = d3.scaleLinear(ext, [0, w]);
    const color = GROUP_COLOR[c.group] || PALETTE.sky;
    const svg = cell.append('svg').attr('class', 'vt-hist').attr('width', w).attr('height', h);
    svg.selectAll('rect').data(bins).join('rect').attr('x', (b) => x(b.x0)).attr('width', (b) => Math.max(1, x(b.x1) - x(b.x0) - 0.5))
      .attr('y', (b) => h - (h - 1) * Math.sqrt(b.length / ymax)).attr('height', (b) => (h - 1) * Math.sqrt(b.length / ymax))
      .attr('fill', color).attr('fill-opacity', 0.55).style('cursor', 'pointer')
      .append('title').text((b) => `${fmtNum(b.x0, c)} … ${fmtNum(b.x1, c)}: ${b.length}\nclick to filter`);
    svg.selectAll('rect').on('click', (ev, b) => this.setFilter(c.key, `${+b.x0.toPrecision(4)}..${+b.x1.toPrecision(4)}`, true));
  }

  setFilter(key, value, rebuildHead = false) {
    if (value) this.filters[key] = value; else delete this.filters[key];
    if (rebuildHead) this.renderHead();
    else this.head.selectAll('.vt-fcell').filter((c) => c.key === key).select('input, select').classed('active', !!value);
    this.apply();
  }
  clearFilters() { this.filters = {}; this.renderHead(); this.apply(); }
  setSelOnly(on) { this.selOnly = on; this.apply(); }
  setHidden(key, hide) { if (hide) this.hidden.add(key); else this.hidden.delete(key); this.renderHead(); this.renderRows(true); }

  predicate(c, f) {
    const arr = this.spec.cols[c.key];
    f = f.trim();
    if (c.type === 'num') {
      if (f === 'blank') return (i) => !Number.isFinite(arr[i]);
      if (f === '!blank' || f === '*') return (i) => Number.isFinite(arr[i]);
      let m = f.match(/^(-?[\d.]+)\s*\.\.\s*(-?[\d.]+)$/);
      if (m) { const a = +m[1], b = +m[2]; return (i) => arr[i] >= a && arr[i] <= b; }
      m = f.match(/^(<=|>=|<|>|=|!=)\s*(-?[\d.]+)$/);
      if (m) {
        const v = +m[2];
        return { '<': (i) => arr[i] < v, '>': (i) => arr[i] > v, '<=': (i) => arr[i] <= v, '>=': (i) => arr[i] >= v,
          '=': (i) => Math.abs(arr[i] - v) < 1e-6 * Math.max(1, Math.abs(v)), '!=': (i) => !(Math.abs(arr[i] - v) < 1e-6) }[m[1]];
      }
      if (/^-?[\d.]+$/.test(f)) { const v = +f; return (i) => Math.abs(arr[i] - v) < 1e-6 * Math.max(1, Math.abs(v)); }
      return null;
    }
    const cats = c.categories;
    if (f.startsWith('=')) { const code = cats.indexOf(f.slice(1)); return (i) => arr[i] === code; }
    const neg = f.startsWith('!'), q = (neg ? f.slice(1) : f).toLowerCase();
    const ok = new Uint8Array(cats.length);
    cats.forEach((s, k) => { ok[k] = s.toLowerCase().includes(q) ? 1 : 0; });
    return neg ? (i) => !ok[arr[i]] : (i) => !!ok[arr[i]];
  }

  apply() {
    const n = this.spec.n, preds = [];
    for (const [key, f] of Object.entries(this.filters)) {
      const c = this.spec.columns.find((x) => x.key === key);
      if (!c || !f) continue;
      const p = this.predicate(c, f);
      if (p) preds.push(p);
    }
    let out = [];
    const selTest = this.selOnly && store.sel.size ? this.spec.involvesSelection : null;
    for (let i = 0; i < n; i++) {
      if (selTest && !selTest(i)) continue;
      let ok = true;
      for (const p of preds) if (!p(i)) { ok = false; break; }
      if (ok) out.push(i);
    }
    if (this.sort) {
      const c = this.spec.columns.find((x) => x.key === this.sort.key), arr = this.spec.cols[this.sort.key], dir = this.sort.dir;
      if (!c) this.sort = null;
    }
    if (this.sort) {
      const c = this.spec.columns.find((x) => x.key === this.sort.key), arr = this.spec.cols[this.sort.key], dir = this.sort.dir;
      if (c.type === 'num') out.sort((a, b) => {
        const va = arr[a], vb = arr[b], fa = Number.isFinite(va), fb = Number.isFinite(vb);
        if (!fa || !fb) return fa === fb ? a - b : fa ? -1 : 1;
        return (va - vb) * dir || a - b;
      });
      else out.sort((a, b) => (c.categories[arr[a]] < c.categories[arr[b]] ? -dir : c.categories[arr[a]] > c.categories[arr[b]] ? dir : a - b));
    }
    // rows involving the selection first (stable: the sort order holds within both parts)
    this.nInvolved = 0;
    if (this.spec.selectionFirst && store.sel.size && !selTest) {
      const inv = [], rest = [];
      for (const i of out) (this.spec.involvesSelection(i) ? inv : rest).push(i);
      this.nInvolved = inv.length;
      out = inv.concat(rest);
    }
    this.order = out;
    this.rowsEl.style('height', `${out.length * ROW_H}px`).style('width', `${d3.sum(this.visibleColumns(), (c) => this.width(c))}px`);
    this.renderRows(true);
    if (this.onChange) this.onChange(this);
  }

  renderRows(force = false) {
    const body = this.body.node(), H = body.clientHeight;
    if (!this.order) return;
    const first = Math.max(0, Math.floor(body.scrollTop / ROW_H) - 8), last = Math.min(this.order.length, Math.ceil((body.scrollTop + H) / ROW_H) + 8);
    if (!force && this.window && this.window[0] === first && this.window[1] === last) { this.markRows(); return; }
    this.window = [first, last];
    const cols = this.visibleColumns(), spec = this.spec;
    const rows = this.rowsEl.selectAll('div.vt-row').data(this.order.slice(first, last), (i) => i);
    rows.exit().remove();
    const enter = rows.enter().append('div').attr('class', 'vt-row')
      .on('click', (ev, i) => spec.onRowClick(i, ev))
      .on('mouseenter', (ev, i) => { const rs = spec.rowResidues(i); store.setHover(rs.length === 1 ? rs[0] : rs[0] ?? null, 'table'); });
    const all = enter.merge(rows);
    const pos = new Map(this.order.slice(first, last).map((r, k) => [r, (first + k) * ROW_H]));
    all.style('top', (i) => `${pos.get(i)}px`);
    if (force || enter.size()) {
      all.each(function (i) {
        const row = d3.select(this);
        row.selectAll('div.vt-cell').data(cols, (c) => c.key).join('div')
          .attr('class', (c) => `vt-cell${c.type === 'num' ? ' num' : ''}`).style('width', (c) => `${c._w}px`)
          .each(function (c) {
            const v = spec.cols[c.key][i];
            if (c.type === 'num') {
              const fin = Number.isFinite(v);
              this.textContent = fin ? fmtNum(v, c) : '–';
              this.classList.toggle('blank', !fin);
            } else {
              const s = c.categories[v];
              this.textContent = s || '';
              this.title = s && s.length > 14 ? s : '';
            }
          });
      });
    }
    this.markRows();
  }

  markRows() {
    const spec = this.spec;
    this.rowsEl.selectAll('div.vt-row').classed('sel', (i) => spec.isSelected(i)).classed('hov', (i) => spec.isHovered(i));
  }

  // first visible-order row that involves the selection (a residue, or a tetrahedron containing one)
  scrollToSelection() { this.scrollToIndex(this.order.findIndex((i) => this.spec.isSelected(i))); }

  scrollToIndex(k) {
    if (k < 0) return;
    const body = this.body.node(), y = k * ROW_H;
    if (y < body.scrollTop || y > body.scrollTop + body.clientHeight - ROW_H) body.scrollTop = Math.max(0, y - body.clientHeight / 3);
  }

  toCSV() {
    const cols = this.visibleColumns(), lines = [cols.map((c) => c.key).join(',')];
    for (const i of this.order) {
      lines.push(cols.map((c) => {
        const v = this.spec.cols[c.key][i];
        if (c.type === 'num') return Number.isFinite(v) ? +v.toPrecision(6) : '';
        const s = c.categories[v] || '';
        return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
      }).join(','));
    }
    return lines.join('\n');
  }
}

export function fmtNum(v, c) {
  const f = c.fmt || '.2f';
  if (f === 'd') return d3.format('d')(v);
  try { return d3.format(f)(v); } catch (e) { return String(v); }
}

function debounce(fn, ms) { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; }
