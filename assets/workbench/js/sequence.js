// Sequence viewer (per seq-viewer.jpeg): chains | residue squares + SS + domains | HETATM checkboxes.
import { store, KIND, SS_COLOR, CLASS_COLOR, KD, PLDDT_COLOR, chainColor, PALETTE } from './store.js';
import { attachScrollbars } from './scrollbars.js';

const d3 = window.d3;
const CELL = 13, GAP = 2, PITCH = CELL + GAP, PAD = 8;
const Y = { ruler: 9, square: 14, letter: 37, ss: 44, dom: 52 };
const LANE = 14, BAR = 11;          // domain lanes: pitch and bar height
const HEIGHT_MIN = 76;
export const DOMAIN_COLOR = { CATH: PALETTE.lavender, Pfam: PALETTE.mint };   // UniProt features: teal
export const domainColor = (d) => DOMAIN_COLOR[d.source] || PALETTE.teal;

// inside the bioGeometry app the domain / ligand / ion / water panel lives in the app's sidebar
// (#bg-selection-panel, same origin); standalone, or from another origin, it stays beside the sequence
function sidebarHost() {
  try { return window.parent !== window ? window.parent.document.getElementById('bg-selection-panel') : null; } catch (e) { return null; }
}

export class SequenceViewer {
  constructor() {
    this.chainsEl = document.getElementById('seq-chains');
    this.hetEl = sidebarHost() || document.getElementById('seq-het');
    this.scroll = document.getElementById('seq-scroll');
    this.spacer = document.getElementById('seq-spacer');
    this.canvas = document.getElementById('seq-canvas');
    this.ctx = this.canvas.getContext('2d');
    this.info = document.getElementById('seq-info');
    this.colorSel = document.getElementById('seq-color');
    this.residues = [];       // residue indices of the focus chain (polymer only)
    this.colors = [];
    this.drag = null;
    this.hoverIdx = null;
    this.hoverDom = null;
    this.domRects = [];
    attachScrollbars(this.scroll, { x: true, y: false });
    attachScrollbars(this.chainsEl, { x: false, y: true });
    if (this.hetEl.id === 'seq-het') attachScrollbars(this.hetEl, { x: false, y: true });
    else {                                  // in the sidebar: the sequence takes the freed column
      document.querySelector('.seq-body').classList.add('no-het');
      const accent = () => this.hetEl.style.setProperty('--sel', PALETTE.sel);
      accent(); store.on('selstyle', accent);
    }

    this.scroll.addEventListener('scroll', () => this.draw());
    new ResizeObserver(() => this.resize()).observe(this.scroll);
    this.colorSel.addEventListener('change', () => { this.recolor(); this.draw(); });
    document.getElementById('seq-clear').addEventListener('click', () => store.clearAll('seq'));

    this.scroll.addEventListener('pointerdown', (e) => this.onDown(e));
    this.scroll.addEventListener('pointermove', (e) => this.onMove(e));
    window.addEventListener('pointerup', (e) => this.onUp(e));
    this.scroll.addEventListener('pointerleave', () => { this.hoverIdx = null; this.hoverDom = null; store.setHover(null, 'seq'); this.draw(); });
    this.scroll.addEventListener('wheel', (e) => {     // vertical wheel scrolls the sequence horizontally
      if (Math.abs(e.deltaY) > Math.abs(e.deltaX)) { this.scroll.scrollLeft += e.deltaY; e.preventDefault(); }
    }, { passive: false });

    store.on('data', () => this.load());
    store.on('features', () => { if (this.colorSel.value === 'rsa') { this.recolor(); this.draw(); } });
    store.on('selection', () => { this.draw(); this.markChains(); this.markHetero(); });
    store.on('hover', ({ source }) => { if (source !== 'seq') this.draw(); });
    store.on('selstyle', () => this.draw());
    store.on('focus', ({ source }) => { this.setChain(store.focusChain, source !== 'seq'); });
  }

  // ---------------------------------------------------------------- data
  load() {
    const S = store.S;
    this.renderChains();
    this.renderHetero();
    const first = S.stats.chain_rows.findIndex((c) => c.n_res > 0);
    store.focusChain = Math.max(0, first);
    this.setChain(store.focusChain, false);
  }

  setChain(ci, scrollToSel) {
    const S = store.S;
    if (!S) return;
    const b = S.buf;
    this.residues = [];
    for (let r = 0; r < b.res_chain.length; r++) if (b.res_chain[r] === ci && b.res_kind[r] <= KIND.nucleic) this.residues.push(r);
    this.pos = new Map(this.residues.map((r, i) => [r, i]));
    const row = S.stats.chain_rows[ci];
    this.domains = (S.stats.domains || []).filter((d) => d.chain === row.id);
    this.info.textContent = `${S.stats.name} · chain ${row.id} · ${row.n_res} ${row.kind === 'nucleic' ? 'nt' : 'res'}` +
      (row.description ? ` · ${row.description}` : '') + (row.uniprot ? ` · UniProt ${row.uniprot}` : '');
    this.spacer.style.width = `${PAD * 2 + this.residues.length * PITCH}px`;
    this.recolor();
    this.markChains();
    if (scrollToSel) this.scrollToSelection();
    this.draw();
  }

  recolor() {
    if (!store.S) return;
    const mode = this.colorSel.value, S = store.S, b = S.buf;
    if (mode === 'rsa' && store.F) { this.colors = store.residueColors('rsa'); return; }
    if (mode === 'ss') { this.colors = store.residueColors('ss'); return; }
    if (mode === 'conf') { this.colors = store.residueColors('conf'); return; }
    if (mode === 'hydropathy') { this.colors = store.residueColors('hydropathy'); return; }
    this.colors = store.residueColors('class');
  }

  // ---------------------------------------------------------------- chains & hetero lists
  renderChains() {
    const S = store.S, n = S.stats.chains.length;
    const rows = d3.select(this.chainsEl).html('').selectAll('div').data(S.stats.chain_rows).join('div')
      .attr('class', 'chain-row').attr('title', (c) => `${c.id}: ${c.description || c.kind} — ${c.n_res} residues` +
        (c.n_hetero ? `, ${c.n_hetero} hetero` : '') + '\nclick: show & select chain · shift: add · ctrl: toggle')
      .on('click', (ev, c) => {
        const rs = store.chainResidues(c.index);
        store.setFocusChain(c.index, 'seq');
        this.setChain(c.index, false);
        store.pick(rs, ev, 'seq');
      });
    rows.append('span').attr('class', 'cid').style('background', (c) => chainColor(c.index, n)).text((c) => c.id);
    rows.append('span').attr('class', 'cdesc').text((c) => c.description || c.kind);
    rows.append('span').attr('class', 'cmeta').text((c) => (c.n_res ? `${c.n_res}` : `${c.n_hetero}h`));
  }

  markChains() {
    if (!store.S) return;
    const selCh = new Set();
    for (const r of store.sel) selCh.add(store.S.buf.res_chain[r]);
    d3.select(this.chainsEl).selectAll('.chain-row')
      .classed('current', (c) => c.index === store.focusChain).classed('sel', (c) => selCh.has(c.index));
  }

  renderHetero() {
    const H = store.S.stats.hetero || [], D = store.S.stats.domains || [];
    const el = d3.select(this.hetEl).html('');
    if (D.length) {
      el.append('div').attr('class', 'het-group').text(`DOMAINS (${D.length})`);
      const it = el.selectAll(null).data(D).enter().append('label').attr('class', 'het-item dom')
        .attr('title', (d) => `${d.source} ${d.id} — ${d.name}\nchain ${d.chain} · residues ${d.start}–${d.end} (${d.n})` +
          (d.detail ? `\n${d.detail}` : '') + '\ntick: add to / remove from the selection');
      it.append('input').attr('type', 'checkbox').on('change', (ev, d) => {
        const rs = store.domainResidues(d);
        if (ev.target.checked) store.add(rs, 'het'); else store.remove(rs, 'het');
        store.setFocusChain(store.S.buf.res_chain[d.r0], 'het');
      });
      it.append('i').attr('class', 'sw').style('background', domainColor);
      it.append('span').text((d) => `${d.source} ${d.name}`);
      it.append('span').attr('class', 'rng').text((d) => `${d.chain} ${d.start}–${d.end}`);
    }
    if (!H.length) { el.append('div').attr('class', 'het-empty').text('no ligands, ions or coordinated waters'); this.markHetero(); return; }
    for (const [kind, title] of [['ligand', 'LIGANDS'], ['ion', 'IONS'], ['water', 'COORDINATED WATERS']]) {
      const items = H.filter((h) => h.kind === kind);
      if (!items.length) continue;
      el.append('div').attr('class', 'het-group').text(`${title} (${items.length})`);
      const it = el.selectAll(null).data(items).enter().append('label').attr('class', 'het-item')
        .attr('title', (h) => `${h.label} · chain ${h.chain} · ${h.name}${h.seq !== undefined ? ' ' + h.seq : ''}` +
          (h.description ? `\n${h.description}` : '') + `\nbinding site: ${(h.site || []).length} residues (◎ selects it)`);
      it.append('input').attr('type', 'checkbox').on('change', (ev, h) => {
        const members = h.members || [h.res];
        if (ev.target.checked) store.add(members, 'het'); else store.remove(members, 'het');
        if (members.length) store.setFocusChain(store.S.buf.res_chain[members[0]], 'het');
      });
      it.append('span').text((h) => h.label);
      it.filter((h) => (h.site || []).length > 0).append('span').attr('class', 'site').text('◎')
        .on('click', (ev, h) => {
          ev.preventDefault(); ev.stopPropagation();
          const members = h.members || [h.res];
          store.setSelection([...members, ...h.site], 'het');
        });
      // ligands and ions: the coordination sphere (binding residues + coordinating waters + what those waters
      // bridge to) selected and measured in one click
      it.filter((h) => h.kind !== 'water').append('span').attr('class', 'site coord').text('⌬')
        .attr('title', 'select this group with its coordinating residues and waters, and measure their geometry')
        .on('click', (ev, h) => {
          ev.preventDefault(); ev.stopPropagation();
          const rs = this.coordinationSphere(h);
          store.setSelection(rs, 'het');
          store.setFocusChain(store.S.buf.res_chain[h.res], 'het');
          store.emit('focus-selection');
          store.emit('measure');
        });
    }
    this.markHetero();
  }

  // the group, every residue bound to it (metal coordination, H-bond, covalent link), the coordinating waters
  // and the residues those waters bridge to; falls back to the 4.5 Å site when nothing binds it polar-wise
  coordinationSphere(h) {
    const b = store.S.buf, P = b.bind_pairs, T = b.bind_type, R = b.atom_res, members = h.members || [h.res];
    const M = new Set(members), out = new Set(members), waters = new Set();
    for (let k = 0; k < T.length; k++) {
      const ra = R[P[2 * k]], rb = R[P[2 * k + 1]];
      if (!M.has(ra) && !M.has(rb)) continue;
      const other = M.has(ra) ? rb : ra;
      out.add(other);
      if (b.res_kind[other] === KIND.water) waters.add(other);
    }
    for (let k = 0; k < T.length; k++) {
      const ra = R[P[2 * k]], rb = R[P[2 * k + 1]];
      if (waters.has(ra) && !M.has(rb)) out.add(rb);
      if (waters.has(rb) && !M.has(ra)) out.add(ra);
    }
    if (out.size === members.length) (h.site || []).forEach((r) => out.add(r));
    return [...out];
  }

  markHetero() {
    d3.select(this.hetEl).selectAll('.het-item').each(function (h) {
      const members = h.r0 !== undefined ? store.domainResidues(h) : h.members || [h.res];
      const on = members.some((m) => store.sel.has(m));
      this.classList.toggle('on', on);
      const box = this.querySelector('input');
      if (box) box.checked = members.every((m) => store.sel.has(m));
    });
  }

  // ---------------------------------------------------------------- drawing
  resize() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = this.scroll.clientWidth, h = Math.max(HEIGHT_MIN, this.scroll.clientHeight);
    this.canvas.width = Math.round(w * dpr); this.canvas.height = Math.round(h * dpr);
    this.canvas.style.width = `${w}px`; this.canvas.style.height = `${h}px`;
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    this.draw();
  }

  xOf(i) { return PAD + i * PITCH - this.scroll.scrollLeft; }
  idxAt(x) {
    const i = Math.floor((x + this.scroll.scrollLeft - PAD) / PITCH);
    return i >= 0 && i < this.residues.length ? i : null;
  }

  draw() {
    const ctx = this.ctx, S = store.S;
    const w = this.scroll.clientWidth, h = this.scroll.clientHeight;
    ctx.clearRect(0, 0, w, h);
    if (!S || !this.residues.length) return;
    const b = S.buf, aa = S.stats.aa;
    const i0 = Math.max(0, Math.floor((this.scroll.scrollLeft - PAD) / PITCH) - 1);
    const i1 = Math.min(this.residues.length, i0 + Math.ceil(w / PITCH) + 3);
    const anySel = store.sel.size > 0;
    ctx.font = '10px "FantasqueSansM Nerd Font", "Fantasque Sans Mono", monospace';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'alphabetic';
    // ruler
    ctx.fillStyle = PALETTE.dim;
    for (let i = i0; i < i1; i++) {
      const r = this.residues[i], n = b.res_seq[r];
      if (n % 10 === 0 || i === 0) {
        const x = this.xOf(i) + CELL / 2;
        ctx.fillText(String(n), x, Y.ruler);
        ctx.fillRect(x - 0.5, Y.ruler + 1, 1, 3);
      }
    }
    // squares + letters
    for (let i = i0; i < i1; i++) {
      const r = this.residues[i], x = this.xOf(i), sel = store.sel.has(r), hov = store.hover === r;
      ctx.globalAlpha = anySel && !sel && store.dim ? 0.45 : 1;
      ctx.fillStyle = this.colors[r] || PALETTE.label;
      ctx.fillRect(x, Y.square, CELL, CELL);
      if (sel) {
        ctx.globalAlpha = 1;
        ctx.strokeStyle = PALETTE.sel; ctx.lineWidth = 1.5;
        ctx.shadowColor = PALETTE.sel; ctx.shadowBlur = 2;
        ctx.strokeRect(x - 1, Y.square - 1, CELL + 2, CELL + 2);
        ctx.shadowBlur = 0;
      }
      if (hov) { ctx.globalAlpha = 1; ctx.strokeStyle = PALETTE.lavender; ctx.lineWidth = 1.5; ctx.strokeRect(x - 2, Y.square - 2, CELL + 4, CELL + 4); }
      ctx.globalAlpha = anySel && !sel && store.dim ? 0.5 : 1;
      ctx.fillStyle = sel ? PALETTE.sel : PALETTE.fg;
      ctx.fillText(aa[r], x + CELL / 2, Y.letter);
    }
    ctx.globalAlpha = 1;
    // secondary-structure track
    let i = i0;
    while (i < i1) {
      const r = this.residues[i], ss = b.res_ss[r];
      let j = i + 1;
      while (j < i1 && b.res_ss[this.residues[j]] === ss) j++;
      const x0 = this.xOf(i), x1 = this.xOf(j - 1) + CELL;
      if (ss === 1) { ctx.fillStyle = SS_COLOR[1]; ctx.fillRect(x0, Y.ss - 2, x1 - x0, 5); }
      else if (ss === 2) {
        ctx.fillStyle = SS_COLOR[2]; ctx.beginPath();
        ctx.moveTo(x0, Y.ss - 1.5); ctx.lineTo(x1 - 5, Y.ss - 1.5); ctx.lineTo(x1 - 5, Y.ss - 4); ctx.lineTo(x1, Y.ss + 0.5);
        ctx.lineTo(x1 - 5, Y.ss + 5); ctx.lineTo(x1 - 5, Y.ss + 2.5); ctx.lineTo(x0, Y.ss + 2.5); ctx.closePath(); ctx.fill();
      } else { ctx.fillStyle = PALETTE.dim; ctx.fillRect(x0, Y.ss, x1 - x0, 1); }
      i = j;
    }
    // domain track: one clickable bar per domain, stacked in lanes where domains overlap
    this.domRects = [];
    const lanes = [];
    for (const d of this.domains) {
      const a = this.pos.get(d.r0), z = this.pos.get(d.r1);
      if (a === undefined || z === undefined) continue;
      let lane = lanes.findIndex((end) => end < a);
      if (lane < 0) { lane = lanes.length; lanes.push(-1); }
      lanes[lane] = z;
      const x0 = this.xOf(a), x1 = this.xOf(z) + CELL, y = Y.dom + lane * LANE;
      if (x1 < 0 || x0 > w || y + BAR > h) continue;
      const col = domainColor(d), rs = store.domainResidues(d);
      let nSel = 0; for (const r of rs) if (store.sel.has(r)) nSel++;
      const full = nSel === rs.length && nSel > 0, part = nSel > 0 && !full, hov = this.hoverDom === d;
      const bx0 = Math.max(x0, -4), bx1 = Math.min(x1, w + 4);
      ctx.globalAlpha = full ? 0.3 : hov ? 0.3 : part ? 0.2 : 0.16;
      ctx.fillStyle = full ? PALETTE.sel : col;
      roundRect(ctx, bx0, y, bx1 - bx0, BAR, 3); ctx.fill();
      ctx.globalAlpha = 1;
      ctx.lineWidth = hov || full ? 1.5 : 1;
      ctx.strokeStyle = full || part ? PALETTE.sel : hov ? PALETTE.fg : col;
      ctx.setLineDash(part ? [3, 2] : []);
      roundRect(ctx, bx0 + 0.5, y + 0.5, bx1 - bx0 - 1, BAR - 1, 3); ctx.stroke();
      ctx.setLineDash([]);
      ctx.save();
      ctx.beginPath(); ctx.rect(bx0 + 2, y, bx1 - bx0 - 4, BAR); ctx.clip();
      ctx.textAlign = 'left';
      ctx.fillStyle = full ? PALETTE.fg : hov ? PALETTE.fg : col;
      ctx.fillText(`${d.source} · ${d.name}${part ? `  (${nSel}/${rs.length})` : ''}`, Math.max(x0 + 5, 5), y + BAR - 2);
      ctx.restore();
      ctx.textAlign = 'center';
      this.domRects.push({ d, x0, x1, y0: y - 1, y1: y + BAR + 1 });
    }
    // drag band
    if (this.drag && this.drag.moved) {
      const a = Math.min(this.drag.i0, this.drag.i1), z = Math.max(this.drag.i0, this.drag.i1);
      ctx.fillStyle = 'rgba(244, 244, 255, 0.10)'; ctx.strokeStyle = PALETTE.sel;
      ctx.fillRect(this.xOf(a) - 1, Y.square - 3, this.xOf(z) - this.xOf(a) + CELL + 2, CELL + 30);
    }
  }

  scrollToSelection() {
    const idx = this.residues.findIndex((r) => store.sel.has(r));
    if (idx < 0) return;
    const x = PAD + idx * PITCH, w = this.scroll.clientWidth;
    if (x < this.scroll.scrollLeft || x > this.scroll.scrollLeft + w - 40) this.scroll.scrollLeft = Math.max(0, x - w / 3);
  }

  // ---------------------------------------------------------------- interaction
  local(e) { const rc = this.scroll.getBoundingClientRect(); return [e.clientX - rc.left, e.clientY - rc.top]; }

  onDown(e) {
    if (e.button !== 0) return;
    const [x, y] = this.local(e);
    const dom = this.domAt(x, y);
    if (dom) { store.pick(store.domainResidues(dom), e, 'seq'); return; }
    const i = this.idxAt(x);
    if (i === null || y > Y.ss + 6) return;
    this.drag = { i0: i, i1: i, moved: false, ev: e };
  }

  domAt(x, y) {
    const q = this.domRects.find((t) => x >= t.x0 && x <= t.x1 && y >= t.y0 && y <= t.y1);
    return q ? q.d : null;
  }

  onMove(e) {
    const [x, y] = this.local(e);
    const i = this.idxAt(x);
    if (this.drag) {
      if (i !== null && i !== this.drag.i1) { this.drag.i1 = i; this.drag.moved = true; this.draw(); }
      return;
    }
    const dom = this.drag ? null : this.domAt(x, y);
    if (dom !== this.hoverDom) { this.hoverDom = dom; this.draw(); }
    this.scroll.style.cursor = dom ? 'pointer' : '';
    if (dom) {
      if (store.hover !== null) store.setHover(null, 'seq');
      this.scroll.title = `${dom.source} ${dom.id} — ${dom.name}\nchain ${dom.chain} · residues ${dom.start}–${dom.end} (${dom.n})` +
        (dom.detail ? `\n${dom.detail}` : '') + '\nclick: select · shift: add · ctrl: toggle · then ⌬ measure';
      return;
    }
    const r = i !== null && y < Y.ss + 6 ? this.residues[i] : null;
    if (r !== store.hover) { store.setHover(r, 'seq'); this.draw(); }
    this.scroll.title = r === null ? '' : this.tooltip(r);
  }

  onUp() {
    if (!this.drag) return;
    const { i0, i1, moved, ev } = this.drag;
    this.drag = null;
    if (moved) {
      const a = Math.min(i0, i1), z = Math.max(i0, i1);
      const rs = this.residues.slice(a, z + 1);
      store.addOrSet(rs, ev, 'seq');
      store.anchor = this.residues[i1];
    } else store.click(this.residues[i0], ev, 'seq');
    this.draw();
  }

  tooltip(r) {
    const S = store.S, b = S.buf;
    const parts = [`${store.label(r)}  ·  ${S.stats.classes[b.res_class[r]]}  ·  ${({ 1: 'helix', 2: 'strand', 0: 'coil' })[b.res_ss[r]] || ''}`];
    if (Number.isFinite(b.res_phi[r]) || Number.isFinite(b.res_psi[r])) parts.push(`φ ${fmtA(b.res_phi[r])}  ψ ${fmtA(b.res_psi[r])}`);
    parts.push(`${S.stats.is_alphafold ? 'pLDDT' : 'B'} ${b.res_b[r].toFixed(1)}`);
    if (store.F) {
      const rsa = store.featureValue('rsa', r), site = store.featureValue('site', r);
      if (Number.isFinite(rsa)) parts.push(`RSA ${rsa.toFixed(2)} (${store.featureValue('exposure', r)})`);
      if (site) parts.push(`ligand site: ${site}`);
    }
    return parts.join('\n');
  }
}

const fmtA = (v) => (Number.isFinite(v) ? `${v.toFixed(0)}°` : '–');
function roundRect(ctx, x, y, w, h, r) {
  ctx.beginPath();
  if (ctx.roundRect) ctx.roundRect(x, y, Math.max(0, w), h, r); else ctx.rect(x, y, Math.max(0, w), h);
}
