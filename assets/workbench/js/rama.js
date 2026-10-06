// Ramachandran plot, residue contact map and AlphaFold PAE — all brushable, all linked to the selection.
import { store, KIND, SS_COLOR, CLASS_COLOR, PLDDT_COLOR, PALETTE, RAMP, toast } from './store.js';

const d3 = window.d3;
// AlphaFold DB convention: dark green = low expected error, white = high
const PAE_SCALE = d3.scaleSequential([0, 30], d3.interpolateRgbBasis(['#0b5d1e', '#2fa84f', '#a8e6a1', '#f4fff2'])).clamp(true);
const REGIONS = [   // matches biochem.rama_region (coarse basins)
  { name: 'αR', color: PALETTE.rose, rects: [[-160, -120, -20, 50]] },
  { name: 'β', color: PALETTE.teal, rects: [[-180, 90, -100, 180], [-100, 90, -45, 100], [-180, -180, -45, -150]] },
  { name: 'PPII', color: PALETTE.sky, rects: [[-100, 100, -50, 180]] },
  { name: 'αL', color: PALETTE.mint, rects: [[30, -20, 100, 80]] },
];
const REGION_COLOR = { 'αR': PALETTE.rose, 'β': PALETTE.teal, 'PPII': PALETTE.sky, 'αL': PALETTE.mint, other: PALETTE.coral };

export class RamaPanel {
  constructor() {
    this.host = document.getElementById('rama-host');
    this.canvas = document.getElementById('rama-canvas');
    this.ctx = this.canvas.getContext('2d');
    this.svg = d3.select('#rama-svg');
    this.tip = document.getElementById('rama-tip');
    this.colorSel = document.getElementById('rama-color');
    this.tab = 'rama';
    this.cmap = null; this.pae = null;
    document.querySelectorAll('#rama-tabs [data-tab]').forEach((b) => b.addEventListener('click', () => this.setTab(b.dataset.tab)));
    this.colorSel.addEventListener('change', () => this.draw());
    new ResizeObserver(() => this.layout()).observe(this.host);
    this.host.addEventListener('pointermove', (e) => this.onMove(e));
    this.host.addEventListener('pointerleave', () => { this.tip.hidden = true; store.setHover(null, 'rama'); });
    store.on('data', () => this.load());
    store.on('selection', () => this.drawSelection());
    store.on('selstyle', () => this.draw());
    store.on('hover', ({ source }) => { if (source !== 'rama') this.drawSelection(); });
  }

  load() {
    const b = store.S.buf;
    this.pts = [];
    for (let r = 0; r < b.res_kind.length; r++) {
      if (b.res_kind[r] === KIND.protein && Number.isFinite(b.res_phi[r]) && Number.isFinite(b.res_psi[r])) this.pts.push(r);
    }
    this.region = this.pts.map((r) => regionOf(b.res_phi[r], b.res_psi[r]));
    this.cmap = null; this.pae = null;
    document.getElementById('tab-pae').hidden = !store.S.stats.is_alphafold;
    if (this.tab === 'pae' && !store.S.stats.is_alphafold) this.tab = 'rama';
    this.setTab(this.tab);
  }

  async setTab(t) {
    this.tab = t;
    document.querySelectorAll('#rama-tabs [data-tab]').forEach((b) => b.classList.toggle('on', b.dataset.tab === t));
    document.getElementById('rama-title').textContent = { rama: 'RAMACHANDRAN', cmap: 'CONTACT MAP', pae: 'PREDICTED ALIGNED ERROR' }[t];
    this.colorSel.hidden = t !== 'rama';
    if (!store.S) return;
    const job = store.job;
    if (t === 'cmap' && !this.cmap) {
      try { const L = await store.requestLayer('contacts', store.params.contacts || {}); if (store.job === job) this.cmap = L; } catch (e) { toast(`contacts: ${e.message}`); }
    }
    if (t === 'pae' && !this.pae) {
      try { const L = await store.requestLayer('pae'); if (store.job === job) { this.pae = L; this.paeImage = null; } } catch (e) { toast(`PAE: ${e.message}`); }
    }
    this.layout();
  }

  layout() {
    const w = this.host.clientWidth, h = this.host.clientHeight;
    if (!w || !h) return;
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    this.canvas.width = Math.round(w * dpr); this.canvas.height = Math.round(h * dpr);
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const m = { l: 34, r: this.tab === 'pae' ? 38 : 8, t: 8, b: 26 };   // PAE keeps room for its colour bar
    const side = Math.max(60, Math.min(w - m.l - m.r, h - m.t - m.b));
    this.box = { x: m.l + (w - m.l - m.r - side) / 2, y: m.t, s: side };
    this.draw();
  }

  // -------------------------------------------------------- scales
  get xs() { const B = this.box; return this.tab === 'rama' ? d3.scaleLinear([-180, 180], [B.x, B.x + B.s]) : d3.scaleLinear([0, this.nAxis], [B.x, B.x + B.s]); }
  get ys() { const B = this.box; return this.tab === 'rama' ? d3.scaleLinear([-180, 180], [B.y + B.s, B.y]) : d3.scaleLinear([0, this.nAxis], [B.y, B.y + B.s]); }

  draw() {
    const ctx = this.ctx;
    ctx.clearRect(0, 0, this.host.clientWidth, this.host.clientHeight);
    this.svg.selectAll('*').remove();
    if (!store.S || !this.box) return;
    if (this.tab === 'rama') this.drawRama();
    else if (this.tab === 'cmap') this.drawMatrix('cmap');
    else this.drawMatrix('pae');
    this.drawSelection();
  }

  drawRama() {
    const ctx = this.ctx, B = this.box, x = this.xs, y = this.ys, b = store.S.buf;
    const g = this.svg.append('g');
    // basins
    for (const reg of REGIONS) {
      for (const [x0, y0, x1, y1] of reg.rects) {
        g.append('rect').attr('class', 'region').attr('x', x(x0)).attr('y', y(y1)).attr('width', x(x1) - x(x0)).attr('height', y(y0) - y(y1))
          .attr('fill', reg.color).attr('stroke', reg.color);
      }
      const [x0, y0, x1, y1] = reg.rects[0];
      g.append('text').attr('class', 'region-label').attr('x', x(x0) + 3).attr('y', y(y1) + 11).attr('fill', reg.color).style('fill', reg.color).text(reg.name);
    }
    // density contours of this structure
    if (this.pts.length > 30) {
      const dens = d3.contourDensity().x((r) => x(b.res_phi[r])).y((r) => y(b.res_psi[r])).size([this.host.clientWidth, this.host.clientHeight])
        .bandwidth(Math.max(6, B.s / 40)).thresholds(6)(this.pts);
      g.append('g').selectAll('path').data(dens).join('path').attr('d', d3.geoPath()).attr('fill', 'none')
        .attr('stroke', PALETTE.lavender).attr('stroke-opacity', (d, i) => 0.12 + 0.05 * i).attr('stroke-width', 1);
    }
    // axes; φ = 0 and ψ = 0 drawn solid so the four quadrants stand out
    g.append('g').attr('class', 'grid').selectAll('line').data([-120, -60, 0, 60, 120]).join('line').attr('class', (v) => (v === 0 ? 'zero' : null))
      .attr('x1', (v) => x(v)).attr('x2', (v) => x(v)).attr('y1', B.y).attr('y2', B.y + B.s);
    g.append('g').attr('class', 'grid').selectAll('line').data([-120, -60, 0, 60, 120]).join('line').attr('class', (v) => (v === 0 ? 'zero' : null))
      .attr('y1', (v) => y(v)).attr('y2', (v) => y(v)).attr('x1', B.x).attr('x2', B.x + B.s);
    g.append('g').attr('class', 'axis').attr('transform', `translate(0,${B.y + B.s})`).call(d3.axisBottom(x).tickValues([-180, -90, 0, 90, 180]).tickSize(3));
    g.append('g').attr('class', 'axis').attr('transform', `translate(${B.x},0)`).call(d3.axisLeft(y).tickValues([-180, -90, 0, 90, 180]).tickSize(3));
    g.append('text').attr('class', 'axis-title').attr('x', B.x + B.s).attr('y', B.y + B.s + 22).attr('text-anchor', 'end').text('φ (°)');
    g.append('text').attr('class', 'axis-title').attr('x', B.x - 28).attr('y', B.y + 8).text('ψ');
    g.append('rect').attr('x', B.x).attr('y', B.y).attr('width', B.s).attr('height', B.s).attr('fill', 'none').attr('stroke', PALETTE.edge);
    this.brushLayer = this.svg.append('g');
    this.selLayer = this.svg.append('g');
    this.attachBrush();
    this.drawPoints();
  }

  pointColor(i, r) {
    const b = store.S.buf, mode = this.colorSel.value;
    if (mode === 'class') return CLASS_COLOR[store.S.stats.classes[b.res_class[r]]] || PALETTE.label;
    if (mode === 'region') return REGION_COLOR[this.region[i]] || PALETTE.label;
    if (mode === 'conf') return store.S.stats.is_alphafold ? PLDDT_COLOR(b.res_b[r]) : RAMP.heat(Math.min(1, b.res_b[r] / 80));
    return SS_COLOR[b.res_ss[r]] || PALETTE.label;
  }

  drawPoints() {
    const ctx = this.ctx, x = this.xs, y = this.ys, b = store.S.buf, any = store.sel.size > 0;
    const big = this.pts.length > 4000 ? 1.6 : 2.4;
    for (let pass = 0; pass < 2; pass++) {
      for (let i = 0; i < this.pts.length; i++) {
        const r = this.pts[i], sel = store.sel.has(r);
        if ((pass === 1) !== sel) continue;
        const px = x(b.res_phi[r]), py = y(b.res_psi[r]);
        ctx.globalAlpha = any && !sel && store.dim ? 0.35 : 0.9;
        ctx.fillStyle = sel ? PALETTE.sel : this.pointColor(i, r);
        const rad = sel ? big + 1.8 : big;
        const name = store.S.stats.res_names[b.res_name[r]];
        ctx.beginPath();
        if (name === 'GLY') { ctx.moveTo(px, py - rad - 1); ctx.lineTo(px + rad + 1, py + rad); ctx.lineTo(px - rad - 1, py + rad); ctx.closePath(); }
        else if (name === 'PRO') ctx.rect(px - rad, py - rad, 2 * rad, 2 * rad);
        else ctx.arc(px, py, rad, 0, 2 * Math.PI);
        ctx.fill();
        if (sel) { ctx.globalAlpha = 1; ctx.strokeStyle = '#fff'; ctx.lineWidth = 0.8; ctx.stroke(); }
      }
    }
    ctx.globalAlpha = 1;
  }

  // vertical gradient bar for a sequential scale, low values at the top like the matrix diagonal
  colorBar(scale, x0, y0, h, ticks) {
    const [lo, hi] = scale.domain(), id = `cbar-${Math.random().toString(36).slice(2, 8)}`;
    const grad = this.svg.append('defs').append('linearGradient').attr('id', id).attr('x1', 0).attr('x2', 0).attr('y1', 0).attr('y2', 1);
    d3.range(0, 1.0001, 0.1).forEach((t) => grad.append('stop').attr('offset', t).attr('stop-color', scale(lo + t * (hi - lo))));
    const g = this.svg.append('g').attr('class', 'cbar');
    g.append('rect').attr('x', x0).attr('y', y0).attr('width', 8).attr('height', h).attr('fill', `url(#${id})`).attr('stroke', PALETTE.edge);
    const ys = d3.scaleLinear([lo, hi], [y0, y0 + h]);
    g.selectAll('text').data(ticks).join('text').attr('x', x0 + 11).attr('y', (v) => ys(v) + 3).text((v) => v);
  }

  drawMatrix(kind) {
    const ctx = this.ctx, B = this.box, b = store.S.buf;
    const poly = [];
    for (let r = 0; r < b.res_kind.length; r++) if (b.res_kind[r] <= KIND.nucleic) poly.push(r);
    this.axisRes = kind === 'pae' && this.pae && this.pae.stats.available ? this.pae.stats.residues : poly;
    this.axisPos = new Map(this.axisRes.map((r, i) => [r, i]));
    this.nAxis = this.axisRes.length;
    const x = this.xs, y = this.ys;
    ctx.fillStyle = '#141420'; ctx.fillRect(B.x, B.y, B.s, B.s);
    if (kind === 'cmap') {
      if (!this.cmap) { this.message('computing contacts…'); return; }
      const e = this.cmap.buf.res_edges, n = e.length / 2, px = Math.max(1, B.s / this.nAxis);
      const ch = b.res_chain;
      for (let k = 0; k < n; k++) {
        const i = this.axisPos.get(e[2 * k]), j = this.axisPos.get(e[2 * k + 1]);
        if (i === undefined || j === undefined) continue;
        ctx.fillStyle = ch[e[2 * k]] === ch[e[2 * k + 1]] ? 'rgba(125,211,252,0.75)' : 'rgba(252,165,165,0.9)';
        ctx.fillRect(x(i), y(j), px, px); ctx.fillRect(x(j), y(i), px, px);
      }
      // chain boundaries
      ctx.strokeStyle = PALETTE.edge; ctx.lineWidth = 1;
      for (let i = 1; i < this.nAxis; i++) if (ch[this.axisRes[i]] !== ch[this.axisRes[i - 1]]) {
        ctx.beginPath(); ctx.moveTo(x(i), B.y); ctx.lineTo(x(i), B.y + B.s); ctx.moveTo(B.x, y(i)); ctx.lineTo(B.x + B.s, y(i)); ctx.stroke();
      }
      this.svg.append('text').attr('x', B.x).attr('y', B.y + B.s + 14).text(`Delaunay contacts ≤ ${this.cmap.stats.cutoff} Å · blue intra / red inter-chain`);
    } else {
      if (!this.pae) { this.message('loading PAE…'); return; }
      if (!this.pae.stats.available) { this.message(this.pae.stats.note || 'no PAE'); return; }
      const n = this.pae.stats.n, q = this.pae.buf.pae, sc = this.pae.stats.scale;
      if (!this.paeImage) {
        const img = new ImageData(n, n), col = PAE_SCALE;
        for (let i = 0; i < n * n; i++) { const c = d3.rgb(col(q[i] * sc)); img.data[4 * i] = c.r; img.data[4 * i + 1] = c.g; img.data[4 * i + 2] = c.b; img.data[4 * i + 3] = 255; }
        const off = document.createElement('canvas'); off.width = n; off.height = n; off.getContext('2d').putImageData(img, 0, 0);
        this.paeImage = off;
      }
      ctx.imageSmoothingEnabled = false;
      ctx.drawImage(this.paeImage, B.x, B.y, B.s, B.s);
      this.svg.append('text').attr('x', B.x).attr('y', B.y + B.s + 14).text(`expected position error (Å): dark = confident relative placement · max ${this.pae.stats.max.toFixed(1)}`);
      this.colorBar(PAE_SCALE, B.x + B.s + 6, B.y, B.s, [0, 5, 10, 15, 20, 25, 30]);
    }
    const g = this.svg.append('g');
    g.append('rect').attr('x', B.x).attr('y', B.y).attr('width', B.s).attr('height', B.s).attr('fill', 'none').attr('stroke', PALETTE.edge);
    const ticks = d3.ticks(0, this.nAxis - 1, 5).filter((t) => t < this.nAxis);
    g.append('g').attr('class', 'axis').attr('transform', `translate(${B.x},0)`).call(d3.axisLeft(y).tickValues(ticks).tickFormat((t) => b.res_seq[this.axisRes[t]]).tickSize(3));
    this.brushLayer = this.svg.append('g');
    this.selLayer = this.svg.append('g');
    this.attachBrush();
  }

  message(t) { this.svg.append('text').attr('x', this.box.x + 10).attr('y', this.box.y + 20).text(t); }

  attachBrush() {
    const B = this.box;
    const brush = d3.brush().extent([[B.x, B.y], [B.x + B.s, B.y + B.s]])
      .on('end', (ev) => {
        if (!ev.selection) return;
        const [[x0, y0], [x1, y1]] = ev.selection, src = ev.sourceEvent;
        let rs = [];
        if (this.tab === 'rama') {
          const b = store.S.buf, x = this.xs, y = this.ys;
          rs = this.pts.filter((r) => { const px = x(b.res_phi[r]), py = y(b.res_psi[r]); return px >= x0 && px <= x1 && py >= y0 && py <= y1; });
        } else {
          const xi = this.xs.invert, yi = this.ys.invert;
          const a = Math.max(0, Math.floor(xi(x0))), z = Math.min(this.nAxis - 1, Math.ceil(xi(x1)));
          const p = Math.max(0, Math.floor(yi(y0))), q = Math.min(this.nAxis - 1, Math.ceil(yi(y1)));
          const s = new Set();
          for (let i = a; i <= z; i++) s.add(this.axisRes[i]);
          for (let i = p; i <= q; i++) s.add(this.axisRes[i]);
          rs = [...s];
        }
        this.brushLayer.call(brush.move, null);
        store.addOrSet(rs, src, 'rama');
      });
    this.brushLayer.call(brush);
    // a plain click (no drag) selects the nearest point
    this.brushLayer.select('.overlay').on('click.pick', (ev) => {
      const r = this.nearest(ev);
      if (r !== null) store.click(r, ev, 'rama');
    });
  }

  nearest(ev) {
    const rc = this.host.getBoundingClientRect(), mx = ev.clientX - rc.left, my = ev.clientY - rc.top;
    if (this.tab === 'rama') {
      const b = store.S.buf, x = this.xs, y = this.ys;
      let best = null, bd = 64;
      for (const r of this.pts) { const d = (x(b.res_phi[r]) - mx) ** 2 + (y(b.res_psi[r]) - my) ** 2; if (d < bd) { bd = d; best = r; } }
      return best;
    }
    if (!this.nAxis || mx < this.box.x || my < this.box.y || mx > this.box.x + this.box.s || my > this.box.y + this.box.s) return null;
    return this.axisRes[Math.min(this.nAxis - 1, Math.floor(this.xs.invert(mx)))];
  }

  drawSelection() {
    if (!this.selLayer || !store.S) return;
    if (this.tab === 'rama') {
      // repaint points (selection changes their emphasis) + labels for small selections
      const ctx = this.ctx, B = this.box;
      ctx.clearRect(B.x - 6, B.y - 6, B.s + 12, B.s + 12);
      this.svg.selectAll('.region, .region-label').raise();
      this.drawPoints();
      const b = store.S.buf, x = this.xs, y = this.ys;
      const sel = this.pts.filter((r) => store.sel.has(r)).slice(0, 28);
      this.selLayer.selectAll('text').data(store.sel.size <= 28 ? sel : []).join('text').attr('class', 'sel-label')
        .attr('x', (r) => x(b.res_phi[r]) + 5).attr('y', (r) => y(b.res_psi[r]) - 4).text((r) => `${store.resName(r)}${b.res_seq[r]}`);
      const h = store.hover;
      this.selLayer.selectAll('circle').data(h !== null && this.pts.includes(h) ? [h] : []).join('circle')
        .attr('cx', (r) => x(b.res_phi[r])).attr('cy', (r) => y(b.res_psi[r])).attr('r', 7).attr('fill', 'none').attr('stroke', PALETTE.lavender).attr('stroke-width', 1.5);
    } else if (this.axisPos) {
      const x = this.xs, y = this.ys, B = this.box, w = Math.max(1, B.s / this.nAxis);
      const idx = [...store.sel].map((r) => this.axisPos.get(r)).filter((i) => i !== undefined).sort((a, b) => a - b);
      const runs = [];
      for (const i of idx) { const last = runs[runs.length - 1]; if (last && i === last[1] + 1) last[1] = i; else runs.push([i, i]); }
      this.selLayer.selectAll('rect.sx').data(runs).join('rect').attr('class', 'sx').attr('x', (d) => x(d[0])).attr('width', (d) => Math.max(w, x(d[1] + 1) - x(d[0])))
        .attr('y', B.y).attr('height', B.s).attr('fill', PALETTE.sel).attr('fill-opacity', 0.12);
      this.selLayer.selectAll('rect.sy').data(runs).join('rect').attr('class', 'sy').attr('y', (d) => y(d[0])).attr('height', (d) => Math.max(w, y(d[1] + 1) - y(d[0])))
        .attr('x', B.x).attr('width', B.s).attr('fill', PALETTE.sel).attr('fill-opacity', 0.12);
    }
  }

  onMove(e) {
    if (!store.S || e.buttons) return;
    const r = this.nearest(e);
    store.setHover(r, 'rama');
    if (r === null) { this.tip.hidden = true; return; }
    const b = store.S.buf;
    let html = `<b>${store.label(r)}</b>`;
    if (this.tab === 'rama') html += `<br><span class="k">φ/ψ</span> ${b.res_phi[r].toFixed(0)}° / ${b.res_psi[r].toFixed(0)}° · ${regionOf(b.res_phi[r], b.res_psi[r])}`;
    else if (this.tab === 'pae' && this.pae && this.pae.stats.available) {
      const rc = this.host.getBoundingClientRect(), j = Math.floor(this.ys.invert(e.clientY - rc.top)), i = this.axisPos.get(r);
      if (j >= 0 && j < this.nAxis) html += `<br><span class="k">PAE vs ${store.label(this.axisRes[j])}</span> ${(this.pae.buf.pae[i * this.pae.stats.n + j] * this.pae.stats.scale).toFixed(1)} Å`;
    }
    this.tip.innerHTML = html;
    const rc = this.host.getBoundingClientRect();
    this.tip.style.left = `${Math.min(e.clientX - rc.left + 12, rc.width - 200)}px`;
    this.tip.style.top = `${e.clientY - rc.top + 12}px`;
    this.tip.hidden = false;
  }
}

function regionOf(phi, psi) {
  if (!Number.isFinite(phi) || !Number.isFinite(psi)) return '';
  if (phi >= -100 && phi <= -50 && psi >= 100) return 'PPII';
  if (phi >= -160 && phi <= -20 && psi >= -120 && psi <= 50) return 'αR';
  if (phi >= -180 && phi <= -45 && (psi >= 90 || psi <= -150)) return 'β';
  if (phi >= 30 && phi <= 100 && psi >= -20 && psi <= 80) return 'αL';
  return 'other';
}
