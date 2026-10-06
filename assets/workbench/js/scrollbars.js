// Scrollbars that live in their own gutters instead of floating over the content.
//
// Firefox on Linux (GTK) and macOS draw overlay scrollbars that appear and widen on hover, covering
// the last column / bottom row and swallowing clicks meant for them. The viewport keeps native
// scrolling (wheel, touchpad, keyboard, momentum) with its own bars hidden; thin tracks are drawn in
// reserved gutters beside it and can be dragged or clicked to page.
export const GUTTER = 9;

export function attachScrollbars(viewport, { x = true, y = true } = {}) {
  const host = viewport.parentElement;
  host.classList.add('cs-host');
  host.classList.toggle('cs-x', x);
  host.classList.toggle('cs-y', y);
  viewport.classList.add('cs-viewport');
  const bars = {};
  for (const axis of ['y', 'x']) {
    if (axis === 'y' ? !y : !x) continue;
    const track = document.createElement('div'), thumb = document.createElement('div');
    track.className = `cs-track cs-${axis}`;
    thumb.className = 'cs-thumb';
    track.append(thumb);
    host.append(track);
    bars[axis] = { track, thumb };
    const pos = axis === 'y' ? 'scrollTop' : 'scrollLeft', client = axis === 'y' ? 'clientY' : 'clientX';
    thumb.addEventListener('pointerdown', (e) => {
      e.preventDefault(); e.stopPropagation();
      thumb.setPointerCapture(e.pointerId);
      track.classList.add('drag');
      const start = viewport[pos], from = e[client], m = metrics(axis);
      const move = (ev) => { viewport[pos] = start + (ev[client] - from) * m.ratio; };
      const up = () => { track.classList.remove('drag'); thumb.removeEventListener('pointermove', move); thumb.removeEventListener('pointerup', up); };
      thumb.addEventListener('pointermove', move);
      thumb.addEventListener('pointerup', up);
    });
    track.addEventListener('pointerdown', (e) => {          // click on the track: page towards the pointer
      if (e.target !== track) return;
      const r = thumb.getBoundingClientRect(), page = axis === 'y' ? viewport.clientHeight : viewport.clientWidth;
      viewport[pos] += (e[client] < (axis === 'y' ? r.top : r.left) ? -1 : 1) * page * 0.9;
    });
  }
  function metrics(axis) {
    const view = axis === 'y' ? viewport.clientHeight : viewport.clientWidth;
    const full = axis === 'y' ? viewport.scrollHeight : viewport.scrollWidth;
    const len = (axis === 'y' ? bars.y.track.clientHeight : bars.x.track.clientWidth) || 1;
    const size = Math.max(18, len * Math.min(1, view / Math.max(1, full)));
    return { view, full, len, size, ratio: (full - view) / Math.max(1, len - size) };
  }
  function update() {
    for (const axis of Object.keys(bars)) {
      const { track, thumb } = bars[axis], m = metrics(axis);
      const scrollable = m.full > m.view + 1;
      track.classList.toggle('idle', !scrollable);
      if (!scrollable) continue;
      const off = (axis === 'y' ? viewport.scrollTop : viewport.scrollLeft) / m.ratio;
      if (axis === 'y') { thumb.style.height = `${m.size}px`; thumb.style.transform = `translateY(${off}px)`; }
      else { thumb.style.width = `${m.size}px`; thumb.style.transform = `translateX(${off}px)`; }
    }
  }
  let queued = false;
  const schedule = () => { if (!queued) { queued = true; requestAnimationFrame(() => { queued = false; update(); }); } };
  viewport.addEventListener('scroll', schedule, { passive: true });
  // content size changes do not resize the viewport: watch the children too
  const ro = new ResizeObserver(schedule);
  ro.observe(viewport);
  const watch = () => { for (const ch of viewport.children) ro.observe(ch); schedule(); };
  new MutationObserver(watch).observe(viewport, { childList: true });
  watch();
  return { update: schedule };
}

// A scrollable box for d3 code: returns the viewport selection to append content into.
export function scrollBox(parent, { x = false, y = true, cls = '' } = {}) {
  const host = parent.append('div').attr('class', `cs-fill ${cls}`);
  const vp = host.append('div');
  attachScrollbars(vp.node(), { x, y });
  return vp;
}
