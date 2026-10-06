"""End-to-end check of the workbench interaction features (standalone server on :8770).

    python scripts/e2e_features.py OUT_DIR [chromium|firefox]

Linked cameras, gutter scrollbars, domain-bar and HETATM selection, selection geometry (⌬ measure),
selection simplices in the geometry view, the selection colour, and the structure downloads.
"""
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

API = "http://127.0.0.1:8770"
# newest Chromium that `playwright install chromium` put in the cache
CHROME = str(max((Path.home() / ".cache/ms-playwright").glob("chromium-*/chrome-linux64/chrome"), default=""))


def main(out, browser="chromium"):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    fails = []

    def check(cond, msg):
        print(("  ok   " if cond else "  FAIL ") + msg)
        if not cond:
            fails.append(msg)

    with sync_playwright() as p:
        if browser == "firefox":
            b = p.firefox.launch(firefox_user_prefs={"webgl.force-enabled": True})
        else:
            b = p.chromium.launch(executable_path=CHROME, args=["--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"])
        page = b.new_page(viewport={"width": 1680, "height": 1000})
        errors = []
        page.on("console", lambda m: errors.append(m.text) if m.type == "error" else None)
        page.on("pageerror", lambda e: errors.append(str(e)))
        page.goto(f"{API}/workbench/index.html?api={API}&pdb=4HHB")
        page.wait_for_function("window.bgStandaloneDone === true", timeout=300_000)
        page.wait_for_function("window.bgStore.F !== null", timeout=300_000)
        page.wait_for_timeout(1500)
        V = "window.bgWorkbench.views"

        print(f"[{browser}] 0. ligands / ions / waters: same colours in both views, binding lines")
        same = page.evaluate(f"""() => {{ const b = window.bgStore.S.buf, mol = {V}.mol, geo = {V}.geo;
            const hemAtom = mol.het.atoms.find((a) => window.bgStore.resName(b.atom_res[a]) === 'HEM' && window.bgStore.S.stats.elements[b.atom_elem[a]] === 'C');
            const i = Array.from(mol.het.atoms).indexOf(hemAtom), j = Array.from(geo.hetExtra.atoms).indexOf(hemAtom);
            const m = mol.het.col.slice(3 * i, 3 * i + 3), g = geo.hetExtra.col.slice(3 * j, 3 * j + 3);
            return {{ j, diff: Math.max(...[0, 1, 2].map((k) => Math.abs(m[k] - g[k]))), hex: '#' + [...m].map((x) => Math.round(x * 255).toString(16).padStart(2, '0')).join('') }}; }}""")
        check(same["j"] >= 0 and same["diff"] < 1e-6, f"a haem carbon is the same colour in both views ({same['hex']}, Δ {same['diff']:.1e})")
        painted = page.evaluate(f"""() => {{ const s = window.bgStore, mol = {V}.mol, geo = {V}.geo, hem = s.S.stats.hetero.find((h) => h.name === 'HEM');
            s.setSelection([hem.res], 't');
            const a = mol.het.pts.geometry.attributes, g = geo.hetExtra.pts.geometry.attributes;
            const at = new Map(Array.from(mol.het.atoms).map((x, i) => [x, i]));
            let dc = 0, dr = 0, n = 0;
            Array.from(geo.hetExtra.atoms).forEach((x, j) => {{ const i = at.get(x); n++;
              for (let k = 0; k < 3; k++) dc = Math.max(dc, Math.abs(a.aColor.array[3 * i + k] - g.aColor.array[3 * j + k]));
              dr = Math.max(dr, Math.abs(a.aRadius.array[i] - g.aRadius.array[j])); }});
            const bc = (m) => (m ? Array.from(m.instanceColor.array).map((v) => v.toFixed(4)).join() : '');
            const lights = (sc) => sc.children.filter((o) => o.isLight).map((l) => `${{l.type}}:${{l.intensity}}:${{l.color.getHexString()}}`).join();
            const out = {{ n, dc, dr, bonds: bc(mol.hetBonds) === bc(geo.topGroup.children.find((o) => o.isInstancedMesh)),
                          lights: lights(mol.scene) === lights(geo.top) }};
            s.setSelection([], 't'); return out; }}""")
        check(painted["n"] > 0 and painted["dc"] < 1e-6 and painted["dr"] < 1e-6 and painted["bonds"] and painted["lights"],
              f"selected haem: same colours, sizes, bonds and lights in both views ({painted})")
        bind = page.evaluate(f"""() => {{ const b = window.bgStore.S.buf, s = window.bgStore; const out = {{n: b.bind_type.length, dots: {V}.mol.bind ? {V}.mol.bind.pair.length : 0, fe: []}};
            for (let k = 0; k < b.bind_type.length; k++) if (b.bind_type[k] === 0) out.fe.push(s.label(b.atom_res[b.bind_pairs[2 * k + 1]]));
            return out; }}""")
        check(bind["dots"] > 0 and sorted(bind["fe"]) == ["A:HIS87", "B:HIS92", "C:HIS87", "D:HIS92"],
              f"dotted binding lines ({bind['n']} contacts, {bind['dots']} dots); Fe coordination: {bind['fe']}")

        print(f"[{browser}] 1. linked cameras")
        mb = page.locator("#mol-host").bounding_box()
        cx, cy = mb["x"] + mb["width"] / 2, mb["y"] + mb["height"] / 2
        page.mouse.move(cx, cy); page.mouse.down(); page.mouse.move(cx + 120, cy + 50, steps=10); page.mouse.up()
        page.mouse.move(cx, cy); page.mouse.wheel(0, -400)
        page.wait_for_timeout(1200)
        d = page.evaluate(f"""() => {{ const m = {V}.mol.camera, g = {V}.geo.camera;
            return [m.position.distanceTo(g.position), m.up.distanceTo(g.up), {V}.mol.controls.target.distanceTo({V}.geo.controls.target)]; }}""")
        check(max(d) < 1e-3, f"structure drag + zoom moved the geometry camera identically (Δ {max(d):.2e})")
        gb = page.locator("#geo-host").bounding_box()
        page.mouse.move(gb["x"] + gb["width"] / 2, gb["y"] + gb["height"] / 2); page.mouse.down()
        page.mouse.move(gb["x"] + gb["width"] / 2 - 90, gb["y"] + gb["height"] / 2 + 40, steps=8); page.mouse.up()
        page.wait_for_timeout(1200)
        d = page.evaluate(f"() => {V}.mol.camera.position.distanceTo({V}.geo.camera.position)")
        check(d < 1e-3, f"…and the other way round (Δ {d:.2e})")

        print(f"[{browser}] 2. gutter scrollbars never cover the table")
        geo = page.evaluate("""() => { const body = document.querySelector('#dock-body .vt-body'), y = document.querySelector('#dock-body .cs-track.cs-y'),
            x = document.querySelector('#dock-body .cs-track.cs-x');
            const b = body.getBoundingClientRect(), ty = y.getBoundingClientRect(), tx = x.getBoundingClientRect();
            return {right: b.right, ty: ty.left, bottom: b.bottom, tx: tx.top, idle: y.classList.contains('idle'),
                    native: body.offsetWidth - body.clientWidth, nativeH: body.offsetHeight - body.clientHeight}; }""")
        check(not geo["idle"] and geo["right"] <= geo["ty"] + 0.5 and geo["bottom"] <= geo["tx"] + 0.5,
              f"table viewport ends before both tracks ({geo['right']:.0f} ≤ {geo['ty']:.0f}, {geo['bottom']:.0f} ≤ {geo['tx']:.0f})")
        check(geo["native"] == 0 and geo["nativeH"] == 0, f"native scrollbars take no space ({geo['native']}, {geo['nativeH']})")
        before = page.evaluate("() => document.querySelector('#dock-body .vt-body').scrollTop")
        th = page.locator("#dock-body .cs-track.cs-y .cs-thumb").bounding_box()
        page.mouse.move(th["x"] + 3, th["y"] + 5); page.mouse.down(); page.mouse.move(th["x"] + 3, th["y"] + 60, steps=5); page.mouse.up()
        after = page.evaluate("() => document.querySelector('#dock-body .vt-body').scrollTop")
        check(after > before + 50, f"dragging the thumb scrolls the rows ({before:.0f} → {after:.0f})")
        # the bottom-most fully visible row, clicked 4 px left of the vertical track: exactly where an overlay
        # scrollbar would sit
        pt = page.evaluate("""() => { const body = document.querySelector('#dock-body .vt-body').getBoundingClientRect();
            const rows = [...document.querySelectorAll('#dock-body .vt-row')].map((e) => e.getBoundingClientRect()).filter((r) => r.bottom <= body.bottom + 0.5);
            const r = rows.reduce((a, q) => (q.bottom > a.bottom ? q : a)); return {x: body.right - 4, y: r.top + r.height / 2}; }""")
        page.mouse.click(pt["x"], pt["y"])
        check(page.evaluate("() => window.bgStore.sel.size") == 1, "the bottom row is clickable right next to the scrollbar (one residue selected)")

        print(f"[{browser}] 3. domain bar + HETATM selection → ⌬ measure")
        page.locator("#seq-chains .chain-row").first.click()       # show chain A in the sequence viewer
        page.keyboard.press("Escape")
        dom = page.evaluate(f"""() => {{ const q = {V}.seq.domRects.find((t) => t.d.source === 'Pfam'); const r = document.getElementById('seq-scroll').getBoundingClientRect();
            return q && {{x: r.left + Math.max(q.x0, 0) + 30, y: r.top + (q.y0 + q.y1) / 2, n: window.bgStore.domainResidues(q.d).length, name: q.d.name}}; }}""")
        check(dom is not None, "a Pfam domain bar is drawn under the sequence")
        page.mouse.click(dom["x"], dom["y"])
        page.wait_for_timeout(500)
        n = page.evaluate("() => window.bgStore.sel.size")
        check(n == dom["n"], f"clicking the {dom['name']} bar selects its {dom['n']} residues ({n})")
        page.locator("#seq-het .het-item", has_text="HEM (i)").locator("input").check()
        page.wait_for_timeout(400)
        n2 = page.evaluate("() => window.bgStore.sel.size")
        check(n2 == dom["n"] + 1, f"ticking HEM (i) adds the haem ({n2})")
        page.locator("#seq-measure").click()
        page.wait_for_function("window.bgStore.selGeo !== null", timeout=60_000)
        page.wait_for_timeout(800)
        res = page.evaluate("""() => { const s = window.bgStore.selGeo.stats; return {groups: s.groups.map((g) => g.label), pair: s.pairs[0],
            rows: document.querySelectorAll('#dock-body table.mini tbody tr').length, tab: document.querySelector('#dock-tabs .tab.on').textContent}; }""")
        check(res["groups"] == ["Pfam Globin · A", "HEM (i)"], f"groups: {res['groups']}")
        check(res["pair"]["min_dist"] < 2.4 and "HIS87" in res["pair"]["min_pair"][0], f"closest contact {res['pair']['min_pair']} {res['pair']['min_dist']:.2f} Å")
        check("Selection" in res["tab"] and res["rows"] >= 3, f"Selection tab shows groups + pairs ({res['rows']} rows)")
        sg = page.evaluate(f"() => {V}.geo.groups.selgeo.children.length")
        check(sg >= 2, f"selection geometry drawn in the geometry view ({sg} objects: group envelopes + contacts)")
        page.screenshot(path=str(out / f"{browser}_01_measure.png"))

        print(f"[{browser}] 4. simplices and hetero in the geometry view")
        st = page.evaluate(f"""() => {{ const g = {V}.geo, st = g.star; return st && {{ tets: st.nTets, nodes: st.nNodes, edges: st.edges.length / 2, inner: st.inner,
            moved: st.moved, extra: st.extra, kinds: g.groups.star.children.map((o) => o.type),
            mesh: (() => {{ const m = g.layerGroups.delaunay.children[0]; return {{ drawn: m.geometry.instanceCount, all: m.geometry.attributes.instanceStart.count,
                                                                                    w: m.material.linewidth, o: m.material.opacity }}; }})(),
            sel: g.groups.star.children.map((o) => ({{ w: o.material.linewidth, o: o.material.opacity, n: o.geometry.attributes.instanceStart.count }})) }}; }}""")
        check(st and st["tets"] > 100 and st["nodes"] > 100 and st["edges"] > st["inner"] > 0, f"tetrahedra touching the selection: {st and {k: st[k] for k in ('tets', 'nodes', 'edges', 'inner')}}")
        check(st and st["kinds"] == ["LineSegments2"] and st["moved"] + st["extra"] == st["edges"] == st["sel"][0]["n"]
              and st["mesh"]["drawn"] == st["mesh"]["all"] - st["moved"],
              f"their edges leave the mesh ({st and st['moved']} of {st and st['mesh']['all']}) and are drawn once, as edges only ({st and st['sel'][0]['n']})")
        check(st and abs(st["sel"][0]["w"] - st["mesh"]["w"] - 0.75) < 1e-6 and st["sel"][0]["o"] == st["mesh"]["o"],
              f"selected edges are 0.75 px wider than the mesh, same opacity ({st and st['mesh']['w']} → {st and st['sel'][0]['w']} px, opacity {st and st['mesh']['o']})")
        glow = page.evaluate(f"() => {V}.geo.groups.het.children.length + {V}.mol.groups.sel.children.filter((o) => o.material && o.material.blending === 2).length")
        check(glow == 0, f"no glow on the selection ({glow})")
        outl = page.evaluate(f"() => [{V}.mol.groups.sel.children.some((o) => o.isMesh && o.material.side === 1), {V}.mol.groups.sel.children.some((o) => o.isLineSegments2)]")
        check(outl[0] and not outl[1], f"the selected cartoon gets an outline, no Cα line across the ribbons ({outl})")
        het = page.evaluate(f"""() => {{ const g = {V}.geo, b = window.bgStore.S.buf; let hidden = 0, n = 0;
            for (let i = 0; i < b.pt_res.length; i++) if (b.res_kind[b.pt_res[i]] >= 2) {{ n++; if (g.pts.geometry.attributes.aRadius.array[i] === 0) hidden++; }}
            return [g.topGroup.children.length, g.top.children.includes(g.topGroup), hidden, n]; }}""")
        check(het[0] >= 2 and het[1] and het[2] == het[3], f"ligands / ions / waters drawn as ball-and-stick on top of the tessellation ({het})")
        txt = page.evaluate("() => document.getElementById('geo-stats').textContent")
        check("tetrahedra" in txt and "nodes" in txt, f"stats overlay: …{txt[-80:]}")
        off = page.evaluate("""() => { const L = window.bgStore.selGeo, P = L.buf.pts, R = L.buf.pt_res, A = window.bgStore.S.buf.res_pos; let d = 0;
            for (let i = 0; i < R.length; i++) d = Math.max(d, Math.hypot(P[3 * i] - A[3 * R[i]], P[3 * i + 1] - A[3 * R[i] + 1], P[3 * i + 2] - A[3 * R[i] + 2]));
            return d; }""")
        check(off < 12, f"the selection overlay sits on the selected residues (max atom–residue distance {off:.1f} Å)")

        print(f"[{browser}] 5. selection colour")
        page.locator("#sel-color").select_option("violet")
        page.wait_for_timeout(500)
        col = page.evaluate("() => [getComputedStyle(document.documentElement).getPropertyValue('--sel').trim(), window.bgStore.constructor && 1]")
        check(col[0].lower() == "#a78bfa", f"--sel follows the picker ({col[0]})")
        page.screenshot(path=str(out / f"{browser}_02_violet.png"))
        page.locator("#sel-color").select_option("magenta")

        print(f"[{browser}] 6. structure downloads")
        codes = page.evaluate(f"""async () => {{ const j = window.bgStore.job, out = [];
            for (const f of ['cif', 'pdb']) {{ const r = await fetch(`{API}/bg/jobs/${{j}}/structure/${{f}}`); out.push([r.status, (await r.text()).slice(0, 6)]); }} return out; }}""")
        check(codes[0][0] == 200 and codes[0][1].startswith("data_") and codes[1][0] == 200, f"mmCIF / PDB download: {codes}")
        opts = page.evaluate("() => [...document.querySelectorAll('#csv-download option')].map((o) => o.value).filter((v) => v.startsWith('structure:'))")
        check(opts == ["structure:cif", "structure:pdb"], f"dock download menu lists the structure files {opts}")
        print(f"[{browser}] 7. rendering, dock tabs, pockets, tetrahedra, coordination sphere, reset")
        px = page.evaluate(f"() => [{V}.geo.pts.geometry.attributes.aRadius.array[0], {V}.geo.layerGroups.delaunay.children[0].isLineSegments2 === true, {V}.mol.renderer.getPixelRatio()]")
        check(abs(px[0] - 4.2 * 0.85) < 0.05 or abs(px[0] - 4.2) < 0.05 or px[0] < 8, f"geometry points are small fixed-size dots ({px[0]:.1f} px)")
        check(px[1] and px[2] >= 1.5, f"thick anti-aliased Delaunay edges, supersampled canvas (pixel ratio {px[2]})")
        page.evaluate("() => window.bgWorkbench.handle({type: 'analysis', kind: 'pockets'})")
        page.wait_for_function("document.querySelector('#dock-tabs .tab.on').textContent.startsWith('Pockets')", timeout=60_000)
        page.locator("#dock-tabs .tab", has_text="Residues").click()
        page.wait_for_timeout(500)
        check(page.evaluate("() => document.querySelectorAll('#dock-body .vt-row').length") > 0, "the Residues table comes back after an analysis tab")
        page.locator("#dock-tabs .tab", has_text="Pockets").click()
        page.wait_for_timeout(400)
        page.locator("#dock-body table.mini tbody tr").first.click(); page.wait_for_timeout(300)
        a = page.evaluate("() => window.bgStore.sel.size")
        page.locator("#dock-body table.mini tbody tr").first.click(); page.wait_for_timeout(300)
        b2 = page.evaluate(f"() => [window.bgStore.sel.size, {V}.geo.pocketSel]")
        check(a > 0 and b2 == [0, None], f"a pocket click selects ({a}) and a second click deselects it ({b2})")
        page.locator("#dock-body .chart-tools .chip").first.click(); page.wait_for_timeout(200)
        check(page.evaluate(f"() => {V}.geo.groups.overlay.visible") is False, "the pocket surfaces can be hidden")
        page.evaluate("() => window.bgStore.setSelection([30, 31, 32, 33], 't')")
        page.locator("#dock-tabs .tab", has_text="Tetrahedra").click()
        page.wait_for_function("() => window.bgStore.T !== null && document.querySelectorAll('#dock-body .vt-row').length > 0", timeout=120_000)
        page.wait_for_timeout(500)
        tt = page.evaluate(f"""() => {{ const t = {V}.dock.tables.tetrahedra, n = t.nInvolved;
            return {{ n, first: t.spec.involvesSelection(t.order[0]), boundary: n < t.order.length ? t.spec.involvesSelection(t.order[n]) : false,
                     heads: [...document.querySelectorAll('#dock-body .vt-hcell')].map((e) => e.textContent) }}; }}""")
        check(tt["n"] > 0 and tt["first"] and not tt["boundary"], f"tetrahedra with the selection come first ({tt['n']} of them)")
        want = ["residues", "residue-class", "tropsha-class", "tetrahedrality", "interaction-type", "voronoi-area-1", "residue-1-depth"]
        check(all(any(h.startswith(w) for h in tt["heads"]) for w in want), f"reference columns present ({[w for w in want if not any(h.startswith(w) for h in tt['heads'])] or 'all'})")
        page.locator("#seq-het .het-item", has_text="HEM (i)").locator(".coord").click()
        page.wait_for_function("window.bgStore.selGeo !== null && window.bgStore.selGeo.stats.groups.some((g) => g.label === 'HEM (i)')", timeout=60_000)
        sph = page.evaluate("() => [...window.bgStore.sel].map((r) => window.bgStore.label(r))")
        check("A:HIS87" in sph and any("HOH" in x for x in sph), f"⌬ selects the haem with its coordinating residues and waters ({len(sph)} residues)")
        page.evaluate("() => window.bgWorkbench.reset()")
        page.wait_for_function("window.bgStore.sel.size === 0 && document.querySelector('#dock-tabs .tab.on').textContent.startsWith('Residues')", timeout=60_000)
        rs = page.evaluate(f"() => [window.bgStore.selGeo, [...{V}.geo.active], {V}.geo.overlay, document.querySelectorAll('#dock-body .vt-row').length > 0]")
        check(rs[0] is None and rs[1] == ["delaunay"] and rs[2] is None and rs[3], f"reset gives a clean slate {rs}")
        popups = []
        page.on("popup", lambda p_: popups.append(p_))
        with page.expect_download(timeout=60_000) as dl:
            page.locator("#csv-download").select_option("tetrahedra")
        check(dl.value.suggested_filename.endswith("_delaunay-tetrahedra.csv") and not popups,
              f"download without a new tab: {dl.value.suggested_filename} ({len(popups)} popups)")
        hdr = page.evaluate(f"""async () => {{ const r = await fetch(`{API}/bg/jobs/${{window.bgStore.job}}/csv/tetrahedra`);
            return [r.headers.get('content-disposition'), (await r.text()).split('\\n')[0].slice(0, 60)]; }}""")
        check("_delaunay-tetrahedra.csv" in hdr[0] and hdr[1].startswith("pdbid-chainid,model,ec-number"), f"named CSV: {hdr[0].split('filename=')[-1]}")
        print(f"[{browser}] 8. multi-select, clear, unchanged shade, node size, depth code, measure")
        page.keyboard.press("Escape")
        box = page.locator("#seq-scroll").bounding_box()
        X = lambda i: box["x"] + 8 + 15 * i + 6  # noqa: E731
        yy = box["y"] + 14 + 6
        page.locator("#seq-multi").click()
        page.mouse.click(X(3), yy); page.mouse.click(X(7), yy); page.mouse.click(X(11), yy); page.wait_for_timeout(200)
        n3 = page.evaluate("() => window.bgStore.sel.size")
        page.mouse.click(X(7), yy); page.wait_for_timeout(200)
        n2 = page.evaluate("() => window.bgStore.sel.size")
        page.locator("#seq-multi").click()
        check(n3 == 3 and n2 == 2, f"multi-select: three clicks select three residues ({n3}); clicking one again removes it ({n2})")
        col0 = page.evaluate(f"() => Array.from({V}.mol.cartoon.mesh.geometry.attributes.color.array.slice(0, 3))")
        page.evaluate("() => window.bgStore.setSelection([60, 61, 62, 63], 't')"); page.wait_for_timeout(300)
        col1 = page.evaluate(f"() => Array.from({V}.mol.cartoon.mesh.geometry.attributes.color.array.slice(0, 3))")
        check(max(abs(a - b) for a, b in zip(col0, col1)) < 1e-6, "a selection leaves the colours of everything else unchanged (dim off)")
        px = page.evaluate(f"() => {V}.geo.basePx()")
        check(abs(px - 6.2) < 1e-6, f"Delaunay nodes are 2 px larger ({px} px all-atom)")
        page.keyboard.press("g")
        page.wait_for_function("window.bgStore.selGeo !== null", timeout=60_000)
        part = page.evaluate("() => window.bgStore.selGeo.stats.summary.participating")
        star = page.evaluate(f"() => {V}.geo.star.nTets")
        check(part["n_tets"] == star > 0, f"measure counts the same tetrahedra the view highlights ({part['n_tets']} = {star})")
        page.locator("#seq-clear").click(); page.wait_for_timeout(300)
        cl = page.evaluate(f"() => [window.bgStore.sel.size, window.bgStore.selGeo, {V}.geo.star, {V}.geo.groups.selgeo.children.length]")
        check(cl == [0, None, None, 0], f"clear empties the selection, its measurement and highlights ({cl})")
        mesh = page.evaluate(f"() => {{ const m = {V}.geo.layerGroups.delaunay.children[0]; return [m.geometry.instanceCount, m.geometry.attributes.instanceStart.count, {V}.geo.groups.star.children.length]; }}")
        check(mesh[0] == mesh[1] and mesh[2] == 0, f"after clear the mesh draws every edge again, no selection overlay ({mesh})")
        page.evaluate("() => window.bgStore.setSelection([86], 't')")
        page.locator("#dock-tabs .tab", has_text="Tetrahedra").click()
        page.wait_for_function("() => window.bgStore.T !== null", timeout=120_000)
        page.wait_for_timeout(500)
        dep = page.evaluate("""() => { const T = window.bgStore.T, c = T.stats.columns.find((x) => x.key === 'depth');
            return c ? c.categories.slice(0, 5) : null; }""")
        check(dep and all(len(x) == 4 and set(x) <= set("SUB") for x in dep), f"depth column with four-letter codes ({dep})")
        check(not errors, f"no console errors ({errors[:3]})")
        b.close()
    print("PASS" if not fails else f"{len(fails)} FAILED")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main(*sys.argv[1:])
