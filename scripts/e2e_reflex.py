"""End-to-end test through the Reflex UI (needs `reflex run` on :3100).

    python scripts/e2e_reflex.py OUT_DIR

Loads 4HHB from the panel, checks the workbench iframe receives it, runs the pocket analysis and a
parameter change, then loads the AlphaFold example and checks the PAE tab.
"""
import sys
import time
from pathlib import Path

from playwright.sync_api import sync_playwright

URL = "http://localhost:3100"
# newest Chromium that `playwright install chromium` put in the cache
CHROME = str(max((Path.home() / ".cache/ms-playwright").glob("chromium-*/chrome-linux64/chrome"), default=""))


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    fails = []

    def check(cond, msg):
        print(("  ok   " if cond else "  FAIL ") + msg)
        if not cond:
            fails.append(msg)

    with sync_playwright() as p:
        b = p.chromium.launch(executable_path=CHROME, args=["--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"])
        page = b.new_page(viewport={"width": 1920, "height": 1080})
        errors = []
        page.on("pageerror", lambda e: errors.append(str(e)))
        page.goto(URL, wait_until="networkidle", timeout=120_000)
        frame = page.frame_locator("#bg-workbench")
        wb = lambda: next(f for f in page.frames if f.url.endswith("/workbench/index.html"))  # noqa: E731
        page.wait_for_function("() => { const f = document.getElementById('bg-workbench'); return f && f.contentWindow && f.contentWindow.bgWorkbench; }", timeout=60_000)

        print("1. load 4HHB from the panel")
        t0 = time.time()
        page.get_by_role("button", name="Load structure").click()
        wb().wait_for_function("() => window.bgStore && window.bgStore.S && window.bgStore.S.stats.name === '4HHB'", timeout=180_000)
        wb().wait_for_function("() => window.bgStore.F !== null", timeout=180_000)
        check(True, f"workbench received 4HHB + features ({time.time() - t0:.1f}s)")
        page.wait_for_timeout(1500)
        page.screenshot(path=str(out / "01_reflex_4hhb.png"))

        print("1b. the domain / ligand / ion / water panel sits in the sidebar")
        items = page.locator("#bg-selection-panel .het-item").all_inner_texts()
        hidden = wb().evaluate("() => getComputedStyle(document.getElementById('seq-het-host')).display === 'none' && !document.getElementById('seq-het').children.length")
        check(any("HEM" in t for t in items) and hidden, f"sidebar lists {len(items)} groups; the sequence viewer's own panel is gone ({hidden})")
        box = lambda loc: loc.first.bounding_box()  # noqa: E731
        panel, frame_box = box(page.locator("#bg-selection-scroll")), box(page.locator("#bg-workbench"))
        order = [box(page.get_by_role("button", name="Load structure"))["y"], panel["y"], box(page.get_by_text("GEOMETRY PARAMETERS", exact=True))["y"]]
        check(panel["x"] + panel["width"] <= frame_box["x"] and order[0] < order[1] < order[2],
              f"in the sidebar, below Load structure and above Geometry parameters (tops {[round(v) for v in order]})")
        page.locator("#bg-selection-panel .het-item", has_text="HEM").first.locator("input").check()
        page.wait_for_timeout(500)
        n = wb().evaluate("() => window.bgStore.sel.size")
        on = page.locator("#bg-selection-panel .het-item.on").count()
        check(n > 0 and on > 0, f"ticking HEM in the sidebar selects it in the workbench ({n} residues, {on} rows marked)")
        page.screenshot(path=str(out / "01b_reflex_sidebar_selection.png"))
        wb().evaluate("() => window.bgStore.clearAll('test')")
        page.wait_for_timeout(300)
        check(page.locator("#bg-selection-panel .het-item.on").count() == 0, "clearing in the workbench unticks the sidebar")

        print("2. run pockets from the Analyses section")
        page.get_by_text("ANALYSES", exact=True).click()
        page.wait_for_timeout(400)
        page.get_by_role("button", name="Pockets/cavities").click()
        wb().wait_for_function("() => document.querySelector('#dock-tabs .tab.on')?.textContent.startsWith('Pockets')", timeout=120_000)
        n = wb().evaluate("() => document.querySelectorAll('#dock-body table.mini tbody tr').length")
        check(n > 0, f"pockets tab shows {n} pockets")
        wb().locator("#dock-body table.mini tbody tr").first.click()
        page.wait_for_timeout(800)
        sel = wb().evaluate("() => window.bgStore.sel.size")
        check(sel > 0, f"clicking the top pocket selects its {sel} lining residues")
        page.screenshot(path=str(out / "02_reflex_pockets.png"))

        print("3. α-shape on, then change α in the panel")
        wb().locator("#geo-toggles [data-layer='alpha']").click()
        wb().wait_for_function("() => document.querySelectorAll('#geo-toggles .busy').length === 0", timeout=60_000)
        page.wait_for_timeout(800)
        before = wb().evaluate("() => document.getElementById('geo-stats').textContent")
        page.get_by_text("GEOMETRY PARAMETERS", exact=True).click() if not page.get_by_text("α radius (Å)").is_visible() else None
        # sliders are hard to drag headlessly, so push the same message the panel's slider sends
        wb().evaluate("() => window.bgWorkbench.handle({type:'params', params:{alpha:{alpha: 6.0}}})")
        wb().wait_for_function("() => document.getElementById('geo-stats').textContent.includes('α 6')", timeout=60_000)
        after = wb().evaluate("() => document.getElementById('geo-stats').textContent")
        check("α 6" in after and before != after, "α-shape recomputed for the new parameter")

        print("4. AlphaFold example")
        page.get_by_text("P69905", exact=True).click()
        wb().wait_for_function("() => window.bgStore.S && window.bgStore.S.stats.is_alphafold", timeout=180_000)
        wb().wait_for_function("() => window.bgStore.F !== null", timeout=180_000)
        wb().locator("#tab-pae").click()
        wb().wait_for_function("() => document.querySelector('#rama-svg').textContent.includes('expected position error')", timeout=60_000)
        check(True, "PAE heatmap rendered")
        wb().locator("#mol-color").select_option("conf")
        page.wait_for_timeout(1200)
        page.screenshot(path=str(out / "03_reflex_alphafold_pae.png"))
        print("5. panel: tetrahedron logo, structure downloads, geometry of the selection")
        circles = page.evaluate("() => document.querySelector('svg circle') ? document.querySelector('svg circle').closest('svg').querySelectorAll('circle').length : 0")
        check(circles == 4, f"the brand logo is a tetrahedron with 4 vertex atoms ({circles})")
        page.get_by_text("DOWNLOADS", exact=True).click()
        page.wait_for_timeout(500)
        popups = []
        page.on("popup", lambda p_: popups.append(p_))
        with page.expect_download(timeout=60_000) as dl:
            page.get_by_text("structure · PDB (.pdb)").click()
        name = dl.value.suggested_filename
        check(name.endswith("_model.pdb") and not popups, f"structure download without a new tab: {name} ({len(popups)} popups)")
        with page.expect_download(timeout=60_000) as dl2:
            page.get_by_text("residues", exact=True).click()
        check(dl2.value.suggested_filename.endswith("_residue-features.csv") and not popups, f"CSV download: {dl2.value.suggested_filename}")
        page.screenshot(path=str(out / "04_reflex_downloads.png"))
        wb().evaluate("() => { const s = window.bgStore; const d = s.S.stats.domains[0]; s.setSelection(s.domainResidues(d), 'test'); }")
        page.get_by_text("ANALYSES", exact=True).click()
        page.wait_for_timeout(400)
        page.get_by_role("button", name="Measure selection").click()
        wb().wait_for_function("() => window.bgStore.selGeo !== null", timeout=60_000)
        tab = wb().evaluate("() => document.querySelector('#dock-tabs .tab.on').textContent")
        check("Selection" in tab, f"the panel button measures the workbench selection ({tab})")
        page.wait_for_timeout(800)
        page.screenshot(path=str(out / "05_reflex_selection.png"))
        print("6. All-atom / Backbone choice before loading")
        page.get_by_text("PDB", exact=True).first.click()          # back to 4HHB (the AlphaFold model has no hetero atoms)
        page.get_by_text("Backbone", exact=True).first.click()
        page.get_by_role("button", name="Load structure").click()
        wb().wait_for_function("() => window.bgStore.S && window.bgStore.S.stats.level === 'ca'", timeout=180_000)
        het_pts = wb().evaluate("() => { const b = window.bgStore.S.buf; let n = 0; for (const r of b.pt_res) if (b.res_kind[r] >= 2) n++; return [n, b.pt_res.length]; }")
        check(het_pts[0] > 0, f"Backbone mode keeps ligand / ion / water atoms as points ({het_pts[0]} of {het_pts[1]})")
        page.get_by_text("All-atom", exact=True).first.click()
        print("7. the reset button gives a clean slate")
        wb().evaluate("() => { window.bgStore.setSelection([10, 11, 12], 't'); window.bgWorkbench.handle({type: 'analysis', kind: 'pockets'}); }")
        page.wait_for_timeout(2500)
        page.locator("button:has(svg.lucide-rotate-ccw)").first.click()
        wb().wait_for_function("() => window.bgStore.sel.size === 0 && document.querySelector('#dock-tabs .tab.on').textContent.startsWith('Residues')", timeout=60_000)
        clean = wb().evaluate("() => [window.bgWorkbench.views.geo.overlay, [...window.bgWorkbench.views.geo.active], document.querySelectorAll('#dock-body .vt-row').length > 0]")
        check(clean[0] is None and clean[1] == ["delaunay"] and clean[2], f"selection, overlays and layers reset; table shown ({clean})")
        toast = wb().evaluate("() => { const t = document.getElementById('toast'); return t && !t.hidden && !t.classList.contains('info') ? t.textContent : ''; }")
        check(not toast, f"no error toast after switching structure ({toast!r})")
        check(not errors, f"no page errors ({errors[:2]})")
        b.close()
    print("PASS" if not fails else f"{len(fails)} FAILED")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main(*sys.argv[1:])
