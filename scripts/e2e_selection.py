"""End-to-end check of the coupled selection in the standalone workbench.

    python scripts/e2e_selection.py OUT_DIR [pdb=4HHB]

Drives the UI like a user (sequence drag, Ramachandran brush, ligand site, table filter → select,
3D click) and asserts that every view reflects the same selection. Needs serve_standalone.py on :8770.
"""
import sys
import time
from pathlib import Path

from playwright.sync_api import sync_playwright

API = "http://127.0.0.1:8770"
# newest Chromium that `playwright install chromium` put in the cache
CHROME = str(max((Path.home() / ".cache/ms-playwright").glob("chromium-*/chrome-linux64/chrome"), default=""))
STATE = """() => {
  const s = window.bgStore;
  const rows = document.querySelectorAll('#dock-body .vt-row.sel').length;
  const labels = document.querySelectorAll('#rama-svg .sel-label').length;
  const het = [...document.querySelectorAll('#seq-het .het-item.on')].map(e => e.textContent.trim());
  return {n: s.sel.size, rows, labels, het, status: document.getElementById('st-sel').textContent,
          chains: [...document.querySelectorAll('#seq-chains .chain-row.sel .cid')].map(e => e.textContent)};
}"""


def main(out, query="pdb=4HHB"):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    fails = []

    def check(cond, msg):
        print(("  ok   " if cond else "  FAIL ") + msg)
        if not cond:
            fails.append(msg)

    with sync_playwright() as p:
        b = p.chromium.launch(executable_path=CHROME, args=["--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"])
        page = b.new_page(viewport={"width": 1680, "height": 1000})
        errors = []
        page.on("console", lambda m: errors.append(m.text) if m.type == "error" else None)
        page.on("pageerror", lambda e: errors.append(str(e)))
        page.goto(f"{API}/workbench/index.html?api={API}&{query}")
        page.wait_for_function("window.bgStandaloneDone === true", timeout=300_000)
        page.wait_for_function("window.bgStore.F !== null", timeout=300_000)
        page.wait_for_timeout(1200)

        print("1. drag across 8 residues in the sequence viewer")
        box = page.locator("#seq-scroll").bounding_box()
        y = box["y"] + 14 + 6
        page.mouse.move(box["x"] + 8 + 15 * 20 + 6, y)
        page.mouse.down()
        page.mouse.move(box["x"] + 8 + 15 * 27 + 6, y, steps=6)
        page.mouse.up()
        page.wait_for_timeout(700)
        st = page.evaluate(STATE)
        check(st["n"] == 8, f"selection has 8 residues ({st['n']})")
        check(st["rows"] >= 8, f"table highlights them ({st['rows']} rows)")
        check(st["labels"] >= 6, f"Ramachandran labels them ({st['labels']})")
        check("8 residues" in st["status"], f"status bar: {st['status'][:70]}")
        page.screenshot(path=str(out / "01_sequence_drag.png"))

        print("2. brush the α-helix basin in the Ramachandran plot")
        rb = page.locator("#rama-host").bounding_box()
        # α region ≈ φ -100..-30, ψ -70..-10 in data space; convert via the plot box
        box_px = page.evaluate("""() => { const svg = document.querySelector('#rama-svg rect[fill="none"]');
                                       const b = svg.getBBox(); return [b.x, b.y, b.width]; }""")
        bx, by, s = box_px
        X = lambda phi: rb["x"] + bx + (phi + 180) / 360 * s
        Y = lambda psi: rb["y"] + by + (180 - psi) / 360 * s
        page.mouse.move(X(-100), Y(-10)); page.mouse.down(); page.mouse.move(X(-30), Y(-70), steps=8); page.mouse.up()
        page.wait_for_timeout(900)
        st = page.evaluate(STATE)
        check(st["n"] > 100, f"helix residues selected from the plot ({st['n']})")
        n_helix = page.evaluate("() => { const s = window.bgStore, b = s.S.buf; let h = 0; for (const r of s.sel) if (b.res_ss[r] === 1) h++; return h; }")
        check(n_helix / max(1, st["n"]) > 0.8, f"…and they are mostly helical by DSSP ({n_helix}/{st['n']})")
        page.screenshot(path=str(out / "02_rama_brush.png"))

        print("3. ligand binding site from the HETATM list")
        site = page.locator("#seq-het .het-item .site").first
        if site.count():
            site.click()
            page.wait_for_timeout(800)
            st = page.evaluate(STATE)
            check(st["n"] > 5 and st["het"], f"site selected: {st['n']} residues, het on: {st['het'][:2]}")
            page.screenshot(path=str(out / "03_ligand_site.png"))
        else:
            print("  (no ligands with a site)")

        print("4. table filter RSA < 0.05 and hydrophobic → select filtered")
        page.locator("#dock-tabs .tab", has_text="Residues").click()
        page.wait_for_timeout(300)
        rsa_idx = page.evaluate("() => [...document.querySelectorAll('#dock-body .vt-hcell')].findIndex(e => e.textContent.startsWith('RSA'))")
        page.locator("#dock-body .vt-fcell").nth(rsa_idx).locator("input").fill("<0.05")
        page.wait_for_timeout(600)
        page.locator("[data-act='selfiltered']").click()
        page.wait_for_timeout(800)
        st = page.evaluate(STATE)
        core = page.evaluate("() => { const s = window.bgStore; let ok = 0; for (const r of s.sel) if (s.featureValue('rsa', r) < 0.05) ok++; return [ok, s.sel.size]; }")
        check(core[0] == core[1] and core[1] > 0, f"selection = buried core residues ({core[0]}/{core[1]})")
        page.screenshot(path=str(out / "04_table_filter_select.png"))
        page.locator("[data-act='clearfilters']").click()

        print("5. click a residue in the 3D structure view")
        page.keyboard.press("Escape")
        mb = page.locator("#mol-host").bounding_box()
        hit = None
        for dx in range(-60, 61, 20):
            for dy in range(-60, 61, 20):
                page.mouse.click(mb["x"] + mb["width"] / 2 + dx, mb["y"] + mb["height"] / 2 + dy)
                page.wait_for_timeout(150)
                st = page.evaluate(STATE)
                if st["n"] == 1:
                    hit = st
                    break
            if hit:
                break
        check(hit is not None, "a 3D click selects exactly one residue")
        if hit:
            check(hit["rows"] == 1, f"table follows the 3D pick ({hit['rows']} row)")
        page.screenshot(path=str(out / "05_3d_pick.png"))

        print("6. geometry layers: α-shape + Voronoi + hull, tetrahedra tab")
        for lay in ("alpha", "voronoi", "hull"):
            page.locator(f"#geo-toggles [data-layer='{lay}']").click()
        page.wait_for_function("window.bgStore && document.querySelectorAll('#geo-toggles .busy').length === 0", timeout=120_000)
        page.wait_for_timeout(1200)
        page.locator("#dock-tabs .tab", has_text="Tetrahedra").click()
        # the residue table stays in the DOM until the tetrahedra arrive: wait for the tetrahedra header itself
        page.wait_for_function("""() => window.bgStore.T !== null
            && [...document.querySelectorAll('#dock-body .vt-hcell')].some(e => e.textContent.startsWith('simplex-type'))
            && document.querySelectorAll('#dock-body .vt-row').length > 0""", timeout=120_000)
        st = page.evaluate(STATE)
        check(st["rows"] >= 1, f"tetrahedra containing the picked residue are highlighted ({st['rows']})")
        page.screenshot(path=str(out / "06_layers_tetrahedra.png"))
        check(not errors, f"no console errors ({errors[:3]})")
        b.close()
    print("PASS" if not fails else f"{len(fails)} FAILED")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main(*sys.argv[1:])
