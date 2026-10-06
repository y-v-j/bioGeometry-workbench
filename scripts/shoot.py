"""Headless screenshots + smoke test of the standalone workbench (Playwright, SwiftShader WebGL).

    python scripts/shoot.py OUT_DIR "pdb=4HHB" "af=P69905&layers=alpha" ...

Prints console errors and the status bar for every URL.
"""
import sys
import time
from pathlib import Path

from playwright.sync_api import sync_playwright

API = "http://127.0.0.1:8770"
# newest Chromium that `playwright install chromium` put in the cache
CHROME = str(max((Path.home() / ".cache/ms-playwright").glob("chromium-*/chrome-linux64/chrome"), default=""))


def main(out: str, queries: list[str], actions: str = ""):
    out_dir = Path(out)
    out_dir.mkdir(parents=True, exist_ok=True)
    failures = 0
    with sync_playwright() as p:
        b = p.chromium.launch(executable_path=CHROME, args=["--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"])
        for i, q in enumerate(queries):
            page = b.new_page(viewport={"width": 1680, "height": 1000})
            errors = []
            page.on("console", lambda m: errors.append(m.text) if m.type == "error" else None)
            page.on("pageerror", lambda e: errors.append(str(e)))
            t0 = time.time()
            page.goto(f"{API}/workbench/index.html?api={API}&{q}")
            page.wait_for_function("window.bgStandaloneDone !== undefined", timeout=300_000)
            try:
                page.wait_for_function("window.bgStore && window.bgStore.F !== null", timeout=180_000)
            except Exception:
                errors.append("features layer did not arrive")
            page.wait_for_timeout(2500)
            status = page.evaluate("window.bgStandaloneDone")
            bar = page.inner_text("#status")
            name = f"{i:02d}_" + "".join(c if c.isalnum() else "_" for c in q)[:60] + ".png"
            page.screenshot(path=str(out_dir / name))
            ok = status is True and not errors
            failures += not ok
            print(f"[{'OK' if ok else 'FAIL'}] {q} ({time.time() - t0:.1f}s) -> {name}")
            print("   ", bar.replace("\n", " | "))
            for e in errors[:8]:
                print("    console:", e[:300])
            page.close()
        b.close()
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
