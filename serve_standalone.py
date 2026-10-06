"""Run the bioGeometry API + workbench without Reflex (testing, embedding, headless screenshots).

    python serve_standalone.py [--port 8770]
    open http://127.0.0.1:8770/workbench/index.html?api=http://127.0.0.1:8770&pdb=4HHB
         ...&af=P69905              (AlphaFold model)
         ...&layers=alpha,voronoi   (geometry layers to switch on)  &analysis=pockets  &tab=tetrahedra
"""
import argparse
from pathlib import Path

import uvicorn
from starlette.routing import Mount
from starlette.staticfiles import StaticFiles

from biogeometry.api import api

ROOT = Path(__file__).resolve().parent
api.router.routes.append(Mount("/workbench", app=StaticFiles(directory=ROOT / "assets" / "workbench", html=True)))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8770)
    a = ap.parse_args()
    print(f"workbench: http://{a.host}:{a.port}/workbench/index.html?api=http://{a.host}:{a.port}&pdb=4HHB")
    uvicorn.run(api, host=a.host, port=a.port, log_level="warning")
