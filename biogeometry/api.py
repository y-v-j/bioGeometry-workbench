"""HTTP API for the workbench (Starlette). Mounted into Reflex via ``api_transformer`` and also
runnable standalone (``python serve_standalone.py``).

GET  /bg/health
POST /bg/jobs                                {spec}          -> summary + structure layer meta
POST /bg/upload  (multipart: file + spec fields)              -> same
GET  /bg/local?dir=...                                       -> structure files in a local directory
POST /bg/jobs/{job}/layers                   {kind, params}  -> layer meta (computes if needed)
GET  /bg/jobs/{job}/layers/{key}/meta                        -> layer meta
GET  /bg/jobs/{job}/layers/{key}/buf/{name}                  -> raw little-endian typed array
GET  /bg/jobs/{job}/csv/{kind}                               -> CSV download (see engine.CSV_KINDS)
GET  /bg/csv-kinds                                           -> available CSV downloads
"""
from __future__ import annotations

import traceback

from starlette.applications import Starlette
from starlette.concurrency import run_in_threadpool
from starlette.middleware import Middleware
from starlette.middleware.gzip import GZipMiddleware
from starlette.requests import Request
from starlette.responses import FileResponse, JSONResponse, Response
from starlette.routing import Route

from . import engine
from .geometry import GeometryError
from .resources import MemoryBudgetError, available_mb
from .structure import StructureError, list_local_structures, store_upload

USER_ERRORS = (engine.EngineError, StructureError, GeometryError, FileNotFoundError, MemoryBudgetError)


def _err(exc: Exception) -> JSONResponse:
    if isinstance(exc, USER_ERRORS):
        return JSONResponse({"error": str(exc)}, status_code=400)
    traceback.print_exc()
    return JSONResponse({"error": f"{type(exc).__name__}: {exc}"}, status_code=500)


def _spec(d: dict) -> engine.JobSpec:
    fields = engine.JobSpec.__dataclass_fields__
    kw = {k: d[k] for k in fields if k in d and d[k] is not None}
    if "ligands" in kw:
        kw["ligands"] = str(kw["ligands"]).lower() in ("1", "true", "yes", "on")
    return engine.JobSpec(**kw)


def _open(spec: engine.JobSpec) -> dict:
    job = engine.open_job(spec)
    meta = engine.build_layer(job, "structure")
    return {"summary": job.summary(), "structure": meta}


def _attachment(name: str) -> dict:
    # the workbench and panel fetch downloads from another port: expose the file name to them (CORS)
    return {"Content-Disposition": f'attachment; filename="{name}"', "Access-Control-Expose-Headers": "Content-Disposition"}


async def health(_: Request):
    mem = available_mb()
    return JSONResponse({"ok": True, "app": "bioGeometry", "mem_available_mb": None if mem is None else round(mem)})


async def create_job(request: Request):
    try:
        spec = _spec(await request.json())
        return JSONResponse(await run_in_threadpool(_open, spec))
    except Exception as exc:  # noqa: BLE001
        return _err(exc)


async def upload(request: Request):
    try:
        form = await request.form()
        f = form.get("file")
        if f is None:
            raise engine.EngineError("no file in upload")
        data = await f.read()
        extra = {k: v for k, v in form.items() if k != "file"}

        def work():
            path = store_upload(data, f.filename)
            return _open(_spec({**extra, "source": "upload", "path": str(path)}))
        return JSONResponse(await run_in_threadpool(work))
    except Exception as exc:  # noqa: BLE001
        return _err(exc)


async def local_files(request: Request):
    try:
        return JSONResponse({"files": list_local_structures(request.query_params.get("dir", "."))})
    except Exception as exc:  # noqa: BLE001
        return _err(exc)


async def make_layer(request: Request):
    try:
        job = engine.get_job(request.path_params["job"])
        if job is None:
            return JSONResponse({"error": "session expired - load the structure again"}, status_code=404)
        body = await request.json()
        meta = await run_in_threadpool(engine.build_layer, job, body.get("kind", ""), body.get("params") or {})
        return JSONResponse(meta)
    except Exception as exc:  # noqa: BLE001
        return _err(exc)


async def layer_meta(request: Request):
    meta = engine.read_meta(request.path_params["job"], request.path_params["key"])
    if meta is None:
        return JSONResponse({"error": "layer not found"}, status_code=404)
    return JSONResponse(meta)


async def layer_buffer(request: Request):
    pp = request.path_params
    path = engine.buffer_path(pp["job"], pp["key"], pp["name"])
    if path is None:
        return JSONResponse({"error": "buffer not found"}, status_code=404)
    return FileResponse(path, media_type="application/octet-stream",
                        headers={"Cache-Control": "public, max-age=86400, immutable"})


async def csv_download(request: Request):
    try:
        job = engine.get_job(request.path_params["job"])
        if job is None:
            return JSONResponse({"error": "session expired - load the structure again"}, status_code=404)
        name, text = await run_in_threadpool(engine.export_csv, job, request.path_params["kind"])
        return Response(text, media_type="text/csv", headers=_attachment(name))
    except Exception as exc:  # noqa: BLE001
        return _err(exc)


async def structure_download(request: Request):
    try:
        job = engine.get_job(request.path_params["job"])
        if job is None:
            return JSONResponse({"error": "session expired - load the structure again"}, status_code=404)
        name, text = await run_in_threadpool(engine.export_structure, job, request.path_params["fmt"])
        media = "chemical/x-mmcif" if name.endswith(".cif") else "chemical/x-pdb"
        return Response(text, media_type=media, headers=_attachment(name))
    except Exception as exc:  # noqa: BLE001
        return _err(exc)


async def csv_kinds(_: Request):
    return JSONResponse(engine.CSV_KINDS)


routes = [
    Route("/bg/health", health),
    Route("/bg/jobs", create_job, methods=["POST"]),
    Route("/bg/upload", upload, methods=["POST"]),
    Route("/bg/local", local_files),
    Route("/bg/csv-kinds", csv_kinds),
    Route("/bg/jobs/{job}/layers", make_layer, methods=["POST"]),
    Route("/bg/jobs/{job}/layers/{key}/meta", layer_meta),
    Route("/bg/jobs/{job}/layers/{key}/buf/{name}", layer_buffer),
    Route("/bg/jobs/{job}/csv/{kind}", csv_download),
    Route("/bg/jobs/{job}/structure/{fmt}", structure_download),
]

api = Starlette(routes=routes, middleware=[Middleware(GZipMiddleware, minimum_size=2048)])
