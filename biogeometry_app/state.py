"""Reflex state: the control plane of the bioGeometry workbench.

Structure loading and analyses run in the engine (worker thread); results are binary layers the
workbench iframe fetches from /bg/... directly. The state only sends small JSON messages to the iframe.
"""
from __future__ import annotations

import asyncio
import json
from pathlib import Path

import reflex as rx

from biogeometry import engine
from biogeometry.structure import STRUCTURE_SUFFIXES, list_local_structures, store_upload

ANALYSES = ["pockets", "contacts", "filtration", "fourbody", "scaling"]
# (source, id, description, recommended mode): the large assemblies default to backbone
EXAMPLES = [
    ("rcsb", "4HHB", "haemoglobin · heme sites", "allatom"),
    ("rcsb", "1CA2", "carbonic anhydrase · Zn site", "allatom"),
    ("rcsb", "1AON", "GroEL–GroES · 21 chains (backbone)", "ca"),
    ("rcsb", "7K00", "70S ribosome · RNA + protein (backbone)", "ca"),
    ("alphafold", "P69905", "AlphaFold · HBA1", "allatom"),
    ("alphafold", "P04637", "AlphaFold · p53 (disorder)", "allatom"),
]
# the two analysis modes: polymer level + ligands, ions and coordinated waters always in the point set
MODES = {"allatom": "heavy", "ca": "ca"}


def _api_base() -> str:
    return rx.config.get_config().api_url.rstrip("/")


# fetch the file and save it through a temporary link with the server's file name
DOWNLOAD_JS = """(async (url) => {
  const r = await fetch(url);
  if (!r.ok) { alert('download failed: HTTP ' + r.status); return; }
  const cd = r.headers.get('content-disposition') || '', m = cd.match(/filename="?([^";]+)"?/);
  const a = document.createElement('a');
  a.href = URL.createObjectURL(await r.blob()); a.download = m ? m[1] : url.split('/').pop();
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 4000);
})(%s)"""


def send_to_workbench(msg: dict):
    """Deliver a message to the workbench iframe, waiting (up to 20 s) until it has booted."""
    payload = json.dumps({"target": "bgworkbench", **msg})
    js = ("(function s(m,n){const f=document.getElementById('bg-workbench');"
          "const w=f&&f.contentWindow;"
          "if(w&&w.bgWorkbench){w.bgWorkbench.handle(m)}else if(n<200){setTimeout(()=>s(m,n+1),100)}})"
          f"({payload},0)")
    return rx.call_script(js)


class State(rx.State):
    # ---- input
    source: str = "rcsb"                 # rcsb | alphafold | upload | local
    pdb_id: str = "4HHB"
    uniprot: str = "P69905"
    chains: str = ""
    mode: str = "allatom"                # allatom | ca  (+ ligands, ions, coordinated waters in both)
    assembly: str = ""
    waters: str = "coordinated"
    local_path: str = str(Path.home())
    local_files: list[str] = []
    uploaded_path: str = ""
    uploaded_name: str = ""

    # ---- job
    job_id: str = ""
    summary: dict = {}
    busy: bool = False
    busy_msg: str = ""
    error: str = ""
    info: str = ""

    # ---- geometry parameters (the workbench toggles which layers are shown)
    d_max_edge: float = 5.0
    d_filter: str = "all"
    d_tets: bool = False
    v_radius: float = 5.0
    v_margin: float = 5.0                  # bounding-box margin that closes the surface cells
    v_subset: str = "bounded"
    v_shrink: float = 0.85
    v_style: str = "cells"
    a_alpha: float = 4.0
    h_scope: str = "global"
    h_layers: bool = True

    # ---- analysis parameters
    p_small: float = 4.0
    p_large: float = 10.0
    p_minvol: float = 200.0
    c_cutoff: float = 5.0
    f_min: float = 1.0
    f_max: float = 12.0
    fb_edge: float = 10.0
    analysis_done: list[str] = []

    # ---- display (view-only)
    mol_style: str = "cartoon"
    mol_color: str = "chain"
    geo_opacity: float = 0.5

    # ================================================================ computed
    @rx.var
    def loaded(self) -> bool:
        return self.job_id != ""

    @rx.var
    def title_line(self) -> str:
        s = self.summary
        if not s:
            return ""
        t = str(s.get("title", "")).strip()
        return f"{s.get('name', '')} — {t[:1].upper() + t[1:].lower() if t and t.isupper() else t}"

    @rx.var
    def summary_line(self) -> str:
        s = self.summary
        if not s:
            return ""
        mode = "all-atom" if s.get("level") in ("heavy", "backbone") else "backbone"
        return (f"{s.get('n_atoms', 0):,} atoms · {s.get('n_protein', 0) + s.get('n_nucleic', 0):,} res · "
                f"{s.get('n_chains', 0)} ch · {mode} + hetero · {s.get('n_points', 0):,} pts")

    @rx.var
    def meta_line(self) -> str:
        m = (self.summary or {}).get("meta", {}) or {}
        parts = [m.get("method", "")]
        if m.get("resolution"):
            parts.append(f"{m['resolution']} Å")
        if m.get("mean_plddt"):
            parts.append(f"mean pLDDT {m['mean_plddt']}")
        if m.get("organism"):
            parts.append(m["organism"])
        return " · ".join(p for p in parts if p)

    @rx.var
    def local_file_names(self) -> list[str]:
        return [Path(p).name for p in self.local_files]

    @rx.var
    def structure_links(self) -> list[list[str]]:
        if not self.job_id:
            return []
        base = _api_base()
        note = "the loaded model: chains, assembly and coordinated waters as shown"
        return [["structure · mmCIF (.cif)", note, f"{base}/bg/jobs/{self.job_id}/structure/cif"],
                ["structure · PDB (.pdb)", note + " (1-letter chain ids; very large assemblies need mmCIF)",
                 f"{base}/bg/jobs/{self.job_id}/structure/pdb"]]

    @rx.var
    def csv_links(self) -> list[list[str]]:
        if not self.job_id:
            return []
        base = _api_base()
        return [[k.replace("_", " "), d, f"{base}/bg/jobs/{self.job_id}/csv/{k}"] for k, d in engine.CSV_KINDS.items()]

    # ================================================================ input handlers
    @rx.event
    def set_text(self, name: str, value: str):
        setattr(self, name, value)

    @rx.event
    def set_source(self, v: str | list[str]):
        self.source = v[0] if isinstance(v, list) else v
        self.error = ""

    @rx.event
    def list_local(self):
        self.error = ""
        p = Path(self.local_path).expanduser()
        try:
            self.local_files = [str(p)] if p.is_file() else list_local_structures(p)
            if not self.local_files:
                self.info = f"no structure files ({', '.join(s for s in STRUCTURE_SUFFIXES if s.count('.') == 1)}) in {p}"
        except Exception as exc:  # noqa: BLE001
            self.local_files = []
            self.error = str(exc)

    @rx.event
    def choose_local(self, name: str):
        match = [p for p in self.local_files if Path(p).name == name]
        if match:
            self.local_path = match[0]
            return State.load_structure

    @rx.event
    async def handle_upload(self, files: list[rx.UploadFile]):
        self.error = ""
        for f in files[:1]:
            data = await f.read()
            try:
                path = store_upload(data, f.name or "upload.cif")
            except Exception as exc:  # noqa: BLE001
                self.error = str(exc)
                return
            self.uploaded_path, self.uploaded_name = str(path), f.name or path.name
        return State.load_structure

    @rx.event
    def load_example(self, source: str, ident: str, mode: str = "allatom"):
        self.source = source
        self.mode = mode
        if source == "alphafold":
            self.uniprot = ident
        else:
            self.pdb_id = ident
        self.chains = ""
        return State.load_structure

    def _spec(self) -> engine.JobSpec:
        common = dict(chains=self.chains, level=MODES.get(self.mode, "heavy"), assembly=self.assembly.strip(),
                      ligands=True, waters=self.waters)
        if self.source == "rcsb":
            return engine.JobSpec(source="rcsb", pdb_id=self.pdb_id.strip(), **common)
        if self.source == "alphafold":
            return engine.JobSpec(source="alphafold", uniprot=self.uniprot.strip(), **common)
        if self.source == "local":
            return engine.JobSpec(source="file", path=self.local_path.strip(), **common)
        return engine.JobSpec(source="upload", path=self.uploaded_path, **common)

    # ================================================================ background work
    @rx.event(background=True)
    async def load_structure(self):
        async with self:
            if self.busy:
                return
            if self.source == "local" and Path(self.local_path).expanduser().is_dir():
                self.error = "choose a file (press ‘List’ and pick one), not a directory"
                return
            if self.source == "upload" and not self.uploaded_path:
                self.error = "drop a PDB / mmCIF file first"
                return
            self.busy, self.busy_msg, self.error, self.info = True, "Loading structure…", "", ""
            spec = self._spec()

        def work():
            job = engine.open_job(spec)
            meta = engine.build_layer(job, "structure")
            return job.summary(), meta

        try:
            summary, meta = await asyncio.to_thread(work)
        except Exception as exc:  # noqa: BLE001
            async with self:
                self.busy, self.error = False, str(exc)
            return
        async with self:
            self.job_id, self.summary = summary["job"], summary
            self._apply_defaults(summary["defaults"])
            self.analysis_done = []
            self.busy = False
            params = self._params()
        yield send_to_workbench({"type": "params", "params": params})
        yield send_to_workbench({"type": "structure", "api": _api_base(), "job": summary["job"], "key": meta["key"]})
        yield self._display_event()

    def _apply_defaults(self, d: dict):
        self.d_max_edge = d["delaunay"]["max_edge"]
        self.v_radius = d["voronoi"]["max_radius"]
        self.v_margin = d["voronoi"].get("box_margin", d["voronoi"]["max_radius"])
        self.a_alpha = d["alpha"]["alpha"]
        self.p_small, self.p_large = d["pockets"]["alpha_small"], d["pockets"]["alpha_large"]
        self.c_cutoff = d["contacts"]["cutoff"]
        self.f_min, self.f_max = d["filtration"]["alpha_min"], d["filtration"]["alpha_max"]

    def _params(self) -> dict:
        return {
            "delaunay": {"max_edge": self.d_max_edge, "edge_filter": self.d_filter, "show_tets": self.d_tets},
            "voronoi": {"max_radius": self.v_radius, "box_margin": self.v_margin, "subset": self.v_subset,
                        "shrink": self.v_shrink, "style": self.v_style},
            "alpha": {"alpha": self.a_alpha},
            "hull": {"scope": self.h_scope, "layers": self.h_layers},
            "contacts": {"cutoff": self.c_cutoff},
        }

    def _analysis_params(self, kind: str) -> dict:
        return {"pockets": {"alpha_small": self.p_small, "alpha_large": self.p_large, "min_volume": self.p_minvol},
                "contacts": {"cutoff": self.c_cutoff}, "filtration": {"alpha_min": self.f_min, "alpha_max": self.f_max},
                "fourbody": {"max_edge": self.fb_edge}}.get(kind, {})

    @rx.event(background=True)
    async def run_analysis(self, kind: str):
        async with self:
            if not self.job_id or kind not in ANALYSES or self.busy:
                return
            job_id, params = self.job_id, self._analysis_params(kind)
            self.busy, self.busy_msg, self.error = True, f"Running {kind}…", ""
        job = engine.get_job(job_id)
        if job is None:
            async with self:
                self.busy, self.error = False, "session expired - press Load again"
            return
        try:
            meta = await asyncio.to_thread(engine.build_layer, job, kind, params)
        except Exception as exc:  # noqa: BLE001
            async with self:
                self.busy, self.error = False, str(exc)
            return
        async with self:
            self.busy = False
            if kind not in self.analysis_done:
                self.analysis_done = self.analysis_done + [kind]
            self.info = _analysis_info(kind, meta["stats"])
        yield send_to_workbench({"type": "analysis", "kind": kind, "key": meta["key"]})

    # ================================================================ parameters + display
    @rx.event
    def set_param(self, name: str, value):
        """Geometry parameter: update and push to the workbench, which recomputes visible layers."""
        if isinstance(value, list):
            value = value[0]
        cur = getattr(self, name)
        try:
            value = bool(value) if isinstance(cur, bool) else float(value) if isinstance(cur, float) else value
        except (TypeError, ValueError):
            return
        setattr(self, name, value)
        if self.job_id:
            return send_to_workbench({"type": "params", "params": self._params()})

    @rx.event
    def set_value(self, name: str, value):
        """Setter that never triggers computation (load options, analysis parameters)."""
        if isinstance(value, list):
            value = value[0]
        cur = getattr(self, name)
        try:
            value = bool(value) if isinstance(cur, bool) else float(value) if isinstance(cur, float) else value
        except (TypeError, ValueError):
            return
        setattr(self, name, value)

    def _display_event(self):
        return send_to_workbench({"type": "display", "opts": {"molStyle": self.mol_style, "molColor": self.mol_color,
                                                              "geoOpacity": self.geo_opacity}})

    @rx.event
    def set_display(self, name: str, value):
        if isinstance(value, list):
            value = value[0]
        setattr(self, name, float(value) if isinstance(getattr(self, name), float) else value)
        return self._display_event()

    @rx.event
    def reset_view(self):
        return send_to_workbench({"type": "reset"})

    @rx.event
    def reset_workbench(self):
        """A clean slate with the loaded structure: no selection, default parameters, layers, views and tables."""
        self.analysis_done, self.info, self.error = [], "", ""
        self.mol_style, self.mol_color, self.geo_opacity = "cartoon", "chain", 0.5
        if self.summary:
            self._apply_defaults(self.summary["defaults"])
        return send_to_workbench({"type": "resetall", "params": self._params(),
                                  "display": {"molStyle": "cartoon", "molColor": "chain", "geoOpacity": 0.5}})

    @rx.event
    def set_mode(self, value):
        self.mode = value[0] if isinstance(value, list) else value

    @rx.event
    def download_file(self, url: str):
        """Fetch a download in the page and save it: no new tab, no flicker."""
        return rx.call_script(DOWNLOAD_JS % json.dumps(url))

    @rx.event
    def measure_selection(self):
        """The selection lives in the workbench, so the workbench requests the layer itself."""
        return send_to_workbench({"type": "analysis", "kind": "selection"})


def _analysis_info(kind: str, s: dict) -> str:
    if kind == "pockets":
        enc = sum(1 for p in s.get("pockets", []) if p["enclosed"])
        return f"{s.get('n_pockets', 0)} pockets ({enc} enclosed) — Pockets tab"
    if kind == "contacts":
        return f"{s.get('n_interchain', 0):,} inter-chain contacts · {s.get('n_interface_residues', 0):,} interface residues"
    if kind == "scaling":
        ex = s.get("exponents", {})
        return ("time ~ nᵏ: " + ", ".join(f"{k} {v:.2f}" for k, v in ex.items())) if ex else "scaling done (size range too small for exponents)"
    if kind == "fourbody":
        return f"{s.get('n_tets', 0):,} residue tetrahedra · {s.get('n_distinct', 0):,} quadruplets"
    return f"α filtration over {len(s.get('alpha', []))} values"
