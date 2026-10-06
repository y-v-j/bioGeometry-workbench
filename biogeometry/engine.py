"""Job cache + layer builders shared by the Reflex app, the HTTP API and scripts.

A *job* is one structure selection (source, chains, assembly, waters, point level). Each view or
analysis is a *layer* computed from a job with a parameter set and written to
``data/cache/jobs/<job>/<layer-key>/`` as ``meta.json`` + raw little-endian ``<name>.bin`` arrays,
so the workbench fetches typed arrays directly (no JSON number bloat).
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import shutil
import threading
import time
from collections import OrderedDict
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from . import analysis as an
from . import annotations, backbone, biochem, features, resources
from . import selection as selgeo
from . import geometry as g
from .sasa import atom_sasa, radii_for
from .structure import (DATA_DIR, ELEMENTS, KIND_ION, KIND_LIGAND, KIND_NAMES, KIND_NUCLEIC, KIND_PROTEIN,
                        KIND_WATER, PointSet, Structure, alphafold_pae, auto_level, dedupe_points,
                        fetch_alphafold, fetch_structure, hetero_groups, load_structure_file,
                        parse_chain_ids, representative_points, residue_points, select_chains,
                        structure_text, chain_tag)

CACHE_DIR = DATA_DIR / "cache" / "jobs"
MAX_POINTS = 400_000
MAX_JOBS_IN_MEMORY = 3
MAX_JOB_DIRS = 40
SCHEMA = 10                      # bump to invalidate cached layers after format changes
MAX_TETS_SENT = 1_000_000        # Delaunay tetrahedra shipped to the viewer (simplex highlighting)

LAYER_KINDS = ("structure", "delaunay", "voronoi", "alpha", "hull", "filtration", "pockets", "contacts",
               "features", "tetrahedra", "fourbody", "scaling", "pae", "selection")


def default_params(level: str, nucleic_fraction: float = 0.0) -> dict:
    coarse = level in ("ca", "centroid")
    k = 1.0 + 0.5 * min(1.0, nucleic_fraction / 0.5) if coarse else 1.0
    return {
        "delaunay": {"max_edge": round((10.0 if coarse else 5.0) * k, 1), "edge_filter": "all",
                     "show_tets": False, "max_tets": 20000},
        "voronoi": {"max_radius": round((10.0 if coarse else 5.0) * k, 1), "box_margin": round((10.0 if coarse else 5.0) * k, 1),
                    "subset": "bounded", "shrink": 0.85,
                    "max_cells": 12000, "style": "cells"},
        "alpha": {"alpha": round((8.0 if coarse else 4.0) * k, 1)},
        "hull": {"scope": "global", "layers": True},
        "filtration": {"alpha_min": 2.0 if coarse else 1.0, "alpha_max": 30.0 if coarse else 12.0, "steps": 40},
        "pockets": {"alpha_small": round((8.0 if coarse else 4.0) * k, 1),
                    "alpha_large": round((25.0 if coarse else 10.0) * k, 1), "min_volume": 200.0, "max_pockets": 25},
        "contacts": {"cutoff": round((10.0 if coarse else 5.0) * k, 1), "min_seq_sep": 3},
        "features": {}, "tetrahedra": {"max_edge": round((10.0 if coarse else 5.0) * k, 1)}, "pae": {},
        "fourbody": {"max_edge": 10.0, "top": 25},
        "scaling": {"steps": 6},
        "selection": {"residues": [], "alpha": None, "max_edge": None},
    }


class EngineError(ValueError):
    pass


@dataclass
class JobSpec:
    source: str = "rcsb"          # rcsb | alphafold | file | upload
    pdb_id: str = ""
    uniprot: str = ""
    path: str = ""
    chains: str = ""
    level: str = "auto"
    assembly: str = ""
    ligands: bool = False         # include ligands, ions and coordinated waters in the geometry point set
    waters: str = "coordinated"   # coordinated | none

    def key(self) -> str:
        if self.source == "rcsb":
            src = f"rcsb:{self.pdb_id.upper()}"
        elif self.source == "alphafold":
            src = f"af:{self.uniprot.upper()}"
        else:
            p = Path(self.path).expanduser().resolve()
            st = p.stat()
            src = f"file:{p}:{st.st_size}:{int(st.st_mtime)}"
        blob = json.dumps([SCHEMA, src, sorted(parse_chain_ids(self.chains) or []), self.level, self.assembly,
                           self.ligands, self.waters])
        return hashlib.sha1(blob.encode()).hexdigest()[:12]


@dataclass
class Job:
    id: str
    spec: JobSpec
    structure: Structure
    points: PointSet
    center: np.ndarray
    timings: dict = field(default_factory=dict)
    n_duplicates: int = 0
    lock: threading.RLock = field(default_factory=threading.RLock)
    _vor: dict = field(default_factory=dict)

    @property
    def dir(self) -> Path:
        return CACHE_DIR / self.id

    def _timed(self, name, fn):
        t0 = time.perf_counter()
        out = fn()
        self.timings[name] = round(time.perf_counter() - t0, 3)
        return out

    # -- geometry, computed lazily and reused across layers
    @cached_property
    def tess(self) -> g.Tess:
        _guard(len(self.points.coords), "Delaunay")
        def run():
            t = g.Tess.build(self.points.coords)
            _ = t.circumradii, t.volumes, t.edges
            return t
        return self._timed("delaunay", run)

    @cached_property
    def pocket_tess(self) -> g.Tess:
        """Pockets are looked for in the polymer alone: a bound ligand would otherwise fill its own pocket.
        Polymer points come first in the point set, so indices stay valid for the full set."""
        poly = np.isin(self.structure.res_kind[self.points.residue], (KIND_PROTEIN, KIND_NUCLEIC))
        n = int(poly.sum())
        if n == len(poly) or not poly[:n].all() or n < 5:
            return self.tess
        return self._timed("delaunay_polymer", lambda: g.Tess.build(self.points.coords[:n]))

    @cached_property
    def hull(self) -> g.Hull:
        return self._timed("hull", lambda: g.convex_hull(self.points.coords))

    @cached_property
    def residue_points(self) -> PointSet:
        ps, _ = dedupe_points(residue_points(self.structure))
        return ps

    @cached_property
    def residue_tess(self) -> g.Tess:
        if self.points.level == "ca" and len(self.points.coords) == len(self.residue_points.coords):
            return self.tess
        return g.Tess.build(self.residue_points.coords)

    def voronoi(self, max_radius: float, box_margin: float | None = None) -> g.VoronoiCells:
        """Cells clipped to the bounding box expanded by box_margin (default: max_radius), so every
        point gets a volume; 'bounded' marks the interior cells (no vertex beyond max_radius)."""
        margin = float(max_radius if box_margin is None else box_margin)
        key = (round(float(max_radius), 3), round(margin, 3))
        if key not in self._vor:
            _guard(len(self.points.coords), "Voronoi")
            self._vor.clear()
            self._vor[key] = self._timed("voronoi", lambda: g.voronoi_cells(self.points.coords, key[0], key[1]))
        return self._vor[key]

    @cached_property
    def layers(self) -> np.ndarray:
        return self._timed("convex_layers", lambda: g.convex_layers(self.points.coords))

    # -- backbone / biochemistry
    @cached_property
    def bb(self) -> dict:
        return backbone.backbone_atoms(self.structure)

    @cached_property
    def dihedrals(self):
        return backbone.dihedrals(self.structure, self.bb)

    @cached_property
    def segments(self) -> list[np.ndarray]:
        return backbone.protein_segments(self.structure, self.bb)

    @cached_property
    def _ss(self):
        return self._timed("secondary_structure", lambda: backbone.secondary_structure(self.structure, self.bb, self.segments))

    @property
    def ss(self) -> np.ndarray:
        return self._ss[0]

    @property
    def hbond_pairs(self) -> np.ndarray:
        return self._ss[1]

    @cached_property
    def hetero(self) -> list[dict]:
        return hetero_groups(self.structure)

    @cached_property
    def sasa(self) -> np.ndarray:
        s = self.structure
        use = (s.res_kind != KIND_WATER)[s.atom_res]
        out = np.zeros(s.n_atoms)
        el = [ELEMENTS[e] if e < len(ELEMENTS) else "C" for e in s.element[use]]
        out[use] = self._timed("sasa", lambda: atom_sasa(s.coords[use], radii_for(el)))
        return out

    @cached_property
    def interactions(self) -> dict:
        return self._timed("interactions", lambda: features.interactions(self.structure, self.hetero, self.hbond_pairs))

    @cached_property
    def default_pockets(self):
        d = self.defaults()["pockets"]
        try:
            return g.find_pockets(self.pocket_tess, d["alpha_small"], d["alpha_large"], d["min_volume"], d["max_pockets"])
        except g.GeometryError:
            return []

    @cached_property
    def _domains(self):
        return self._timed("annotations", lambda: annotations.domains(self.structure))

    @property
    def domains(self) -> list[dict]:
        return self._domains[0]

    @property
    def nucleic_fraction(self) -> float:
        poly = np.isin(self.structure.res_kind, (KIND_PROTEIN, KIND_NUCLEIC))
        return float((self.structure.res_kind[poly] == KIND_NUCLEIC).mean()) if poly.any() else 0.0

    def defaults(self) -> dict:
        return default_params(self.points.level, self.nucleic_fraction)

    def summary(self) -> dict:
        s = self.structure
        kinds = np.bincount(s.res_kind, minlength=5)
        return {
            "job": self.id, "name": s.name, "title": s.title, "source": s.source,
            "source_kind": self.spec.source, "assembly": s.assembly or "", "chains": s.chains,
            "n_chains": len(s.chains), "n_atoms": int(s.n_atoms), "n_residues": int(s.n_residues),
            "n_protein": int(kinds[KIND_PROTEIN]), "n_nucleic": int(kinds[KIND_NUCLEIC]),
            "n_ligand": int(kinds[KIND_LIGAND]), "n_ion": int(kinds[KIND_ION]), "n_water": int(kinds[KIND_WATER]),
            "level": self.points.level, "n_points": int(len(self.points.coords)), "n_duplicates": self.n_duplicates,
            "is_alphafold": s.is_alphafold, "meta": s.meta, "defaults": self.defaults(),
            "timings": dict(self.timings),
        }


# ------------------------------------------------------------------ job cache

_jobs: "OrderedDict[str, Job]" = OrderedDict()
_jobs_lock = threading.Lock()


def _guard(n: int, what: str):
    if n > MAX_POINTS:
        raise EngineError(f"{what} on {n:,} points exceeds the {MAX_POINTS:,}-point limit; "
                          f"use a coarser level (ca / centroid) or select chains")


def open_job(spec: JobSpec) -> Job:
    if spec.source == "rcsb" and not spec.pdb_id.strip():
        raise EngineError("enter a PDB id")
    if spec.source == "alphafold" and not spec.uniprot.strip():
        raise EngineError("enter a UniProt accession (e.g. P69905)")
    if spec.source in ("file", "upload") and not spec.path:
        raise EngineError("no file given")
    jid = spec.key()
    with _jobs_lock:
        if jid in _jobs:
            _jobs.move_to_end(jid)
            return _jobs[jid]
    t0 = time.perf_counter()
    kw = {"waters": spec.waters}
    if spec.source == "rcsb":
        st = fetch_structure(spec.pdb_id, assembly=spec.assembly.strip() or None, **kw)
    elif spec.source == "alphafold":
        st = fetch_alphafold(spec.uniprot, **kw)
    else:
        st = load_structure_file(spec.path, assembly=spec.assembly.strip() or None, **kw)
    t_parse = time.perf_counter() - t0
    st = select_chains(st, parse_chain_ids(spec.chains))
    level = auto_level(st) if spec.level == "auto" else spec.level
    ps, dup = dedupe_points(representative_points(st, level, include_hetero=spec.ligands))
    job = Job(jid, spec, st, ps, center=st.coords.mean(axis=0), n_duplicates=dup)
    job.timings["parse"] = round(t_parse, 3)
    with _jobs_lock:
        _jobs[jid] = job
        while len(_jobs) > MAX_JOBS_IN_MEMORY:
            _jobs.popitem(last=False)
    _prune_cache()
    return job


def get_job(job_id: str) -> Job | None:
    with _jobs_lock:
        return _jobs.get(job_id)


def _prune_cache():
    if not CACHE_DIR.exists():
        return
    dirs = sorted((d for d in CACHE_DIR.iterdir() if d.is_dir()), key=lambda d: d.stat().st_mtime)
    for d in dirs[:-MAX_JOB_DIRS]:
        shutil.rmtree(d, ignore_errors=True)


# ------------------------------------------------------------------ layer storage

def layer_key(kind: str, params: dict) -> str:
    h = hashlib.sha1(json.dumps(params, sort_keys=True, default=str).encode()).hexdigest()[:10]
    return f"{kind}-{h}"


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    if isinstance(x, (np.floating, float)):
        v = float(x)
        return v if np.isfinite(v) else None
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def _save(job: Job, key: str, kind: str, params: dict, stats: dict, buffers: dict[str, np.ndarray]):
    d = job.dir / key
    tmp = job.dir / (key + ".tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    meta_bufs = {}
    for name, arr in buffers.items():
        a = np.ascontiguousarray(arr)
        if a.dtype == np.float64:
            a = a.astype(np.float32)
        elif a.dtype == np.int64:
            a = a.astype(np.uint32 if (a.size == 0 or a.min() >= 0) else np.int32)
        elif a.dtype == bool:
            a = a.astype(np.uint8)
        a = a.astype(a.dtype.newbyteorder("<"), copy=False)
        (tmp / f"{name}.bin").write_bytes(a.tobytes())
        meta_bufs[name] = {"dtype": a.dtype.name, "shape": list(a.shape)}
    meta = {"job": job.id, "key": key, "kind": kind, "params": params, "stats": _jsonable(stats),
            "buffers": meta_bufs}
    (tmp / "meta.json").write_text(json.dumps(meta, allow_nan=False))
    shutil.rmtree(d, ignore_errors=True)
    tmp.rename(d)
    return meta


def _safe(*parts: str) -> bool:
    """Path components from URLs: letters, digits, '-' and '_' only (no '..', no separators)."""
    return all(p and all(c.isalnum() or c in "-_" for c in p) for p in parts)


def read_meta(job_id: str, key: str) -> dict | None:
    if not _safe(job_id, key):
        return None
    p = CACHE_DIR / job_id / key / "meta.json"
    return json.loads(p.read_text()) if p.exists() else None


def buffer_path(job_id: str, key: str, name: str) -> Path | None:
    if not _safe(job_id, key, name):
        return None
    p = CACHE_DIR / job_id / key / f"{name}.bin"
    return p if p.exists() else None


def build_layer(job: Job, kind: str, params: dict | None = None) -> dict:
    if kind not in LAYER_KINDS:
        raise EngineError(f"unknown layer kind '{kind}'")
    p = dict(job.defaults().get(kind, {}))
    p.update({k: v for k, v in (params or {}).items() if v is not None})
    key = layer_key(kind, p)
    meta = read_meta(job.id, key)
    if meta is not None:
        return meta
    with job.lock:
        meta = read_meta(job.id, key)
        if meta is not None:
            return meta
        if kind == "selection":
            n_sel = int(np.isin(job.structure.atom_res, np.asarray(p.get("residues") or [], dtype=np.int64)).sum())
            need = resources.selection_estimate_mb(min(n_sel, selgeo.MAX_SELECTION_ATOMS))
        else:
            need = resources.layer_estimate_mb(kind, len(job.points.coords), job.structure.n_atoms,
                                               job.points.level in ("ca", "centroid"))
        with resources.budget(need, f"{kind} on {len(job.points.coords):,} {job.points.level} points"):
            t0 = time.perf_counter()
            stats, bufs = BUILDERS[kind](job, p)
            stats["compute_s"] = round(time.perf_counter() - t0, 3)
        stats["bytes"] = int(sum(np.asarray(b).nbytes for b in bufs.values()))
        return _save(job, key, kind, p, stats, bufs)


# ------------------------------------------------------------------ builders

def _codes(values) -> tuple[list, np.ndarray]:
    uniq, code = np.unique(np.asarray(values, dtype=object).astype(str), return_inverse=True)
    return uniq.tolist(), code.astype(np.uint16)


def _bonds(s: Structure) -> np.ndarray:
    """Covalent bonds by distance (≤ 1.9 Å, ≤ 2.2 Å involving S/Se/P), never to metals or waters."""
    elem = np.array([ELEMENTS[e] if e < len(ELEMENTS) else "X" for e in s.element])
    ok = ~np.isin(elem, ["Mg", "Zn", "Fe", "Ca", "Mn", "Na", "K", "Cu", "Co", "Ni"]) & (s.res_kind[s.atom_res] != KIND_WATER)
    idx = np.flatnonzero(ok)
    if len(idx) < 2:
        return np.zeros((0, 2), np.uint32)
    pairs = cKDTree(s.coords[idx]).query_pairs(2.2, output_type="ndarray")
    if not len(pairs):
        return np.zeros((0, 2), np.uint32)
    a, b = idx[pairs[:, 0]], idx[pairs[:, 1]]
    d = np.linalg.norm(s.coords[a] - s.coords[b], axis=1)
    heavy = np.isin(elem[a], ["S", "Se", "P"]) | np.isin(elem[b], ["S", "Se", "P"])
    keep = (d > 0.4) & ((d <= 1.9) | (heavy & (d <= 2.2)))
    return np.column_stack([a[keep], b[keep]]).astype(np.uint32)


def _structure(job: Job, p: dict):
    s, ps = job.structure, job.points
    c = job.center
    R = s.n_residues
    res_names, name_code = _codes(s.res_name)
    atom_names, atom_name_code = _codes(s.atom_name)
    classes, class_code = _codes(biochem.residue_class(s.res_name.astype(str), s.res_kind))
    icodes, icode_code = _codes(s.res_icode)
    ch = s.chain_index()
    # residue anchor: CA / P for polymers, centroid otherwise
    cnt = np.maximum(np.bincount(s.atom_res, minlength=R), 1)
    cen = np.stack([np.bincount(s.atom_res, weights=s.coords[:, k], minlength=R) for k in range(3)], axis=1) / cnt[:, None]
    anchor = cen.copy()
    rp = job.residue_points
    anchor[rp.residue] = rp.coords
    bb = job.bb
    # trace for coarse styles: consecutive polymer residues in a chain that are bonded
    same = ch[:-1] == ch[1:]
    dd = np.linalg.norm(anchor[1:] - anchor[:-1], axis=1)
    poly = np.isin(s.res_kind, (KIND_PROTEIN, KIND_NUCLEIC))
    lim = np.where(s.res_kind[:-1] == KIND_NUCLEIC, 8.5, 4.4)
    link = same & poly[:-1] & poly[1:] & (dd < lim)
    i = np.flatnonzero(link)
    phi, psi, omega = job.dihedrals
    cart = job._timed("cartoon", lambda: backbone.cartoon_mesh(s, bb, job.ss, job.segments))
    hetero = [dict(h) for h in job.hetero]
    poly_atoms = np.flatnonzero(poly[s.atom_res])
    if len(poly_atoms) and hetero:
        ptree = cKDTree(s.coords[poly_atoms])
        for h in hetero:   # binding-site residues: polymer atoms within 4.5 Å of the group
            members = h.get("members", [h["res"]])
            atoms = np.flatnonzero(np.isin(s.atom_res, members))
            hits = ptree.query_ball_point(s.coords[atoms], features.CONTACT)
            site = np.unique(s.atom_res[poly_atoms[np.fromiter((j for x in hits for j in x), dtype=np.int64)]]) \
                if any(len(x) for x in hits) else np.zeros(0, int)
            h["site"] = site.tolist()
    doms, dom_note = job._domains
    bind_pairs, bind_type, bind_dist = features.binding_contacts(s)
    aa = biochem.one_letter(s.res_name.astype(str))
    chain_rows = []
    for ci, cid in enumerate(s.chains):
        m = ch == ci
        kinds = np.bincount(s.res_kind[m], minlength=5)
        polyres = np.flatnonzero(m & poly)
        info = s.chain_info.get(cid, {})
        chain_rows.append({
            "id": cid, "index": ci, "n_res": int(len(polyres)),
            "kind": "nucleic" if kinds[KIND_NUCLEIC] > kinds[KIND_PROTEIN] else ("protein" if len(polyres) else "hetero"),
            "first": int(polyres[0]) if len(polyres) else int(np.flatnonzero(m)[0]),
            "last": int(polyres[-1]) if len(polyres) else int(np.flatnonzero(m)[-1]),
            "seq_from": int(s.res_seq[polyres[0]]) if len(polyres) else None,
            "seq_to": int(s.res_seq[polyres[-1]]) if len(polyres) else None,
            "description": info.get("description", ""), "uniprot": info.get("uniprot", ""),
            "n_hetero": int(kinds[KIND_LIGAND] + kinds[KIND_ION]), "n_water": int(kinds[KIND_WATER]),
        })
    stats = job.summary()
    stats.update(res_names=res_names, atom_names=atom_names, classes=classes, icodes=icodes,
                 elements=ELEMENTS + ["X"], kinds=list(KIND_NAMES), center=c.tolist(),
                 extent=float(np.linalg.norm(s.coords - c, axis=1).max()), chain_rows=chain_rows,
                 hetero=hetero, domains=doms, domain_note=dom_note,
                 cartoon={"sides": cart["sides"], "tube_sides": cart.get("tube_sides", cart["sides"])},
                 aa="".join(aa), bind_types=list(features.BIND_NAMES),
                 bind_counts=np.bincount(bind_type, minlength=len(features.BIND_NAMES)).tolist(),
                 mode="all-atom" if ps.level in ("heavy", "backbone") else "Cα", hetero_points=bool(job.spec.ligands))
    bufs = {
        "atom_pos": (s.coords - c).astype(np.float32), "atom_res": s.atom_res.astype(np.uint32),
        "atom_elem": s.element, "atom_b": s.bfactor, "atom_name": atom_name_code,
        "res_pos": (anchor - c).astype(np.float32), "res_chain": ch.astype(np.uint16), "res_kind": s.res_kind,
        "res_name": name_code, "res_seq": s.res_seq, "res_icode": icode_code, "res_class": class_code,
        "res_ss": job.ss, "res_phi": phi.astype(np.float32), "res_psi": psi.astype(np.float32),
        "res_omega": omega.astype(np.float32),
        "res_b": an.per_residue(s.bfactor, s.atom_res, R, "mean").astype(np.float32),
        "trace": np.stack([i, i + 1], axis=1).astype(np.uint32),
        "bonds": _bonds(s), "ladder": backbone.base_ladder(s, bb).astype(np.uint32),
        "cart_pos": cart["pos"] - c.astype(np.float32) if len(cart["pos"]) else cart["pos"],
        "cart_res": cart["res"], "cart_rings": cart["rings"],
        "pt_pos": (ps.coords - c).astype(np.float32), "pt_res": ps.residue.astype(np.uint32),
        "pt_atom": ps.atom.astype(np.int32),
        "bind_pairs": bind_pairs.astype(np.uint32), "bind_type": bind_type, "bind_dist": bind_dist,
    }
    return stats, bufs


def _sample(n: int, k: int, seed: int = 0) -> np.ndarray:
    if n <= k:
        return np.arange(n)
    return np.sort(np.random.default_rng(seed).choice(n, size=k, replace=False))


def _delaunay(job: Job, p: dict):
    t = job.tess
    ps = job.points
    ch = job.structure.chain_index()[ps.residue]
    e = t.edges
    L = t.edge_lengths()
    inter = ch[e[:, 0]] != ch[e[:, 1]]
    keep = L <= float(p["max_edge"])
    if p["edge_filter"] == "inter":
        keep &= inter
    elif p["edge_filter"] == "intra":
        keep &= ~inter
    bufs = {"edges": e[keep], "edge_len": L[keep], "edge_inter": inter[keep]}
    R = t.circumradii
    finite_r = R[np.isfinite(R)]
    # the compact simplices (every edge within max_edge) go to the viewer: it highlights the ones
    # touching the selection and draws them all, shrunken, when "show tetrahedra" is on
    compact = np.nonzero(t.max_edge_per_tet() <= float(p["max_edge"]))[0]
    pick = compact[_sample(len(compact), MAX_TETS_SENT)]
    bufs["tets"] = t.simplices[pick]
    bufs["tet_r"] = np.minimum(R[pick], 1e3)
    bufs["hist_edge_len"] = L[_sample(len(L), 200_000)]
    bufs["hist_circumradius"] = np.minimum(finite_r[_sample(len(finite_r), 200_000)], 100.0)
    stats = {
        "n_points": len(t.points), "n_tets": t.n_tets, "n_edges": len(e), "n_shown": int(keep.sum()),
        "n_interchain": int((inter & (L <= p["max_edge"])).sum()), "mean_edge": float(L.mean()),
        "median_edge": float(np.median(L)), "mean_degree": float(2 * len(e) / len(t.points)),
        "tets_per_point": float(t.n_tets / len(t.points)), "delaunay_s": job.timings.get("delaunay"),
        "n_compact_tets": int(len(compact)), "tets_sampled": bool(len(pick) < len(compact)),
    }
    return stats, bufs


def _voronoi(job: Job, p: dict):
    ps = job.points
    cells = job.voronoi(float(p["max_radius"]), float(p.get("box_margin") or p["max_radius"]))
    subset = str(p["subset"])
    sel = np.ones(len(cells.volume), dtype=bool) if subset == "all" else cells.bounded.copy()
    if subset.startswith("chain:"):
        want = set(parse_chain_ids(subset[6:]) or [])
        sel &= np.isin(job.structure.res_chain[ps.residue], list(want))
    elif subset == "buried":
        sel &= g.hull_depth(job.hull, ps.coords) > 2 * float(p["max_radius"])
    elif subset == "surface":
        rp = cells.ridge_points
        open_nb = np.zeros(len(cells.vor.points), dtype=bool)
        m = ~cells.padded(cells.bounded)[rp]
        open_nb[rp[:, 0][m[:, 1]]] = True
        open_nb[rp[:, 1][m[:, 0]]] = True
        sel &= open_nb[: len(sel)]
    sites = np.nonzero(sel)[0]
    n_sel = len(sites)
    if n_sel > int(p["max_cells"]):
        d = np.linalg.norm(ps.coords[sites] - ps.coords.mean(axis=0), axis=1)
        sites = np.sort(sites[np.argsort(d)[: int(p["max_cells"])]])
    vol = cells.volume
    bufs = {"site_volume": vol.astype(np.float32), "cell_sites": sites}
    if p.get("style") == "skeleton":
        bufs["vertices"] = (cells.vor.vertices - job.center).astype(np.float32)
        bufs["skeleton"] = g.voronoi_ridge_edges(cells, float(p["max_radius"]))
    else:
        pos, site = g.voronoi_cell_mesh(cells, sites, float(p["shrink"]))
        bufs["cell_pos"] = pos - job.center.astype(np.float32)
        bufs["cell_site"] = site
    s = job.structure
    rn = s.res_name[ps.residue].astype(str)
    fin = np.isfinite(vol)
    by_type = []
    for name in np.unique(rn[fin]):
        v = vol[fin & (rn == name)]
        by_type.append({"name": name, "mean": float(v.mean()), "sd": float(v.std()), "n": int(len(v))})
    by_type.sort(key=lambda r: -r["n"])
    inner = fin & cells.bounded
    stats = {"n_sites": len(vol), "n_bounded": int(cells.bounded.sum()), "bounded_fraction": float(cells.bounded.mean()),
             "n_closed": int(np.isfinite(vol).sum()), "box_margin": cells.box_margin,
             "n_selected": n_sel, "n_drawn": len(sites),
             "mean_volume": float(vol[inner].mean()) if inner.any() else None,
             "median_volume": float(np.median(vol[inner])) if inner.any() else None,
             "total_volume": float(np.nansum(vol)), "mean_faces": float(cells.n_faces.mean()),
             "by_residue_type": by_type[:40], "voronoi_s": job.timings.get("voronoi")}
    return stats, bufs


def _alpha(job: Job, p: dict):
    t = job.tess
    a = g.alpha_shape(t, float(p["alpha"]))
    hull_v = job.hull.volume
    comp_sizes = np.bincount(a.face_component) if len(a.face_component) else np.zeros(0)
    surf = np.zeros(len(t.points), dtype=bool)
    surf[a.surface_points] = True
    res_surface = np.zeros(job.structure.n_residues, dtype=bool)
    res_surface[job.points.residue[surf]] = True
    stats = {"alpha": a.alpha, "volume": a.volume, "area": a.area, "n_faces": len(a.faces),
             "n_components": a.n_components, "n_cavities": a.n_cavities, "n_surface_points": int(surf.sum()),
             "n_isolated_points": int(len(a.isolated_points)), "surface_residues": int(res_surface.sum()),
             "hull_volume": hull_v, "volume_over_hull": a.volume / hull_v if hull_v else None,
             "tets_kept": int(a.keep.sum()), "n_tets": t.n_tets,
             "largest_component_faces": int(comp_sizes.max()) if len(comp_sizes) else 0}
    return stats, {"faces": a.faces, "face_comp": a.face_component, "surface": surf}


def _hull(job: Job, p: dict):
    ps, s = job.points, job.structure
    h = job.hull
    depth = g.hull_depth(h, ps.coords)
    bufs = {"faces": h.faces, "depth": depth}
    stats = {"volume": h.volume, "area": h.area, "sphericity": h.sphericity, "n_vertices": int(len(h.vertices)),
             "n_faces": int(len(h.faces)), "max_depth": float(depth.max()), "mean_depth": float(depth.mean()),
             "shape": g.gyration_descriptors(ps.coords)}
    if p.get("layers"):
        t0 = time.perf_counter()
        bufs["layer"] = job.layers
        stats["n_layers"] = int(job.layers.max() + 1)
        stats["layers_s"] = time.perf_counter() - t0
    if p.get("scope") == "chains":
        ch = s.chain_index()[ps.residue]
        faces, fch, per = [], [], []
        for ci, name in enumerate(s.chains):
            idx = np.nonzero(ch == ci)[0]
            if len(idx) < 4:
                continue
            try:
                hc = g.convex_hull(ps.coords[idx])
            except g.GeometryError:
                continue
            faces.append(idx[hc.faces])
            fch.append(np.full(len(hc.faces), ci, dtype=np.uint16))
            per.append({"chain": name, "volume": hc.volume, "area": hc.area, "sphericity": hc.sphericity,
                        "points": int(len(idx))})
        if faces:
            bufs["chain_faces"] = np.concatenate(faces)
            bufs["chain_face_chain"] = np.concatenate(fch)
        stats["per_chain"] = per
        stats["sum_chain_volume"] = float(sum(r["volume"] for r in per))
    return stats, bufs


def _filtration(job: Job, p: dict):
    f = g.alpha_filtration(job.tess, np.linspace(float(p["alpha_min"]), float(p["alpha_max"]), int(p["steps"])))
    f["hull_volume"] = job.hull.volume
    f["hull_area"] = job.hull.area
    return f, {}


def _pockets(job: Job, p: dict):
    s, ps = job.structure, job.points
    pk = g.find_pockets(job.pocket_tess, float(p["alpha_small"]), float(p["alpha_large"]), float(p["min_volume"]),
                        int(p["max_pockets"]))
    rows, faces, fid, members = [], [], [], []
    for k in pk:
        res = np.unique(ps.residue[k.points])
        rows.append({"id": k.id, "volume": k.volume, "mouth_area": k.mouth_area, "enclosed": k.enclosed,
                     "n_tets": k.n_tets, "n_residues": int(len(res)),
                     "chains": sorted(set(s.res_chain[res].tolist())), "residues": [s.residue_label(r) for r in res[:60]],
                     "res_index": res.tolist()})
        faces.append(k.faces)
        fid.append(np.full(len(k.faces), k.id, dtype=np.uint16))
    bufs = {"faces": np.concatenate(faces), "face_pocket": np.concatenate(fid)} if faces else {}
    return {"pockets": rows, "n_pockets": len(rows), "alpha_small": p["alpha_small"],
            "alpha_large": p["alpha_large"], "polymer_only": job.pocket_tess is not job.tess}, bufs


def _contacts(job: Job, p: dict):
    s = job.structure
    net = an.contact_network(s, job.tess, job.points, float(p["cutoff"]), int(p["min_seq_sep"]))
    stats = {k: net[k] for k in ("cutoff", "n_contacts", "n_interchain", "mean_degree", "pairs")}
    stats["pairs"] = stats["pairs"][:200]
    stats["chain_matrix"] = net["chain_matrix"]
    stats["chains"] = s.chains
    stats["n_interface_residues"] = int(net["interface"].sum())
    return stats, {"res_edges": net["edges"], "res_degree": net["degree"].astype(np.float32),
                   "res_interface": net["interface"], "res_interchain": net["interchain_edges"]}


def _features(job: Job, p: dict):
    t = features.residue_table(job)
    vals = {c["key"]: t.data[c["key"]] for c in t.columns if c["type"] == "num"}
    summary = {k: {"mean": float(np.nanmean(v)) if np.isfinite(v).any() else None,
                   "min": float(np.nanmin(v)) if np.isfinite(v).any() else None,
                   "max": float(np.nanmax(v)) if np.isfinite(v).any() else None} for k, v in vals.items()}
    corr_keys = [k for k in ("rsa", "hull_depth", "vor_volume", "del_degree", "cn10", "bfactor", "hydropathy",
                             "convex_layer", "alpha_surface") if k in vals]
    stats = {"columns": t.columns, "n": t.n, "summary": summary, "correlations": an.spearman(vals, corr_keys),
             "timings": dict(job.timings)}
    return stats, t.buffers("c_")


def _tetrahedra(job: Job, p: dict):
    t, T = features.tetra_table(job, p.get("max_edge"))
    stats = {"columns": t.columns, "n": t.n, "n_compact": T["n_compact"], "n_tets": T["n_tets"], "sampled": T["sampled"],
             "max_edge": T["max_edge"], "consts": T["consts"]}
    return stats, {**t.buffers("c_"), "tets": T["tets"].astype(np.uint32), "res4": T["res4"].astype(np.uint32)}


def _pae(job: Job, p: dict):
    s = job.structure
    if not s.is_alphafold:
        return {"available": False, "note": "PAE is available for AlphaFold models only"}, {}
    m = alphafold_pae(s.meta.get("uniprot", ""))
    if m is None:
        return {"available": False, "note": "PAE file not available"}, {}
    # quantise to 0.25 Å steps in uint8 (0–63.75 Å): an N² matrix stays small
    q = np.clip(np.rint(m * 4), 0, 255).astype(np.uint8)
    prot = np.flatnonzero(s.res_kind == KIND_PROTEIN)
    return {"available": True, "n": int(m.shape[0]), "max": float(m.max()), "scale": 0.25,
            "residues": prot[: m.shape[0]].tolist()}, {"pae": q}


def _fourbody(job: Job, p: dict):
    return an.four_body_statistics(job.structure, job.residue_tess, job.residue_points, float(p["max_edge"]),
                                   int(p["top"])), {}


def _scaling(job: Job, p: dict):
    pts = job.points.coords
    n = len(pts)
    sizes = sorted({int(x) for x in np.geomspace(min(500, n), n, int(p["steps"]))})
    d = job.defaults()
    rows = an.scaling_benchmark(pts, sizes, d["voronoi"]["max_radius"], d["alpha"]["alpha"])
    fits = {}
    for alg in ("delaunay", "alpha", "voronoi", "hull"):
        x = np.log([r["n"] for r in rows])
        y = np.log([max(r[alg], 1e-6) for r in rows])
        if len(rows) >= 3 and x[-1] - x[0] > np.log(4):   # meaningless on a narrow size range
            fits[alg] = float(np.polyfit(x[len(x) // 2 - 1:], y[len(y) // 2 - 1:], 1)[0])
    return {"rows": rows, "exponents": fits, "level": job.points.level}, {}


def _selection(job: Job, p: dict):
    res = np.asarray(p.get("residues") or [], dtype=np.int64)
    if not len(res):
        raise EngineError("select residues, a domain or a ligand first")
    if res.min() < 0 or res.max() >= job.structure.n_residues:
        raise EngineError("selection does not match the loaded structure")
    alpha, me = p.get("alpha"), p.get("max_edge")
    return selgeo.selection_geometry(job, res, None if alpha in (None, "") else float(alpha),
                                     None if me in (None, "") else float(me))


BUILDERS = {"structure": _structure, "delaunay": _delaunay, "voronoi": _voronoi, "alpha": _alpha, "hull": _hull,
            "filtration": _filtration, "pockets": _pockets, "contacts": _contacts, "features": _features,
            "tetrahedra": _tetrahedra, "pae": _pae, "fourbody": _fourbody, "scaling": _scaling,
            "selection": _selection}


# ------------------------------------------------------------------ CSV exports

CSV_KINDS = {
    "residues": "per-residue features (geometry + biochemistry)",
    "tetrahedra": "Delaunay tetrahedra of the point set, in the reference column format",
    "hull_residue": "convex-hull triangles of the residue-level point set",
    "hull_points": "convex-hull triangles of the current point set",
    "pockets": "two-probe pockets (default parameters)",
    "contacts": "residue contacts from Delaunay edges (default cutoff)",
}
# file-name suffix per kind: what the rows are
CSV_SUFFIX = {"residues": "residue-features", "tetrahedra": "delaunay-tetrahedra", "hull_residue": "hull-triangles",
              "hull_points": "hull-triangles", "pockets": "pockets", "contacts": "residue-contacts"}


def export_name(job: Job, content: str, ext: str = "csv", atom_set: str | None = None) -> str:
    """bioGeometry_<id>-<chains>_<atom set>_<content>.<ext>, e.g. bioGeometry_4HHB-ABCD_all-atom_delaunay-tetrahedra.csv
    (atom_set='' leaves that part out)."""
    s = job.structure
    tag = chain_tag(s).replace("/", "+")
    if atom_set is None:
        atom_set = "c-alpha" if job.points.level in ("ca", "centroid") else "all-atom"
    return f"bioGeometry_{s.name}-{tag}_{atom_set + '_' if atom_set else ''}{content}.{ext}"


def _table_rows(columns, data, n):
    header = [c["key"] for c in columns]
    def cell(c, i):
        v = data[c["key"]][i]
        if c["type"] == "cat":
            return c["categories"][int(v)]
        v = float(v)
        return "" if not np.isfinite(v) else f"{v:.4g}"
    yield header
    for i in range(n):
        yield [cell(c, i) for c in columns]


def export_structure(job: Job, fmt: str) -> tuple[str, str]:
    """(filename, text) of the loaded model as mmCIF or PDB."""
    return export_name(job, "model", ext=fmt, atom_set=""), structure_text(job.structure, fmt)


def export_csv(job: Job, kind: str) -> tuple[str, str]:
    """(filename, csv text) for one of CSV_KINDS."""
    s = job.structure
    if kind not in CSV_KINDS:
        raise EngineError(f"unknown CSV kind '{kind}'")
    name = export_name(job, CSV_SUFFIX[kind], atom_set="c-alpha" if kind == "hull_residue" else None)
    buf = io.StringIO()
    w = csv.writer(buf)
    if kind == "residues":
        t = features.residue_table(job)
        w.writerow(["residue"] + [c["key"] for c in t.columns])
        rows = _table_rows(t.columns, t.data, t.n)
        next(rows)
        for i, row in enumerate(rows):
            w.writerow([s.residue_label(i)] + row)
    elif kind == "tetrahedra":
        t, T = features.tetra_table(job, job.defaults()["tetrahedra"]["max_edge"])
        cols = {c["key"]: c for c in t.columns}
        # the reference column order, with the text columns built here
        order = ["pdbid-chainid", "model", "ec-number", "ec-category", "atom-set", "simplex-indices", "residues",
                 "residue-class", "residue-type", "participating-chains"] + [c["key"] for c in t.columns
                 if c["key"] not in ("pdbid-chainid", "model", "ec-number", "ec-category", "atom-set",
                                     "residue-class", "residue-type", "participating-chains")]
        w.writerow(order)
        tets, r4, lab = T["tets"], T["res4"], T["res_label"]
        def cell(key, i):
            if key == "simplex-indices":
                return "-".join(map(str, tets[i]))
            if key == "residues":
                return "-".join(lab[r4[i]])
            c, v = cols[key], t.data[key][i]
            if c["type"] == "cat":
                return c["categories"][int(v)]
            v = float(v)
            return "" if not np.isfinite(v) else (str(int(v)) if c.get("fmt") == "d" else f"{v:.4g}")
        for i in range(t.n):
            w.writerow([cell(k, i) for k in order])
    elif kind in ("hull_residue", "hull_points"):
        pts, idx = ((job.residue_points.coords, job.residue_points.residue) if kind == "hull_residue"
                    else (job.points.coords, job.points.residue))
        h = g.convex_hull(pts)
        bres = an.per_residue(s.bfactor, s.atom_res, s.n_residues, "mean")
        w.writerow(["amino_acid_residues", "residue_numbers", "side_a", "side_b", "side_c", "triangle_area",
                    "tf_1", "tf_2", "tf_3", "avg_tf"])
        for f in h.faces:
            q = idx[f]
            a, b, c = pts[f]
            sides = [np.linalg.norm(a - b), np.linalg.norm(b - c), np.linalg.norm(c - a)]
            area = 0.5 * np.linalg.norm(np.cross(b - a, c - a))
            tf = bres[q]
            w.writerow(["-".join(s.res_name[q].astype(str)), "-".join(map(str, s.res_seq[q]))]
                       + [f"{x:.4g}" for x in sides] + [f"{area:.4g}"] + [f"{x:.4g}" for x in tf] + [f"{tf.mean():.4g}"])
        w.writerow([])
        w.writerow(["# hull_area", f"{h.area:.4g}", "hull_volume", f"{h.volume:.4g}"])
    elif kind == "pockets":
        w.writerow(["pocket", "volume_A3", "mouth_area_A2", "enclosed", "n_tetrahedra", "chains", "lining_residues"])
        for k in job.default_pockets:
            res = np.unique(job.points.residue[k.points])
            w.writerow([k.id + 1, f"{k.volume:.1f}", f"{k.mouth_area:.1f}", int(k.enclosed), k.n_tets,
                        ",".join(sorted(set(s.res_chain[res]))), " ".join(s.residue_label(r) for r in res)])
    elif kind == "contacts":
        d = job.defaults()["contacts"]
        net = an.contact_network(s, job.tess, job.points, d["cutoff"], d["min_seq_sep"])
        w.writerow(["residue_a", "residue_b", "interchain"])
        for a, b in net["edges"]:
            w.writerow([s.residue_label(a), s.residue_label(b), int(s.res_chain[a] != s.res_chain[b])])
    else:
        raise EngineError(f"unknown CSV kind '{kind}'")
    return name, buf.getvalue()
