"""Residue-level analyses built on the geometry primitives.

* contact network   - Delaunay edges (<= cutoff) lifted to residues; degree, chain-chain matrix,
                      interface residues. Delaunay contacts are a parameter-light alternative
                      to distance-cutoff contacts: only *non-occluded* neighbours are connected.
* four-body stats   - residue composition of Delaunay tetrahedra (Tropsha/Vaisman "SNAPP" style
                      log-likelihood of each quadruplet vs a multinomial expectation).
* residue profiles  - per-residue Voronoi volume, hull depth, convex layer, contact degree.
* scaling benchmark - wall time of each algorithm vs number of points (the point of this app:
                      how far does each method go on GroEL / ribosome-sized inputs).
"""
from __future__ import annotations

import time
from collections import Counter
from math import lgamma, log

import numpy as np

from . import geometry as g
from .structure import PointSet, Structure

AA3 = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
       "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]


def residue_edges(tess: g.Tess, ps: PointSet, cutoff: float) -> np.ndarray:
    """Unique residue pairs (r1<r2) joined by a Delaunay edge no longer than cutoff."""
    e = tess.edges
    keep = tess.edge_lengths() <= cutoff
    r = ps.residue[e[keep]]
    r = r[r[:, 0] != r[:, 1]]
    r.sort(axis=1)
    return np.unique(r, axis=0) if len(r) else np.zeros((0, 2), dtype=np.int64)


def contact_network(s: Structure, tess: g.Tess, ps: PointSet, cutoff: float,
                    min_seq_sep: int = 2) -> dict:
    re = residue_edges(tess, ps, cutoff)
    ch = s.chain_index()
    same = ch[re[:, 0]] == ch[re[:, 1]]
    sep = np.abs(s.res_seq[re[:, 0]] - s.res_seq[re[:, 1]])
    long_range = ~same | (sep >= min_seq_sep)
    re = re[long_range]
    same = same[long_range]
    nres = s.n_residues
    degree = np.bincount(re.ravel(), minlength=nres)
    nch = len(s.chains)
    mat = np.zeros((nch, nch), dtype=np.int64)
    inter = re[~same]
    np.add.at(mat, (ch[inter[:, 0]], ch[inter[:, 1]]), 1)
    mat = mat + mat.T
    interface = np.zeros(nres, dtype=bool)
    interface[inter.ravel()] = True
    pairs = [
        {"a": s.chains[i], "b": s.chains[j], "contacts": int(mat[i, j]),
         "residues": int(len(np.unique(inter[((ch[inter[:, 0]] == i) & (ch[inter[:, 1]] == j))
                                             | ((ch[inter[:, 0]] == j) & (ch[inter[:, 1]] == i))])))}
        for i in range(nch) for j in range(i + 1, nch) if mat[i, j] > 0
    ]
    pairs.sort(key=lambda d: -d["contacts"])
    return {
        "cutoff": cutoff,
        "n_contacts": int(len(re)),
        "n_interchain": int(len(inter)),
        "mean_degree": float(degree.mean()) if nres else 0.0,
        "degree": degree,
        "interface": interface,
        "chain_matrix": mat,
        "pairs": pairs,
        "edges": re,
        "interchain_edges": inter,
    }


def four_body_statistics(s: Structure, tess: g.Tess, ps: PointSet, max_edge: float,
                         top: int = 25) -> dict:
    """Composition of residue quadruplets in compact Delaunay tetrahedra.

    Only tetrahedra whose four points belong to four *different* standard amino acids and
    whose longest edge is <= max_edge are counted. The log-likelihood of quadruplet q is
    log(f_obs(q) / p_exp(q)), with p_exp the multinomial probability of drawing that
    composition from the residue-type frequencies in the same tetrahedra.
    """
    res = ps.residue[tess.simplices]
    distinct = np.ones(len(res), dtype=bool)
    for i in range(4):
        for j in range(i + 1, 4):
            distinct &= res[:, i] != res[:, j]
    ok = distinct & (tess.max_edge_per_tet() <= max_edge)
    names = s.res_name[res[ok]]
    std = np.isin(names, AA3).all(axis=1)
    names = np.sort(names[std].astype(str), axis=1)
    n = len(names)
    if n == 0:
        return {"n_tets": 0, "top": [], "max_edge": max_edge}
    quads = Counter(map(tuple, names))
    freq = Counter(names.ravel())
    total = 4 * n
    a = {k: v / total for k, v in freq.items()}
    rows = []
    for q, cnt in quads.items():
        comp = Counter(q)
        log_mult = lgamma(5) - sum(lgamma(c + 1) for c in comp.values())
        log_p = log_mult + sum(c * log(a[t]) for t, c in comp.items())
        rows.append({"quad": "-".join(q), "count": cnt, "log_likelihood": log(cnt / n) - log_p})
    rows.sort(key=lambda r: -r["count"])
    common = [r for r in rows if r["count"] >= 5]
    common.sort(key=lambda r: -r["log_likelihood"])
    return {
        "n_tets": int(n), "max_edge": max_edge, "n_distinct": len(quads),
        "most_common": rows[:top],
        "most_enriched": common[:top],
        "most_depleted": common[::-1][:top],
    }


def per_residue(values: np.ndarray, residue: np.ndarray, n_res: int, how: str = "sum") -> np.ndarray:
    """Aggregate a per-point quantity to residues (NaN-aware)."""
    v = np.asarray(values, dtype=float)
    ok = np.isfinite(v)
    if how == "sum":
        tot = np.bincount(residue[ok], weights=v[ok], minlength=n_res)
        any_nan = np.bincount(residue[~ok], minlength=n_res) > 0
        tot[any_nan] = np.nan
        return tot
    if how == "mean":
        tot = np.bincount(residue[ok], weights=v[ok], minlength=n_res)
        cnt = np.bincount(residue[ok], minlength=n_res)
        with np.errstate(invalid="ignore"):
            return np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan)
    if how == "min":
        out = np.full(n_res, np.inf)
        np.minimum.at(out, residue[ok], v[ok])
        out[np.isinf(out)] = np.nan
        return out
    raise ValueError(how)


def scaling_benchmark(points: np.ndarray, sizes: list[int] | None = None, max_radius: float = 8.0,
                      alpha: float = 8.0, seed: int = 0) -> list[dict]:
    """Time Delaunay / alpha / Voronoi / hull on random subsets of increasing size."""
    n = len(points)
    if sizes is None:
        sizes = sorted({int(x) for x in np.geomspace(min(500, n), n, 7)})
    rng = np.random.default_rng(seed)
    rows = []
    for m in sizes:
        sub = points[rng.choice(n, size=m, replace=False)] if m < n else points
        row = {"n": int(m)}
        t = time.perf_counter(); tess = g.Tess.build(sub); _ = tess.circumradii
        row["delaunay"] = time.perf_counter() - t
        t = time.perf_counter(); g.alpha_shape(tess, alpha); row["alpha"] = time.perf_counter() - t
        t = time.perf_counter(); g.voronoi_cells(sub, max_radius); row["voronoi"] = time.perf_counter() - t
        t = time.perf_counter(); g.convex_hull(sub); row["hull"] = time.perf_counter() - t
        row["tets"] = int(tess.n_tets)
        rows.append(row)
    return rows


def spearman(metrics: dict, keys: list[str]) -> list[dict]:
    """Spearman ρ for every pair of metrics, with tied values given their average rank."""
    from scipy.stats import rankdata
    out = []
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            x, y = np.asarray(metrics[a], float), np.asarray(metrics[b], float)
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() > 3 and np.ptp(x[ok]) > 0 and np.ptp(y[ok]) > 0:
                rx, ry = rankdata(x[ok]), rankdata(y[ok])
                out.append({"a": a, "b": b, "spearman": float(np.corrcoef(rx, ry)[0, 1]), "n": int(ok.sum())})
    return out
