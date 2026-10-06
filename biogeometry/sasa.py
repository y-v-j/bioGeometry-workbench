"""Solvent-accessible surface area, Shrake & Rupley (1973), vectorised.

Each atom's sphere (van der Waals radius + 1.4 Å water probe) is sampled with a Fibonacci spiral;
a sample point is buried when it lies inside any other expanded sphere. Instead of a Python loop
per atom, all sample points of a chunk of atoms are tested at once against their K nearest atom
centres with one KD-tree query. Points whose K-th neighbour is still close enough to bury them are
re-tested with a larger K, so the result is exact for the chosen sampling.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

PROBE = 1.4
# Bondi radii (Mantina for metals), Å
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "SE": 1.90, "MG": 1.73, "ZN": 1.39,
       "FE": 1.94, "CA": 2.31, "MN": 1.97, "CL": 1.75, "NA": 2.27, "K": 2.75, "CU": 1.40, "CO": 1.92,
       "NI": 1.63, "BR": 1.85, "I": 1.98, "F": 1.47}


def sphere_points(n: int) -> np.ndarray:
    i = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])


def radii_for(elements: list[str]) -> np.ndarray:
    return np.array([VDW.get(str(e).upper(), 1.80) for e in elements], dtype=float)


def atom_sasa(coords: np.ndarray, radii: np.ndarray, n_points: int = 96, probe: float = PROBE,
              k: int = 4, chunk: int = 20000, workers: int = -1) -> np.ndarray:
    """Per-atom SASA in Å² (exact for the sampling up to floating point)."""
    coords = np.asarray(coords, float)
    n = len(coords)
    if n == 0:
        return np.zeros(0)
    R = np.asarray(radii, float) + probe
    rmax = float(R.max())
    Rpad = np.append(R, 0.0)
    unit = sphere_points(n_points)
    tree = cKDTree(coords)
    area = np.zeros(n)
    k = min(k, n)
    for s in range(0, n, chunk):
        idx = np.arange(s, min(n, s + chunk))
        pts = (coords[idx, None, :] + unit[None, :, :] * R[idx, None, None]).reshape(-1, 3)
        owner = np.repeat(idx, n_points)
        # atoms farther than rmax cannot bury a point, so the search is bounded (missing neighbours: d=inf, j=n)
        d, j = tree.query(pts, k=k, workers=workers, distance_upper_bound=rmax)
        if k == 1:
            d, j = d[:, None], j[:, None]
        buried = ((d < Rpad[j]) & (j != owner[:, None])).any(axis=1)
        # points whose K-th neighbour is still within reach may be buried by a farther atom: exact search
        unsure = ~buried & np.isfinite(d[:, -1]) if k < n else np.zeros(len(pts), bool)
        if unsure.any():
            hits = tree.query_ball_point(pts[unsure], r=rmax)
            ui = np.flatnonzero(unsure)
            for q, js in zip(ui, hits):
                js = np.asarray(js)
                js = js[js != owner[q]]
                if len(js) and (np.linalg.norm(coords[js] - pts[q], axis=1) < R[js]).any():
                    buried[q] = True
        exposed = (~buried).reshape(len(idx), n_points).sum(axis=1)
        area[idx] = 4.0 * np.pi * R[idx] ** 2 * exposed / n_points
    return area
