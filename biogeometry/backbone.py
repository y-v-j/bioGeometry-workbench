"""Backbone geometry: dihedrals, secondary structure (Kabsch & Sander) and cartoon meshes.

* ``backbone_atoms``   per-residue N, CA, C, O (protein) and P / C4' (nucleic) coordinates
* ``dihedrals``        φ, ψ, ω in degrees, NaN where a neighbour is missing or not bonded
* ``secondary_structure``  DSSP-equivalent H / E / C from backbone H-bond energies, plus the H-bond list
* ``cartoon_mesh``     PyMOL-style swept ribbons (helix, strand arrows, coil tubes) for proteins and
                       smoothed tubes for nucleic acids, as one vertex buffer with a residue index per vertex

The secondary-structure and ribbon code is ported from mayaviBiomolecules (``mbio.geom.cartoon``) and
adapted to the residue-table layout used here.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import splev, splprep
from scipy.spatial import cKDTree

from .structure import KIND_NUCLEIC, KIND_PROTEIN, Structure

COIL, HELIX, SHEET = 0, 1, 2
SS_CHARS = np.array(["C", "H", "E"])
PROFILES = {COIL: (0.30, 0.30), HELIX: (1.45, 0.28), SHEET: (1.25, 0.25)}
ARROW_HALF_WIDTH, ARROW_TIP_WIDTH = 1.85, 0.04
HBOND_ENERGY_CUTOFF = -0.5          # kcal/mol (DSSP)
_HBOND_Q = 0.084 * 332.0
MAX_PEPTIDE_BOND = 2.0              # Å, C(i)-N(i+1)
NUCLEIC_TUBE_RADIUS = 1.1


# ------------------------------------------------------------------ atoms and dihedrals

def backbone_atoms(s: Structure) -> dict[str, np.ndarray]:
    """{'N','CA','C','O','P','C4'} -> (R, 3) arrays, NaN where the residue lacks the atom."""
    out = {k: np.full((s.n_residues, 3), np.nan) for k in ("N", "CA", "C", "O", "P", "C4'", "C1'")}
    kind = s.res_kind[s.atom_res]
    for name, kinds in (("N", (KIND_PROTEIN,)), ("CA", (KIND_PROTEIN,)), ("C", (KIND_PROTEIN,)),
                        ("O", (KIND_PROTEIN,)), ("P", (KIND_NUCLEIC,)), ("C4'", (KIND_NUCLEIC,)),
                        ("C1'", (KIND_NUCLEIC,))):
        hit = np.flatnonzero((s.atom_name == name) & np.isin(kind, kinds))
        out[name][s.atom_res[hit]] = s.coords[hit]
    return out


def _dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1, axis=1, keepdims=True)
    v = b0 - np.sum(b0 * b1n, axis=1, keepdims=True) * b1n
    w = b2 - np.sum(b2 * b1n, axis=1, keepdims=True) * b1n
    x = np.sum(v * w, axis=1)
    y = np.sum(np.cross(b1n, v) * w, axis=1)
    return np.degrees(np.arctan2(y, x))


def peptide_links(s: Structure, bb: dict) -> np.ndarray:
    """(R,) bool: residue i is peptide-bonded to residue i+1 (same chain, C–N ≤ 2 Å)."""
    link = np.zeros(s.n_residues, dtype=bool)
    if s.n_residues < 2:
        return link
    same = (s.res_chain[:-1] == s.res_chain[1:]) & (s.res_kind[:-1] == KIND_PROTEIN) & (s.res_kind[1:] == KIND_PROTEIN)
    with np.errstate(invalid="ignore"):
        d = np.linalg.norm(bb["C"][:-1] - bb["N"][1:], axis=1)
    link[:-1] = same & (d <= MAX_PEPTIDE_BOND)
    return link


def dihedrals(s: Structure, bb: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """φ, ψ, ω (degrees) per residue; ω is the peptide bond preceding residue i (≈180 trans, ≈0 cis)."""
    R = s.n_residues
    phi = np.full(R, np.nan)
    psi = np.full(R, np.nan)
    omega = np.full(R, np.nan)
    link = peptide_links(s, bb)
    i = np.flatnonzero(link[:-1]) if R > 1 else np.array([], int)   # bond i -> i+1
    if len(i):
        with np.errstate(invalid="ignore"):
            phi[i + 1] = _dihedral(bb["C"][i], bb["N"][i + 1], bb["CA"][i + 1], bb["C"][i + 1])
            psi[i] = _dihedral(bb["N"][i], bb["CA"][i], bb["C"][i], bb["N"][i + 1])
            omega[i + 1] = _dihedral(bb["CA"][i], bb["C"][i], bb["N"][i + 1], bb["CA"][i + 1])
    return phi, psi, omega


# ------------------------------------------------------------------ segments

def protein_segments(s: Structure, bb: dict, min_length: int = 3) -> list[np.ndarray]:
    """Runs of consecutive, peptide-bonded protein residues that have N, CA, C and O."""
    full = (s.res_kind == KIND_PROTEIN) & np.all(np.isfinite(np.concatenate(
        [bb["N"], bb["CA"], bb["C"], bb["O"]], axis=1)), axis=1)
    link = peptide_links(s, bb)
    segs, run = [], []
    for r in range(s.n_residues):
        if not full[r]:
            if len(run) >= min_length:
                segs.append(np.array(run))
            run = []
            continue
        if run and not link[run[-1]]:
            if len(run) >= min_length:
                segs.append(np.array(run))
            run = []
        run.append(r)
    if len(run) >= min_length:
        segs.append(np.array(run))
    return segs


# ------------------------------------------------------------------ secondary structure (Kabsch & Sander)

def _hbond_pairs(n, h, c, o, seg_id, radius=6.0):
    neighbours = cKDTree(o).query_ball_point(n, r=radius)
    counts = np.fromiter((len(x) for x in neighbours), dtype=int, count=len(neighbours))
    donor = np.repeat(np.arange(len(n)), counts)
    acceptor = np.fromiter((j for hit in neighbours for j in hit), dtype=int, count=int(counts.sum()))
    if donor.size == 0:
        return np.empty((0, 2), int), np.empty(0)
    keep = (np.abs(acceptor - donor) >= 2) & (seg_id[acceptor] == seg_id[donor])
    acceptor, donor = acceptor[keep], donor[keep]
    if acceptor.size == 0:
        return np.empty((0, 2), int), np.empty(0)

    def inv(a, b):
        with np.errstate(divide="ignore", invalid="ignore"):
            return 1.0 / np.linalg.norm(a - b, axis=1)

    e = _HBOND_Q * (inv(o[acceptor], n[donor]) + inv(c[acceptor], h[donor])
                    - inv(o[acceptor], h[donor]) - inv(c[acceptor], n[donor]))
    e = np.nan_to_num(e, nan=0.0, posinf=0.0, neginf=0.0)
    hit = e < HBOND_ENERGY_CUTOFF
    return np.column_stack([acceptor[hit], donor[hit]]), e[hit]


class _Bonded:
    def __init__(self, pairs, total):
        self.total = total
        self.keys = np.sort(pairs[:, 0].astype(np.int64) * total + pairs[:, 1]) if len(pairs) else np.empty(0, np.int64)

    def __call__(self, i, j):
        i, j = np.asarray(i), np.asarray(j)
        ok = (i >= 0) & (i < self.total) & (j >= 0) & (j < self.total)
        if self.keys.size == 0:
            return np.zeros(np.broadcast(i, j).shape, dtype=bool)
        key = np.where(ok, i.astype(np.int64) * self.total + j, -1)
        pos = np.clip(np.searchsorted(self.keys, key), 0, self.keys.size - 1)
        return ok & (self.keys[pos] == key)


def _runs(values, seg_id):
    edges = np.flatnonzero((np.diff(values) != 0) | (np.diff(seg_id) != 0)) + 1
    bounds = np.concatenate([[0], edges, [len(values)]])
    for a, b in zip(bounds[:-1], bounds[1:]):
        yield a, b, values[a]


def secondary_structure(s: Structure, bb: dict, segs: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """(ss_code per residue: 0 coil / 1 helix / 2 sheet, -1 for non-protein;
        H-bond partners as (acceptor_residue, donor_residue) global indices)."""
    ss_out = np.full(s.n_residues, -1, dtype=np.int8)
    if not segs:
        return ss_out, np.empty((0, 2), int)
    idx = np.concatenate(segs)
    lengths = [len(x) for x in segs]
    seg_id = np.repeat(np.arange(len(segs)), lengths)
    n, ca, c, o = bb["N"][idx], bb["CA"][idx], bb["C"][idx], bb["O"][idx]
    starts = np.cumsum([0] + lengths)[:-1]
    h = np.full_like(n, np.nan)
    co = c[:-1] - o[:-1]
    co /= np.linalg.norm(co, axis=1, keepdims=True)
    h[1:] = n[1:] + co
    h[starts] = np.nan
    total = len(idx)
    pairs, _ = _hbond_pairs(n, h, c, o, seg_id)
    bonded = _Bonded(pairs, total)
    ss = np.full(total, COIL, dtype=int)
    k = np.arange(total)
    for turn in (4, 3, 5):
        t = bonded(k, k + turn)
        run = t[:-1] & t[1:]
        for i in np.flatnonzero(run):
            ss[i + 1:i + turn + 1] = HELIX
    if len(pairs):
        shifts = (-1, 0, 1)
        cand = np.concatenate([np.column_stack([pairs[:, 0] + a, pairs[:, 1] + b]) for a in shifts for b in shifts])
        cand = np.unique(np.concatenate([cand, cand[:, ::-1]]), axis=0)
        bi, bj = cand[:, 0], cand[:, 1]
        inside = (bi >= 0) & (bi < total) & (bj >= 0) & (bj < total) & (np.abs(bi - bj) > 2)
        bi, bj = bi[inside], bj[inside]
        parallel = (bonded(bi - 1, bj) & bonded(bj, bi + 1)) | (bonded(bj - 1, bi) & bonded(bi, bj + 1))
        anti = (bonded(bi, bj) & bonded(bj, bi)) | (bonded(bi - 1, bj + 1) & bonded(bj - 1, bi + 1))
        partnered = np.zeros(total, dtype=bool)
        hit = parallel | anti
        partnered[np.unique(np.concatenate([bi[hit], bj[hit]]))] = True
        ss[partnered & (ss != HELIX)] = SHEET
    # tidy: drop runs too short to be real, fill single-residue holes
    out = ss.copy()
    for a, b, v in _runs(ss, seg_id):
        if v in (HELIX, SHEET) and b - a < {HELIX: 4, SHEET: 2}[v]:
            out[a:b] = COIL
    ss = out.copy()
    for a, b, v in _runs(out, seg_id):
        if v == COIL and b - a == 1 and 0 < a and b < total and out[a - 1] == out[b] != COIL \
                and seg_id[a - 1] == seg_id[b]:
            ss[a:b] = out[a - 1]
    ss_out[idx] = ss
    return ss_out, idx[pairs] if len(pairs) else np.empty((0, 2), int)


# ------------------------------------------------------------------ cartoon mesh

def _normalize(v):
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.where(n == 0, 1.0, n)


def _ribbon(ca, c, o, ss, samples, sides):
    """Swept elliptical cross-section along one protein segment -> (points (rings, sides+1, 3), u)."""
    m = len(ca)
    if m < 4:
        return None
    ca = ca.copy()
    for _ in range(2):  # flatten strands (PyMOL cartoon_flat_sheets)
        sm = ca.copy()
        sm[1:-1] = 0.25 * ca[:-2] + 0.5 * ca[1:-1] + 0.25 * ca[2:]
        inner = ss == SHEET
        inner[0] = inner[-1] = False
        ca[inner] = sm[inner]
    tangent = np.empty_like(ca)
    tangent[1:-1] = ca[2:] - ca[:-2]
    tangent[0], tangent[-1] = ca[1] - ca[0], ca[-1] - ca[-2]
    tangent = _normalize(tangent)
    carbonyl = _normalize(o - c)
    side = carbonyl - np.sum(carbonyl * tangent, axis=1, keepdims=True) * tangent
    weak = np.linalg.norm(side, axis=1) < 1e-3
    if weak.any():
        side[weak] = np.cross(np.roll(ca, 1, axis=0) - ca, np.roll(ca, -1, axis=0) - ca)[weak]
    side = _normalize(side)
    for i in range(1, m):
        if np.dot(side[i], side[i - 1]) < 0:
            side[i] = -side[i]
    keep = np.concatenate([[True], np.linalg.norm(np.diff(ca, axis=0), axis=1) > 1e-6])
    if keep.sum() < 4:
        return None
    try:
        tck, _ = splprep(ca[keep].T, s=0, k=3)
    except Exception:  # noqa: BLE001 - degenerate geometry
        return None
    rings = max(int((m - 1) * samples) + 1, 4)
    t = np.linspace(0.0, 1.0, rings)
    centres = np.asarray(splev(t, tck)).T
    tan = _normalize(np.asarray(splev(t, tck, der=1)).T)
    u = t * (m - 1)
    grid = np.arange(m)
    ring_side = _normalize(np.column_stack([np.interp(u, grid, side[:, k]) for k in range(3)]))
    curv = np.asarray(splev(t, tck, der=2)).T
    radial = curv - np.sum(curv * tan, axis=1, keepdims=True) * tan
    helical = _normalize(np.cross(_normalize(radial), tan))
    wgt = np.interp(u, grid, (ss == HELIX).astype(float))[:, None]
    opposed = np.sum(helical * ring_side, axis=1) < 0
    helical[opposed] *= -1.0
    ring_side = _normalize((1.0 - wgt) * ring_side + wgt * helical)
    ring_side = _normalize(ring_side - np.sum(ring_side * tan, axis=1, keepdims=True) * tan)
    ring_normal = _normalize(np.cross(tan, ring_side))
    width = np.interp(u, grid, np.array([PROFILES[int(x)][0] for x in ss]))
    thick = np.interp(u, grid, np.array([PROFILES[int(x)][1] for x in ss]))
    for a, b, v in _runs(ss, np.zeros(m, int)):
        if v != SHEET:
            continue
        head, tip = b - 2.0, b - 1.0
        inside = (u >= head) & (u <= tip)
        if inside.any():
            frac = (u[inside] - head) / max(tip - head, 1e-9)
            width[inside] = ARROW_HALF_WIDTH + frac * (ARROW_TIP_WIDTH - ARROW_HALF_WIDTH)
    theta = np.linspace(0.0, 2.0 * np.pi, sides + 1)
    pts = (centres[:, None, :]
           + (width[:, None] * np.cos(theta)[None, :])[:, :, None] * ring_side[:, None, :]
           + (thick[:, None] * np.sin(theta)[None, :])[:, :, None] * ring_normal[:, None, :])
    pts = np.concatenate([np.repeat(centres[:1][:, None, :], sides + 1, axis=1), pts,
                          np.repeat(centres[-1:][:, None, :], sides + 1, axis=1)])
    u = np.concatenate([u[:1], u, u[-1:]])
    return pts, u


def _tube(points, samples, sides, radius):
    """Smooth circular tube through points (nucleic backbone)."""
    m = len(points)
    keep = np.concatenate([[True], np.linalg.norm(np.diff(points, axis=0), axis=1) > 1e-6])
    p = points[keep]
    if len(p) < 2:
        return None
    k = min(3, len(p) - 1)
    try:
        tck, _ = splprep(p.T, s=0, k=k)
    except Exception:  # noqa: BLE001
        return None
    rings = max(int((m - 1) * samples) + 1, 4)
    t = np.linspace(0, 1, rings)
    centres = np.asarray(splev(t, tck)).T
    tan = _normalize(np.asarray(splev(t, tck, der=1)).T)
    ref = np.where(np.abs(tan[:, :1]) < 0.9, np.array([[1.0, 0, 0]]), np.array([[0, 1.0, 0]]))
    side = _normalize(np.cross(tan, ref))
    for i in range(1, rings):  # parallel transport to avoid twisting
        s_i = side[i - 1] - np.dot(side[i - 1], tan[i]) * tan[i]
        n = np.linalg.norm(s_i)
        side[i] = s_i / n if n > 1e-6 else side[i]
    normal = _normalize(np.cross(tan, side))
    theta = np.linspace(0, 2 * np.pi, sides + 1)
    pts = (centres[:, None, :] + radius * (np.cos(theta)[None, :, None] * side[:, None, :]
                                           + np.sin(theta)[None, :, None] * normal[:, None, :]))
    pts = np.concatenate([np.repeat(centres[:1][:, None, :], sides + 1, axis=1), pts,
                          np.repeat(centres[-1:][:, None, :], sides + 1, axis=1)])
    u = np.concatenate([[0.0], t * (m - 1), [m - 1.0]])
    return pts, u


def mesh_budget(n_residues: int) -> tuple[int, int]:
    if n_residues <= 3_000:
        return 8, 12
    if n_residues <= 10_000:
        return 5, 10
    return 3, 8


def cartoon_mesh(s: Structure, bb: dict, ss: np.ndarray, segs: list[np.ndarray]) -> dict:
    """One vertex buffer for every ribbon/tube; JS rebuilds the triangle strips from ``rings``.

    Returns {"pos": (V,3) f32, "res": (V,) u32, "rings": (S,) u32, "sides": int}.
    """
    n_poly = int(sum(len(x) for x in segs)) + int((s.res_kind == KIND_NUCLEIC).sum())
    samples, sides = mesh_budget(n_poly)
    pos, res, rings = [], [], []
    for seg in segs:
        out = _ribbon(bb["CA"][seg], bb["C"][seg], bb["O"][seg], np.maximum(ss[seg], 0), samples, sides)
        if out is None:
            continue
        pts, u = out
        r_idx = seg[np.clip(np.rint(u).astype(int), 0, len(seg) - 1)]
        pos.append(pts.reshape(-1, 3))
        res.append(np.repeat(r_idx, sides + 1))
        rings.append(pts.shape[0])
    # nucleic acid tubes through P (fallback C4')
    nuc = np.flatnonzero(s.res_kind == KIND_NUCLEIC)
    if len(nuc):
        anchor = np.where(np.isfinite(bb["P"][:, :1]), bb["P"], bb["C4'"])
        run = []
        for r in list(nuc) + [-1]:
            ok = r >= 0 and np.isfinite(anchor[r]).all()
            if ok and run and (s.res_chain[r] != s.res_chain[run[-1]]
                               or np.linalg.norm(anchor[r] - anchor[run[-1]]) > 8.5):
                ok_break = True
            else:
                ok_break = False
            if (not ok or ok_break) and len(run) >= 2:
                out = _tube(anchor[run], samples, max(sides - 2, 6), NUCLEIC_TUBE_RADIUS)
                if out is not None:
                    pts, u = out
                    seg = np.array(run)
                    pos.append(pts.reshape(-1, 3))
                    res.append(np.repeat(seg[np.clip(np.rint(u).astype(int), 0, len(seg) - 1)], pts.shape[1]))
                    rings.append(-pts.shape[0])  # negative = tube (different side count)
            if not ok or ok_break:
                run = [r] if ok else []
            else:
                run.append(r)
    if not pos:
        return {"pos": np.zeros((0, 3), np.float32), "res": np.zeros(0, np.uint32),
                "rings": np.zeros(0, np.int32), "sides": sides}
    return {"pos": np.concatenate(pos).astype(np.float32), "res": np.concatenate(res).astype(np.uint32),
            "rings": np.asarray(rings, dtype=np.int32), "sides": sides, "tube_sides": max(sides - 2, 6)}


def base_ladder(s: Structure, bb: dict) -> np.ndarray:
    """(B, 2) atom index pairs C1' -> N1/N9 for drawing nucleobase 'rungs' on the tube."""
    nuc = s.res_kind[s.atom_res] == KIND_NUCLEIC
    if not nuc.any():
        return np.zeros((0, 2), np.int64)
    out = []
    names = s.atom_name
    first = {}
    for i in np.flatnonzero(nuc):
        r = s.atom_res[i]
        first.setdefault(r, {})[names[i]] = i
    for r, d in first.items():
        a = d.get("C1'")
        b = d.get("N9") or d.get("N1")
        tip = d.get("N1") if d.get("N9") is not None else d.get("N3")
        if a is not None and b is not None:
            out.append((a, b))
            if tip is not None and tip != b:
                out.append((b, tip))
    return np.asarray(out, dtype=np.int64).reshape(-1, 2)
