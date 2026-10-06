"""Computational-geometry primitives on 3-D point sets (numpy + scipy/Qhull).

Everything here is vectorised so that it scales to ribosome-sized inputs
(~10^5 points, ~7x10^5 tetrahedra). Conventions:

* A Delaunay tessellation is used through :class:`Tess`, a thin cache around
  ``scipy.spatial.Delaunay`` holding per-tetrahedron volume and circumradius.
* Alpha shapes are the *regularised* alpha shape: the union of Delaunay
  tetrahedra whose circumradius is <= alpha (the interior of the alpha complex).
  Its boundary is the set of triangles shared by one kept and one rejected (or
  missing) tetrahedron. As alpha -> infinity this is the convex hull.
* Voronoi cell volumes are computed from the Voronoi ridges (each ridge is a
  polygon on the perpendicular bisector of two sites, so it contributes a pyramid
  of height d/2 to each of them). Cells that are unbounded, or that reach further
  than ``max_radius`` from their site, are flagged as *open* and get NaN volume:
  such surface cells extend into solvent and their volume is meaningless.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from itertools import chain as ichain

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import ConvexHull, Delaunay, QhullError, Voronoi

# face k of a tetrahedron is the triangle opposite vertex k
TET_FACES = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]])
TET_EDGES = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])


class GeometryError(ValueError):
    pass


# ------------------------------------------------------------------ Delaunay

@dataclass
class Tess:
    points: np.ndarray       # (P,3)
    simplices: np.ndarray    # (T,4) int
    neighbors: np.ndarray    # (T,4) int, -1 = outside the hull

    @classmethod
    def build(cls, points: np.ndarray) -> "Tess":
        pts = np.asarray(points, dtype=np.float64)
        if len(pts) < 5:
            raise GeometryError(f"need at least 5 points for a 3-D tessellation (got {len(pts)})")
        try:
            dt = Delaunay(pts)
        except QhullError as exc:
            raise GeometryError(f"Qhull failed (points coplanar?): {str(exc).splitlines()[0]}") from exc
        return cls(pts, dt.simplices.astype(np.int64), dt.neighbors.astype(np.int64))

    @property
    def n_tets(self) -> int:
        return len(self.simplices)

    @cached_property
    def volumes(self) -> np.ndarray:
        a, b, c, d = (self.points[self.simplices[:, k]] for k in range(4))
        return np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))) / 6.0

    @cached_property
    def circumradii(self) -> np.ndarray:
        return tet_circumradius(self.points[self.simplices])

    @cached_property
    def edges(self) -> np.ndarray:
        """Unique Delaunay edges (E,2), i<j."""
        return unique_edges(self.simplices)

    def edge_lengths(self, edges: np.ndarray | None = None) -> np.ndarray:
        e = self.edges if edges is None else edges
        return np.linalg.norm(self.points[e[:, 0]] - self.points[e[:, 1]], axis=1)

    def max_edge_per_tet(self) -> np.ndarray:
        p = self.points[self.simplices]
        lens = np.stack([np.linalg.norm(p[:, i] - p[:, j], axis=1) for i, j in TET_EDGES], axis=1)
        return lens.max(axis=1)


def tet_circumradius(tets: np.ndarray) -> np.ndarray:
    """Circumradius of tetrahedra given as (T,4,3) coordinates. Degenerate -> inf."""
    a = tets[:, 0]
    u, v, w = tets[:, 1] - a, tets[:, 2] - a, tets[:, 3] - a
    vw, wu, uv = np.cross(v, w), np.cross(w, u), np.cross(u, v)
    det = np.einsum("ij,ij->i", u, vw)
    num = (np.einsum("ij,ij->i", u, u)[:, None] * vw
           + np.einsum("ij,ij->i", v, v)[:, None] * wu
           + np.einsum("ij,ij->i", w, w)[:, None] * uv)
    scale = np.maximum(np.abs(u).max(axis=1), 1e-12) ** 3
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.linalg.norm(num, axis=1) / (2.0 * np.abs(det))
    r[np.abs(det) < 1e-10 * scale] = np.inf
    return r


def unique_edges(simplices: np.ndarray) -> np.ndarray:
    e = simplices[:, TET_EDGES].reshape(-1, 2)
    e.sort(axis=1)
    n = int(simplices.max()) + 1
    key = np.unique(e[:, 0] * n + e[:, 1])
    return np.stack([key // n, key % n], axis=1)


def tet_components(neighbors: np.ndarray, mask: np.ndarray) -> tuple[int, np.ndarray]:
    """Connected components of the tetrahedra selected by ``mask`` (face adjacency).

    Returns (n_components, labels) with labels = -1 for unselected tetrahedra.
    """
    idx = np.nonzero(mask)[0]
    labels = np.full(len(mask), -1, dtype=np.int64)
    if len(idx) == 0:
        return 0, labels
    local = np.full(len(mask), -1, dtype=np.int64)
    local[idx] = np.arange(len(idx))
    nb = neighbors[idx]
    rows = np.repeat(np.arange(len(idx)), 4)
    cols = local[np.where(nb >= 0, nb, 0)].ravel()
    ok = (nb.ravel() >= 0) & (cols >= 0)
    g = coo_matrix((np.ones(ok.sum(), dtype=np.int8), (rows[ok], cols[ok])), shape=(len(idx), len(idx)))
    n, lab = connected_components(g, directed=False)
    labels[idx] = lab
    return n, labels


def _boundary_faces(tess: Tess, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Triangles separating masked tets from unmasked/outside, oriented outward.

    Returns (faces (F,3), owning tet index (F,)).
    """
    nb = tess.neighbors
    outside = nb < 0
    nb_in = np.where(outside, False, mask[np.where(outside, 0, nb)])
    is_bnd = mask[:, None] & ~nb_in
    t_idx, k_idx = np.nonzero(is_bnd)
    faces = tess.simplices[t_idx[:, None], TET_FACES[k_idx]]
    # orient so the normal points away from the opposite vertex
    p = tess.points
    a, b, c = p[faces[:, 0]], p[faces[:, 1]], p[faces[:, 2]]
    n = np.cross(b - a, c - a)
    opp = p[tess.simplices[t_idx, k_idx]]
    flip = np.einsum("ij,ij->i", n, opp - a) > 0
    faces[flip] = faces[flip][:, [0, 2, 1]]
    return faces, t_idx


def triangle_areas(points: np.ndarray, faces: np.ndarray) -> np.ndarray:
    a, b, c = points[faces[:, 0]], points[faces[:, 1]], points[faces[:, 2]]
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)


# ------------------------------------------------------------------ alpha shape

@dataclass
class AlphaShape:
    alpha: float
    keep: np.ndarray          # (T,) bool  tetrahedra in the shape
    faces: np.ndarray         # (F,3) boundary triangles, outward oriented
    face_component: np.ndarray  # (F,) component id of the owning tet
    volume: float
    area: float
    n_components: int
    n_cavities: int           # enclosed voids (empty regions not reaching the outside)
    surface_points: np.ndarray  # (S,) indices of points on the boundary
    isolated_points: np.ndarray  # points in no kept tetrahedron


def alpha_shape(tess: Tess, alpha: float) -> AlphaShape:
    keep = tess.circumradii <= alpha
    faces, owner = _boundary_faces(tess, keep)
    ncomp, comp = tet_components(tess.neighbors, keep)
    ncav = count_cavities(tess, keep)
    in_shape = np.zeros(len(tess.points), dtype=bool)
    in_shape[tess.simplices[keep].ravel()] = True
    return AlphaShape(
        alpha=float(alpha), keep=keep, faces=faces, face_component=comp[owner],
        volume=float(tess.volumes[keep].sum()),
        area=float(triangle_areas(tess.points, faces).sum()) if len(faces) else 0.0,
        n_components=int(ncomp), n_cavities=ncav,
        surface_points=np.unique(faces) if len(faces) else np.zeros(0, dtype=np.int64),
        isolated_points=np.nonzero(~in_shape)[0],
    )


def count_cavities(tess: Tess, keep: np.ndarray) -> int:
    """Number of connected empty regions (within the hull) that do not touch the hull boundary."""
    empty = ~keep
    n, lab = tet_components(tess.neighbors, empty)
    if n == 0:
        return 0
    touches = np.zeros(n, dtype=bool)
    on_hull = (tess.neighbors < 0).any(axis=1) & empty
    touches[lab[on_hull]] = True
    return int((~touches).sum())


def alpha_filtration(tess: Tess, alphas: np.ndarray) -> dict[str, list[float]]:
    """Volume, area, #components and #cavities of the alpha shape over an alpha grid."""
    out = {"alpha": [], "volume": [], "area": [], "components": [], "cavities": [], "surface_fraction": []}
    order = np.argsort(tess.circumradii)
    r_sorted = tess.circumradii[order]
    cumvol = np.concatenate([[0.0], np.cumsum(tess.volumes[order])])
    n_pts = len(tess.points)
    for a in alphas:
        keep = tess.circumradii <= a
        faces, _ = _boundary_faces(tess, keep)
        ncomp, _ = tet_components(tess.neighbors, keep)
        out["alpha"].append(float(a))
        out["volume"].append(float(cumvol[np.searchsorted(r_sorted, a, side="right")]))
        out["area"].append(float(triangle_areas(tess.points, faces).sum()) if len(faces) else 0.0)
        out["components"].append(int(ncomp))
        out["cavities"].append(count_cavities(tess, keep))
        out["surface_fraction"].append(float(len(np.unique(faces)) / n_pts) if len(faces) else 0.0)
    return out


@dataclass
class Pocket:
    id: int
    volume: float
    n_tets: int
    mouth_area: float          # area of the opening toward the large-probe exterior
    enclosed: bool             # True = cavity with no opening
    points: np.ndarray         # lining point indices
    faces: np.ndarray          # boundary triangles of the pocket region


def find_pockets(tess: Tess, alpha_small: float, alpha_large: float, min_volume: float = 0.0,
                 max_pockets: int = 50) -> list[Pocket]:
    """Two-probe pocket detection (Kawabata & Go 2007 style) on the Delaunay complex.

    Pocket space = tetrahedra *outside* the small-alpha shape but *inside* the
    large-alpha shape: regions a large probe cannot reach but a small one can.
    Each connected region is one pocket; its mouth is the part of its boundary
    facing tetrahedra outside the large shape (or the hull exterior). Pockets
    without a mouth are enclosed cavities.
    """
    if alpha_large <= alpha_small:
        raise GeometryError("alpha_large must be larger than alpha_small")
    R = tess.circumradii
    space = (R > alpha_small) & (R <= alpha_large)
    n, lab = tet_components(tess.neighbors, space)
    if n == 0:
        return []
    vol = np.bincount(lab[space], weights=tess.volumes[space], minlength=n)
    order = np.argsort(-vol)
    order = order[vol[order] >= min_volume][:max_pockets]
    large = R <= alpha_large
    pockets = []
    for new_id, c in enumerate(order):
        member = lab == c
        faces, owner = _boundary_faces(tess, member)
        # classify each boundary face: mouth if the far side is outside the large shape
        t_idx, k_idx = np.nonzero(member[:, None] & ~_nb_mask(tess, member))
        far = tess.neighbors[t_idx, k_idx]
        mouth = (far < 0) | ~large[np.maximum(far, 0)]
        mouth_area = float(triangle_areas(tess.points, faces[mouth]).sum()) if mouth.any() else 0.0
        pockets.append(Pocket(
            id=new_id, volume=float(vol[c]), n_tets=int(member.sum()), mouth_area=mouth_area,
            enclosed=not bool(mouth.any()),
            points=np.unique(tess.simplices[member]), faces=faces[~mouth] if mouth.any() else faces,
        ))
    return pockets


def _nb_mask(tess: Tess, mask: np.ndarray) -> np.ndarray:
    nb = tess.neighbors
    return np.where(nb < 0, False, mask[np.where(nb < 0, 0, nb)])


# ------------------------------------------------------------------ Voronoi

@dataclass
class VoronoiCells:
    vor: Voronoi              # built on the points plus their mirror images when a box is used
    volume: np.ndarray        # (P,) cell volume: every cell closed by the box; NaN for open cells without one
    area: np.ndarray          # (P,) cell surface area, same convention
    bounded: np.ndarray       # (P,) bool  'interior': closed and no vertex farther than max_radius (Richards test)
    n_faces: np.ndarray       # (P,) int number of ridges (= Delaunay degree)
    # flattened ridge polygons, triangulated as fans
    tri: np.ndarray           # (Q,3) indices into vor.vertices
    tri_ridge: np.ndarray     # (Q,) ridge index of each triangle
    ridge_points: np.ndarray  # (Rg,2) may index mirror images (>= n_real)
    ridge_finite: np.ndarray  # (Rg,) bool
    closed: np.ndarray = None  # (P,) bool  finite cell (always True with a box)
    n_real: int = 0
    box_margin: float | None = None

    def padded(self, mask: np.ndarray) -> np.ndarray:
        """A per-point mask extended with False for the mirror images (indexable by ridge_points)."""
        out = np.zeros(len(self.vor.points), dtype=bool)
        out[: len(mask)] = mask
        return out


def mirror_points(points: np.ndarray, margin: float, slab: float | None = None) -> np.ndarray:
    """Reflect the cloud across the six faces of its bounding box expanded by ``margin`` (after
    py_della ``_mirror``). The images close every cell, i.e. the diagram is clipped to the box.
    Only points within ``slab`` of a face are reflected (None: all of them)."""
    lo, hi = points.min(axis=0) - margin, points.max(axis=0) + margin
    out = [points]
    for axis in range(3):
        for plane in (lo[axis], hi[axis]):
            near = np.ones(len(points), dtype=bool) if slab is None else np.abs(points[:, axis] - plane) <= slab
            img = points[near].copy()
            img[:, axis] = 2 * plane - img[:, axis]
            out.append(img)
    return np.vstack(out)


def voronoi_cells(points: np.ndarray, max_radius: float, box_margin: float | None = None) -> VoronoiCells:
    """Voronoi cells of the points with volume, area and the Richards 'interior' test.

    With ``box_margin`` the cloud is mirrored across the faces of its bounding box expanded by that
    margin, so every cell is closed and gets a volume (surface cells are clipped by the box and are
    margin-dependent; ``bounded`` still tells interior cells apart). Without it, open cells get NaN.
    """
    pts = np.asarray(points, dtype=np.float64)
    n = len(pts)
    # with a box: reflect the points near each face first; if a real cell is still open, reflect them all
    slabs = (None,) if box_margin is None else (2 * float(box_margin) + 15.0, None)
    for slab in slabs:
        all_pts = pts if box_margin is None else mirror_points(pts, float(box_margin), slab)
        try:
            vor = Voronoi(all_pts)
        except QhullError as exc:
            raise GeometryError(f"Qhull Voronoi failed: {str(exc).splitlines()[0]}") from exc
        rv = vor.ridge_vertices
        rp = np.asarray(vor.ridge_points, dtype=np.int64)
        lens = np.fromiter(map(len, rv), dtype=np.int64, count=len(rv))
        flat = np.fromiter(ichain.from_iterable(rv), dtype=np.int64, count=int(lens.sum()))
        has_inf = np.zeros(len(rv), dtype=bool)
        np.logical_or.at(has_inf, np.repeat(np.arange(len(rv)), lens), flat < 0)
        if box_margin is None or not (has_inf & (rp < n).any(axis=1)).any():
            break
    N = len(all_pts)
    starts = np.concatenate([[0], np.cumsum(lens)[:-1]])

    closed = np.ones(N, dtype=bool)
    closed[rp[has_inf].ravel()] = False
    # Richards test: a vertex farther than max_radius from its site -> not an interior cell
    V = vor.vertices
    owner_a = np.repeat(rp[:, 0], lens)
    owner_b = np.repeat(rp[:, 1], lens)
    fin = flat >= 0
    d = np.full(len(flat), np.inf)
    d[fin] = np.linalg.norm(V[flat[fin]] - all_pts[owner_a[fin]], axis=1)  # equidistant to both sites
    too_far = d > max_radius
    interior = closed.copy()
    interior[owner_a[too_far]] = False
    interior[owner_b[too_far]] = False

    # fan triangulation of finite ridges (Qhull lists 3-D ridge vertices in polygon order)
    fr = ~has_inf & (lens >= 3)
    m = lens[fr] - 2
    ridge_ids = np.repeat(np.nonzero(fr)[0], m)
    k = np.arange(m.sum()) - np.repeat(np.cumsum(m) - m, m) + 1
    s0 = starts[ridge_ids]
    tri = np.stack([flat[s0], flat[s0 + k], flat[s0 + k + 1]], axis=1)
    area = 0.5 * np.linalg.norm(np.cross(V[tri[:, 1]] - V[tri[:, 0]], V[tri[:, 2]] - V[tri[:, 0]]), axis=1)
    ridge_area = np.bincount(ridge_ids, weights=area, minlength=len(rv))
    h = 0.5 * np.linalg.norm(all_pts[rp[:, 0]] - all_pts[rp[:, 1]], axis=1)
    pyr = ridge_area * h / 3.0
    volume = np.bincount(rp[:, 0], weights=pyr, minlength=N) + np.bincount(rp[:, 1], weights=pyr, minlength=N)
    cell_area = np.bincount(rp[:, 0], weights=ridge_area, minlength=N) + np.bincount(rp[:, 1], weights=ridge_area, minlength=N)
    volume = np.where(closed, volume, np.nan)[:n]
    cell_area = np.where(closed, cell_area, np.nan)[:n]
    n_faces = np.bincount(rp.ravel(), minlength=N)[:n]
    return VoronoiCells(vor, volume, cell_area, interior[:n], n_faces, tri, ridge_ids, rp, ~has_inf,
                        closed=closed[:n], n_real=n, box_margin=box_margin)


def voronoi_cell_mesh(cells: VoronoiCells, sites: np.ndarray, shrink: float = 0.9
                      ) -> tuple[np.ndarray, np.ndarray]:
    """Triangle soup of the selected (bounded) cells, each shrunk toward its site.

    Returns (positions (Q*3,3) float32, site index per vertex (Q*3,) int32).
    Each ridge is emitted once per selected owner, so neighbouring cells stay
    visually separate when ``shrink < 1``.
    """
    pts = cells.vor.points
    V = cells.vor.vertices
    sel = np.zeros(len(pts), dtype=bool)
    sel[sites] = True                               # the caller picks the cells (interior or all closed)
    rp = cells.ridge_points[cells.tri_ridge]
    out_pos, out_site = [], []
    for side in (0, 1):
        owner = rp[:, side]
        m = sel[owner]
        if not m.any():
            continue
        tris = V[cells.tri[m]]                       # (q,3,3)
        c = pts[owner[m]][:, None, :]
        tris = c + shrink * (tris - c)
        out_pos.append(tris.reshape(-1, 3))
        out_site.append(np.repeat(owner[m], 3))
    if not out_pos:
        return np.zeros((0, 3), np.float32), np.zeros(0, np.int32)
    return np.concatenate(out_pos).astype(np.float32), np.concatenate(out_site).astype(np.int32)


def voronoi_ridge_edges(cells: VoronoiCells, max_radius: float) -> np.ndarray:
    """Unique finite Voronoi edges (pairs of Voronoi vertex indices) of bounded ridges."""
    rv = cells.vor.ridge_vertices
    pairs = []
    ok = cells.ridge_finite & cells.padded(cells.bounded)[cells.ridge_points].any(axis=1)
    for r in np.nonzero(ok)[0]:
        poly = rv[r]
        pairs.extend(zip(poly, poly[1:] + poly[:1]))
    if not pairs:
        return np.zeros((0, 2), np.int64)
    e = np.sort(np.asarray(pairs, dtype=np.int64), axis=1)
    return np.unique(e, axis=0)


# ------------------------------------------------------------------ convex hull

@dataclass
class Hull:
    faces: np.ndarray         # (F,3) outward oriented
    vertices: np.ndarray      # hull vertex indices
    volume: float
    area: float
    equations: np.ndarray     # (F,4) outward unit normals + offset

    @property
    def sphericity(self) -> float:
        """psi = pi^(1/3) (6V)^(2/3) / A   (1 for a sphere)."""
        return float(np.pi ** (1 / 3) * (6 * self.volume) ** (2 / 3) / self.area) if self.area else 0.0


def convex_hull(points: np.ndarray) -> Hull:
    pts = np.asarray(points, dtype=np.float64)
    if len(pts) < 4:
        raise GeometryError("need at least 4 points for a 3-D hull")
    try:
        h = ConvexHull(pts)
    except QhullError as exc:
        raise GeometryError(f"Qhull hull failed: {str(exc).splitlines()[0]}") from exc
    faces = h.simplices.copy()
    # scipy does not guarantee winding; orient with the facet normals
    a, b, c = pts[faces[:, 0]], pts[faces[:, 1]], pts[faces[:, 2]]
    flip = np.einsum("ij,ij->i", np.cross(b - a, c - a), h.equations[:, :3]) < 0
    faces[flip] = faces[flip][:, [0, 2, 1]]
    return Hull(faces, h.vertices, float(h.volume), float(h.area), h.equations)


def hull_depth(hull: Hull, points: np.ndarray, chunk: int = 20000) -> np.ndarray:
    """Distance from each (interior) point to the hull boundary = min over facet planes."""
    eq = hull.equations
    out = np.empty(len(points))
    for s in range(0, len(points), chunk):
        p = points[s:s + chunk]
        out[s:s + chunk] = -(p @ eq[:, :3].T + eq[:, 3]).max(axis=1)
    return np.maximum(out, 0.0)


def convex_layers(points: np.ndarray, max_layers: int = 10_000) -> np.ndarray:
    """Onion peeling: layer index of each point (0 = outermost hull)."""
    layer = np.full(len(points), -1, dtype=np.int32)
    remaining = np.arange(len(points))
    k = 0
    while len(remaining) >= 4 and k < max_layers:
        try:
            h = ConvexHull(points[remaining])
        except QhullError:
            break
        layer[remaining[h.vertices]] = k
        remaining = np.delete(remaining, h.vertices)
        k += 1
    layer[remaining] = k
    return layer


# ------------------------------------------------------------------ shape

def gyration_descriptors(points: np.ndarray) -> dict[str, float]:
    """Radius of gyration and shape descriptors from the gyration tensor."""
    c = points - points.mean(axis=0)
    S = c.T @ c / len(points)
    l3, l2, l1 = np.sort(np.linalg.eigvalsh(S))  # l1 >= l2 >= l3
    rg2 = l1 + l2 + l3
    return {
        "rg": float(np.sqrt(rg2)),
        "asphericity": float(l1 - 0.5 * (l2 + l3)),
        "acylindricity": float(l2 - l3),
        "kappa2": float(1.5 * (l1**2 + l2**2 + l3**2) / rg2**2 - 0.5) if rg2 else 0.0,
        "principal_lengths": [float(np.sqrt(x)) for x in (l1, l2, l3)],
    }
