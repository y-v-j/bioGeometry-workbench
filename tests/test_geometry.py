"""Geometry checked against independent oracles (closed forms, scipy, brute force)."""
import numpy as np
import pytest
from scipy.spatial import ConvexHull

from biogeometry import geometry as g

RNG = np.random.default_rng(7)


def fib_sphere(n, r, center=(0, 0, 0)):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)] * r + center


@pytest.fixture(scope="module")
def cloud():
    return RNG.normal(size=(400, 3)) * [10, 6, 4]


def test_tet_volumes_tile_hull(cloud):
    t = g.Tess.build(cloud)
    assert t.volumes.sum() == pytest.approx(ConvexHull(cloud).volume, rel=1e-9)


def test_circumradius_matches_linear_solve(cloud):
    t = g.Tess.build(cloud)
    P = cloud[t.simplices[:200]]
    for tet, r in zip(P, t.circumradii[:200]):
        A = 2 * (tet[1:] - tet[0])
        b = (tet[1:] ** 2).sum(1) - (tet[0] ** 2).sum()
        c = np.linalg.solve(A, b)
        dists = np.linalg.norm(tet - c, axis=1)
        assert np.allclose(dists, r, rtol=1e-7)


def test_unique_edges_brute_force(cloud):
    t = g.Tess.build(cloud)
    brute = {tuple(sorted((s[i], s[j]))) for s in t.simplices for i in range(4) for j in range(i + 1, 4)}
    assert {tuple(e) for e in t.edges} == brute


def test_alpha_infinite_is_hull(cloud):
    t = g.Tess.build(cloud)
    a = g.alpha_shape(t, np.inf)
    h = ConvexHull(cloud)
    assert a.volume == pytest.approx(h.volume, rel=1e-9)
    assert a.area == pytest.approx(h.area, rel=1e-9)
    assert a.n_components == 1 and a.n_cavities == 0
    # outward orientation: signed volume from faces equals +volume (divergence theorem)
    p = cloud[a.faces]
    signed = np.einsum("ij,ij->i", p[:, 0], np.cross(p[:, 1], p[:, 2])).sum() / 6
    assert signed == pytest.approx(h.volume, rel=1e-9)


def test_alpha_tiny_is_empty(cloud):
    a = g.alpha_shape(g.Tess.build(cloud), 1e-6)
    assert a.volume == 0 and len(a.faces) == 0 and len(a.isolated_points) == len(cloud)


def test_alpha_volume_monotone_and_filtration_consistent(cloud):
    t = g.Tess.build(cloud)
    alphas = np.linspace(1, 30, 12)
    f = g.alpha_filtration(t, alphas)
    assert np.all(np.diff(f["volume"]) >= -1e-9)
    for a, v, ar in zip(alphas, f["volume"], f["area"]):
        s = g.alpha_shape(t, a)
        assert v == pytest.approx(s.volume) and ar == pytest.approx(s.area)


def test_hollow_shell_has_one_cavity():
    pts = np.vstack([fib_sphere(700, 10), fib_sphere(900, 11.5)])
    t = g.Tess.build(pts)
    a = g.alpha_shape(t, 2.0)            # shell closes, interior (radius ~10) stays empty
    assert a.n_cavities == 1
    pockets = g.find_pockets(t, 2.0, np.inf)
    big = pockets[0]
    assert big.enclosed
    assert big.volume == pytest.approx(4 / 3 * np.pi * 10 ** 3, rel=0.05)


def test_open_pocket_has_mouth():
    # shell with a polar cap removed -> the interior is a pocket with an opening
    pts = np.vstack([fib_sphere(700, 10), fib_sphere(900, 11.5)])
    pts = pts[pts[:, 2] < 8.5]
    t = g.Tess.build(pts)
    pk = g.find_pockets(t, 2.0, 30.0)[0]
    assert not pk.enclosed and pk.mouth_area > 0
    assert pk.volume > 0.5 * 4 / 3 * np.pi * 10 ** 3


def test_voronoi_volume_matches_convex_hull_of_cell():
    # jittered lattice: interior cells are bounded; ConvexHull of the cell's vertices is an
    # independent volume oracle (and checks that the ridge polygons are traversed in order)
    grid = np.stack(np.meshgrid(*[np.arange(7.0)] * 3, indexing="ij"), -1).reshape(-1, 3)
    pts = grid + RNG.normal(scale=0.08, size=grid.shape)
    cells = g.voronoi_cells(pts, max_radius=2.0)
    vor = cells.vor
    interior = np.nonzero(cells.bounded)[0]
    assert len(interior) >= 125  # the 5x5x5 core at least
    for i in interior:
        verts = vor.vertices[vor.regions[vor.point_region[i]]]
        assert cells.volume[i] == pytest.approx(ConvexHull(verts).volume, rel=1e-8)
    core = np.all((grid >= 1) & (grid <= 5), axis=1)
    assert np.nanmean(cells.volume[core]) == pytest.approx(1.0, rel=0.02)


def test_voronoi_open_cells_on_surface():
    pts = RNG.uniform(-10, 10, size=(500, 3))
    cells = g.voronoi_cells(pts, max_radius=6.0)
    hull_v = ConvexHull(pts).vertices
    assert not cells.bounded[hull_v].any()          # hull points have unbounded cells
    assert np.isnan(cells.volume[hull_v]).all()


def test_voronoi_mesh_shrink():
    grid = np.stack(np.meshgrid(*[np.arange(5.0)] * 3, indexing="ij"), -1).reshape(-1, 3)
    pts = grid + RNG.normal(scale=0.05, size=grid.shape)
    cells = g.voronoi_cells(pts, 2.0)
    site = int(np.nonzero(cells.bounded)[0][0])
    pos, sid = g.voronoi_cell_mesh(cells, np.array([site]), shrink=0.5)
    assert set(sid) == {site}
    # all shrunk vertices within half the max cell reach of the site
    assert np.linalg.norm(pos - pts[site], axis=1).max() <= 0.5 * 2.0 + 1e-6


def test_hull_depth_cube_center():
    cube = np.array([[x, y, z] for x in (0, 1) for y in (0, 1) for z in (0, 1)], float)
    h = g.convex_hull(cube)
    assert h.volume == pytest.approx(1) and h.area == pytest.approx(6)
    d = g.hull_depth(h, np.array([[0.5, 0.5, 0.5], [0.5, 0.5, 0.1], [0, 0, 0]]))
    assert d == pytest.approx([0.5, 0.1, 0.0])


def test_hull_faces_outward():
    pts = RNG.normal(size=(300, 3))
    h = g.convex_hull(pts)
    p = pts[h.faces]
    signed = np.einsum("ij,ij->i", p[:, 0], np.cross(p[:, 1], p[:, 2])).sum() / 6
    assert signed == pytest.approx(h.volume)


def test_sphericity_of_sphere():
    h = g.convex_hull(fib_sphere(4000, 3.0))
    assert h.sphericity == pytest.approx(1.0, abs=0.01)


def test_convex_layers_nested_cubes():
    cube = np.array([[x, y, z] for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)], float)
    pts = np.vstack([cube * 3, cube * 2, cube * 1])
    assert list(g.convex_layers(pts)) == [0] * 8 + [1] * 8 + [2] * 8


def test_gyration_line_and_sphere():
    line = np.c_[np.linspace(-1, 1, 101), np.zeros(101), np.zeros(101)]
    assert g.gyration_descriptors(line)["kappa2"] == pytest.approx(1.0)
    sph = g.gyration_descriptors(fib_sphere(5000, 1.0))
    assert sph["kappa2"] == pytest.approx(0.0, abs=1e-3)
    assert sph["rg"] == pytest.approx(1.0, rel=1e-3)
