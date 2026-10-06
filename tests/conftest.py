"""Shared fixtures: a synthetic structure built from internal coordinates (exact, offline) and the
cached 4HHB entry for integration tests (skipped when it is neither cached nor downloadable)."""
from __future__ import annotations

import numpy as np
import pytest


# ideal backbone geometry (Engh & Huber)
B_NCA, B_CAC, B_CN, B_CO = 1.458, 1.525, 1.329, 1.231
A_NCAC, A_CACN, A_CNCA, A_CACO = 111.2, 116.2, 121.7, 120.5


def place(a, b, c, bond, angle, torsion):
    """NeRF: the atom d with |cd| = bond, angle bcd and dihedral abcd (degrees)."""
    angle, torsion = np.radians(angle), np.radians(torsion)
    bc = c - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d2 = np.array([-bond * np.cos(angle), bond * np.sin(angle) * np.cos(torsion), bond * np.sin(angle) * np.sin(torsion)])
    return c + d2[0] * bc + d2[1] * m + d2[2] * n


def backbone(n_res, phi, psi, omega=180.0, origin=(0.0, 0.0, 0.0)):
    """[(N, CA, C, O, CB)] for a chain with constant φ/ψ/ω."""
    o = np.asarray(origin, float)
    N = o + np.array([0.0, 0.0, 0.0])
    CA = o + np.array([B_NCA, 0.0, 0.0])
    C = CA + B_CAC * np.array([-np.cos(np.radians(A_NCAC)), np.sin(np.radians(A_NCAC)), 0.0])
    out = []
    for i in range(n_res):
        if i:
            N_new = place(prev[0], prev[1], prev[2], B_CN, A_CACN, psi)
            CA_new = place(prev[1], prev[2], N_new, B_NCA, A_CNCA, omega)
            C_new = place(prev[2], N_new, CA_new, B_CAC, A_NCAC, phi)
            N, CA, C = N_new, CA_new, C_new
        O = place(N, CA, C, B_CO, A_CACO, psi + 180.0)
        CB = place(C, N, CA, 1.53, 110.5, -122.5)
        prev = (N, CA, C)
        out.append((N, CA, C, O, CB))
    return out


def pdb_atom(serial, name, resn, chain, resi, xyz, element, het=False, b=20.0):
    rec = "HETATM" if het else "ATOM  "
    nm = f" {name:<3}" if len(name) < 4 else name
    return (f"{rec}{serial:5d} {nm:<4} {resn:>3} {chain}{resi:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00{b:6.2f}          {element:>2}")


def synthetic_pdb(helix_res=16):
    """Chain A: ideal α-helix (ALA, GLY at 8); chain a: same helix shifted 30 Å (tests case-sensitive ids);
    hetero groups far from the protein: two ZN, an acetate ligand, and four waters
    (ion shell, ligand-bound, bulk, protein-only)."""
    lines, k = [], 1
    for chain, shift in (("A", 0.0), ("a", 30.0)):
        for i, (N, CA, C, O, CB) in enumerate(backbone(helix_res, -57.0, -47.0, origin=(shift, 0, 0))):
            resn = "GLY" if i == 7 else "ALA"
            atoms = [("N", N, "N"), ("CA", CA, "C"), ("C", C, "C"), ("O", O, "O")]
            if resn == "ALA":
                atoms.append(("CB", CB, "C"))
            for name, xyz, el in atoms:
                lines.append(pdb_atom(k, name, resn, chain, i + 1, xyz, el))
                k += 1
        lines.append(f"TER   {k:5d}      {resn} {chain}{helix_res:4d}")
        k += 1
    zn1, zn2 = np.array([0.0, 40.0, 0.0]), np.array([0.0, 60.0, 0.0])
    act = np.array([20.0, 40.0, 0.0])
    het = [("ZN", "ZN", "B", 101, zn1, "ZN"), ("ZN", "ZN", "B", 102, zn2, "ZN"),
           ("C", "ACT", "B", 103, act, "C"), ("O", "ACT", "B", 103, act + [1.25, 0, 0], "O"),
           ("OXT", "ACT", "B", 103, act + [-0.6, 1.1, 0], "O"), ("CH3", "ACT", "B", 103, act + [-0.75, -1.25, 0], "C"),
           ("O", "HOH", "W", 201, zn1 + [2.1, 0, 0], "O"),                 # ion shell
           ("O", "HOH", "W", 202, act + [1.25 + 2.9, 0, 0], "O"),         # H-bond to the acetate O
           ("O", "HOH", "W", 203, np.array([-40.0, -40.0, -40.0]), "O"),  # bulk solvent
           ("O", "HOH", "W", 204, np.array([0.0, 0.0, -2.9]), "O")]       # touches protein only
    for name, resn, chain, resi, xyz, el in het:
        lines.append(pdb_atom(k, name, resn, chain, resi, xyz, el, het=True))
        k += 1
    lines.append("END")
    return "\n".join(lines) + "\n"


@pytest.fixture(scope="session")
def synthetic_path(tmp_path_factory):
    p = tmp_path_factory.mktemp("struct") / "synthetic.pdb"
    p.write_text(synthetic_pdb())
    return p


@pytest.fixture(scope="session")
def synthetic(synthetic_path):
    from biogeometry.structure import load_structure_file
    return load_structure_file(synthetic_path)


@pytest.fixture(scope="session")
def hhb():
    """4HHB from the local cache (downloaded once if the network is available)."""
    from biogeometry.structure import StructureError, fetch_structure
    try:
        return fetch_structure("4HHB")
    except StructureError as exc:
        pytest.skip(f"4HHB not cached and not downloadable: {exc}")


@pytest.fixture(scope="session")
def hhb_job():
    from biogeometry import engine
    from biogeometry.structure import StructureError
    try:
        return engine.open_job(engine.JobSpec(source="rcsb", pdb_id="4HHB"))
    except StructureError as exc:
        pytest.skip(f"4HHB not available: {exc}")


@pytest.fixture(scope="session", autouse=True)
def isolated_cache(tmp_path_factory):
    """Keep layer caches and uploads written by tests out of the project's data/ directory."""
    from biogeometry import engine, structure
    old = engine.CACHE_DIR, structure.DATA_DIR
    engine.CACHE_DIR = tmp_path_factory.mktemp("jobs")
    structure.DATA_DIR = tmp_path_factory.mktemp("data")
    yield
    engine.CACHE_DIR, structure.DATA_DIR = old
