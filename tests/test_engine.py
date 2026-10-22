"""Structure handling, backbone, biochemistry, SASA, feature tables, layer cache and HTTP API."""
from __future__ import annotations

import gzip
import json

import numpy as np
import pytest

from biogeometry import backbone, biochem, engine, features
from biogeometry.sasa import atom_sasa, radii_for
from biogeometry.structure import (ELEMENTS, KIND_ION, KIND_LIGAND, KIND_PROTEIN, KIND_WATER, StructureError,
                                   coordinated_water_roles, hetero_groups, load_structure_file,
                                   normalize_pdb_id, normalize_uniprot, representative_points, select_chains)

from conftest import backbone as build_backbone


# --------------------------------------------------------------------------- structure

def test_only_coordinated_waters_are_kept(synthetic, synthetic_path):
    waters = np.flatnonzero(synthetic.res_kind == KIND_WATER)
    assert sorted(synthetic.res_seq[waters]) == [201, 202]          # bulk (203) and protein-only (204) dropped
    roles = coordinated_water_roles(synthetic)
    by_seq = {int(synthetic.res_seq[r]): role for r, role in roles.items()}
    assert by_seq == {201: "ion shell", 202: "ligand"}
    none = load_structure_file(synthetic_path, waters="none")
    assert not (none.res_kind == KIND_WATER).any()


def test_residue_kinds(synthetic):
    kinds = dict(zip(synthetic.res_name, synthetic.res_kind))
    assert kinds["ALA"] == KIND_PROTEIN and kinds["GLY"] == KIND_PROTEIN
    assert kinds["ZN"] == KIND_ION and kinds["ACT"] == KIND_LIGAND


def test_hetero_numbering_and_water_groups(synthetic):
    groups = hetero_groups(synthetic)
    labels = [g["label"] for g in groups]
    assert labels[:3] == ["ZN (i)", "ZN (ii)", "ACT (i)"]
    water = {g["label"]: g for g in groups if g["kind"] == "water"}
    assert set(water) == {"H₂O ion shell (1)", "H₂O ligand (1)"}
    assert all(len(g["members"]) == 1 for g in water.values())


def test_roman():
    assert [biochem.roman(n) for n in (1, 4, 9, 14, 40, 90, 400, 1994)] == \
        ["i", "iv", "ix", "xiv", "xl", "xc", "cd", "mcmxciv"]


def test_chain_selection_is_case_sensitive(synthetic):
    assert select_chains(synthetic, ["a"]).chains == ["a"]
    assert select_chains(synthetic, ["A"]).chains == ["A"]
    assert select_chains(synthetic, ["w"]).chains == ["W"]     # wrong case accepted when unambiguous
    with pytest.raises(StructureError):
        select_chains(synthetic, ["Z"])


def test_chain_info_handles_hetero_only_chains(synthetic, hhb):
    # PDB-format file whose ions/waters sit in their own chains (regression: get_entity_of(Chain) crashed)
    assert set(synthetic.chain_info) == {"A", "a", "B", "W"}
    assert hhb.chain_info["A"]["description"] == "Hemoglobin subunit alpha"
    assert hhb.chain_info["A"]["uniprot"] == "P69905"


@pytest.mark.parametrize("level", ["heavy", "ca"])
def test_hetero_atoms_join_both_analysis_modes(synthetic, level):
    """All-atom and Cα modes both add every ligand, ion and coordinated-water atom; polymer points come first."""
    ps = representative_points(synthetic, level, include_hetero=True)
    kinds = synthetic.res_kind[ps.residue]
    het_atoms = int((synthetic.res_kind[synthetic.atom_res] >= KIND_LIGAND).sum())
    assert (kinds >= KIND_LIGAND).sum() == het_atoms == 2 + 4 + 2        # 2 Zn, acetate, 2 coordinated waters
    n_poly = int((kinds < KIND_LIGAND).sum())
    assert (kinds[:n_poly] < KIND_LIGAND).all()
    assert n_poly == (32 if level == "ca" else int((synthetic.res_kind[synthetic.atom_res] < KIND_LIGAND).sum()))
    assert np.allclose(ps.coords[n_poly:], synthetic.coords[ps.atom[n_poly:]])
    plain = representative_points(synthetic, level)
    assert set(synthetic.res_kind[plain.residue]) == {KIND_PROTEIN}


def test_id_normalisation():
    assert normalize_pdb_id(" 4hhb ") == "4HHB"
    assert normalize_uniprot("AF-P69905-F1-model_v6") == "P69905"
    with pytest.raises(StructureError):
        normalize_pdb_id("hemoglobin")
    with pytest.raises(StructureError):
        normalize_uniprot("not-an-accession")


def test_4hhb_coordinated_waters(hhb):
    kinds = np.bincount(hhb.res_kind, minlength=5)
    assert kinds[KIND_PROTEIN] == 574
    roles = coordinated_water_roles(hhb)
    assert kinds[KIND_WATER] == len(roles) > 0           # every kept water has a role
    het = [g["label"] for g in hetero_groups(hhb) if g["kind"] != "water"]
    assert het.count("HEM (i)") == 1 and "HEM (iv)" in het


# --------------------------------------------------------------------------- backbone

def test_dihedrals_of_ideal_helix(synthetic):
    bb = backbone.backbone_atoms(synthetic)
    phi, psi, omega = backbone.dihedrals(synthetic, bb)
    a = slice(1, 15)                                      # interior residues of chain A
    assert np.allclose(phi[a], -57.0, atol=0.3)
    assert np.allclose(psi[a], -47.0, atol=0.3)
    assert np.allclose(np.abs(omega[a]), 180.0, atol=0.3)
    assert np.isnan(phi[0]) and np.isnan(psi[15])        # chain termini


def test_dihedrals_of_strand(tmp_path):
    from conftest import pdb_atom
    lines, k = [], 1
    for i, (N, CA, C, O, CB) in enumerate(build_backbone(8, -120.0, 130.0)):
        for name, xyz, el in (("N", N, "N"), ("CA", CA, "C"), ("C", C, "C"), ("O", O, "O"), ("CB", CB, "C")):
            lines.append(pdb_atom(k, name, "VAL", "A", i + 1, xyz, el))
            k += 1
    p = tmp_path / "strand.pdb"
    p.write_text("\n".join(lines) + "\nEND\n")
    s = load_structure_file(p)
    phi, psi, _ = backbone.dihedrals(s, backbone.backbone_atoms(s))
    assert np.allclose(phi[1:], -120.0, atol=0.3) and np.allclose(psi[:-1], 130.0, atol=0.3)
    assert set(biochem.rama_region(phi[1:-1], psi[1:-1])) == {"β"}


def test_secondary_structure_of_ideal_helix(synthetic):
    bb = backbone.backbone_atoms(synthetic)
    ss, hb = backbone.secondary_structure(synthetic, bb, backbone.protein_segments(synthetic, bb))
    assert (ss[1:15] == 1).all() and ss[0] == 0 and ss[15] == 0
    assert (ss[synthetic.res_kind != KIND_PROTEIN] == -1).all()
    assert len(hb) >= 2 * 11                              # i -> i+4 bonds in both helices


def test_4hhb_secondary_structure(hhb_job):
    ss = hhb_job.ss[hhb_job.structure.res_kind == KIND_PROTEIN]
    assert 0.70 < (ss == 1).mean() < 0.90                 # haemoglobin is ~75-80 % helix by DSSP
    assert (ss == 2).mean() < 0.02                        # and has no β-sheet


def test_cartoon_mesh_is_well_formed(hhb_job):
    s, bb = hhb_job.structure, hhb_job.bb
    mesh = backbone.cartoon_mesh(s, bb, hhb_job.ss, hhb_job.segments)
    pos = mesh["pos"].reshape(-1, 3)
    assert np.isfinite(pos).all() and len(pos) == len(mesh["res"])
    assert mesh["res"].min() >= 0 and mesh["res"].max() < s.n_residues
    lo, hi = s.coords.min(0) - 5, s.coords.max(0) + 5
    assert ((pos >= lo) & (pos <= hi)).all()


def test_mesh_budget_shrinks_with_size():
    budgets = [backbone.mesh_budget(n) for n in (500, 5000, 50000)]
    assert budgets[0] > budgets[1] > budgets[2]


# --------------------------------------------------------------------------- biochemistry

def test_rama_regions():
    phi = np.array([-63, -120, -75, 60, 60, np.nan])
    psi = np.array([-43, 130, 145, 45, -150, 0])
    assert list(biochem.rama_region(phi, psi)) == ["αR", "β", "PPII", "αL", "other", ""]
    out = biochem.rama_outlier(np.array(["other", "other", "αR"], dtype=object), ["GLY", "ALA", "ALA"])
    assert list(out) == [False, True, False]


def test_bands():
    assert list(biochem.exposure_class(np.array([0.01, 0.1, 0.3, 0.8, np.nan]))) == \
        ["core", "buried", "intermediate", "exposed", ""]
    assert list(biochem.plddt_band(np.array([95, 80, 60, 30]))) == ["very high", "confident", "low", "very low"]
    assert list(biochem.residue_class(["LYS", "ASP", "PHE", "GLY", "HEM"], [0, 0, 0, 0, KIND_LIGAND])) == \
        ["basic", "acidic", "aromatic", "special", "ligand"]


# --------------------------------------------------------------------------- SASA

def test_sasa_isolated_sphere():
    a = atom_sasa(np.zeros((1, 3)), np.array([1.7]), n_points=200)
    assert a[0] == pytest.approx(4 * np.pi * (1.7 + 1.4) ** 2, rel=1e-9)


def test_sasa_two_overlapping_spheres():
    r, d = 1.7, 3.0
    R = r + 1.4
    exact = 4 * np.pi * R ** 2 - 2 * np.pi * R * (R - d / 2)      # sphere minus the buried cap
    a = atom_sasa(np.array([[0, 0, 0], [d, 0, 0]], float), np.array([r, r]), n_points=4000)
    assert a == pytest.approx([exact, exact], rel=0.01)


def test_sasa_neighbour_fallback_is_exact():
    rng = np.random.default_rng(1)
    xyz = rng.uniform(0, 12, (300, 3))
    rad = rng.uniform(1.5, 2.0, 300)
    assert np.allclose(atom_sasa(xyz, rad, k=2), atom_sasa(xyz, rad, k=300))


def test_sasa_matches_biopython(hhb):
    pytest.importorskip("Bio")
    from Bio.PDB import MMCIFParser, ShrakeRupley
    from biogeometry.structure import STRUCT_DIR
    import io
    text = gzip.decompress((STRUCT_DIR / "4HHB.cif.gz").read_bytes()).decode()
    bio = MMCIFParser(QUIET=True).get_structure("4HHB", io.StringIO(text))[0]
    for chain in bio:
        for res in [r for r in chain if r.id[0] != " "]:
            chain.detach_child(res.id)                   # polymer only, like the comparison below
        for res in chain:
            for atom in [a for a in res if a.element == "H"]:
                res.detach_child(atom.id)
    ShrakeRupley().compute(bio, level="A")
    ref = sum(a.sasa for a in bio.get_atoms())
    poly = hhb.polymer_mask()[hhb.atom_res]
    elements = [ELEMENTS[e] if e < len(ELEMENTS) else "X" for e in hhb.element[poly]]
    ours = atom_sasa(hhb.coords[poly], radii_for(elements)).sum()
    assert ours == pytest.approx(ref, rel=0.03)


# --------------------------------------------------------------------------- tetrahedra

def test_tetra_metrics_regular_tetrahedron():
    a = 2.0
    pts = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float) * a / (2 * np.sqrt(2))
    m = features.tetra_metrics(pts, np.array([[0, 1, 2, 3]]))
    assert all(m[f"L{i}"][0] == pytest.approx(a) for i in range(1, 7))
    assert m["Tvol"][0] == pytest.approx(a ** 3 / (6 * np.sqrt(2)))
    assert m["Tsurface_area"][0] == pytest.approx(np.sqrt(3) * a ** 2)
    assert all(m[f"O{i}"][0] == pytest.approx(np.arccos(23 / 27)) for i in range(1, 5))
    assert m["SumL"][0] == pytest.approx(6 * a) and m["AvgL"][0] == pytest.approx(a)


def test_tetra_metrics_corner_tetrahedron():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    m = features.tetra_metrics(pts, np.array([[0, 1, 2, 3]]))
    assert m["Tvol"][0] == pytest.approx(1 / 6)
    assert m["O1"][0] == pytest.approx(np.pi / 2)          # an octant
    assert m["T1_area"][0] == pytest.approx(0.5) and m["T4_area"][0] == pytest.approx(np.sqrt(3) / 2)
    assert m["L1"][0] == pytest.approx(1) and m["L6"][0] == pytest.approx(np.sqrt(2))


def test_simplex_type():
    chain = np.array([0] * 10 + [1] * 10)
    quads = np.array([[3, 4, 5, 6], [6, 3, 5, 4], [3, 4, 5, 9], [3, 4, 7, 15], [1, 3, 5, 7], [8, 9, 10, 11]])
    assert list(features.simplex_type(quads, chain)) == [4, 4, 3, 2, 0, 2]


# --------------------------------------------------------------------------- feature tables

def test_residue_table_columns_and_values(hhb_job):
    t = features.residue_table(hhb_job)
    s = hhb_job.structure
    keys = {c["key"] for c in t.columns}
    for k in ("chain", "resi", "resn", "class", "hydropathy", "ss", "phi", "psi", "rama", "sasa", "rsa",
              "exposure", "cn10", "vor_volume", "hull_depth", "convex_layer", "alpha_surface", "del_degree",
              "hbonds", "salt_bridges", "interface", "site", "pocket", "domain"):
        assert k in keys, k
    assert all(len(v) == t.n == s.n_residues for v in t.data.values())
    col = {c["key"]: c for c in t.columns}

    def cat(key, r):
        return col[key]["categories"][t.data[key][r]]
    prot = np.flatnonzero(s.res_kind == KIND_PROTEIN)
    rsa = t.data["rsa"][prot]
    assert np.nanmin(rsa) >= 0 and np.nanmedian(rsa) < 0.5
    # proximal histidine HIS F8 (α87) coordinates the haem iron
    his87 = next(r for r in prot if s.res_chain[r] == "A" and s.res_seq[r] == 87)
    assert "HEM" in cat("site", his87)
    helix = prot[hhb_job.ss[prot] == 1]
    assert np.mean(t.data["hbonds"][helix] > 0) > 0.9
    assert (t.data["salt_bridges"][prot] > 0).sum() > 10
    assert sum(cat("interface", r) != "" for r in prot) > 50  # α1β1 / α1β2 interfaces


REFERENCE_COLUMNS = (
    "pdbid-chainid model ec-number ec-category atom-set simplex-indices residues residue-class residue-type "
    "participating-chains het-id edge-len-1 edge-len-2 edge-len-3 edge-len-4 edge-len-5 edge-len-6 sum-edge-len avg-edge-len "
    "triangle-1-area triangle-2-area triangle-3-area triangle-4-area surface-area-tetrahedron tetrahedron-volume "
    "omega-1 omega-2 omega-3 omega-4 voronoi-volume-1 voronoi-volume-2 voronoi-volume-3 voronoi-volume-4 "
    "voronoi-area-1 voronoi-area-2 voronoi-area-3 voronoi-area-4 sasa-1 sasa-2 sasa-3 sasa-4 residue-1-depth "
    "residue-2-depth residue-3-depth residue-4-depth depth circumradius circumsphere-volume tetrahedrality tropsha-class simplex-type "
    "alpha-value alpha-shape alpha-shape-volume alpha-shape-surface-area hull-volume hull-surface-area "
    "bfactor-residue-1 bfactor-residue-2 bfactor-residue-3 bfactor-residue-4 sum-bfactor avg-bfactor interaction-type").split()


def test_tetra_table_has_the_reference_columns(hhb_job):
    t, T = features.tetra_table(hhb_job)
    keys = [c["key"] for c in t.columns]
    assert keys == [k for k in REFERENCE_COLUMNS if k not in ("simplex-indices", "residues")]   # built from 'tets'
    assert t.n == T["n_compact"] == len(T["tets"]) > 1000 and not T["sampled"]
    col = {c["key"]: c for c in t.columns}
    assert col["atom-set"]["categories"] == ["all-atom"] and col["pdbid-chainid"]["categories"] == ["4HHB-ABCD"]
    assert set(col["tropsha-class"]["categories"]) == {"{1,1,1,1}", "{2,1,1}", "{2,2}", "{3,1}", "{4}"}
    P = hhb_job.points.coords
    edges = [np.linalg.norm(P[T["tets"][:, i]] - P[T["tets"][:, j]], axis=1) for i in range(4) for j in range(i + 1, 4)]
    assert np.max(edges) <= T["max_edge"] + 1e-4
    assert np.isfinite(t.data["sasa-1"]).all() and np.isfinite(t.data["voronoi-volume-1"]).mean() > 0.3


def test_tetra_csv_header_is_the_reference(hhb_job):
    name, text = engine.export_csv(hhb_job, "tetrahedra")
    assert name == "bioGeometry_4HHB-ABCD_all-atom_delaunay-tetrahedra.csv"
    rows = text.splitlines()
    assert rows[0].split(",") == REFERENCE_COLUMNS
    first = rows[1].split(",")
    assert first[0] == "4HHB-ABCD" and first[1] == "1" and first[4] == "all-atom"
    assert len(first[5].split("-")) == 4 and len(first[6].split("-")) == 4    # simplex indices, residues


def test_tetrahedrality():
    a = 2.0
    reg = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float) * a / (2 * np.sqrt(2))
    flat = np.array([[0, 0, 0], [4, 0, 0], [0, 1, 0], [0, 0, 0.5]], float)
    L = []
    for P in (reg, flat):
        m = features.tetra_metrics(P, np.array([[0, 1, 2, 3]]))
        L.append([m[f"L{i + 1}"][0] for i in range(6)])
    t = features.tetrahedrality(np.array(L))
    assert t[0] == pytest.approx(0, abs=1e-12) and t[1] > 0.1


def test_tropsha_class_and_simplex_type_follow_the_sequence(synthetic_path):
    job = engine.open_job(engine.JobSpec(source="file", path=str(synthetic_path), level="ca", ligands=True))
    t, T = features.tetra_table(job, 12.0)
    col = {c["key"]: c for c in t.columns}
    s, chain = job.structure, job.structure.chain_index()
    poly = np.isin(s.res_kind, (0, 1))
    for i in range(t.n):
        r = T["res4"][i]
        runs, cur = [], 1
        for k in range(3):                               # independent re-derivation of the contiguous runs
            a, b = r[k], r[k + 1]
            if a == b or (b - a == 1 and chain[a] == chain[b] and poly[a] and poly[b]):
                cur += 1
            else:
                runs.append(cur); cur = 1
        runs.append(cur)
        label = "{" + ",".join(map(str, sorted(runs, reverse=True))) + "}"
        assert col["tropsha-class"]["categories"][t.data["tropsha-class"][i]] == label
        assert t.data["simplex-type"][i] == (max(runs) if max(runs) >= 2 else 0)


def test_zinc_site_tetrahedra():
    from biogeometry.structure import StructureError
    try:
        job = engine.open_job(engine.JobSpec(source="rcsb", pdb_id="1CA2", level="heavy", ligands=True))
    except StructureError as exc:
        pytest.skip(str(exc))
    t, T = features.tetra_table(job)
    col = {c["key"]: c for c in t.columns}
    cat = lambda k, i: col[k]["categories"][t.data[k][i]]  # noqa: E731
    assert col["ec-number"]["categories"] == ["4.2.1.1"] and col["ec-category"]["categories"] == ["Lyases"]
    zn = [i for i in range(t.n) if "ZN" in cat("residue-type", i).split("-")]
    assert zn and all("ion" in cat("residue-class", i) and cat("interaction-type", i).startswith("protein-ion") for i in zn)
    assert any("HOH" in cat("residue-type", i) and cat("interaction-type", i) == "protein-ion-water" for i in zn)
    # het-id: each non-polymer name once, in vertex order; blank for polymer-only tetrahedra
    for i in range(t.n):
        pairs = zip(cat("residue-type", i).split("-"), cat("residue-class", i).split("-"))
        assert cat("het-id", i) == "-".join(dict.fromkeys(n for n, c in pairs if c in ("lig", "ion", "hoh")))
    assert {cat("het-id", i) for i in zn} >= {"ZN", "ZN-HOH"}
    r = np.minimum(t.data["circumradius"].astype(float), 1e4)
    assert np.allclose(t.data["circumsphere-volume"], 4 / 3 * np.pi * r ** 3, rtol=1e-5)


def test_vertex_classes_use_the_reference_name_sets():
    names = ["ALA", "MSE", "DA", "U", "HOH", "WAT", "ZN", "HEM", "SEP", "PSU"]
    kinds = [0, 0, 1, 1, 4, 4, 3, 2, 0, 1]
    o2 = [False] * 9 + [True]
    got = [biochem.VERTEX_CLASSES[c] for c in biochem.vertex_class(names, kinds, o2)]
    assert got == ["prt", "prt", "dna", "rna", "hoh", "hoh", "ion", "lig", "prt", "rna"]


def test_residue_depth_classes(hhb):
    sasa = atom_sasa(hhb.coords, radii_for([ELEMENTS[e] if e < len(ELEMENTS) else "C" for e in hhb.element]))
    depth, cls = features.residue_depth(hhb, sasa)
    prot = hhb.res_kind == KIND_PROTEIN
    frac = np.bincount(cls[prot], minlength=3) / prot.sum()
    assert 0.3 < frac[0] < 0.6 and frac[1] > 0.15 and frac[2] > 0.1   # surface, just under, buried all occur
    assert (depth[cls == 2] > features.UNDER_SURFACE).all() and depth[cls == 2].mean() > depth[cls == 1].mean()


# --------------------------------------------------------------------------- engine

def test_job_keys():
    k = engine.JobSpec(source="rcsb", pdb_id="4hhb").key()
    assert k == engine.JobSpec(source="rcsb", pdb_id="4HHB").key()
    assert k != engine.JobSpec(source="rcsb", pdb_id="4HHB", chains="A").key()
    assert k != engine.JobSpec(source="rcsb", pdb_id="4HHB", level="ca").key()
    assert engine.JobSpec(source="rcsb", pdb_id="4HHB", chains="A,B").key() == \
        engine.JobSpec(source="rcsb", pdb_id="4HHB", chains="B A").key()


def test_layers_are_cached_by_params(synthetic_path):
    job = engine.open_job(engine.JobSpec(source="file", path=str(synthetic_path)))
    a = engine.build_layer(job, "alpha", {"alpha": 4.0})
    b = engine.build_layer(job, "alpha", {"alpha": 4.0})
    c = engine.build_layer(job, "alpha", {"alpha": 6.0})
    assert a["key"] == b["key"] != c["key"]
    assert engine.read_meta(job.id, a["key"]) == a
    with pytest.raises(engine.EngineError):
        engine.build_layer(job, "nonsense")


def test_every_layer_kind_builds(synthetic_path):
    job = engine.open_job(engine.JobSpec(source="file", path=str(synthetic_path)))
    for kind in engine.LAYER_KINDS:
        params = {"residues": list(range(10))} if kind == "selection" else None
        meta = engine.build_layer(job, kind, params)
        for name, b in meta["buffers"].items():
            path = engine.buffer_path(job.id, meta["key"], name)
            assert path is not None and path.stat().st_size == np.prod(b["shape"], dtype=int) * np.dtype(b["dtype"]).itemsize
    pae = engine.read_meta(job.id, engine.build_layer(job, "pae")["key"])
    assert pae["stats"]["available"] is False             # not an AlphaFold model


def test_structure_layer_lists_hetero_and_sites(hhb_job):
    st = engine.build_layer(hhb_job, "structure")["stats"]
    hem = [h for h in st["hetero"] if h["name"] == "HEM"]
    assert len(hem) == 4 and all(len(h["site"]) > 10 for h in hem)
    rows = {r["id"]: r for r in st["chain_rows"]}
    assert rows["A"]["description"] == "Hemoglobin subunit alpha" and rows["A"]["n_res"] == 141


def test_csv_exports(hhb_job):
    for kind in engine.CSV_KINDS:
        name, text = engine.export_csv(hhb_job, kind)
        assert name.startswith("bioGeometry_4HHB-ABCD_") and name.endswith(f"_{engine.CSV_SUFFIX[kind]}.csv")
        lines = text.splitlines()
        assert len(lines) > 1, kind
    header = engine.export_csv(hhb_job, "residues")[1].splitlines()[0].split(",")
    assert header[:4] == ["residue", "chain", "resi", "resn"]


def test_path_components_are_validated():
    assert engine.read_meta("..", "x") is None
    assert engine.buffer_path("abc", "..", "x") is None
    assert engine.buffer_path("abc", "key", "") is None


# --------------------------------------------------------------------------- HTTP API

@pytest.fixture(scope="module")
def client():
    from starlette.testclient import TestClient
    from biogeometry.api import api
    return TestClient(api)


def test_api_roundtrip(client, synthetic_path):
    assert client.get("/bg/health").json()["ok"] is True
    r = client.post("/bg/jobs", json={"source": "file", "path": str(synthetic_path)})
    assert r.status_code == 200, r.text
    body = r.json()
    job, s = body["summary"]["job"], body["structure"]
    assert body["summary"]["n_water"] == 2 and body["summary"]["n_ion"] == 2
    meta = client.post(f"/bg/jobs/{job}/layers", json={"kind": "features"}).json()
    assert meta["kind"] == "features"
    name, info = next(iter(meta["buffers"].items()))
    buf = client.get(f"/bg/jobs/{job}/layers/{meta['key']}/buf/{name}")
    assert len(buf.content) == np.prod(info["shape"], dtype=int) * np.dtype(info["dtype"]).itemsize
    assert client.get(f"/bg/jobs/{job}/layers/{s['key']}/meta").json()["kind"] == "structure"
    csv = client.get(f"/bg/jobs/{job}/csv/residues")
    assert csv.status_code == 200 and "attachment" in csv.headers["content-disposition"]


def test_api_upload(client, synthetic_path):
    r = client.post("/bg/upload", files={"file": ("mini.pdb", synthetic_path.read_bytes())}, data={"chains": "A"})
    assert r.status_code == 200, r.text
    assert r.json()["summary"]["chains"] == ["A"]


def test_api_errors(client):
    assert client.post("/bg/jobs", json={"source": "rcsb", "pdb_id": "nope!"}).status_code == 400
    assert client.post("/bg/jobs/doesnotexist/layers", json={"kind": "alpha"}).status_code == 404
    assert client.get("/bg/jobs/x/layers/y/meta").status_code == 404
    assert client.get("/bg/jobs/x/layers/y/buf/z").status_code == 404
    bad = client.post("/bg/upload", files={"file": ("notes.txt", b"hello")})
    assert bad.status_code == 400 and "expected one of" in bad.json()["error"]
    assert set(client.get("/bg/csv-kinds").json()) == set(engine.CSV_KINDS)
    assert json.loads(client.get("/bg/health").text)["app"] == "bioGeometry"


# --------------------------------------------------------------------------- memory budget

def test_layer_refused_when_memory_is_short(synthetic_path, monkeypatch):
    from biogeometry import resources
    job = engine.open_job(engine.JobSpec(source="file", path=str(synthetic_path), level="heavy"))
    cached = engine.build_layer(job, "delaunay")
    monkeypatch.setattr(resources, "available_mb", lambda: 100.0)      # less than the reserve
    assert engine.build_layer(job, "delaunay") == cached              # cache hits need no memory
    with pytest.raises(resources.MemoryBudgetError, match="coarser point level"):
        engine.build_layer(job, "voronoi", {"max_radius": 7.5})


def test_parse_refused_when_memory_is_short(synthetic_path, monkeypatch):
    from biogeometry import resources
    monkeypatch.setattr(resources, "available_mb", lambda: 10.0)
    with pytest.raises(resources.MemoryBudgetError, match="loading synthetic"):
        load_structure_file(synthetic_path)


def test_large_computations_run_one_at_a_time(monkeypatch):
    import threading
    import time
    from biogeometry import resources
    monkeypatch.setattr(resources, "available_mb", lambda: 1e6)
    active, peak = [0], [0]

    def work():
        with resources.budget(resources.HEAVY_MB + 1, "test"):
            active[0] += 1
            peak[0] = max(peak[0], active[0])
            time.sleep(0.05)
            active[0] -= 1
    threads = [threading.Thread(target=work) for _ in range(4)]
    [t.start() for t in threads]
    [t.join() for t in threads]
    assert peak[0] == 1


def test_estimates_scale_with_size():
    from biogeometry import resources
    small = resources.layer_estimate_mb("voronoi", 10_000, 10_000, False)
    big = resources.layer_estimate_mb("voronoi", 400_000, 400_000, False)
    assert big == pytest.approx(40 * small) and 3000 < big < 5000      # ~4 GB at the point limit
    assert resources.layer_estimate_mb("features", 8_000, 59_000, True) < 300   # 1AON at Cα level


def test_api_reports_memory_budget_errors(client, synthetic_path, monkeypatch):
    from biogeometry import resources
    job = client.post("/bg/jobs", json={"source": "file", "path": str(synthetic_path)}).json()["summary"]["job"]
    monkeypatch.setattr(resources, "available_mb", lambda: 100.0)
    r = client.post(f"/bg/jobs/{job}/layers", json={"kind": "voronoi", "params": {"max_radius": 6.5}})
    assert r.status_code == 400 and "memory" in r.json()["error"]
    assert client.get("/bg/health").json()["mem_available_mb"] is not None


def test_cgroup_limit_caps_available_memory(tmp_path, monkeypatch):
    from pathlib import Path
    from biogeometry import resources
    rel = next((ln[3:].strip() for ln in Path("/proc/self/cgroup").read_text().splitlines() if ln.startswith("0::")),
               None) if Path("/proc/self/cgroup").exists() else None
    if rel is None:
        pytest.skip("no cgroup v2 on this system")
    own = tmp_path / rel.lstrip("/")
    own.mkdir(parents=True)
    (own / "memory.max").write_text("max\n")
    (own / "memory.current").write_text(str(50 * 2 ** 20))
    (own.parent / "memory.max").write_text(str(1024 * 2 ** 20))          # a 1 GB limit one level up
    (own.parent / "memory.current").write_text(str(400 * 2 ** 20))
    monkeypatch.setattr(resources, "CGROUP_ROOT", tmp_path)
    monkeypatch.setattr(resources, "_meminfo_mb", lambda: 10_000.0)
    assert resources.available_mb() == pytest.approx(624)                # 1024 − 400, not the 10 GB free
    with pytest.raises(resources.MemoryBudgetError):
        resources.require(200, "test")                                   # 200 + 768 reserve > 624


# --------------------------------------------------------------------------- selection geometry

def test_selection_of_domain_and_haem(hhb_job):
    job = hhb_job
    dom = next(d for d in job.domains if d["chain"] == "A" and d["source"] == "Pfam")
    hem = next(h for h in job.hetero if h["label"] == "HEM (i)")
    res = list(range(dom["r0"], dom["r1"] + 1)) + [hem["res"]]
    meta = engine.build_layer(job, "selection", {"residues": res})
    st = meta["stats"]
    labels = [g["label"] for g in st["groups"]]
    assert labels == ["Pfam Globin · A", "HEM (i)"]
    domain, haem = st["groups"]
    assert domain["kind"] == "domain" and haem["kind"] == "ligand" and haem["n_atoms"] == 43
    assert 0.6 < haem["buried_in_selection"] < 0.95                     # the haem sits deep in its pocket
    pair = st["pairs"][0]
    assert pair["min_dist"] < 2.4 and "HIS87 NE2" in pair["min_pair"][0]  # proximal His - Fe bond
    assert "A:HIS58" in pair["labels_a"] and "A:HIS87" in pair["labels_a"]
    assert 400 < pair["interface_area"] < 900 and pair["n_polar"] > 0
    S = st["summary"]
    assert S["hull_volume"] > S["alpha_volume"] > 0 and S["n_groups"] == 2
    assert set(meta["buffers"]) >= {"pts", "pt_group", "pt_res", "hull_faces", "alpha_faces", "contacts"}


def test_selection_ion_coordination(synthetic_path, tmp_path):
    """A Zn with four N/O ligands at 2.1 Å is reported as metal coordination with those atoms."""
    from conftest import pdb_atom
    lines, k = [], 1
    zn = np.array([0.0, 0.0, 0.0])
    ligands = [("NE2", "HIS", 1, [2.1, 0, 0], "N"), ("ND1", "HIS", 2, [-2.1, 0, 0], "N"),
               ("SG", "CYS", 3, [0, 2.3, 0], "S"), ("OD1", "ASP", 4, [0, -2.0, 0.5], "O")]
    for name, resn, resi, xyz, el in ligands:
        base = zn + np.array(xyz)
        for an, off, e in (("N", [0, 0, 3.5], "N"), ("CA", [0, 0, 2.5], "C"), ("C", [0.8, 0, 3.0], "C"), (name, [0, 0, 0], el)):
            lines.append(pdb_atom(k, an, resn, "A", resi, base + np.array(off) * np.sign(xyz[0] + xyz[1] + 0.1), e)); k += 1
    lines.append(pdb_atom(k, "ZN", "ZN", "A", 101, zn, "ZN", het=True))
    p = tmp_path / "zinc.pdb"
    p.write_text("\n".join(lines) + "\nEND\n")
    job = engine.open_job(engine.JobSpec(source="file", path=str(p), level="heavy"))
    res = list(range(job.structure.n_residues))
    st = engine.build_layer(job, "selection", {"residues": res})["stats"]
    pair = next(q for q in st["pairs"] if st["groups"][q["b"]]["kind"] == "ion")
    assert len(pair["metal"]) == 4 and max(m["dist"] for m in pair["metal"]) < 2.4
    assert {m["atom"].split()[-1] for m in pair["metal"]} == {"NE2", "ND1", "SG", "OD1"}


def test_selection_requires_residues(synthetic_path):
    job = engine.open_job(engine.JobSpec(source="file", path=str(synthetic_path)))
    with pytest.raises(engine.EngineError):
        engine.build_layer(job, "selection", {"residues": []})
    with pytest.raises(engine.EngineError):
        engine.build_layer(job, "selection", {"residues": [10 ** 6]})


def test_delaunay_layer_ships_compact_tetrahedra(hhb_job):
    meta = engine.build_layer(hhb_job, "delaunay")
    path = engine.buffer_path(hhb_job.id, meta["key"], "tets")
    tets = np.frombuffer(path.read_bytes(), dtype="<u4").reshape(-1, 4)
    assert len(tets) == meta["stats"]["n_compact_tets"] > 1000 and not meta["stats"]["tets_sampled"]
    P = hhb_job.points.coords
    edges = [np.linalg.norm(P[tets[:, i]] - P[tets[:, j]], axis=1) for i in range(4) for j in range(i + 1, 4)]
    assert np.max(edges) <= meta["params"]["max_edge"] + 1e-4


# --------------------------------------------------------------------------- structure export

@pytest.mark.parametrize("fmt", ["cif", "pdb"])
def test_structure_export_roundtrip(hhb_job, tmp_path, fmt):
    import gemmi
    name, text = engine.export_structure(hhb_job, fmt)
    assert name == f"bioGeometry_4HHB-ABCD_model.{fmt}"
    p = tmp_path / name
    p.write_text(text)
    back = gemmi.read_structure(str(p))
    s = hhb_job.structure
    assert back[0].count_atom_sites() == s.n_atoms
    assert [c.name for c in back[0]] == s.chains
    hem = [r for c in back[0] for r in c if r.name == "HEM"]
    assert len(hem) == 4 and all(r.het_flag == "H" for r in hem)
    a = back[0]["A"][0][0]
    assert a.pos.dist(gemmi.Position(*s.coords[0])) < 1e-3


def test_pdb_export_shortens_long_chain_ids(tmp_path):
    from biogeometry.structure import load_structure_file, structure_text
    import gemmi
    st = gemmi.Structure()
    model = gemmi.Model(1)
    for name in ("AA", "BB"):
        ch = gemmi.Chain(name)
        for i in range(6):
            res = gemmi.Residue(); res.name = "ALA"; res.seqid = gemmi.SeqId(i + 1, " "); res.entity_type = gemmi.EntityType.Polymer
            at = gemmi.Atom(); at.name = "CA"; at.element = gemmi.Element("C")
            at.pos = gemmi.Position(3.8 * i, 10.0 * (name == "BB"), 0.5 * (i % 2)); res.add_atom(at)
            ch.add_residue(res)
        model.add_chain(ch)
    st.add_model(model)
    p = tmp_path / "long.cif"
    st.setup_entities()
    st.make_mmcif_document().write_file(str(p))
    s = load_structure_file(p)
    assert s.chains == ["AA", "BB"]
    pdb = structure_text(s, "pdb")
    chains = {ln[21] for ln in pdb.splitlines() if ln.startswith("ATOM")}
    assert len(chains) == 2 and all(len(c) == 1 for c in chains)


def test_api_structure_download(client, synthetic_path):
    job = client.post("/bg/jobs", json={"source": "file", "path": str(synthetic_path)}).json()["summary"]["job"]
    r = client.get(f"/bg/jobs/{job}/structure/cif")
    assert r.status_code == 200 and r.text.startswith("data_") and ".cif" in r.headers["content-disposition"]
    r = client.get(f"/bg/jobs/{job}/structure/pdb")
    assert r.status_code == 200 and "ATOM" in r.text and "HETATM" in r.text
    assert client.get(f"/bg/jobs/{job}/structure/xyz").status_code == 400


# --------------------------------------------------------------------------- binding contacts, modes, centring

def test_binding_contacts(synthetic, hhb):
    pairs, typ, d = features.binding_contacts(synthetic)
    got = {(synthetic.res_name[synthetic.atom_res[a]], synthetic.res_name[synthetic.atom_res[b]], int(t)) for (a, b), t in zip(pairs, typ)}
    assert ("ZN", "HOH", features.BIND_METAL) in got               # ion-shell water, 2.1 Å
    assert ("ACT", "HOH", features.BIND_WATER) in got               # water H-bonded to the acetate O, 2.9 Å
    assert all(synthetic.atom_res[a] != synthetic.atom_res[b] for a, b in pairs)
    pairs, typ, d = features.binding_contacts(hhb)
    fe = [(hhb.residue_label(hhb.atom_res[b]), hhb.atom_name[b]) for (a, b), t in zip(pairs, typ)
          if t == features.BIND_METAL and hhb.atom_name[a] == "FE"]
    assert sorted(fe) == [("A:HIS87", "NE2"), ("B:HIS92", "NE2"), ("C:HIS87", "NE2"), ("D:HIS92", "NE2")]


def test_structure_layer_ships_binding_and_point_atoms(synthetic_path):
    job = engine.open_job(engine.JobSpec(source="file", path=str(synthetic_path), level="ca", ligands=True))
    meta = engine.build_layer(job, "structure")
    assert {"bind_pairs", "bind_type", "bind_dist", "pt_atom"} <= set(meta["buffers"])
    assert meta["stats"]["mode"] == "Cα" and meta["stats"]["hetero_points"] is True
    assert sum(meta["stats"]["bind_counts"]) == meta["buffers"]["bind_type"]["shape"][0] > 0


def test_selection_points_are_centred_like_the_structure(hhb_job):
    """Regression: the selection overlay was drawn one centroid away from the selection."""
    hem = next(h for h in hhb_job.hetero if h["label"] == "HEM (i)")
    meta = engine.build_layer(hhb_job, "selection", {"residues": [hem["res"]] + list(range(80, 95))})
    pts = np.frombuffer(engine.buffer_path(hhb_job.id, meta["key"], "pts").read_bytes(), dtype="<f4").reshape(-1, 3)
    res = np.frombuffer(engine.buffer_path(hhb_job.id, meta["key"], "pt_res").read_bytes(), dtype="<u4")
    s = hhb_job.structure
    atoms = np.flatnonzero(np.isin(s.atom_res, res))
    assert np.allclose(np.sort(pts, axis=0), np.sort(s.coords[atoms] - hhb_job.center, axis=0), atol=1e-3)


def test_pockets_ignore_hetero_points():
    """With ligands in the point set the haem would fill its own pocket: pockets use the polymer alone, so
    they match the pockets of the same structure loaded without hetero points."""
    from biogeometry.structure import StructureError
    try:
        with_het = engine.open_job(engine.JobSpec(source="rcsb", pdb_id="4HHB", level="heavy", ligands=True))
        polymer = engine.open_job(engine.JobSpec(source="rcsb", pdb_id="4HHB", level="heavy", ligands=False))
    except StructureError as exc:
        pytest.skip(str(exc))
    a = engine.build_layer(with_het, "pockets")["stats"]
    b = engine.build_layer(polymer, "pockets")["stats"]
    assert a["polymer_only"] is True and b["polymer_only"] is False
    assert len(with_het.pocket_tess.points) == len(polymer.points.coords) < len(with_het.points.coords)
    assert [round(p["volume"], 3) for p in a["pockets"]] == [round(p["volume"], 3) for p in b["pockets"]]
    assert a["n_pockets"] > 0


# --------------------------------------------------------------------------- bounding-box Voronoi, depth labels

def test_box_voronoi_tiles_the_expanded_box():
    """Mirroring across the expanded bounding box clips every cell to the box: all cells close, they tile
    the box exactly, and interior cells are unchanged."""
    from biogeometry import geometry as g
    rng = np.random.default_rng(3)
    P = rng.uniform(0, 20, (600, 3))
    margin = 4.0
    boxed = g.voronoi_cells(P, 3.0, box_margin=margin)
    assert np.isfinite(boxed.volume).all() and np.isfinite(boxed.area).all() and boxed.closed.all()
    side = P.max(axis=0) - P.min(axis=0) + 2 * margin
    assert boxed.volume.sum() == pytest.approx(np.prod(side), rel=1e-6)
    raw = g.voronoi_cells(P, 3.0)
    both = raw.bounded & boxed.bounded
    assert both.sum() > 50 and np.allclose(raw.volume[both], boxed.volume[both])
    assert np.isnan(raw.volume).any()                       # without the box, surface cells stay open


def test_every_residue_gets_a_voronoi_volume(hhb_job):
    t = features.residue_table(hhb_job)
    col = {c["key"]: c for c in t.columns}
    prot = hhb_job.structure.res_kind == KIND_PROTEIN
    assert np.isfinite(t.data["vor_volume"][prot]).all()
    inter = np.array(col["vor_interior"]["categories"])[t.data["vor_interior"]]
    assert {"interior", "surface"} <= set(inter[prot])
    assert np.isnan(t.data["packing"][prot][inter[prot] == "surface"]).all()   # packing only for interior residues
    assert col["depth_class"]["categories"] == ["S", "U", "B"]
    T, _ = features.tetra_table(hhb_job)
    assert np.isfinite(T.data["voronoi-volume-1"]).all() and np.isfinite(T.data["voronoi-area-4"]).all()


def test_depth_code_and_participating_tetrahedra(hhb_job):
    t, T = features.tetra_table(hhb_job)
    col = {c["key"]: c for c in t.columns}
    cat = lambda k, i: col[k]["categories"][t.data[k][i]]  # noqa: E731
    for i in range(0, t.n, max(1, t.n // 50)):
        assert cat("depth", i) == "".join(cat(f"residue-{k}-depth", i) for k in (1, 2, 3, 4))
    hem = next(h for h in hhb_job.hetero if h["label"] == "HEM (i)")
    res = [hem["res"], 86, 87]
    st = engine.build_layer(hhb_job, "selection", {"residues": res, "max_edge": 5.0})["stats"]["summary"]["participating"]
    # the same tetrahedra as the table's rows that involve those residues
    involved = np.isin(T["res4"], res).any(axis=1).sum()
    assert st["n_tets"] == involved > 0 and st["n_nodes"] > 4 and st["n_edges"] > st["n_nodes"]
    assert st["n_neighbour_residues"] > 0 and st["volume"] > 0
