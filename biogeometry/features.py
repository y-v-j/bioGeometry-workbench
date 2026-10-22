"""Feature tables for the workbench: one row per residue and one row per residue-level tetrahedron.

Every column carries a descriptor (label, group, unit, type, explanation) that the workbench uses
for headers, tooltips and filters. Numeric columns travel as float32 buffers; categorical columns as
uint16 codes plus their category list, so even a ribosome's table is a few megabytes.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from . import analysis as an
from . import biochem
from . import geometry as g
from .sasa import atom_sasa, radii_for
from .structure import ELEMENTS, KIND_ION, KIND_LIGAND, KIND_NUCLEIC, KIND_PROTEIN, KIND_WATER, KIND_NAMES, METALS, chain_tag

ANIONIC = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}
CATIONIC = {("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"), ("HIS", "ND1"), ("HIS", "NE2")}
SALT_BRIDGE = 4.0
DISULFIDE = 2.5
CONTACT = 4.5          # heavy-atom contact for interfaces and ligand sites


@dataclass
class Table:
    columns: list[dict] = field(default_factory=list)   # descriptors, in display order
    data: dict[str, np.ndarray] = field(default_factory=dict)
    n: int = 0

    def num(self, key, values, label, group, unit="", desc="", fmt=".2f", hidden=False):
        self.columns.append({"key": key, "label": label, "group": group, "unit": unit, "type": "num",
                             "desc": desc, "fmt": fmt, "hidden": hidden})
        self.data[key] = np.asarray(values, dtype=np.float32)

    def cat(self, key, values, label, group, desc="", order=None, hidden=False):
        vals = np.asarray(values, dtype=object).astype(str)
        cats = list(order) if order else []
        for v in dict.fromkeys(vals):
            if v not in cats:
                cats.append(v)
        lookup = {c: i for i, c in enumerate(cats)}
        self.columns.append({"key": key, "label": label, "group": group, "unit": "", "type": "cat",
                             "desc": desc, "categories": cats, "hidden": hidden})
        self.data[key] = np.fromiter((lookup[v] for v in vals), dtype=np.uint16, count=len(vals))

    def buffers(self, prefix: str) -> dict[str, np.ndarray]:
        return {f"{prefix}{k}": v for k, v in self.data.items()}


# ------------------------------------------------------------------ interactions

def _pairs_within(coords, mask_a, mask_b, r):
    """Atom index pairs (a in A, b in B) within r Å."""
    ia, ib = np.flatnonzero(mask_a), np.flatnonzero(mask_b)
    if not len(ia) or not len(ib):
        return np.zeros((0, 2), dtype=np.int64)
    ta, tb = cKDTree(coords[ia]), cKDTree(coords[ib])
    hits = ta.query_ball_tree(tb, r)
    cnt = np.fromiter((len(h) for h in hits), dtype=np.int64, count=len(hits))
    a = np.repeat(ia, cnt)
    b = ib[np.fromiter((j for h in hits for j in h), dtype=np.int64, count=int(cnt.sum()))]
    return np.column_stack([a, b])


def interactions(s, hetero, bb_hbonds):
    """Per-residue interaction descriptors (lists/strings) from heavy-atom geometry."""
    R = s.n_residues
    names = s.atom_name.astype(str)
    resn = s.res_name.astype(str)[s.atom_res]
    keys = list(zip(resn, names))
    anion = np.fromiter((k in ANIONIC for k in keys), dtype=bool, count=len(keys))
    cation = np.fromiter((k in CATIONIC for k in keys), dtype=bool, count=len(keys))
    salt = np.zeros(R, dtype=np.int32)
    pairs = _pairs_within(s.coords, anion, cation, SALT_BRIDGE)
    if len(pairs):
        rr = np.unique(np.sort(s.atom_res[pairs], axis=1), axis=0)
        rr = rr[rr[:, 0] != rr[:, 1]]
        np.add.at(salt, rr.ravel(), 1)

    sg = (names == "SG") & (resn == "CYS")
    disulf = np.full(R, "", dtype=object)
    pairs = _pairs_within(s.coords, sg, sg, DISULFIDE)
    for a, b in pairs:
        ra, rb = s.atom_res[a], s.atom_res[b]
        if ra != rb:
            disulf[ra] = f"{s.res_chain[rb]}:CYS{s.res_seq[rb]}"

    poly = np.isin(s.res_kind, (KIND_PROTEIN, KIND_NUCLEIC))
    poly_atom = poly[s.atom_res]
    iface = [set() for _ in range(R)]
    tree = cKDTree(s.coords[poly_atom])
    pidx = np.flatnonzero(poly_atom)
    pp = tree.query_pairs(CONTACT, output_type="ndarray")
    if len(pp):
        a, b = pidx[pp[:, 0]], pidx[pp[:, 1]]
        ca, cb = s.res_chain[s.atom_res[a]], s.res_chain[s.atom_res[b]]
        diff = ca != cb
        ra, rb = s.atom_res[a[diff]], s.atom_res[b[diff]]
        for x, y, cx, cy in zip(ra, rb, ca[diff], cb[diff]):
            iface[x].add(cy)
            iface[y].add(cx)

    site = [set() for _ in range(R)]
    label_of = {}
    for h in hetero:
        for m in h.get("members", [h["res"]]):
            label_of[m] = h["label"] if h["kind"] != "water" else "H₂O"
    het_atom = np.isin(s.res_kind, (KIND_LIGAND, KIND_ION))[s.atom_res]
    pairs = _pairs_within(s.coords, het_atom, poly_atom, CONTACT)
    for a, b in pairs:
        hr, pr = s.atom_res[a], s.atom_res[b]
        site[pr].add(label_of.get(hr, str(s.res_name[hr])))
        site[hr].add(f"{s.res_chain[pr]}")
    hb = np.zeros(R, dtype=np.int32)
    if len(bb_hbonds):
        np.add.at(hb, bb_hbonds.ravel(), 1)
    return {
        "salt": salt, "disulfide": disulf, "hbonds": hb,
        "interface": np.array([",".join(sorted(x)) for x in iface], dtype=object),
        "site": np.array([", ".join(sorted(x)) for x in site], dtype=object),
    }


# ------------------------------------------------------------------ binding contacts (structure view)

BIND_METAL, BIND_HBOND, BIND_WATER, BIND_COVALENT = 0, 1, 2, 3
BIND_NAMES = ("metal coordination", "hydrogen bond", "water bridge", "covalent link")
BIND_METAL_DIST = 3.0       # ion / metal atom ··· N, O, S (Mg-O 2.1, Zn-N 2.1, Na-O 2.4, K-O 2.8)
BIND_HBOND_DIST = 3.5       # N/O ··· N/O
BIND_COVALENT_DIST = 2.0    # ligand ··· polymer heavy atoms closer than this are bonded (glycans, PLP-Lys, adducts)


def binding_contacts(s) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Atom pairs that bind ligands, ions and coordinated waters to the rest of the structure.

    Returns (pairs (K,2) [hetero atom, partner atom], type (K,) BIND_*, distance (K,)). Pairs inside one
    residue, water-water and polymer-polymer pairs are left out.
    """
    kind_atom = s.res_kind[s.atom_res]
    het = np.flatnonzero(kind_atom >= KIND_LIGAND)
    empty = (np.zeros((0, 2), np.int64), np.zeros(0, np.uint8), np.zeros(0, np.float32))
    if not len(het):
        return empty
    sym = np.array([ELEMENTS[e].upper() if e < len(ELEMENTS) else "X" for e in s.element], dtype=object)
    metal = np.isin(sym, list(METALS))
    polar = np.isin(sym, ("N", "O"))
    donor_acc = polar | (sym == "S")
    hits = cKDTree(s.coords).query_ball_point(s.coords[het], r=BIND_HBOND_DIST)
    a = np.repeat(het, [len(h) for h in hits])
    b = np.fromiter((j for h in hits for j in h), dtype=np.int64, count=len(a))
    keep = s.atom_res[a] != s.atom_res[b]
    both_het = kind_atom[b] >= KIND_LIGAND
    keep &= ~both_het | (a < b)                                   # hetero-hetero pairs once
    keep &= ~((kind_atom[a] == KIND_WATER) & (kind_atom[b] == KIND_WATER))
    a, b = a[keep], b[keep]
    d = np.linalg.norm(s.coords[a] - s.coords[b], axis=1)
    typ = np.full(len(a), 255, dtype=np.uint8)
    m = (metal[a] & donor_acc[b]) | (metal[b] & donor_acc[a])
    typ[m & (d <= BIND_METAL_DIST)] = BIND_METAL
    cov = ~m & ~metal[a] & ~metal[b] & (d <= BIND_COVALENT_DIST) & (kind_atom[a] == KIND_LIGAND) & np.isin(kind_atom[b], (KIND_PROTEIN, KIND_NUCLEIC))
    typ[cov] = BIND_COVALENT
    hb = (typ == 255) & ~m & polar[a] & polar[b]
    water = (kind_atom[a] == KIND_WATER) | (kind_atom[b] == KIND_WATER)
    typ[hb & water] = BIND_WATER
    typ[hb & ~water] = BIND_HBOND
    ok = typ != 255
    return np.stack([a[ok], b[ok]], axis=1), typ[ok], d[ok].astype(np.float32)


# ------------------------------------------------------------------ residue table

def residue_table(job) -> Table:
    s = job.structure
    R = s.n_residues
    names = s.res_name.astype(str)
    kinds = s.res_kind
    t = Table(n=R)
    is_af = s.is_alphafold

    # ---- identity
    t.cat("chain", s.res_chain, "chain", "identity", "author chain id")
    t.num("resi", s.res_seq, "res #", "identity", desc="author residue number", fmt="d")
    t.cat("resn", names, "residue", "identity", "residue / ligand name (3-letter code)")
    t.cat("aa", biochem.one_letter(names), "aa", "identity", "one-letter code", hidden=True)
    t.cat("kind", [KIND_NAMES[k] for k in kinds], "kind", "identity",
          "protein, nucleic acid, ligand, ion or coordinated water", order=KIND_NAMES)

    # ---- biochemistry
    t.cat("class", biochem.residue_class(names, kinds), "class", "biochemistry",
          "side-chain chemistry: aliphatic, aromatic, polar, acidic, basic, special (G/P/C); purine/pyrimidine for nucleotides",
          order=biochem.CLASS_ORDER)
    t.num("charge", biochem.lookup(biochem.CHARGE, names, 0.0) * (kinds == KIND_PROTEIN), "charge", "biochemistry",
          "e", "formal side-chain charge at pH 7 (His counted neutral)", fmt="+.0f")
    t.num("hydropathy", biochem.lookup(biochem.KYTE_DOOLITTLE, names), "hydropathy", "biochemistry", "",
          "Kyte–Doolittle hydropathy: > 0 hydrophobic, < 0 hydrophilic", fmt="+.1f")
    t.num("ref_volume", biochem.lookup(biochem.VOLUME, names), "ref. volume", "biochemistry", "Å³",
          "residue volume in solution (Zamyatnin 1972) - the reference for the packing ratio", fmt=".0f", hidden=True)

    # ---- backbone and secondary structure
    ss = job.ss
    ss_txt = np.where(ss >= 0, np.array(["C", "H", "E"])[np.maximum(ss, 0)], "")
    phi, psi, omega = job.dihedrals
    region = biochem.rama_region(phi, psi)
    t.cat("ss", ss_txt, "SS", "backbone", "secondary structure from backbone H-bonds (Kabsch–Sander): "
          "H helix (α/3₁₀/π), E strand, C coil", order=["H", "E", "C", ""])
    t.num("phi", phi, "φ", "backbone", "°", "backbone dihedral C(i-1)–N–Cα–C", fmt=".0f")
    t.num("psi", psi, "ψ", "backbone", "°", "backbone dihedral N–Cα–C–N(i+1)", fmt=".0f")
    t.num("omega", omega, "ω", "backbone", "°", "peptide-bond dihedral preceding the residue (≈180 trans, ≈0 cis)",
          fmt=".0f", hidden=True)
    t.cat("rama", region, "Rama region", "backbone",
          "Ramachandran basin: αR right-handed helix, β strand, PPII polyproline II, αL left-handed helix, other",
          order=list(biochem.RAMA_REGIONS) + [""])
    t.cat("rama_outlier", np.where(biochem.rama_outlier(region, names), "outlier", ""), "Rama outlier", "backbone",
          "non-glycine residue outside every populated basin - strained or mis-modelled backbone")
    t.cat("cis", np.where(np.abs(omega) < 30, "cis", ""), "cis peptide", "backbone",
          "cis peptide bond (|ω| < 30°); common before Pro, rare and often functional otherwise")

    # ---- confidence
    bres = an.per_residue(s.bfactor, s.atom_res, R, "mean")
    if is_af:
        t.num("bfactor", bres, "pLDDT", "confidence", "",
              "AlphaFold per-residue confidence (0–100): > 90 very high, 70–90 confident, 50–70 low, < 50 very low",
              fmt=".1f")
        t.cat("plddt_band", biochem.plddt_band(bres), "pLDDT band", "confidence",
              "AlphaFold confidence band", order=["very high", "confident", "low", "very low", ""])
    else:
        t.num("bfactor", bres, "B-factor", "confidence", "Å²",
              "mean atomic displacement (B) factor: high = mobile or poorly ordered", fmt=".1f")

    # ---- surface and burial
    sasa_atom = job.sasa
    sasa_res = np.bincount(s.atom_res, weights=sasa_atom, minlength=R)
    rsa = sasa_res / biochem.lookup(biochem.MAX_ASA, names)
    rsa[kinds >= KIND_LIGAND] = np.nan
    t.num("sasa", sasa_res, "SASA", "surface", "Å²", "solvent-accessible surface (Shrake–Rupley, 1.4 Å probe)",
          fmt=".1f")
    t.num("rsa", rsa, "RSA", "surface", "", "relative accessibility = SASA / maximum (Tien 2013); < 0.25 buried",
          fmt=".2f")
    t.cat("exposure", biochem.exposure_class(rsa), "exposure", "surface",
          "core (< 5 %), buried (< 25 %), intermediate (< 50 %), exposed", order=["core", "buried", "intermediate",
                                                                               "exposed", ""])
    depth, depth_cls = residue_depth(s, job.sasa)
    t.num("depth", depth, "depth", "surface", "Å", "residue depth: mean distance of its atoms to the nearest "
          "solvent-exposed atom (exposed atoms count 0); defined for ligands and waters too", fmt=".1f")
    t.cat("depth_class", np.asarray(DEPTH_CLASSES, dtype=object)[depth_cls], "depth class", "surface",
          f"S surface (RSA ≥ {SURFACE_RSA:.0%}), U just under the surface (depth ≤ {UNDER_SURFACE} Å), B buried",
          order=list(DEPTH_CLASSES))
    rp = job.residue_points
    cn = np.full(R, np.nan)
    if len(rp.coords) > 1:
        tree = cKDTree(rp.coords)
        cn[rp.residue] = np.array([len(x) - 1 for x in tree.query_ball_point(rp.coords, 10.0)])
    t.num("cn10", cn, "neighbours 10 Å", "surface", "", "residues with Cα/P within 10 Å (contact number): "
          "a coarse burial measure independent of the surface", fmt=".0f")

    # ---- computational geometry (at the job's point level)
    d = job.defaults()
    ps = job.points
    cells = job.voronoi(d["voronoi"]["max_radius"])
    depth = g.hull_depth(job.hull, ps.coords)
    alpha = g.alpha_shape(job.tess, d["alpha"]["alpha"])
    surf = np.zeros(len(ps.coords))
    surf[alpha.surface_points] = 1.0
    net = an.contact_network(s, job.tess, ps, d["contacts"]["cutoff"], d["contacts"]["min_seq_sep"])
    vol = an.per_residue(cells.volume, ps.residue, R, "sum")
    has_pt = np.bincount(ps.residue, minlength=R) > 0
    interior = (an.per_residue((~cells.bounded).astype(float), ps.residue, R, "sum") == 0) & has_pt
    t.num("vor_volume", np.where(has_pt, vol, np.nan), "Voronoi vol.", "geometry", "Å³",
          f"volume of the residue's Voronoi cell(s) at the {ps.level} level; every cell is closed by the bounding box "
          f"(expanded by {cells.box_margin:g} Å), so surface cells depend on that margin - see 'Voronoi interior'",
          fmt=".1f")
    t.cat("vor_interior", np.where(~has_pt, "", np.where(interior, "interior", "surface")), "Voronoi interior", "geometry",
          "interior: every cell closed by the structure itself (no vertex beyond the max cell reach); "
          "surface: at least one cell reaches out to the bounding box", order=["interior", "surface", ""])
    if ps.level == "heavy":
        t.num("packing", np.where(interior, vol / biochem.lookup(biochem.VOLUME, names), np.nan), "packing ratio", "geometry", "",
              "Voronoi volume / reference volume for interior residues: < 1 tightly packed, > 1 loosely packed "
              "or near a cavity", fmt=".2f")
    t.num("hull_depth", an.per_residue(depth, ps.residue, R, "min"), "hull depth", "geometry", "Å",
          "distance below the convex hull - how deep inside the molecule the residue sits", fmt=".1f")
    t.num("convex_layer", an.per_residue(job.layers.astype(float), ps.residue, R, "min"), "convex layer",
          "geometry", "", "onion-peeling layer of the convex hull (0 = outermost)", fmt=".0f")
    t.num("alpha_surface", an.per_residue(surf, ps.residue, R, "mean"), "α-surface", "geometry", "",
          f"fraction of the residue's points on the α-shape boundary (α = {d['alpha']['alpha']} Å)", fmt=".2f")
    t.num("del_degree", net["degree"].astype(float), "Delaunay contacts", "geometry", "",
          f"residues joined by Delaunay edges ≤ {d['contacts']['cutoff']} Å (|Δseq| ≥ {d['contacts']['min_seq_sep']}): "
          "non-occluded neighbours", fmt=".0f")

    # ---- interactions and annotations
    it = job.interactions
    t.num("hbonds", it["hbonds"], "backbone H-bonds", "interactions", "",
          "backbone N–H···O=C hydrogen bonds (DSSP energy < −0.5 kcal/mol), as donor or acceptor", fmt=".0f")
    t.num("salt_bridges", it["salt"], "salt bridges", "interactions", "",
          "Asp/Glu carboxylate ↔ Lys/Arg/His nitrogen pairs within 4 Å", fmt=".0f")
    t.cat("disulfide", it["disulfide"], "disulfide", "interactions", "partner of a Cys–Cys disulfide (SG–SG ≤ 2.5 Å)")
    t.cat("interface", it["interface"], "interface", "interactions",
          "other chains with a heavy atom within 4.5 Å - protein–protein / protein–nucleic interface")
    t.cat("site", it["site"], "ligand site", "interactions",
          "ligands / ions within 4.5 Å (binding-site residues); on a ligand row, the chains it touches")
    pk = np.full(R, np.nan)
    for p in reversed(job.default_pockets):  # largest pocket wins (listed first)
        pk[np.unique(ps.residue[p.points])] = p.id + 1
    t.num("pocket", pk, "pocket", "interactions", "",
          "largest two-probe pocket lined by the residue (see the Pockets analysis); blank if none", fmt=".0f")
    dom = np.full(R, "", dtype=object)
    for dm in job.domains:
        sl = slice(dm["r0"], dm["r1"] + 1)
        cur = dom[sl]
        dom[sl] = np.where(cur == "", f"{dm['source']}:{dm['name']}", cur + f" · {dm['source']}:{dm['name']}")
    t.cat("domain", dom, "domain", "annotation", "Pfam / CATH (PDB entries) or UniProt (AlphaFold) domain")
    return t


# ------------------------------------------------------------------ tetrahedron table

def tetra_metrics(points: np.ndarray, simplices: np.ndarray) -> dict[str, np.ndarray]:
    """The 20 classical per-tetrahedron metrics (edges, volume, faces, solid angles), vectorised.

    Order follows the original della CSV: L1..L6 = ab, ac, ad, bc, bd, cd; T1..T4 = faces abc, abd, acd, bcd;
    O1..O4 = solid angles at a, b, c, d.
    """
    p = points[simplices]
    a, b, c, dd = p[:, 0], p[:, 1], p[:, 2], p[:, 3]
    nrm = lambda v: np.linalg.norm(v, axis=1)  # noqa: E731
    L = np.stack([nrm(a - b), nrm(a - c), nrm(a - dd), nrm(b - c), nrm(b - dd), nrm(c - dd)], axis=1)
    vol = np.abs(np.einsum("ij,ij->i", a - dd, np.cross(b - dd, c - dd))) / 6.0
    tri = lambda x, y, z: 0.5 * nrm(np.cross(y - x, z - x))  # noqa: E731
    T = np.stack([tri(a, b, c), tri(a, b, dd), tri(a, c, dd), tri(b, c, dd)], axis=1)
    O = np.empty((len(p), 4))
    for k, (v, i, j, l) in enumerate([(0, 1, 2, 3), (1, 0, 2, 3), (2, 0, 1, 3), (3, 0, 1, 2)]):
        x, y, z = p[:, i] - p[:, v], p[:, j] - p[:, v], p[:, l] - p[:, v]
        nx, ny, nz = nrm(x), nrm(y), nrm(z)
        num = np.abs(np.einsum("ij,ij->i", x, np.cross(y, z)))
        den = nx * ny * nz + np.einsum("ij,ij->i", x, y) * nz + np.einsum("ij,ij->i", x, z) * ny \
            + np.einsum("ij,ij->i", y, z) * nx
        O[:, k] = 2.0 * np.arctan2(num, den)
    out = {f"L{i + 1}": L[:, i] for i in range(6)}
    out.update(SumL=L.sum(1), AvgL=L.mean(1), Tvol=vol)
    out.update({f"T{i + 1}_area": T[:, i] for i in range(4)})
    out["Tsurface_area"] = T.sum(1)
    out.update({f"O{i + 1}": O[:, i] for i in range(4)})
    return out


def simplex_type(res_idx: np.ndarray, chain: np.ndarray) -> np.ndarray:
    """Longest run of sequence-adjacent residues (same chain, consecutive indices) among the four: 0, 2, 3, 4."""
    srt = np.sort(res_idx, axis=1)
    ch = chain[srt]
    step = (np.diff(srt, axis=1) == 1) & (ch[:, 1:] == ch[:, :-1])
    best = np.zeros(len(srt), dtype=int)
    run = np.ones(len(srt), dtype=int)
    for k in range(3):
        run = np.where(step[:, k], run + 1, 1)
        best = np.maximum(best, run)
    return np.where(best >= 2, best, 0)


# ------------------------------------------------------------------ residue depth

DEPTH_CLASSES = ("S", "U", "B")    # surface, just under the surface, buried
SURFACE_RSA = 0.20          # polymer residues with this relative accessibility or more are on the surface
SURFACE_SASA = 5.0          # Å²: the same for ligands and ions (no reference maximum to scale by)
EXPOSED_ATOM = 1.0          # Å² for an atom to count as solvent-exposed
UNDER_SURFACE = 2.0         # Å of mean depth: deeper residues are buried


def residue_depth(s, sasa: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(depth Å, class code into DEPTH_CLASSES) per residue.

    Depth is the mean, over the residue's atoms, of the distance to the nearest solvent-exposed atom
    (an exposed atom counts 0), in the spirit of residue depth (Chakravarty & Varadarajan 1999). A residue
    is on the surface when its relative accessibility is ≥ 20 % (≥ 5 Å² of SASA for ligands and ions);
    otherwise it is just under the surface up to 2 Å of mean depth, and buried beyond. Coordinated waters
    are classed by depth alone.
    """
    R = s.n_residues
    res_sasa = np.bincount(s.atom_res, weights=sasa, minlength=R)
    exposed = np.flatnonzero(sasa > EXPOSED_ATOM)
    counts = np.maximum(np.bincount(s.atom_res, minlength=R), 1)
    if len(exposed):
        d, _ = cKDTree(s.coords[exposed]).query(s.coords, k=1)
        d[sasa > EXPOSED_ATOM] = 0.0
    else:
        d = np.full(s.n_atoms, 99.0)
    depth = np.bincount(s.atom_res, weights=d, minlength=R) / counts
    max_asa = biochem.lookup(biochem.MAX_ASA, s.res_name.astype(str))
    rsa = np.where(np.isfinite(max_asa) & (max_asa > 0), res_sasa / np.where(max_asa > 0, max_asa, 1), np.nan)
    surface = np.where(np.isfinite(rsa), rsa >= SURFACE_RSA, res_sasa >= SURFACE_SASA) & (s.res_kind != KIND_WATER)
    cls = np.where(surface, 0, np.where(depth <= UNDER_SURFACE, 1, 2)).astype(np.uint8)
    return depth, cls


# ------------------------------------------------------------------ tetrahedra table (the point set of the job)

MAX_TETRA_ROWS = 400_000
EC_CLASSES = {"1": "Oxidoreductases", "2": "Transferases", "3": "Hydrolases", "4": "Lyases",
              "5": "Isomerases", "6": "Ligases", "7": "Translocases"}
# runs of sequence-contiguous vertices -> (Tropsha class, simplex type), indexed by the three 'link' bits
# between the sorted vertices (bit k set when vertex k and k+1 are the same residue or i, i+1 in one chain)
_RUNS = {0: ("{1,1,1,1}", 0), 1: ("{2,1,1}", 2), 2: ("{2,1,1}", 2), 4: ("{2,1,1}", 2),
         3: ("{3,1}", 3), 6: ("{3,1}", 3), 5: ("{2,2}", 2), 7: ("{4}", 4)}
_INTERACTION_ORDER = ("protein", "DNA", "RNA", "ligand", "ion", "water")


def tetrahedrality(L: np.ndarray) -> np.ndarray:
    """Σ_{i<j} (l_i − l_j)² / (15 l̄²) over the six edges: 0 for a regular tetrahedron."""
    mean = L.mean(axis=1)
    s = np.zeros(len(L))
    for i in range(6):
        for j in range(i + 1, 6):
            s += (L[:, i] - L[:, j]) ** 2
    return s / (15.0 * mean ** 2)


def _interaction_type(cls4: np.ndarray, chain4: np.ndarray) -> np.ndarray:
    """protein-only, protein-protein, protein-ligand, protein-ion, protein-DNA, protein-RNA, RNA-ion …"""
    names = {0: "protein", 1: "DNA", 2: "RNA", 3: "ligand", 4: "ion", 5: "water"}
    cache, out = {}, np.empty(len(cls4), dtype=object)
    for i in range(len(cls4)):
        key = (tuple(cls4[i]), tuple(chain4[i]))
        lab = cache.get(key)
        if lab is None:
            parts = []
            for c in (0, 1, 2):                                  # polymers: once per distinct chain (max twice)
                n = len({chain4[i][k] for k in range(4) if cls4[i][k] == c})
                parts += [names[c]] * min(n, 2)
            parts += [names[c] for c in (3, 4, 5) if c in key[0]]
            lab = parts[0] + "-only" if len(parts) == 1 else "-".join(parts)
            cache[key] = lab
        out[i] = lab
    return out


def tetra_rows(job, max_edge: float | None = None) -> dict:
    """Everything the tetrahedra table and its CSV need, for the compact Delaunay tetrahedra of the job's
    point set (every edge ≤ max_edge; ligand, ion and coordinated-water atoms are vertices like any other)."""
    from .annotations import chain_ec
    s, ps, tess = job.structure, job.points, job.tess
    d = job.defaults()
    max_edge = float(max_edge or d["delaunay"]["max_edge"])
    keep = np.flatnonzero(tess.max_edge_per_tet() <= max_edge)
    sampled = len(keep) > MAX_TETRA_ROWS
    if sampled:
        keep = np.sort(np.random.default_rng(0).choice(keep, MAX_TETRA_ROWS, replace=False))
    tets = np.sort(tess.simplices[keep], axis=1)                  # vertices in point (≈ sequence) order
    res4 = ps.residue[tets]
    coarse = ps.level in ("ca", "centroid")
    m = tetra_metrics(tess.points, tets)
    L = np.stack([m[f"L{i + 1}"] for i in range(6)], axis=1)
    circ = tess.circumradii[keep]

    # per point: Voronoi cell, SASA (atom, or the residue's for Cα/centroid polymer points), B (residue mean)
    vor = job.voronoi(d["voronoi"]["max_radius"])
    atom = ps.atom
    res_sasa = np.bincount(s.atom_res, weights=job.sasa, minlength=s.n_residues)
    poly_pt = np.isin(s.res_kind[ps.residue], (KIND_PROTEIN, KIND_NUCLEIC))
    pt_sasa = np.where(coarse & poly_pt, res_sasa[ps.residue], job.sasa[np.maximum(atom, 0)])
    res_b = an.per_residue(s.bfactor, s.atom_res, s.n_residues, "mean")
    depth, depth_cls = residue_depth(s, job.sasa)

    # per residue: class code, chain, contiguity
    names = s.res_name.astype(str)
    o2 = np.zeros(s.n_residues, dtype=bool)
    o2[s.atom_res[s.atom_name.astype(str) == "O2'"]] = True
    vclass = biochem.vertex_class(names, s.res_kind, o2)
    chain = s.chain_index()
    poly = np.isin(s.res_kind, (KIND_PROTEIN, KIND_NUCLEIC))
    r = res4
    links = np.zeros((len(r), 3), dtype=bool)
    for k in range(3):
        a, b = r[:, k], r[:, k + 1]
        links[:, k] = (a == b) | ((b - a == 1) & (chain[a] == chain[b]) & poly[a] & poly[b])
    code = links[:, 0] * 1 + links[:, 1] * 2 + links[:, 2] * 4
    tropsha = np.array([_RUNS[c][0] for c in range(8)], dtype=object)[code]
    stype = np.array([_RUNS[c][1] for c in range(8)])[code]

    cls4, ch4 = vclass[r], chain[r]
    vnames = np.where(vclass == 5, "HOH", names)                  # every water is written as HOH
    rclass = np.array(["-".join(biochem.VERTEX_CLASSES[c] for c in row) for row in cls4], dtype=object)
    rtype = np.array(["-".join(row) for row in vnames[r]], dtype=object)
    chains = np.asarray(s.chains, dtype=object)
    sep = "" if all(len(c) == 1 for c in s.chains) else "/"
    pchains = np.array([sep.join(row) for row in chains[ch4]], dtype=object)
    # het ids: the distinct non-polymer residue names among the vertices, in vertex order (HEM, NA-MG, NA-HEM-HOH)
    het = np.isin(s.res_kind, (KIND_LIGAND, KIND_ION, KIND_WATER))
    hetid = np.full(len(r), "", dtype=object)
    for i in np.flatnonzero(het[r].any(axis=1)):
        hetid[i] = "-".join(dict.fromkeys(vnames[x] for x in r[i] if het[x]))
    ec_of = chain_ec(s)
    ec_chain = np.array([ec_of.get(c, "") for c in s.chains], dtype=object)
    ec = np.array([";".join(dict.fromkeys(x for x in ec_chain[row] if x)) for row in ch4], dtype=object)
    ec_cat = np.array([";".join(dict.fromkeys(EC_CLASSES.get(e.split(".")[0], "") for e in x.split(";") if e)) if x else ""
                       for x in ec], dtype=object)

    alpha = float(d["alpha"]["alpha"])
    ash = g.alpha_shape(tess, alpha)
    hull = job.hull
    tag = chain_tag(s)
    return {
        "n": len(keep), "n_compact": int((tess.max_edge_per_tet() <= max_edge).sum()), "n_tets": tess.n_tets,
        "sampled": sampled, "max_edge": max_edge, "tets": tets, "res4": r,
        "consts": {"pdbid-chainid": f"{s.name}-{tag}", "model": int(s.meta.get("model", 1)),
                   "atom-set": "c-alpha" if coarse else "all-atom", "alpha-value": alpha,
                   "alpha-shape-volume": ash.volume, "alpha-shape-surface-area": ash.area,
                   "hull-volume": hull.volume, "hull-surface-area": hull.area},
        "ec": ec, "ec_cat": ec_cat, "rclass": rclass, "rtype": rtype, "pchains": pchains, "hetid": hetid,
        "L": L, "m": m, "circ": circ, "tetrahedrality": tetrahedrality(L), "tropsha": tropsha, "stype": stype,
        "vor_vol": vor.volume[tets], "vor_area": vor.area[tets], "sasa": pt_sasa[tets],
        "depth_cls": depth_cls[r], "depth": depth, "b": res_b[r],
        "inside": np.where(circ <= alpha, "inside", "outside"),
        "interaction": _interaction_type(cls4, ch4),
        "res_label": np.array([f"{n}{q}{(i or '').strip()}" for n, q, i in zip(vnames, s.res_seq, s.res_icode)], dtype=object),
    }


def tetra_table(job, max_edge: float | None = None) -> tuple[Table, dict]:
    """The tetrahedra table in the column order of the reference output (see README)."""
    T = tetra_rows(job, max_edge)
    n, c = T["n"], T["consts"]
    t = Table(n=n)
    const = lambda v: np.full(n, v, dtype=object)  # noqa: E731
    t.cat("pdbid-chainid", const(c["pdbid-chainid"]), "pdbid-chainid", "identity", "PDB / AlphaFold id and the chains analysed", hidden=True)
    t.num("model", np.full(n, c["model"]), "model", "identity", desc="model number in the file", fmt="d", hidden=True)
    t.cat("ec-number", T["ec"], "ec-number", "identity", "EC number of the enzyme chain(s) (mmCIF entity or UniProt); blank if none")
    t.cat("ec-category", T["ec_cat"], "ec-category", "identity", "enzyme class from the first EC digit (1 oxidoreductases … 7 translocases)")
    t.cat("atom-set", const(c["atom-set"]), "atom-set", "identity", "c-alpha or all-atom point set (+ ligands, ions, coordinated waters)", hidden=True)
    t.cat("residue-class", T["rclass"], "residue-class", "identity", "prt / dna / rna / lig / ion / hoh for each vertex")
    t.cat("residue-type", T["rtype"], "residue-type", "identity", "residue names of the four vertices (waters as HOH)")
    t.cat("participating-chains", T["pchains"], "participating-chains", "identity", "chain of each vertex")
    t.cat("het-id", T["hetid"], "het-id", "identity",
          "names of the ligands, ions and waters among the vertices, each once, e.g. MG, NA-MG, NA-HEM-HOH; blank if none")
    for i in range(6):
        t.num(f"edge-len-{i + 1}", T["L"][:, i], f"edge-len-{i + 1}", "edges", "Å", "edge length " + ["ab", "ac", "ad", "bc", "bd", "cd"][i])
    t.num("sum-edge-len", T["m"]["SumL"], "sum-edge-len", "edges", "Å", "sum of the six edge lengths")
    t.num("avg-edge-len", T["m"]["AvgL"], "avg-edge-len", "edges", "Å", "mean edge length")
    for i in range(4):
        t.num(f"triangle-{i + 1}-area", T["m"][f"T{i + 1}_area"], f"triangle-{i + 1}-area", "faces", "Å²",
              "area of face " + ["abc", "abd", "acd", "bcd"][i])
    t.num("surface-area-tetrahedron", T["m"]["Tsurface_area"], "surface-area-tetrahedron", "faces", "Å²", "total face area")
    t.num("tetrahedron-volume", T["m"]["Tvol"], "tetrahedron-volume", "volume", "Å³", "volume (scalar triple product)")
    for i in range(4):
        t.num(f"omega-{i + 1}", T["m"][f"O{i + 1}"], f"omega-{i + 1}", "angles", "sr", f"solid angle at vertex {'abcd'[i]}", fmt=".3f")
    for i in range(4):
        t.num(f"voronoi-volume-{i + 1}", T["vor_vol"][:, i], f"voronoi-volume-{i + 1}", "voronoi", "Å³",
              f"Voronoi cell volume of vertex {'abcd'[i]} (cells closed by the bounding box; surface cells depend on its margin)")
    for i in range(4):
        t.num(f"voronoi-area-{i + 1}", T["vor_area"][:, i], f"voronoi-area-{i + 1}", "voronoi", "Å²",
              f"Voronoi cell surface area of vertex {'abcd'[i]} (cells closed by the bounding box)")
    for i in range(4):
        t.num(f"sasa-{i + 1}", T["sasa"][:, i], f"sasa-{i + 1}", "surface", "Å²",
              f"solvent-accessible surface of vertex {'abcd'[i]} (the atom; the whole residue for Cα points)")
    for i in range(4):
        t.cat(f"residue-{i + 1}-depth", np.asarray(DEPTH_CLASSES, dtype=object)[T["depth_cls"][:, i]], f"residue-{i + 1}-depth", "surface",
              f"residue of vertex {'abcd'[i]}: S surface (RSA ≥ {SURFACE_RSA:.0%}), U just under the surface (depth ≤ {UNDER_SURFACE} Å), B buried",
              order=list(DEPTH_CLASSES))
    codes = np.asarray(DEPTH_CLASSES, dtype=object)[T["depth_cls"]]
    t.cat("depth", np.array(["".join(row) for row in codes], dtype=object), "depth", "surface",
          "the four residue depths in vertex order, e.g. SSUU (S surface, U just under the surface, B buried)")
    t.num("circumradius", np.minimum(T["circ"], 1e4), "circumradius", "shape", "Å",
          "radius of the sphere through the four vertices (the α at which the tetrahedron enters the α-shape)")
    t.num("circumsphere-volume", 4.0 / 3.0 * np.pi * np.minimum(T["circ"], 1e4) ** 3, "circumsphere-volume", "shape", "Å³",
          "volume of the sphere through the four vertices, 4/3 π r³ with r the circumradius")
    t.num("tetrahedrality", T["tetrahedrality"], "tetrahedrality", "shape", "",
          "Σ(lᵢ − lⱼ)² / 15 l̄² over the six edges: 0 for a regular tetrahedron, larger when distorted", fmt=".3f")
    t.cat("tropsha-class", T["tropsha"], "tropsha-class", "shape",
          "Singh-Tropsha-Vaisman class: sizes of the sequence-contiguous runs among the four residues",
          order=["{1,1,1,1}", "{2,1,1}", "{2,2}", "{3,1}", "{4}"])
    t.num("simplex-type", T["stype"], "simplex-type", "shape", "",
          "longest run of sequence-contiguous residues (0 = none, 2, 3, 4 = a backbone stretch)", fmt="d")
    t.num("alpha-value", np.full(n, c["alpha-value"]), "alpha-value", "alpha", "Å", "α of the α-shape", hidden=True)
    t.cat("alpha-shape", T["inside"], "alpha-shape", "alpha", "inside the α-shape (circumradius ≤ α) or outside", order=["inside", "outside"])
    t.num("alpha-shape-volume", np.full(n, c["alpha-shape-volume"]), "alpha-shape-volume", "alpha", "Å³", "volume of the whole α-shape", hidden=True)
    t.num("alpha-shape-surface-area", np.full(n, c["alpha-shape-surface-area"]), "alpha-shape-surface-area", "alpha", "Å²", "area of the whole α-shape", hidden=True)
    t.num("hull-volume", np.full(n, c["hull-volume"]), "hull-volume", "alpha", "Å³", "convex-hull volume of the point set", hidden=True)
    t.num("hull-surface-area", np.full(n, c["hull-surface-area"]), "hull-surface-area", "alpha", "Å²", "convex-hull area of the point set", hidden=True)
    for i in range(4):
        t.num(f"bfactor-residue-{i + 1}", T["b"][:, i], f"bfactor-residue-{i + 1}", "confidence", "Å²",
              f"mean B-factor (pLDDT for AlphaFold) of the residue of vertex {'abcd'[i]}")
    t.num("sum-bfactor", T["b"].sum(axis=1), "sum-bfactor", "confidence", "Å²", "sum of the four residue B-factors")
    t.num("avg-bfactor", T["b"].mean(axis=1), "avg-bfactor", "confidence", "Å²", "mean of the four residue B-factors")
    t.cat("interaction-type", T["interaction"], "interaction-type", "interactions",
          "protein-only, protein-protein, protein-ligand, protein-ion, protein-DNA, protein-RNA, …")
    return t, T
