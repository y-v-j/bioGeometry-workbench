"""Geometry of a user selection: domains, chains, ligands and ions picked in the workbench.

The selection is split into groups (a domain or chain segment, each ligand / ion instance, the
coordinated waters) and measured three ways:

* as a whole: convex hull, α-shape, Delaunay complex, gyration shape, solvent accessibility;
* per group: size, hull, shape, SASA alone / inside the selection / inside the whole structure;
* per pair of groups: minimum and centroid distance, atom contacts, polar (H-bond candidate) and
  metal-coordination contacts, Delaunay contacts and the surface buried on association.
"""
from __future__ import annotations

import re

import numpy as np
from scipy.spatial import cKDTree

from . import geometry as g
from .sasa import atom_sasa, radii_for
from .structure import ELEMENTS, KIND_NAMES, KIND_NUCLEIC, KIND_PROTEIN, KIND_WATER

MAX_SELECTION_ATOMS = 60_000    # larger selections represent the polymer by its Cα / P atoms
MAX_PAIR_GROUPS = 16            # pairwise metrics for the first groups only
CONTACT = 4.0                   # heavy-atom contact (Å)
POLAR = 3.5                     # N/O ··· N/O: hydrogen-bond candidate
METAL = 2.8                     # ion ··· N/O/S: first-shell coordination
NEAR = 8.0                      # pairs farther apart than this skip contact and buried-surface work
MAX_CONTACT_LINES = 20_000
MAX_BSA_ATOMS = 40_000          # buried surface is skipped for pairs larger than this
POLYMER = (KIND_PROTEIN, KIND_NUCLEIC)


def _elements(s, atoms):
    return np.array([ELEMENTS[e] if e < len(ELEMENTS) else "X" for e in s.element[atoms]], dtype=object)


def selection_groups(job, residues: np.ndarray) -> list[dict]:
    """Split residue indices into labelled groups: domains (≥ half of the domain selected), chain
    segments, one group per ligand / ion instance and one per coordinated-water role."""
    s = job.structure
    sel = np.zeros(s.n_residues, dtype=bool)
    sel[residues] = True
    het_label = {}
    for h in job.hetero:
        for m in h.get("members", [h["res"]]):
            het_label[int(m)] = h["label"]
    dom_of = {}
    for d in job.domains:
        rr = np.arange(d["r0"], d["r1"] + 1)
        rr = rr[np.isin(s.res_kind[rr], POLYMER)]
        if len(rr) and sel[rr].mean() >= 0.5:
            name = f"{d['source']} {d['name']} · {d['chain']}"
            for r in rr[sel[rr]]:
                dom_of.setdefault(int(r), name)
    groups: dict[str, dict] = {}
    for r in map(int, residues):
        k = int(s.res_kind[r])
        if k in POLYMER:
            lab, kind = (dom_of[r], "domain") if r in dom_of else (f"chain {s.res_chain[r]}", "chain")
        else:
            lab = het_label.get(r, f"{s.res_name[r]} {s.res_chain[r]}{s.res_seq[r]}")
            kind = KIND_NAMES[k]
        groups.setdefault(lab, {"label": lab, "kind": kind, "residues": []})["residues"].append(r)
    for gr in groups.values():          # waters: count the selected ones, not the structure-wide group size
        if gr["kind"] == "water":
            gr["label"] = re.sub(r"\s*\(\d+\)$", "", gr["label"]) + f" ×{len(gr['residues'])}"
    order = {"domain": 0, "chain": 1, "ligand": 2, "ion": 3, "water": 4}
    return sorted(groups.values(), key=lambda gr: (order.get(gr["kind"], 5), gr["residues"][0]))


def _hull_stats(P):
    try:
        h = g.convex_hull(P)
    except g.GeometryError:
        return None, {}
    return h, {"hull_volume": h.volume, "hull_area": h.area, "sphericity": h.sphericity}


def participating_tetrahedra(job, residues: np.ndarray, max_edge: float) -> dict:
    """The Delaunay tetrahedra of the job's point set that the selection takes part in (any vertex in a
    selected residue, every edge ≤ max_edge): the same set the geometry view highlights."""
    s, ps, t = job.structure, job.points, job.tess
    sel = np.zeros(s.n_residues, dtype=bool)
    sel[residues] = True
    pres = ps.residue[t.simplices]
    star = np.flatnonzero(sel[pres].any(axis=1) & (t.max_edge_per_tet() <= max_edge))
    if not len(star):
        return {"n_tets": 0, "max_edge": max_edge}
    tets = t.simplices[star]
    inside = sel[pres[star]].all(axis=1)
    e = np.unique(np.sort(np.concatenate([tets[:, [i, j]] for i in range(4) for j in range(i + 1, 4)]), axis=1), axis=0)
    L = np.linalg.norm(t.points[e[:, 0]] - t.points[e[:, 1]], axis=1)
    res = np.unique(pres[star])
    nb = res[~sel[res]]
    return {"n_tets": int(len(star)), "n_inside": int(inside.sum()), "n_nodes": int(len(np.unique(tets))),
            "n_edges": int(len(e)), "volume": float(t.volumes[star].sum()), "mean_edge": float(L.mean()),
            "max_edge": max_edge, "n_neighbour_residues": int(len(nb)),
            "neighbours": [s.residue_label(int(r)) for r in nb[:80]]}


def selection_geometry(job, residues, alpha: float | None = None, max_edge: float | None = None) -> tuple[dict, dict]:
    s = job.structure
    residues = np.unique(np.asarray(residues, dtype=np.int64))
    groups = selection_groups(job, residues)
    sel = np.zeros(s.n_residues, dtype=bool)
    sel[residues] = True
    atoms = np.flatnonzero(sel[s.atom_res])
    level = "heavy"
    if len(atoms) > MAX_SELECTION_ATOMS:
        poly = np.isin(s.res_kind[s.atom_res[atoms]], POLYMER)
        rep = np.isin(s.atom_name[atoms].astype(str), ("CA", "P"))
        atoms = atoms[~poly | rep]
        level = "ca"
    if alpha is None:
        alpha = 4.0 if level == "heavy" else 8.0
    P = s.coords[atoms]
    res_of = s.atom_res[atoms]
    el = _elements(s, atoms)
    gid_of_res = {r: gi for gi, gr in enumerate(groups) for r in gr["residues"]}
    gid = np.fromiter((gid_of_res[int(r)] for r in res_of), dtype=np.int32, count=len(res_of))
    not_water = s.res_kind[res_of] != KIND_WATER
    radii = radii_for(list(el))
    sasa_sel = atom_sasa(P, radii)                       # every atom with only the selection present
    sasa_ctx = job.sasa[atoms]                           # with the whole structure present

    summary = {
        "n_residues": int(len(residues)), "n_atoms": int(len(atoms)), "level": level, "alpha": float(alpha),
        "chains": sorted(set(s.res_chain[residues].tolist())), "n_groups": len(groups),
        "sasa_selection": float(sasa_sel[not_water].sum()), "sasa_in_structure": float(sasa_ctx[not_water].sum()),
        "shape": g.gyration_descriptors(P) if len(P) >= 2 else None,
    }
    summary["sasa_buried_by_rest"] = summary["sasa_selection"] - summary["sasa_in_structure"]
    summary["participating"] = participating_tetrahedra(job, residues, float(max_edge or job.defaults()["delaunay"]["max_edge"]))
    # the viewer works in coordinates centred on the structure (like every other layer)
    bufs = {"pts": P - job.center, "pt_group": gid.astype(np.uint16), "pt_res": res_of.astype(np.uint32)}
    hull, hs = _hull_stats(P)
    summary.update(hs)
    if hull is not None:
        bufs["hull_faces"] = hull.faces
    tess = None
    try:
        tess = g.Tess.build(P)
    except g.GeometryError as exc:
        summary["note"] = f"no tessellation: {exc}"
    if tess is not None:
        a = g.alpha_shape(tess, float(alpha))
        e = tess.edges
        L = tess.edge_lengths(e)
        summary.update({"n_tets": tess.n_tets, "n_edges": int(len(e)), "mean_edge": float(L.mean()),
                        "alpha_volume": a.volume, "alpha_area": a.area, "alpha_components": a.n_components,
                        "alpha_cavities": a.n_cavities})
        if hull is not None and hull.volume:
            summary["alpha_over_hull"] = a.volume / hull.volume
        bufs["alpha_faces"] = a.faces
        cut = 5.0 if level == "heavy" else 10.0
        cross = (gid[e[:, 0]] != gid[e[:, 1]]) & (L <= cut)
        del_pairs = np.sort(gid[e[cross]], axis=1)
    else:
        del_pairs = np.zeros((0, 2), dtype=np.int32)

    # ---------------------------------------------------------------- per group
    idx_of = [np.flatnonzero(gid == k) for k in range(len(groups))]
    trees = [cKDTree(P[ix]) if len(ix) else None for ix in idx_of]
    iso_area = []
    out_groups, gh_faces, gh_gid = [], [], []
    for k, gr in enumerate(groups):
        ix = idx_of[k]
        Pk = P[ix]
        iso = atom_sasa(Pk, radii[ix]) if len(ix) else np.zeros(0)
        nw = not_water[ix]
        iso_area.append(iso)
        row = {"label": gr["label"], "kind": gr["kind"], "residues": gr["residues"], "n_residues": len(gr["residues"]),
               "n_atoms": int(len(ix)), "chains": sorted(set(s.res_chain[gr["residues"]].tolist())),
               "centroid": Pk.mean(axis=0).tolist() if len(ix) else None,
               "rg": g.gyration_descriptors(Pk)["rg"] if len(ix) >= 2 else 0.0,
               "sasa_alone": float(iso[nw].sum()), "sasa_in_selection": float(sasa_sel[ix][nw].sum()),
               "sasa_in_structure": float(sasa_ctx[ix][nw].sum())}
        if row["sasa_alone"] > 0:
            row["buried_in_selection"] = 1 - row["sasa_in_selection"] / row["sasa_alone"]
            row["buried_in_structure"] = 1 - row["sasa_in_structure"] / row["sasa_alone"]
        if len(ix) >= 4:
            h, hs = _hull_stats(Pk)
            row.update(hs)
            if h is not None:
                gh_faces.append(ix[h.faces])
                gh_gid.append(np.full(len(h.faces), k, dtype=np.uint16))
        out_groups.append(row)
    if gh_faces:
        bufs["group_hull_faces"] = np.concatenate(gh_faces)
        bufs["group_hull_group"] = np.concatenate(gh_gid)

    # ---------------------------------------------------------------- per pair
    pairs, lines = [], []
    polar_el = np.isin(el, ("N", "O"))
    coord_el = np.isin(el, ("N", "O", "S"))
    n_pair = min(len(groups), MAX_PAIR_GROUPS)
    for a in range(n_pair):
        for b in range(a + 1, n_pair):
            ia, ib = idx_of[a], idx_of[b]
            if not len(ia) or not len(ib):
                continue
            d, j = trees[b].query(P[ia], k=1)
            m = int(np.argmin(d))
            row = {"a": a, "b": b, "min_dist": float(d[m]),
                   "min_pair": [s.residue_label(int(res_of[ia[m]])) + f" {s.atom_name[atoms[ia[m]]]}",
                                s.residue_label(int(res_of[ib[j[m]]])) + f" {s.atom_name[atoms[ib[j[m]]]]}"],
                   "centroid_dist": float(np.linalg.norm(P[ia].mean(0) - P[ib].mean(0))),
                   "delaunay_contacts": int(((del_pairs[:, 0] == a) & (del_pairs[:, 1] == b)).sum())}
            if row["min_dist"] <= NEAR:
                hits = trees[a].query_ball_tree(trees[b], CONTACT)
                cp = np.array([(ia[x], ib[y]) for x, ys in enumerate(hits) for y in ys], dtype=np.int64).reshape(-1, 2)
                dist = np.linalg.norm(P[cp[:, 0]] - P[cp[:, 1]], axis=1) if len(cp) else np.zeros(0)
                polar = polar_el[cp[:, 0]] & polar_el[cp[:, 1]] & (dist <= POLAR) if len(cp) else np.zeros(0, bool)
                ra, rb = np.unique(res_of[cp[:, 0]]), np.unique(res_of[cp[:, 1]])
                row.update({"n_contacts": int(len(cp)), "n_polar": int(polar.sum()),
                            "residues_a": ra.tolist(), "residues_b": rb.tolist(),
                            "labels_a": [s.residue_label(int(r)) for r in ra[:40]],
                            "labels_b": [s.residue_label(int(r)) for r in rb[:40]]})
                metal = []
                for ion_side, other in ((a, b), (b, a)):
                    if groups[ion_side]["kind"] != "ion":
                        continue
                    for x in idx_of[ion_side]:
                        near = trees[other].query_ball_point(P[x], METAL)
                        for y in near:
                            q = idx_of[other][y]
                            if coord_el[q]:
                                metal.append({"atom": s.residue_label(int(res_of[q])) + f" {s.atom_name[atoms[q]]}",
                                              "dist": float(np.linalg.norm(P[x] - P[q]))})
                if metal:
                    row["metal"] = sorted(metal, key=lambda z: z["dist"])
                if len(ia) + len(ib) <= MAX_BSA_ATOMS:
                    union = np.concatenate([ia, ib])
                    u = atom_sasa(P[union], radii[union])
                    ua, ub = u[:len(ia)][not_water[ia]].sum(), u[len(ia):][not_water[ib]].sum()
                    sa, sb = iso_area[a][not_water[ia]].sum(), iso_area[b][not_water[ib]].sum()
                    row["bsa"] = float(sa + sb - ua - ub)
                    row["interface_area"] = row["bsa"] / 2
                    row["buried_a"] = float(1 - ua / sa) if sa else None
                    row["buried_b"] = float(1 - ub / sb) if sb else None
                if len(lines) < MAX_CONTACT_LINES and len(cp):
                    lines.append(cp[: MAX_CONTACT_LINES - len(lines)])
            pairs.append(row)
    if lines:
        bufs["contacts"] = np.concatenate(lines)
    stats = {"summary": summary, "groups": out_groups, "pairs": pairs,
             "pairs_truncated": len(groups) > MAX_PAIR_GROUPS}
    return stats, bufs
