"""Fetch, parse and select macromolecular structures (PDB, mmCIF, AlphaFold DB).

Sources:
  * ``fetch_structure("1AON")``          RCSB entry (mmCIF, cached on disk; biological assemblies)
  * ``fetch_alphafold("P69905")``        AlphaFold DB model for a UniProt accession (pLDDT in the B column)
  * ``load_structure_file(path)``        a local PDB/mmCIF file (.gz accepted)
  * ``load_structure_bytes(data, name)`` an upload from the browser

Everything returns a :class:`Structure`: flat numpy arrays, one row per atom (first model, no
hydrogens, first altloc) plus a residue table. Polymer residues, ligands, ions and *coordinated*
waters are all kept and told apart by ``res_kind``; geometry point sets choose from them.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import re
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from . import biochem, resources

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
STRUCT_DIR = DATA_DIR / "structures"
AF_DIR = DATA_DIR / "alphafold"

RCSB_URL = "https://files.rcsb.org/download/{id}.{ext}.gz"
RCSB_ASSEMBLY_URL = "https://files.rcsb.org/download/{id}-assembly{n}.cif.gz"
AF_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
STRUCTURE_SUFFIXES = (".cif", ".mmcif", ".pdb", ".ent", ".cif.gz", ".mmcif.gz", ".pdb.gz", ".ent.gz")
USER_AGENT = {"User-Agent": "bioGeometry/0.1 (structure workbench)"}

# Residue kinds, shared with the workbench (keep in sync with workbench.js KIND).
KIND_PROTEIN, KIND_NUCLEIC, KIND_LIGAND, KIND_ION, KIND_WATER = 0, 1, 2, 3, 4
KIND_NAMES = ("protein", "nucleic", "ligand", "ion", "water")
POLYMER_KINDS = (KIND_PROTEIN, KIND_NUCLEIC)

REP_ATOMS = {KIND_PROTEIN: ("CA",), KIND_NUCLEIC: ("P", "C4'", "C4*")}
BACKBONE_ATOMS = {
    "N", "CA", "C", "O", "OXT",
    "P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'",
}
LEVELS = ("ca", "centroid", "backbone", "heavy")

ELEMENTS = ["C", "N", "O", "S", "P", "Se", "Mg", "Zn", "Fe", "Ca", "Mn", "Cl", "Na", "K", "Cu", "Co", "Ni", "Br", "I", "F",
            # appended (codes above stay stable): less common metals, so they are still recognised as metals
            "Cd", "Hg", "Sr", "Ba", "Li", "Rb", "Cs", "Mo", "W", "V", "Cr", "Pt", "Au", "Ag", "Al", "Ga", "Pb", "Tl",
            "Yb", "Gd", "Tb", "Sm", "Eu", "Lu", "B", "As", "Si"]
ELEMENT_CODE = {e.upper(): i for i, e in enumerate(ELEMENTS)}  # unknown -> len(ELEMENTS)
METALS = {"MG", "ZN", "FE", "CA", "MN", "NA", "K", "CU", "CO", "NI", "CD", "HG", "SR", "BA", "LI", "RB", "CS",
          "AL", "GA", "PT", "AU", "AG", "W", "MO", "V", "CR", "YB", "LU", "TB", "GD", "SM", "EU", "PB", "TL"}
HALIDES = {"CL", "BR", "I", "F", "IOD"}
WATER_NAMES = {"HOH", "WAT", "DOD", "H2O"}
# Only waters taking part in ion / ligand binding are kept ("coordinated"):
ION_SHELL_CUTOFF = 3.0      # water O in the first coordination shell of an ion (Mg-O 2.1, Na-O 2.4, K-O 2.8)
LIGAND_WATER_CUTOFF = 3.5   # water O hydrogen-bonded to (or bridging) a ligand heavy atom


class StructureError(ValueError):
    """User-facing error: bad id, missing chain, unreadable file."""


@dataclass
class Structure:
    """Atoms of one model as parallel arrays, plus a residue table and entry metadata."""

    name: str
    title: str
    source: str                  # "rcsb:mmcif", "alphafold:AF-P69905-F1", "file:/path", "upload:name"
    coords: np.ndarray           # (N, 3) float64
    element: np.ndarray          # (N,) uint8 code into ELEMENTS
    atom_name: np.ndarray        # (N,) object
    bfactor: np.ndarray          # (N,) float32 (pLDDT for AlphaFold models)
    atom_res: np.ndarray         # (N,) int32 index into the residue table
    res_chain: np.ndarray        # (R,) object  auth chain id
    res_name: np.ndarray         # (R,) object
    res_seq: np.ndarray          # (R,) int32
    res_icode: np.ndarray        # (R,) object
    res_kind: np.ndarray         # (R,) uint8 KIND_*
    chains: list[str] = field(default_factory=list)
    assembly: str | None = None
    meta: dict = field(default_factory=dict)        # method, resolution, organism, uniprot, af metrics...
    chain_info: dict = field(default_factory=dict)  # chain -> {"description", "entity", "uniprot"}

    @property
    def n_atoms(self) -> int:
        return len(self.coords)

    @property
    def n_residues(self) -> int:
        return len(self.res_name)

    @property
    def is_alphafold(self) -> bool:
        return self.source.startswith("alphafold:")

    def residue_label(self, r: int) -> str:
        ic = self.res_icode[r].strip()
        return f"{self.res_chain[r]}:{self.res_name[r]}{self.res_seq[r]}{ic}"

    def chain_index(self) -> np.ndarray:
        lookup = {c: i for i, c in enumerate(self.chains)}
        return np.fromiter((lookup[c] for c in self.res_chain), dtype=np.int32, count=self.n_residues)

    def polymer_mask(self) -> np.ndarray:
        return np.isin(self.res_kind, POLYMER_KINDS)


# --------------------------------------------------------------------------- download helpers

def _get(url: str, timeout: float = 60) -> bytes:
    req = urllib.request.Request(url, headers=USER_AGENT)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def _atomic_write(path: Path, blob: bytes):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.write_bytes(blob)
    tmp.replace(path)


def normalize_pdb_id(pdb_id: str) -> str:
    pid = (pdb_id or "").strip().upper()
    if not re.fullmatch(r"[0-9][A-Z0-9]{3}|PDB_[0-9A-Z]{8}", pid):
        raise StructureError(f"'{pdb_id}' is not a valid PDB id (expected e.g. 1AON)")
    return pid


def normalize_uniprot(acc: str) -> str:
    """Accept 'P69905', 'AF-P69905-F1', 'AF-P69905-F1-model_v6' and return the accession."""
    s = (acc or "").strip().upper()
    m = re.fullmatch(r"(?:AF-)?([A-Z0-9]{6,10})(?:-F\d+)?(?:-MODEL_V\d+)?", s)
    if not m or not re.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}", m.group(1)):
        raise StructureError(f"'{acc}' is not a UniProt accession (expected e.g. P69905 or AF-P69905-F1)")
    return m.group(1)


def fetch_structure(pdb_id: str, assembly: str | None = None, timeout: float = 60, **kw) -> Structure:
    """Download (or reuse the cached copy of) an RCSB entry as mmCIF and parse it.

    ``assembly`` = "1", "2", ... downloads the biological assembly file instead of building it.
    """
    pid = normalize_pdb_id(pdb_id)
    if assembly and str(assembly).isdigit():
        path = STRUCT_DIR / f"{pid}-assembly{assembly}.cif.gz"
        url = RCSB_ASSEMBLY_URL.format(id=pid, n=assembly)
    else:
        path = STRUCT_DIR / f"{pid}.cif.gz"
        url = RCSB_URL.format(id=pid, ext="cif")
    if not path.exists() or path.stat().st_size == 0:
        try:
            _atomic_write(path, _get(url, timeout))
        except Exception as exc:
            raise StructureError(f"could not download {pid} from RCSB ({exc})") from exc
    st = _parse_file(path)
    return _to_structure(st, pid, "rcsb:mmcif", None, **kw)


def alphafold_metadata(acc: str, timeout: float = 30) -> dict:
    """The AlphaFold DB prediction record (cached), for the first fragment of the accession."""
    acc = normalize_uniprot(acc)
    cache = AF_DIR / f"{acc}.json"
    if cache.exists():
        recs = json.loads(cache.read_text())
    else:
        try:
            recs = json.loads(_get(AF_API.format(acc=acc), timeout))
        except Exception as exc:
            raise StructureError(f"no AlphaFold DB prediction for {acc} ({exc})") from exc
        _atomic_write(cache, json.dumps(recs).encode())
    if not recs:
        raise StructureError(f"no AlphaFold DB prediction for {acc}")
    return recs[0]


def fetch_alphafold(acc: str, timeout: float = 60, **kw) -> Structure:
    """Download the AlphaFold DB model for a UniProt accession (pLDDT is in the B-factor column)."""
    rec = alphafold_metadata(acc, timeout)
    entry = rec["entryId"]
    version = rec.get("latestVersion", 4)
    path = AF_DIR / f"{entry}-model_v{version}.cif"
    if not path.exists() or path.stat().st_size == 0:
        try:
            _atomic_write(path, _get(rec["cifUrl"], timeout))
        except Exception as exc:
            raise StructureError(f"could not download {entry} from AlphaFold DB ({exc})") from exc
    st = _parse_file(path)
    s = _to_structure(st, entry, f"alphafold:{entry}", None, **kw)
    s.title = f"{rec.get('uniprotDescription', '')} ({rec.get('gene', '') or rec['uniprotAccession']}) · " \
              f"{rec.get('organismScientificName', '')}".strip(" ·")
    s.meta.update(method="AlphaFold prediction", uniprot=rec["uniprotAccession"], af_entry=entry,
                  af_version=version, organism=rec.get("organismScientificName", ""),
                  gene=rec.get("gene", ""), mean_plddt=rec.get("globalMetricValue"),
                  pae_url=rec.get("paeDocUrl", ""), uniprot_start=rec.get("uniprotStart"),
                  uniprot_end=rec.get("uniprotEnd"))
    for c in s.chain_info.values():
        c["uniprot"] = rec["uniprotAccession"]
        c["description"] = rec.get("uniprotDescription", c.get("description", ""))
    return s


def alphafold_pae(acc_or_entry: str, timeout: float = 120) -> np.ndarray | None:
    """Predicted aligned error matrix (N×N, Å) for an AlphaFold model, or None if unavailable."""
    rec = alphafold_metadata(acc_or_entry, timeout)
    url = rec.get("paeDocUrl")
    if not url:
        return None
    path = AF_DIR / Path(url).name
    if not path.exists():
        try:
            _atomic_write(path, _get(url, timeout))
        except Exception:
            return None
    doc = json.loads(path.read_text())
    doc = doc[0] if isinstance(doc, list) else doc
    if "predicted_aligned_error" in doc:
        return np.asarray(doc["predicted_aligned_error"], dtype=np.float32)
    if "distance" in doc:  # very old format: flat lists
        n = int(max(doc["residue1"]))
        m = np.zeros((n, n), np.float32)
        m[np.asarray(doc["residue1"]) - 1, np.asarray(doc["residue2"]) - 1] = doc["distance"]
        return m
    return None


def load_structure_file(path: str | Path, assembly: str | None = None, **kw) -> Structure:
    p = Path(path).expanduser()
    if not p.is_file():
        raise StructureError(f"file not found: {p}")
    if not p.name.lower().endswith(STRUCTURE_SUFFIXES):
        raise StructureError(f"{p.name}: expected one of {', '.join(STRUCTURE_SUFFIXES)}")
    st = _parse_file(p)
    src = f"upload:{_display_name(p)}" if p.parent.name == "uploads" else f"file:{p}"
    return _to_structure(st, _display_name(p), src, assembly, **kw)


def store_upload(data: bytes, filename: str) -> Path:
    name = Path(filename).name
    if not name.lower().endswith(STRUCTURE_SUFFIXES):
        raise StructureError(f"{name}: expected one of {', '.join(STRUCTURE_SUFFIXES)}")
    digest = hashlib.sha1(data).hexdigest()[:12]
    path = DATA_DIR / "uploads" / f"{digest}_{name}"
    if not path.exists():
        _atomic_write(path, data)
    return path


def load_structure_bytes(data: bytes, filename: str, **kw) -> tuple[Structure, Path]:
    path = store_upload(data, filename)
    return load_structure_file(path, **kw), path


def list_local_structures(directory: str | Path) -> list[str]:
    d = Path(directory).expanduser()
    if not d.is_dir():
        raise StructureError(f"not a directory: {d}")
    return sorted(str(p) for p in d.iterdir() if p.is_file() and p.name.lower().endswith(STRUCTURE_SUFFIXES))


def _display_name(path: Path) -> str:
    m = re.fullmatch(r"[0-9a-f]{12}_(.+)", path.name)
    return _stem(m.group(1) if m and path.parent.name == "uploads" else path.name)


def _stem(name: str) -> str:
    n = name
    for suf in (".gz", ".cif", ".mmcif", ".pdb", ".ent"):
        if n.lower().endswith(suf):
            n = n[: -len(suf)]
    return n


_NAMES: dict[int, tuple[dict, dict]] = {}   # id(gemmi.Structure) -> names, consumed by _to_structure


def _cif_names(path: Path) -> tuple[dict, dict]:
    """(entity id -> description, ligand comp id -> name) from an mmCIF file; empty for PDB format."""
    try:
        block = gemmi.cif.read(str(path)).sole_block()
    except Exception:  # noqa: BLE001 - PDB format or unreadable
        return {}, {}
    ent = {row[0]: gemmi.cif.as_string(row[1]) for row in block.find("_entity.", ["id", "pdbx_description"])}
    het = {row[0]: gemmi.cif.as_string(row[1]) for row in block.find("_pdbx_entity_nonpoly.", ["comp_id", "name"])}
    for row in block.find("_entity.", ["id", "?pdbx_ec"]):            # EC number(s) of enzymes, e.g. 4.2.1.1
        if row.has(1) and gemmi.cif.as_string(row[1]).strip() not in ("", "?", "."):
            het[f"__ec_{row[0]}"] = gemmi.cif.as_string(row[1]).strip()
    for row in block.find("_struct_ref.", ["entity_id", "db_name", "pdbx_db_accession"]):
        if gemmi.cif.as_string(row[1]).upper() == "UNP":
            het[f"__unp_{row[0]}"] = gemmi.cif.as_string(row[2])
    return ent, het


def _parse_file(path: Path) -> gemmi.Structure:
    try:
        st = gemmi.read_structure(str(path), merge_chain_parts=True)
    except Exception as exc:
        try:
            raw = path.read_bytes()
            if raw[:2] == b"\x1f\x8b":
                raw = gzip.decompress(raw)
            text = raw.decode("utf-8", "replace")
            if text.lstrip().startswith("data_"):
                st = gemmi.make_structure_from_block(gemmi.cif.read_string(text).sole_block())
            else:
                st = gemmi.read_pdb_string(text)
        except Exception:
            raise StructureError(f"cannot parse {path.name}: {exc}") from exc
    if len(st) == 0:
        raise StructureError(f"{path.name}: no models found")
    st.setup_entities()
    _NAMES[id(st)] = _cif_names(path)
    return st


# --------------------------------------------------------------------------- residue kinds

def _residue_kind(res: gemmi.Residue) -> int:
    name = res.name.upper()
    if name in WATER_NAMES or res.is_water():
        return KIND_WATER
    info = gemmi.find_tabulated_residue(res.name)
    if res.entity_type == gemmi.EntityType.Polymer or res.het_flag == "A":
        if info is not None and info.is_nucleic_acid():
            return KIND_NUCLEIC
        if info is not None and info.is_amino_acid():
            return KIND_PROTEIN
    if info is not None and info.is_amino_acid() and res.find_atom("CA", "*") is not None \
            and res.entity_type != gemmi.EntityType.NonPolymer:
        return KIND_PROTEIN
    heavy = [a for a in res if not a.is_hydrogen()]
    if len(heavy) == 1 and (heavy[0].element.name.upper() in METALS | HALIDES or name in METALS | HALIDES):
        return KIND_ION
    return KIND_LIGAND


def _metadata(st: gemmi.Structure) -> dict:
    meta = {}
    try:
        meta["resolution"] = float(st.resolution) if st.resolution else None
    except Exception:  # noqa: BLE001
        meta["resolution"] = None
    for key, tag in (("method", "_exptl.method"), ("keywords", "_struct_keywords.pdbx_keywords"),
                     ("deposited", "_pdbx_database_status.recvd_initial_deposition_date")):
        if tag in st.info:
            meta[key] = st.info[tag]
    return meta


def _chain_info(st: gemmi.Structure, model: gemmi.Model) -> dict:
    ent_names = _NAMES.get(id(st), ({}, {}))[0]
    out = {}
    for chain in model:
        # polymer subchain first; chains holding only ions / ligands / waters fall back to their subchains
        spans = [chain.get_polymer()] + list(chain.subchains())
        ent = next((e for e in (st.get_entity_of(sp) for sp in spans if len(sp)) if e is not None), None)
        desc, uniprot, ec = "", "", ""
        if ent is not None:
            ec = _NAMES.get(id(st), ({}, {}))[1].get(f"__ec_{ent.name}", "")
            desc = ent_names.get(ent.name, "") or getattr(ent, "full_name", "") or ""
            acc = list(getattr(ent, "sifts_unp_acc", []) or [])
            uniprot = acc[0].split("-")[0] if acc else _NAMES.get(id(st), ({}, {}))[1].get(f"__unp_{ent.name}", "")
        out[chain.name] = {"description": desc.strip(), "entity": ent.name if ent is not None else "", "uniprot": uniprot,
                           "ec": ec}
    return out


def _to_structure(st: gemmi.Structure, name: str, source: str, assembly: str | None,
                  waters: str = "coordinated", **_ignored) -> Structure:
    title = st.info["_struct.title"] if "_struct.title" in st.info else ""
    model = st[0]
    if assembly and not str(assembly).isdigit():
        assemblies = {a.name: a for a in st.assemblies}
        if assembly not in assemblies:
            raise StructureError(f"assembly '{assembly}' not in file (have: {', '.join(assemblies) or 'none'})")
        model = gemmi.make_assembly(assemblies[assembly], model, gemmi.HowToNameCopiedChain.Short)
    n_sites = model.count_atom_sites()
    resources.require(resources.parse_estimate_mb(n_sites), f"loading {name} ({n_sites:,} atoms)")
    chain_info = _chain_info(st, model)

    coords, elem, aname, bfac, ares = [], [], [], [], []
    rchain, rname, rseq, ricode, rkind = [], [], [], [], []
    unknown = len(ELEMENTS)
    for chain in model:
        for res in chain:
            kind = _residue_kind(res)
            if kind == KIND_WATER and waters == "none":
                continue
            ridx = len(rname)
            n_before = len(coords)
            seen: set[str] = set()
            for atom in res:
                if atom.is_hydrogen() or atom.name in seen:
                    continue
                seen.add(atom.name)
                p = atom.pos
                coords.append((p.x, p.y, p.z))
                elem.append(ELEMENT_CODE.get(atom.element.name.upper(), unknown))
                aname.append(atom.name)
                bfac.append(atom.b_iso)
                ares.append(ridx)
            if len(coords) == n_before:
                continue
            rchain.append(chain.name)
            rname.append(res.name)
            rseq.append(res.seqid.num)
            ricode.append(res.seqid.icode.strip())
            rkind.append(kind)
    if not coords:
        raise StructureError(f"{name}: no atoms found")
    s = Structure(
        name=name, title=title, source=source,
        coords=np.asarray(coords, dtype=np.float64),
        element=np.asarray(elem, dtype=np.uint8),
        atom_name=np.asarray(aname, dtype=object),
        bfactor=np.asarray(bfac, dtype=np.float32),
        atom_res=np.asarray(ares, dtype=np.int32),
        res_chain=np.asarray(rchain, dtype=object),
        res_name=np.asarray(rname, dtype=object),
        res_seq=np.asarray(rseq, dtype=np.int32),
        res_icode=np.asarray(ricode, dtype=object),
        res_kind=np.asarray(rkind, dtype=np.uint8),
        chains=list(dict.fromkeys(rchain)), assembly=assembly,
        meta=_metadata(st), chain_info=chain_info,
    )
    s.meta["het_names"] = _NAMES.pop(id(st), ({}, {}))[1]
    if waters == "coordinated":
        s = _keep_coordinated_waters(s)
    if not s.polymer_mask().any():
        raise StructureError(f"{name}: no protein or nucleic-acid chains found")
    return s


def coordinated_water_roles(s: Structure) -> dict[int, str]:
    """{water residue index: role} for waters that take part in ion / ligand binding.

    role is "ion shell" (O within 3.0 Å of an ion: first coordination sphere) or "ligand"
    (O within 3.5 Å of a ligand heavy atom: hydrogen bond / water bridge to the ligand).
    Waters touching only protein or nucleic acid, and bulk solvent, are not coordinated.
    """
    water_res = s.res_kind == KIND_WATER
    if not water_res.any():
        return {}
    w_atoms = np.flatnonzero(water_res[s.atom_res])
    roles: dict[int, str] = {}
    for kind, cutoff, role in ((KIND_LIGAND, LIGAND_WATER_CUTOFF, "ligand"), (KIND_ION, ION_SHELL_CUTOFF, "ion shell")):
        anchors = np.flatnonzero(s.res_kind[s.atom_res] == kind)
        if not len(anchors):
            continue
        d, _ = cKDTree(s.coords[anchors]).query(s.coords[w_atoms], k=1, distance_upper_bound=cutoff)
        for r in np.unique(s.atom_res[w_atoms[np.isfinite(d)]]):
            roles[int(r)] = role          # ion shell overrides ligand (assigned second)
    return roles


def _keep_coordinated_waters(s: Structure) -> Structure:
    """Drop every water except those coordinating an ion or bound to a ligand."""
    water_res = s.res_kind == KIND_WATER
    if not water_res.any():
        return s
    keep_res = ~water_res
    keep_res[list(coordinated_water_roles(s))] = True
    return subset_residues(s, keep_res)


# ----------------------------------------------------------------------- export

def chain_tag(s: Structure) -> str:
    """'ABCD' for single-letter chain ids, 'A1/B1' otherwise, '56ch' for very many chains."""
    if len(s.chains) > 26:
        return f"{len(s.chains)}ch"
    return ("" if all(len(c) == 1 for c in s.chains) else "/").join(s.chains)


def to_gemmi(s: Structure) -> gemmi.Structure:
    """The loaded model (chains, assembly and coordinated waters as shown) as a gemmi structure."""
    st = gemmi.Structure()
    st.name = s.name
    if s.title:
        st.info["_struct.title"] = s.title
    model = gemmi.Model(1)
    starts = np.searchsorted(s.atom_res, np.arange(s.n_residues + 1))     # atoms are stored per residue
    chain = None
    for r in range(s.n_residues):
        if chain is None or chain.name != s.res_chain[r]:
            if chain is not None:
                model.add_chain(chain)
            chain = gemmi.Chain(str(s.res_chain[r]))
        res = gemmi.Residue()
        res.name = str(s.res_name[r])
        res.seqid = gemmi.SeqId(int(s.res_seq[r]), (s.res_icode[r] or " ")[0])
        res.het_flag = "A" if s.res_kind[r] in POLYMER_KINDS else "H"
        for a in range(starts[r], starts[r + 1]):
            at = gemmi.Atom()
            at.name = str(s.atom_name[a])
            e = int(s.element[a])
            at.element = gemmi.Element(ELEMENTS[e] if e < len(ELEMENTS) else "X")
            x, y, z = s.coords[a]
            at.pos = gemmi.Position(float(x), float(y), float(z))
            at.b_iso = float(s.bfactor[a])
            at.occ = 1.0
            res.add_atom(at)
        chain.add_residue(res)
    if chain is not None:
        model.add_chain(chain)
    st.add_model(model)
    st.setup_entities()
    return st


def structure_text(s: Structure, fmt: str) -> str:
    """mmCIF or PDB text of the loaded model. PDB needs 1-character chain ids: longer ids are
    shortened when there are few enough chains, otherwise PDB is refused (use mmCIF)."""
    st = to_gemmi(s)
    if fmt == "cif":
        return st.make_mmcif_document().as_string()
    if fmt != "pdb":
        raise StructureError(f"unknown structure format '{fmt}' (cif or pdb)")
    if any(len(c) > 1 for c in s.chains):
        if len(s.chains) > 62:
            raise StructureError(f"{s.name} has {len(s.chains)} chains - too many for PDB format; download mmCIF")
        st.shorten_chain_names()
    return st.make_pdb_string()


# ----------------------------------------------------------------------- selection

def parse_chain_ids(text: str | None) -> list[str] | None:
    if not text or not text.strip():
        return None
    return [c for c in re.split(r"[\s,;]+", text.strip()) if c]


def select_chains(s: Structure, chains: list[str] | None) -> Structure:
    """Restrict to the given auth chain ids. Matching is case-sensitive (ribosomes use A and a);
    a chain given in the wrong case is accepted only when the match is unique."""
    if not chains:
        return s
    have = set(s.chains)
    wanted = [c for c in chains if c in have]
    for c in chains:
        if c in have:
            continue
        alts = [x for x in s.chains if x.upper() == c.upper()]
        if len(alts) == 1:
            wanted.append(alts[0])
        else:
            raise StructureError(f"chain '{c}' not in {s.name} (chains: {', '.join(s.chains[:60])}"
                                 f"{' ...' if len(s.chains) > 60 else ''})")
    return subset_residues(s, np.isin(s.res_chain, np.asarray(wanted, dtype=object)))


def subset_residues(s: Structure, keep_res: np.ndarray) -> Structure:
    new_index = np.full(s.n_residues, -1, dtype=np.int32)
    new_index[keep_res] = np.arange(int(keep_res.sum()), dtype=np.int32)
    keep_atom = keep_res[s.atom_res]
    kept_chains = set(s.res_chain[keep_res])
    chains = [c for c in s.chains if c in kept_chains]
    return Structure(
        name=s.name, title=s.title, source=s.source,
        coords=s.coords[keep_atom], element=s.element[keep_atom], atom_name=s.atom_name[keep_atom],
        bfactor=s.bfactor[keep_atom], atom_res=new_index[s.atom_res[keep_atom]],
        res_chain=s.res_chain[keep_res], res_name=s.res_name[keep_res], res_seq=s.res_seq[keep_res],
        res_icode=s.res_icode[keep_res], res_kind=s.res_kind[keep_res], chains=chains, assembly=s.assembly,
        meta=dict(s.meta), chain_info={c: v for c, v in s.chain_info.items() if c in kept_chains},
    )


# ----------------------------------------------------------------------- hetero groups

def hetero_groups(s: Structure) -> list[dict]:
    """Ligand / ion / water instances, numbered per name in order of appearance: NA (i), NA (ii)…

    Coordinated waters are summarised as one group per role (ion shell, ligand) rather than listed one by one.
    """
    out, counter = [], {}
    het = np.flatnonzero(np.isin(s.res_kind, (KIND_LIGAND, KIND_ION)))
    for r in het:
        nm = str(s.res_name[r])
        counter[nm] = counter.get(nm, 0) + 1
        out.append({"res": int(r), "name": nm, "label": f"{nm} ({biochem.roman(counter[nm])})",
                    "kind": KIND_NAMES[int(s.res_kind[r])], "chain": str(s.res_chain[r]),
                    "seq": int(s.res_seq[r]),
                    "description": biochem.hetero_description(nm) or s.meta.get("het_names", {}).get(nm, "").lower()})
    waters = np.flatnonzero(s.res_kind == KIND_WATER)
    if len(waters):
        roles = coordinated_water_roles(s)
        for role, desc in (("ion shell", "waters in an ion's first coordination shell (O ≤ 3.0 Å)"),
                           ("ligand", "waters hydrogen-bonded to / bridging a ligand (O ≤ 3.5 Å)")):
            members = np.array([w for w in waters if roles.get(int(w)) == role], dtype=int)
            if len(members):
                out.append({"res": int(members[0]), "members": members.tolist(), "name": "HOH",
                            "label": f"H₂O {role} ({len(members)})", "kind": "water",
                            "chain": ",".join(dict.fromkeys(s.res_chain[members])),
                            "seq": int(s.res_seq[members[0]]), "description": desc})
    return out


# ----------------------------------------------------------------------- point sets

@dataclass
class PointSet:
    level: str
    coords: np.ndarray        # (P, 3) float64
    residue: np.ndarray       # (P,) int32 residue index
    atom: np.ndarray          # (P,) int32 atom index (-1 for centroids)


def auto_level(s: Structure) -> str:
    n = int(s.polymer_mask()[s.atom_res].sum())
    if n <= 12_000:
        return "heavy"
    if n <= 60_000:
        return "backbone" if int(s.polymer_mask().sum()) <= 4_000 else "ca"
    return "ca"


def representative_points(s: Structure, level: str, include_hetero: bool = False) -> PointSet:
    """Reduce the structure to the point set used for tessellation.

    ca       one point per residue: CA (protein), P or C4' (nucleic), centroid (ligand)
    centroid one point per residue: side-chain centroid (residue centroid for Gly/nucleic/ligand)
    backbone backbone heavy atoms
    heavy    all heavy atoms
    The level applies to the polymer. With ``include_hetero`` every heavy atom of the ligands, ions and
    (coordinated) waters is added at any level, so "Cα + hetero" keeps the full shape of a ligand.
    """
    if level == "auto":
        level = auto_level(s)
    if level not in LEVELS:
        raise StructureError(f"unknown level '{level}' (choose from {', '.join(LEVELS)})")
    poly = _polymer_points(s, level)
    if not include_hetero:
        return poly
    het = np.flatnonzero(np.isin(s.res_kind, (KIND_LIGAND, KIND_ION, KIND_WATER))[s.atom_res])
    if not len(het):
        return poly
    return PointSet(level, np.concatenate([poly.coords, s.coords[het]]),
                    np.concatenate([poly.residue, s.atom_res[het]]).astype(np.int32),
                    np.concatenate([poly.atom, het]).astype(np.int32))


def _polymer_points(s: Structure, level: str) -> PointSet:
    use_res = np.isin(s.res_kind, POLYMER_KINDS)
    use_atom = use_res[s.atom_res]
    atom_idx = np.arange(s.n_atoms, dtype=np.int32)

    if level in ("heavy", "backbone"):
        sel = use_atom.copy()
        if level == "backbone":
            is_bb = np.fromiter((n in BACKBONE_ATOMS for n in s.atom_name), dtype=bool, count=s.n_atoms)
            sel &= is_bb
        return PointSet(level, s.coords[sel], s.atom_res[sel], atom_idx[sel])

    res_ids = np.flatnonzero(use_res)
    nres = s.n_residues
    if level == "centroid":
        side = np.fromiter((n not in BACKBONE_ATOMS for n in s.atom_name), dtype=bool, count=s.n_atoms)
        sc = np.bincount(s.atom_res, weights=(side & use_atom).astype(float), minlength=nres) > 0
        w = (np.where(sc[s.atom_res], side, True) & use_atom).astype(float)
    else:
        w = use_atom.astype(float)
    denom = np.maximum(np.bincount(s.atom_res, weights=w, minlength=nres), 1e-12)
    cen = np.stack([np.bincount(s.atom_res, weights=w * s.coords[:, k], minlength=nres) for k in range(3)],
                   axis=1) / denom[:, None]
    if level == "centroid":
        return PointSet(level, cen[res_ids], res_ids.astype(np.int32), np.full(len(res_ids), -1, np.int32))

    rep = np.full(nres, -1, dtype=np.int64)
    for kind, names in REP_ATOMS.items():
        for nm in reversed(names):
            hit = np.nonzero((s.atom_name == nm) & (s.res_kind[s.atom_res] == kind))[0]
            rep[s.atom_res[hit]] = hit
    coords = np.where((rep >= 0)[:, None], s.coords[np.maximum(rep, 0)], cen)
    return PointSet(level, coords[res_ids], res_ids.astype(np.int32), rep[res_ids].astype(np.int32))


def residue_points(s: Structure) -> PointSet:
    """One point per polymer residue (CA / P), for residue-level analyses at any level."""
    return representative_points(s, "ca")


def dedupe_points(ps: PointSet, tol: float = 1e-3) -> tuple[PointSet, int]:
    key = np.round(ps.coords / tol).astype(np.int64)
    _, first = np.unique(key, axis=0, return_index=True)
    first.sort()
    dropped = len(ps.coords) - len(first)
    if dropped == 0:
        return ps, 0
    return PointSet(ps.level, ps.coords[first], ps.residue[first], ps.atom[first]), dropped
