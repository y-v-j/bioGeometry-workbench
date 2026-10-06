"""Sequence-domain annotations for the sequence viewer's domain track.

* PDB entries: Pfam and CATH domain ranges per chain from the PDBe API (SIFTS mappings), in author
  residue numbering.
* AlphaFold models: UniProt features (Domain, Repeat, Zinc finger, Transmembrane, Coiled coil, Motif,
  Region) from the UniProt REST API; AlphaFold residue numbers are UniProt positions.

Responses are cached under data/annotations/. Network failures return an empty list plus a note, so
the workbench works offline without the domain track.
"""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import numpy as np

from .structure import DATA_DIR, Structure

ANNOT_DIR = DATA_DIR / "annotations"
PDBE = "https://www.ebi.ac.uk/pdbe/api/mappings/{db}/{pdb}"
UNIPROT = "https://rest.uniprot.org/uniprotkb/{acc}.json?fields=ft_domain,ft_region,ft_repeat,ft_zn_fing,ft_transmem,ft_coiled,ft_motif"
UNIPROT_TYPES = ("Domain", "Repeat", "Zinc finger", "Transmembrane", "Coiled coil", "Motif", "Region")


def _cached_json(key: str, url: str, timeout: float = 20):
    path = ANNOT_DIR / f"{key}.json"
    if path.exists():
        return json.loads(path.read_text())
    req = urllib.request.Request(url, headers={"User-Agent": "bioGeometry/0.1", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        data = json.loads(resp.read())
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    return data


UNIPROT_EC = "https://rest.uniprot.org/uniprotkb/{acc}.json?fields=ec"


def uniprot_ec(acc: str) -> str:
    """EC number(s) of a UniProt entry ('4.2.1.1', or several joined by ';'); '' if none or offline."""
    if not acc:
        return ""
    try:
        data = _cached_json(f"uniprot_ec_{acc}", UNIPROT_EC.format(acc=acc))
    except Exception:  # noqa: BLE001
        return ""
    out = []
    def walk(x):
        if isinstance(x, dict):
            if "ecNumbers" in x:
                out.extend(e.get("value", "") for e in x["ecNumbers"] if isinstance(e, dict))
            for v in x.values():
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
    walk(data.get("proteinDescription", {}))
    return ";".join(dict.fromkeys(e for e in out if e))


def chain_ec(s: Structure) -> dict[str, str]:
    """chain -> EC number(s): from the mmCIF entity (PDB entries), else from the chain's UniProt entry."""
    out = {}
    for ch, info in s.chain_info.items():
        ec = info.get("ec", "")
        if not ec and (info.get("uniprot") or (s.is_alphafold and s.meta.get("uniprot"))):
            ec = uniprot_ec(info.get("uniprot") or s.meta.get("uniprot", ""))
        out[ch] = ec
    if s.is_alphafold and not any(out.values()) and s.meta.get("uniprot"):
        ec = uniprot_ec(s.meta["uniprot"])
        out = {ch: ec for ch in s.chains}
    return out


def _pdbe_domains(pdb_id: str) -> list[dict]:
    pid = pdb_id.lower()
    out = []
    for db, key in (("pfam", "Pfam"), ("cath", "CATH")):
        try:
            data = _cached_json(f"pdbe_{db}_{pid}", PDBE.format(db=db, pdb=pid))
        except Exception:  # noqa: BLE001 - 404 means "no domains in this database"
            continue
        entries = (data.get(pid) or {}).get(key) or {}
        for acc, rec in entries.items():
            name = rec.get("identifier") or rec.get("description") or acc
            if db == "cath":
                name = rec.get("topology") or rec.get("homology") or name
            for m in rec.get("mappings", []):
                out.append({"source": key, "id": m.get("domain", acc) if db == "cath" else acc, "name": name,
                            "chain": m.get("chain_id"), "start": m["start"].get("author_residue_number"),
                            "end": m["end"].get("author_residue_number"),
                            "detail": rec.get("description") or rec.get("class", "")})
    return out


def _uniprot_domains(acc: str, chains: list[str]) -> list[dict]:
    try:
        data = _cached_json(f"uniprot_{acc}", UNIPROT.format(acc=acc))
    except Exception:  # noqa: BLE001
        return []
    out = []
    for f in data.get("features", []):
        if f.get("type") not in UNIPROT_TYPES:
            continue
        loc = f.get("location", {})
        a, b = loc.get("start", {}).get("value"), loc.get("end", {}).get("value")
        if a is None or b is None:
            continue
        for ch in chains:
            out.append({"source": "UniProt", "id": f.get("type"), "name": f.get("description") or f.get("type"),
                        "chain": ch, "start": int(a), "end": int(b), "detail": f.get("type")})
    return out


def domains(s: Structure) -> tuple[list[dict], str]:
    """Domain ranges mapped onto residue indices: [{source, id, name, chain, start, end, r0, r1, detail}].

    Returns (domains, note). r0/r1 are inclusive residue-table indices of the first/last residue in range.
    """
    note = ""
    try:
        if s.is_alphafold:
            raw = _uniprot_domains(s.meta.get("uniprot", ""), s.chains)
        elif s.source.startswith("rcsb:"):
            raw = _pdbe_domains(s.name)
        else:
            # local files: use UniProt references from the file if any
            raw = []
            for ch, info in s.chain_info.items():
                if info.get("uniprot"):
                    raw += _uniprot_domains(info["uniprot"], [ch])
            if not raw:
                note = "no domain source for local files without UniProt references"
    except Exception as exc:  # noqa: BLE001
        return [], f"domain lookup failed ({exc})"
    out = []
    seq = s.res_seq
    for d in raw:
        if d.get("chain") is None or d.get("start") is None or d.get("end") is None:
            continue
        in_chain = np.flatnonzero((s.res_chain == d["chain"]) & np.isin(s.res_kind, (0, 1)))
        if not len(in_chain):
            continue
        hit = in_chain[(seq[in_chain] >= d["start"]) & (seq[in_chain] <= d["end"])]
        if not len(hit):
            continue
        out.append({**d, "r0": int(hit.min()), "r1": int(hit.max()), "n": int(len(hit))})
    if not out and not note:
        note = "no domain annotations found"
    return out, note
