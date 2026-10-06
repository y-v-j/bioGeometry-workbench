"""Residue biochemistry: reference tables and classifiers used by the feature table.

Sources of the constants:
  * Kyte & Doolittle (1982) hydropathy
  * Zamyatnin (1972) residue volumes in solution (Å³) - the reference for packing ratios
  * Tien et al. (2013) theoretical maximum accessible surface area (Å²) - for RSA
  * Ramachandran regions: coarse rectangles around the populated basins (Lovell et al. 2003 contours
    are not reproduced); "outlier" only for non-Gly residues outside every basin.
"""
from __future__ import annotations

import numpy as np

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O", "MSE": "M", "ASX": "B",
    "GLX": "Z", "UNK": "X",
    # nucleotides
    "A": "A", "C": "C", "G": "G", "U": "U", "T": "T", "I": "I",
    "DA": "A", "DC": "C", "DG": "G", "DT": "T", "DU": "U", "DI": "I",
}

KYTE_DOOLITTLE = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9, "ALA": 1.8, "GLY": -0.4,
    "THR": -0.7, "SER": -0.8, "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5,
    "GLN": -3.5, "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5, "MSE": 1.9,
}

# formal side-chain charge at pH 7 (His treated as neutral; it is flagged as titratable instead)
CHARGE = {"ASP": -1, "GLU": -1, "LYS": 1, "ARG": 1}
TITRATABLE = {"ASP", "GLU", "LYS", "ARG", "HIS", "CYS", "TYR"}

VOLUME = {  # Zamyatnin 1972, Å³
    "ALA": 88.6, "ARG": 173.4, "ASN": 114.1, "ASP": 111.1, "CYS": 108.5, "GLN": 143.8, "GLU": 138.4,
    "GLY": 60.1, "HIS": 153.2, "ILE": 166.7, "LEU": 166.7, "LYS": 168.6, "MET": 162.9, "PHE": 189.9,
    "PRO": 112.7, "SER": 89.0, "THR": 116.1, "TRP": 227.8, "TYR": 193.6, "VAL": 140.0, "MSE": 162.9,
}

MAX_ASA = {  # Tien et al. 2013, theoretical, Å²
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0, "GLN": 225.0, "GLU": 223.0,
    "GLY": 104.0, "HIS": 224.0, "ILE": 197.0, "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0,
    "PRO": 159.0, "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0, "MSE": 224.0,
    # nucleotides (approximate maxima for an exposed nucleotide in a strand)
    "A": 400.0, "G": 420.0, "C": 370.0, "U": 370.0, "DA": 390.0, "DG": 410.0, "DC": 360.0, "DT": 380.0,
}

CLASSES = {  # biochemical class shown in the sequence viewer and table
    "ALA": "aliphatic", "VAL": "aliphatic", "LEU": "aliphatic", "ILE": "aliphatic", "MET": "aliphatic",
    "MSE": "aliphatic", "PHE": "aromatic", "TRP": "aromatic", "TYR": "aromatic",
    "SER": "polar", "THR": "polar", "ASN": "polar", "GLN": "polar",
    "ASP": "acidic", "GLU": "acidic", "LYS": "basic", "ARG": "basic", "HIS": "basic",
    "GLY": "special", "PRO": "special", "CYS": "special",
    "A": "purine", "G": "purine", "DA": "purine", "DG": "purine", "I": "purine", "DI": "purine",
    "C": "pyrimidine", "U": "pyrimidine", "T": "pyrimidine", "DC": "pyrimidine", "DT": "pyrimidine", "DU": "pyrimidine",
}
CLASS_ORDER = ["aliphatic", "aromatic", "polar", "acidic", "basic", "special", "purine", "pyrimidine",
               "ligand", "ion", "water", "other"]

HETERO_DESCRIPTIONS = {
    "HEM": "heme (protoporphyrin IX + Fe)", "HEC": "heme C", "FAD": "flavin adenine dinucleotide",
    "FMN": "flavin mononucleotide", "NAD": "NAD⁺", "NAP": "NADP⁺", "NDP": "NADPH", "ATP": "ATP",
    "ADP": "ADP", "ANP": "AMP-PNP (ATP analogue)", "GTP": "GTP", "GDP": "GDP", "SAM": "S-adenosylmethionine",
    "SAH": "S-adenosylhomocysteine", "COA": "coenzyme A", "PLP": "pyridoxal 5′-phosphate", "BCL": "bacteriochlorophyll",
    "CLA": "chlorophyll a", "SF4": "[4Fe-4S] cluster", "FES": "[2Fe-2S] cluster", "MG": "magnesium ion",
    "ZN": "zinc ion", "CA": "calcium ion", "NA": "sodium ion", "K": "potassium ion", "MN": "manganese ion",
    "FE": "iron ion", "FE2": "iron(II) ion", "CU": "copper ion", "CO": "cobalt ion", "NI": "nickel ion",
    "CL": "chloride ion", "BR": "bromide ion", "IOD": "iodide ion", "SO4": "sulfate", "PO4": "phosphate",
    "GOL": "glycerol (cryoprotectant)", "EDO": "ethylene glycol (additive)", "PEG": "PEG fragment",
    "ACT": "acetate", "DMS": "DMSO", "MPD": "MPD (additive)", "TRS": "Tris buffer", "EPE": "HEPES buffer",
    "MES": "MES buffer", "NAG": "N-acetylglucosamine (glycan)", "MAN": "mannose", "BMA": "β-mannose",
    "FUC": "fucose", "GAL": "galactose", "GLC": "glucose", "BGC": "β-glucose", "OXY": "dioxygen",
    "CMO": "carbon monoxide", "XE": "xenon", "HOH": "water",
}


# Residue-name sets for the vertex classes of the tetrahedra table (prt / dna / rna / hoh); any water is
# written as HOH. Names outside these sets fall back to the residue kind read from the structure.
_WATER = {"HOH", "WAT", "DOD", "H2O"}
_AMINO = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU",
          "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE", "SEC", "PYL"}
_DNA = {"DA", "DC", "DG", "DT", "DU", "DI"}
_RNA = {"A", "C", "G", "U", "I", "ADE", "CYT", "GUA", "THY", "URA"}
_NUCLEIC = _DNA | _RNA
VERTEX_CLASSES = ("prt", "dna", "rna", "lig", "ion", "hoh")


def vertex_class(names, kinds, has_o2prime=None) -> np.ndarray:
    """Per-residue class code into VERTEX_CLASSES. Unlisted polymer residues (modified amino acids or
    nucleotides) follow the residue kind; a nucleotide with an O2' atom is RNA."""
    out = np.empty(len(names), dtype=np.uint8)
    for i, (n, k) in enumerate(zip(names, kinds)):
        n = str(n).upper()
        if n in _WATER or k == 4:
            c = 5
        elif n in _AMINO:
            c = 0
        elif n in _DNA:
            c = 1
        elif n in _RNA:
            c = 2
        elif k == 0:
            c = 0
        elif k == 1:
            c = 2 if has_o2prime is not None and has_o2prime[i] else 1
        elif k == 3:
            c = 4
        else:
            c = 3
        out[i] = c
    return out


def roman(n: int) -> str:
    vals = [(1000, "m"), (900, "cm"), (500, "d"), (400, "cd"), (100, "c"), (90, "xc"), (50, "l"), (40, "xl"),
            (10, "x"), (9, "ix"), (5, "v"), (4, "iv"), (1, "i")]
    out = ""
    for v, s in vals:
        while n >= v:
            out += s
            n -= v
    return out or "0"


def hetero_description(name: str) -> str:
    return HETERO_DESCRIPTIONS.get(name.upper(), "")


def one_letter(names) -> np.ndarray:
    return np.array([THREE_TO_ONE.get(str(n).upper(), "X") for n in names], dtype=object)


def lookup(table: dict, names, default=np.nan) -> np.ndarray:
    return np.array([table.get(str(n).upper(), default) for n in names], dtype=float)


def residue_class(names, kinds) -> np.ndarray:
    kind_class = {2: "ligand", 3: "ion", 4: "water"}
    out = []
    for n, k in zip(names, kinds):
        out.append(kind_class.get(int(k)) or CLASSES.get(str(n).upper(), "other"))
    return np.array(out, dtype=object)


# --------------------------------------------------------------- Ramachandran regions
RAMA_REGIONS = ("αR", "β", "PPII", "αL", "other")


def rama_region(phi: np.ndarray, psi: np.ndarray) -> np.ndarray:
    """Coarse backbone-conformation basin for each (φ, ψ) pair in degrees ('' where undefined)."""
    phi = np.asarray(phi, float)
    psi = np.asarray(psi, float)
    out = np.full(len(phi), "", dtype=object)
    ok = np.isfinite(phi) & np.isfinite(psi)
    out[ok] = "other"
    alpha_r = ok & (phi >= -160) & (phi <= -20) & (psi >= -120) & (psi <= 50)
    ppii = ok & (phi >= -100) & (phi <= -50) & (psi >= 100)
    beta = ok & (phi >= -180) & (phi <= -45) & ((psi >= 90) | (psi <= -150)) & ~ppii
    alpha_l = ok & (phi >= 30) & (phi <= 100) & (psi >= -20) & (psi <= 80)
    out[alpha_r] = "αR"
    out[beta] = "β"
    out[ppii] = "PPII"
    out[alpha_l] = "αL"
    return out


def rama_outlier(region: np.ndarray, names) -> np.ndarray:
    """True for residues outside every basin, ignoring Gly (which legitimately populates them)."""
    gly = np.array([str(n).upper() == "GLY" for n in names])
    return (region == "other") & ~gly


def plddt_band(plddt: np.ndarray) -> np.ndarray:
    """AlphaFold confidence band: very high (>90), confident (70–90), low (50–70), very low (<50)."""
    p = np.asarray(plddt, float)
    out = np.full(len(p), "", dtype=object)
    out[p > 90] = "very high"
    out[(p > 70) & (p <= 90)] = "confident"
    out[(p > 50) & (p <= 70)] = "low"
    out[np.isfinite(p) & (p <= 50)] = "very low"
    return out


def exposure_class(rsa: np.ndarray) -> np.ndarray:
    """core (< 5 %), buried (< 25 %), intermediate (< 50 %), exposed (≥ 50 %) relative accessibility."""
    r = np.asarray(rsa, float)
    out = np.full(len(r), "", dtype=object)
    out[np.isfinite(r)] = "exposed"
    out[r < 0.5] = "intermediate"
    out[r < 0.25] = "buried"
    out[r < 0.05] = "core"
    return out
