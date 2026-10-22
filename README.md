# bioGeometry

**Computational geometry workbench for biomolecules.**

bioGeometry loads a protein, nucleic-acid or complex structure from the **PDB**, the **AlphaFold DB** or a local file. It describes the structure through the geometry of its atoms: convex hull, Delaunay tessellation, α-shape, Voronoi cells, pockets and interfaces. Five linked views share one selection, so a residue, ligand, ion or domain picked in one view is highlighted in all of them. Every residue gets geometric and biochemical features, and any selection (a binding site, a metal centre, two chains) can be measured.

![bioGeometry with haemoglobin (4HHB) loaded](docs/screenshots/guide/02-app-overview-4hhb.png)

| | |
|---|---|
| ![Zinc site of carbonic anhydrase II](docs/screenshots/guide/05-1ca2-zinc-site-tetrahedra.png) | ![Pockets of thermolysin](docs/screenshots/guide/06-1tlp-pockets.png) |
| Carbonic anhydrase II (1CA2): the zinc with its coordination sphere, and the Delaunay tetrahedra it takes part in | Thermolysin (1TLP): two-probe pockets; the active-site cleft is selected |

**How to use it:** see the [User guide](docs/USER_GUIDE.md). **How it was built, and why:** see [development-notes.md](development-notes.md).

---

## What it does

- **Five linked views:**
  - a sequence viewer, with chains, secondary structure and Pfam / CATH / UniProt domain bars
  - a 3D structure view
  - a 3D geometry view: hull, Delaunay, α-shape, Voronoi
  - a Ramachandran / Delaunay contact-map / AlphaFold PAE panel
  - a data dock of tables and charts
- **One selection everywhere:**
  - Click, ctrl-click or shift-click residues, domains, chains, ligands and ions in any view.
  - The geometry view highlights every Delaunay tetrahedron the selection takes part in.
- **Ligands, ions and coordinated waters** are drawn identically in both 3D views. Dotted lines show metal coordination, hydrogen bonds, water bridges and covalent links.
- **Residue features:**
  - Voronoi volume and packing, hull depth, convex layer, α-surface, Delaunay degree
  - SASA / RSA and depth class (surface, under the surface, buried)
  - secondary structure, φ/ψ/ω and Ramachandran basin
  - class, charge and hydropathy
  - H-bonds, salt bridges, interfaces, ligand sites and domains
- **Tetrahedra table:** every compact Delaunay tetrahedron, with:
  - het ids of the ligands, ions and waters among its vertices (e.g. `ZN`, `NA-HEM-HOH`)
  - edge lengths, face areas, volume, solid angles and circumsphere volume
  - Voronoi volume, SASA and depth (S / U / B) of each vertex
  - tetrahedrality, Tropsha class, α-shape membership, EC class and interaction type
- **Selection measurement:** for any selection, split into groups (domains, chains, ligands, ions, waters):
  - hull and α-shape, radius of gyration
  - buried surface, contacts, metal coordination
  - the tetrahedra it takes part in
- **Analyses:** pockets and cavities, chain interfaces, α-filtration, four-body statistics, and a scaling benchmark.
- **Downloads:** the structure as mmCIF or PDB, and every table as CSV, saved without opening a tab.
- **Large assemblies:** the 70S ribosome (143 k atoms) loads in seconds in *Backbone* mode, within a memory budget that refuses work rather than crashing.

## Algorithms

| Method | Implementation |
|---|---|
| **Delaunay tessellation** | `scipy.spatial.Delaunay` (Qhull) on the analysis point set, with per-tetrahedron volume and circumradius. Edges are kept up to a *max edge* (5 Å all-atom, 10 Å backbone). |
| **Convex hull and onion layers** | `scipy.spatial.ConvexHull`; the layers come from peeling successive hulls. Hull depth is a residue's distance to the hull. |
| **α-shape** | The regularised α-shape: the union of Delaunay tetrahedra with circumradius ≤ α. Its boundary is the set of triangles between a kept and a rejected tetrahedron. Reported: volume, area, connected components and cavities. |
| **α-filtration** | The α-shape's volume, area, components (β₀) and enclosed cavities as α grows. |
| **Voronoi cells** | `scipy.spatial.Voronoi` with the points mirrored across the faces of their bounding box, expanded by a margin, following py_della. Every residue's cell is closed and the cells tile the box. Interior cells pass the Richards test (no vertex farther than the max reach). The packing ratio is the reference residue volume (Zamyatnin) / Voronoi volume. |
| **Pockets / cavities** | Two-probe detection on the Delaunay complex (after Kawabata & Go 2007): tetrahedra outside the small-α shape but inside the large-α shape, split into connected components. A pocket with a mouth is open; one without is an enclosed cavity. Pockets are computed on the polymer alone, so a bound ligand does not fill its own pocket. |
| **Chain interfaces** | Inter-chain Delaunay edges within the contact cutoff, summarised as a chain × chain matrix. |
| **Four-body statistics** | Residue quadruplets of compact tetrahedra, scored by log(f_obs / p_multinomial), in the tradition of Singh, Tropsha & Vaisman (1996). Tetrahedra are also classed by sequence contiguity: {1,1,1,1} … {4}. |
| **Tetrahedrality** | Σ(lᵢ − lⱼ)² / 15 l̄² over the six edges (0 for a regular tetrahedron). |
| **SASA / RSA / depth** | A vectorised Shrake–Rupley SASA (1.4 Å probe), and RSA after Tien et al. (2013). Depth is the mean distance of a residue's atoms to the nearest solvent-exposed atom. A residue is S (surface) at RSA ≥ 20 %, U (just under the surface) at depth ≤ 2 Å, otherwise B (buried). |
| **Secondary structure** | Kabsch–Sander (DSSP) H-bond energy, φ/ψ/ω, and Ramachandran basins (αR, β, PPII, αL). |
| **Binding contacts** | Geometric rules (no hydrogens): metal coordination (N/O/S ≤ 3.0 Å), hydrogen bond (N/O ↔ N/O ≤ 3.5 Å), water bridge, covalent link (≤ 2.0 Å). Only coordinated waters are kept: ≤ 3.0 Å from an ion or ≤ 3.5 Å from a ligand. |
| **Selection geometry** | The hull, α-shape, radius of gyration and κ² of each group. Also SASA alone vs in the structure, buried surface on association, closest atoms, and contacts. |

Structures are read with **gemmi**. All geometry is vectorised **numpy / scipy**. Results travel to the browser as binary buffers, which are drawn with **three.js** (3D) and **d3.js** (sequence, plots, tables).

## Installation

Requirements: Python ≥ 3.10 (developed on 3.14) and an internet connection the first time (frontend packages, RCSB / AlphaFold downloads).

```bash
git clone <repository-url> bioGeometry-workbench
cd bioGeometry-workbench
conda create -n biogeometry python=3.12 -y      # or any virtualenv
conda activate biogeometry
pip install -r requirements.txt                 # reflex, gemmi, numpy, scipy, starlette, uvicorn (+ pytest, biopython)
```

Optional, for the browser checks in `scripts/`: `pip install playwright && playwright install chromium firefox`.

## Basic usage

### The app (sidebar + workbench)

```bash
reflex run                     # first run installs the frontend (~1 min)
# open http://localhost:3100   (backend and API on :8100)
```

1. In the left sidebar, pick a source:
   - **PDB** id
   - **AlphaFold** UniProt accession
   - **Upload** (drop a file)
   - **Local** (a path on this machine; **List**, then pick a file)
   - or click one of the examples (4HHB, 1CA2, 1AON, 7K00, P69905, P04637)
2. Choose **All-atom** (every heavy atom) or **Backbone** (one point per residue: Cα for protein, P for nucleic acids). Ligands, ions and coordinated waters are included in both.
3. Click **Load structure** (or **Load selected file**). The views fill in within seconds; the residue feature table follows.
4. Select something: residues in the sequence, a domain or ligand in the sidebar's **SELECTION** list, an atom in 3D. Then press **⌬ measure**, open the dock tabs, or run an **Analysis** from the sidebar.

Files: `.pdb`, `.ent`, `.cif`, `.mmcif`, also gzipped. Downloads are cached in `data/`.

### The standalone workbench (no sidebar)

The same workbench also runs from the API server alone, configured through the URL. This is handy for linking to a structure, embedding, or headless checks:

```bash
python serve_standalone.py     # :8770
```

```
http://localhost:8770/workbench/index.html?api=http://localhost:8770&pdb=8UCW
  &af=P69905                   an AlphaFold model instead of a PDB entry
  &path=/abs/path/file.cif     a local file (~ allowed; spaces as %20)
  &chains=C,D                  a chain subset (case-sensitive)
  &mode=ca                     Backbone mode (default allatom); &assembly=1 for a biological assembly
  &layers=alpha,voronoi  &analysis=pockets  &tab=tetrahedra  &select=87,92
```

**DNA in the standalone workbench.** The DNA screenshots below were taken this way, so the workbench appears **without the sidebar**: the domain / ligand / ion list sits beside the sequence instead. A structure that contains only DNA or RNA shows an **empty Ramachandran panel**, because nucleotides have no φ/ψ. 8UCW is a DNA polymerase α–DNA complex, so loaded whole, its plot shows the polymerase residues; with `&chains=C,D` (the DNA template and the RNA–DNA primer) the plot is empty. DNA can also be loaded in the full app, through **PDB**, **Upload** or **Local**.

| | |
|---|---|
| ![8UCW in the standalone workbench](docs/screenshots/guide/01-dna-8ucw-standalone.png) | ![8UCW DNA chains only](docs/screenshots/guide/01b-dna-8ucw-nucleic-chains-standalone.png) |
| 8UCW: DNA polymerase α with its DNA template and RNA–DNA primer (`&pdb=8UCW`) | Only the nucleic-acid chains (`&chains=C,D`): an empty Ramachandran panel |

> **Note on nucleic-acid structures.** The workbench loads DNA and RNA and computes their computational geometry as it does for proteins: tessellation, hull, α-shape, Voronoi cells, surface, depth and contacts. Many residue features and analysis labels, however, are still defined for proteins, for example charge, secondary structure and Ramachandran basins. For nucleic acids these features and their nomenclature are therefore not yet accurate, and should be read with care. Nucleic-acid-specific analysis is work in progress.

## Project layout

```
biogeometry/          engine: structure loading, backbone & SS, SASA, features, geometry, analyses,
                      selection geometry, layer cache, memory budget, Starlette API (/bg/*)
biogeometry_app/      Reflex app: the sidebar (sources, selection list, parameters, analyses, downloads)
assets/workbench/     the workbench page (ES modules: store, sequence, view3d, rama, table, dock, main)
serve_standalone.py   API + workbench without Reflex
scripts/              Playwright end-to-end checks, screenshot helper, memory probe
tests/                pytest suite
docs/                 user guide and screenshots
```

## Tests

```bash
pytest                                              # 79 tests, ~21 s
python serve_standalone.py &                        # then, with Playwright:
python scripts/e2e_features.py  OUT_DIR firefox     # views, selection, colours, measure, pockets, downloads …
python scripts/e2e_selection.py OUT_DIR             # coupled selection across all views
python scripts/e2e_reflex.py    OUT_DIR             # through the app's sidebar (needs `reflex run`)
```

## Acknowledgements

bioGeometry combines the geometry engine of **compGeoD3** with the residue-tetrahedra and hull-triangle features of **protein-geometry-dashboard**. The bounding-box Voronoi follows **py_della**. Data come from the RCSB PDB, the AlphaFold Protein Structure Database, PDBe (Pfam / CATH mappings) and UniProt.

## License

bioGeometry is released under the [MIT License](LICENSE). The bundled third-party code in `assets/workbench/vendor/` keeps its own licences: three.js r170 and its add-ons (MIT) and D3 v7 (ISC). See [THIRD-PARTY-LICENSES.md](assets/workbench/vendor/THIRD-PARTY-LICENSES.md). Runtime dependencies installed with pip (Reflex, gemmi, numpy, scipy, Starlette, uvicorn) are under their own licences.
