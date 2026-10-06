# bioGeometry: development notes

A record of how bioGeometry was built, the decisions taken along the way and why, the problems found and how they were fixed, and what is left open. The overview, algorithms and installation are in [README.md](README.md); usage is in the [user guide](docs/USER_GUIDE.md). The work below took place on 2026-10-05 and 2026-10-06, in the `biopyc` conda environment. The project lives in `github-new/gemini/bioGeometry-workbench/` (renamed from `protein-geometry-dashboard-2/` in section 3.18).

---

## 1. Origin

### 1.1 Groundwork: analysing the predecessor projects

Before anything was built, four existing projects were run and analysed. Each analysis answered the same five questions:

- which algorithms are implemented
- how the results are rendered
- what is calculated and output
- which CSV features are generated
- how time and memory behave on large structures

The notes are in the Obsidian vault, under `/home/kestrel/Documents/kestrel_obsidian_vault/bioGeometry packages analysis/`:

| Note | Project | Notable content |
|---|---|---|
| `compGeoD3.md` | **A**: `github-new/compGeoD3/` (Python engine + d3/three.js viewer) | the geometry engine that bioGeometry reuses |
| `protein-geometry-dashboard.md` | **B**: `github-new/gemini/protein-geometry-dashboard/` (Dash) | residue tetrahedra with 20 metrics, hull triangles, CSV exports; commentary on its many AI-generated patches |
| `pyDella.md` | `computational_geometry_projects/py_della/` | includes a compute-resource estimate for a ribosome |
| `mayaviBiomolecules.md` | `computational_geometry_projects/mayaviBiomolecules/` | |

### 1.2 The request

> Examine the code for projects A and B. Integrate the functionality and create a new project at
> `github-new/gemini/protein-geometry-dashboard-2/` as an interactive computational geometry workbench
> for the analysis of three-dimensional biological macromolecules, named **bioGeometry**.

The folder was later renamed `bioGeometry-workbench/` at the user's request (section 3.18).

The requirements were:

- use the `biopyc` conda environment
- ask questions during development, each with a recommended solution
- use a Reflex + d3.js + Python stack (suited to larger structures)
- be creative with the visualisation and the interface
- keep the **selection tightly coupled** across the sequence viewer, molecule viewer and Ramachandran chart, reflected in a data table with its own filtering
- add **table features** that explain structure geometrically and biochemically
- take the running Dash dashboard (project B, `http://localhost:3000`) as a reference
- **fetch, load and analyse from the PDB and AlphaFold**

### 1.3 Decisions taken with the user

| Question | Answer |
|---|---|
| Reference material | Use the files the user saved: `seq-viewer.jpeg` and `seq-viewer-instructions.md` for the sequence viewer, and `bio_dashboard_image_sample*.png` for the layout. They were kept in `docs/reference/` until the user removed it as no longer needed (3.19). |
| Domain annotations | PDBe / SIFTS (Pfam, CATH) plus UniProt features (recommended option) |
| Code base | Copy and extend compGeoD3 (recommended option) |
| `biopyc` environment | Install packages freely |

The main points of `seq-viewer-instructions.md`:

- **Layout:** three scrolling regions: a chain list, a horizontal track of residue squares, and a HETATM list (ligands, ions, coordinated waters). Domains are drawn below the squares.
- **Coupling:** selecting chains, residues, ligands, ions or domains selects them in every other view, and back.
- **Geometry panel:** toggles for hull, Delaunay, α-shape and Voronoi, with their parameters in a separate left panel.
- **Style:** tight spacing and the "Midnight Ink" theme shared with PySysMon, pyQuotes, pyReader and kBar.

Later instructions from the user:

- **During the first build:** keep only the coordinated waters, the ones taking part in ion or ligand binding, not all waters.
- **After an out-of-memory incident:** monitor system resources so processes do not end abruptly (section 6).

---

## 2. Architecture and the reasons for it

| Decision | Why |
|---|---|
| **The Reflex app is only the control panel**; the workbench is a static page (three.js 0.170, d3 v7, ES modules) in an iframe | Reflex state is serialised to the browser on every change, which suits forms and buttons, not arrays of 10⁵–10⁶ coordinates. The workbench talks to the engine directly and receives messages from the panel (`bgWorkbench.handle`). |
| **Binary layer protocol** (from compGeoD3): `meta.json` plus little-endian `.bin` buffers per layer, cached in `data/cache/jobs/<job>/<layer-key>/` | Typed arrays load without JSON parsing; layers are keyed by job and parameters, so repeats are instant. `engine.SCHEMA` (now 10) invalidates the cache when the format changes. |
| **One selection bus** (`store.js`): a `Set` of residue indices with `selection`, `hover`, `focus`, `data`, `features`, `params` and `selstyle` events | Every component speaks residue indices and the same click rules, so coupling is automatic and new views join with a few listeners. |
| The engine API is a Starlette app mounted into Reflex (`api_transformer`) and also served by `serve_standalone.py` | The same engine can run without Reflex, for headless tests and embedding. |
| Own implementations of Kabsch–Sander secondary structure (numpy), Shrake–Rupley SASA (KD-tree), cartoon meshes and interactions | No external DSSP binary or Biopython at runtime; everything vectorised for large assemblies. Each is checked against references (section 5). |
| gemmi for parsing; author chain ids, **case-sensitive** | Ribosomes use both `A` and `a`. A wrong-case id is accepted only when it matches a single chain. |
| Hetero instances numbered with roman numerals (`HEM (i)`, `HEM (ii)`) | Distinguishable in lists, tables and labels. |
| Coordinated waters: O ≤ 3.0 Å from an ion ("ion shell") or ≤ 3.5 Å from a ligand heavy atom ("ligand") | The user's instruction. These are the waters that matter for binding; bulk solvent only adds noise. |
| Midnight Ink theme, in `theme.css` tokens and `assets/biogeometry.css` for the Radix overrides | Consistency with the user's other projects. |

### Layout of the code

```
biogeometry/          engine: structure, backbone, sasa, biochem, annotations, features, selection,
                      geometry + analysis (from compGeoD3), engine (jobs/layers), resources, api
biogeometry_app/      Reflex app: the sidebar (state.py, biogeometry_app.py)
assets/workbench/     workbench page: store, sequence, view3d, rama, table, dock, scrollbars, main
serve_standalone.py   API + workbench without Reflex (port 8770)
scripts/              Playwright end-to-end checks, screenshots, memory probe
tests/                pytest suite (synthetic structures built from internal coordinates + cached 4HHB)
docs/                 USER_GUIDE.md; screenshots/ (guide, seaborn-palette)
```

---

## 3. How the work went

| Section | Round | Topics |
|---|---|---|
| 3.1 | First build | engine, workbench, Reflex panel |
| 3.2 | Hardening | job pinning, tests, large structures, memory budget |
| 3.3 | Feedback 1 | gutter scrollbars, linked cameras, selection colour, selection geometry, logo, structure download |
| 3.4 | Feedback 2 | binding lines, all-atom / Cα modes, shared hetero palette, the selection's tetrahedra |
| 3.5 | Clean-up | caches, references, documentation |
| 3.6 | Feedback 3 | crisp rendering, reset, tetrahedra table in the user's columns, coordination sphere |
| 3.7 | Feedback 4 | S / U / B depth, selection look, bounding-box Voronoi, flicker-free downloads, hetero on top |
| 3.8 | Feedback 5 | multi-select, participating tetrahedra, `depth` column, no shade change, larger nodes |
| 3.9 | Feedback 6 | selected-edge colour, no gold anywhere |
| 3.10 | Feedback 7 | sky-blue mesh, lighter selected edges, ochre ligands |
| 3.11 | Feedback 8 | seaborn palette, identical ligands in both views, sRGB shader fix |
| 3.12 | Feedback 9 | app opened in the browser; selected edges recoloured, not thickened |
| 3.13 | Feedback 10 | selected edges 0.75 px wider |
| 3.14 | Feedback 11 | selection list moved to the sidebar |
| 3.15 | Feedback 12 | readable Ramachandran axes, dotted quadrant lines |
| 3.16 | Feedback 13 | reference dashboard, tagline, DNA test, local DNA files |
| 3.17 | Feedback 14 | *Backbone* instead of *Cα* |
| 3.18 | Feedback 15 | folder rename, README, user guide, screenshots, nucleic-acid note, these notes |
| 3.19 | Feedback 16 | `docs/reference/` removed |
| 3.20 | Feedback 17 | published on GitHub; MIT licence |

### 3.1 First build

**Engine**

- **`structure.py`:** RCSB, AlphaFold DB (with the PAE matrix) and local files; entity names and UniProt accessions read from the mmCIF; coordinated-water filtering; hetero groups with binding sites.
- **`backbone.py`:** φ/ψ/ω, Kabsch–Sander secondary structure, and cartoon ribbon and nucleic-tube meshes sized to the structure.
- **`sasa.py`:** a vectorised Shrake–Rupley solvent-accessible surface.
- **`biochem.py`:** residue tables and Ramachandran basins.
- **`annotations.py`:** PDBe Pfam/CATH and UniProt domains.
- **`features.py`:** the residue table, in eight column groups with descriptions, and the tetrahedra table (20 metrics, vectorised).
- **`engine.py`:** 13 layer kinds and the CSV exports.

**Workbench**

- **Sequence viewer:** virtualised canvas, domain lanes, HETATM list with ◎ site selection.
- **Structure and geometry views:** two three.js views with impostor spheres, render-on-demand and screen-space picking.
- **Ramachandran / contact / PAE panel:** brushable.
- **Data table:** virtual rows, histograms and filters.
- **Analysis dock.**
- **Reflex panel:** sources, examples, geometry parameters, analyses, display, load options and downloads.

**First verification**

- **4HHB:** 4 HEM, 2 PO4 and 8 coordinated waters; Pfam/CATH domains; secondary structure about 81% helix.
- **AlphaFold P69905:** pLDDT 98 and a 142 × 142 PAE matrix.
- **Standalone coupled-selection test:** passed. A drag across the sequence, a Ramachandran brush, ligand-site selection, a table filter followed by *select filtered*, and a 3D pick each kept every view in step.
- **Reflex panel test:** passed.

**Problems fixed during the first build**

- gemmi objects have no `__dict__`, so cif names were kept in a module-level map instead.
- Chain descriptions showed "protein" until `_entity.pdbx_description` was read.
- A favicon 404 appeared in the console.
- Stale cached layers needed a schema bump (3 → 4).
- The Playwright and Chromium versions did not match.

### 3.2 Hardening

| Change | Detail |
|---|---|
| Switching structures could show "layer not found" | A `params` message started a layer refresh for the old job, which then fetched under the new job id. The job is now pinned through every request, and late results are dropped. |
| PAE plot looked uniformly green | The data were correct (median 2 Å for P69905); a 0–30 Å colour bar was added. |
| Test suite | 35 engine and API tests were added to the 16 geometry tests copied from compGeoD3. Synthetic structures are built from internal coordinates (ideal helix, strand, zinc site), so expected values are exact. |
| `_chain_info` crashed on PDB files with ion- or water-only chains | Found by the new tests. Subchains are now walked. |
| `read_meta` accepted unvalidated path parts | Now shares the `_safe()` guard with `buffer_path`. |
| The e2e tetrahedra check passed on leftover residue rows | It now waits for the tetrahedra header itself. |
| Large structure (7K00, 143 k atoms) | Opens in 4 s, with every geometry layer ≤ 1 s. Features took 13.5 s, 10 s of it SASA. Bounding the KD-tree search by the largest inflated radius made SASA 1.7× faster with identical results (7K00 10.2 → 6.1 s, 1AON 4.2 → 2.3 s). |
| Out-of-memory incident | See section 6: resource monitor, memory caps, and an engine memory budget measured with `scripts/memprobe.py`. |
| Reflex served stale workbench files | `assets/` had been excluded from hot reload, and Reflex copies assets only when it recompiles. `assets/` is watched again; only runtime data, scripts, tests and docs are excluded. |

### 3.3 First round of user feedback: scrollbars, cameras, selection geometry

| Request | What was done |
|---|---|
| The data-table scrollbar overlaps rows on hover | The cause was Firefox/GTK overlay scrollbars, which CSS cannot turn off. `scrollbars.js` hides the native bars and draws its own in reserved gutters, keeping native scrolling. Applied to the tables, dock panels, sequence strip and side lists. Verified in Firefox with Playwright's Firefox build (installed for this). |
| Move and zoom both 3D views together | The camera link is on by default and works both ways; tests check the two cameras stay identical to within 1×10⁻³. |
| Nicer selection colours; highlight the simplices of a selection | A single selection-colour token (`--sel`, ◉ picker), cylinder bonds and ligand glow were added. Two real defects were found: *show tetrahedra* never rendered, and a selected ligand was invisible in the geometry view. Both fixed; the Delaunay layer now ships its compact tetrahedra. |
| Domain bars and HETATM entries as selection tools; compute the geometry of the selection | Domain bars became clickable, with a DOMAINS list beside them. **Selection geometry** (`selection.py`, the ⌬ Selection tab) reports groups, hulls, α-shape, contacts, metal coordination and buried surface. Validated on the haemoglobin haem site (HIS87 NE2–Fe 2.14 Å, 83% buried) and the carbonic anhydrase zinc (His94/96/119 at 1.9–2.1 Å). |
| A tetrahedron logo | Two lit faces, a dashed hidden edge, and atoms at the vertices; used in the sidebar, the favicon and the empty screen. |
| Download the structure as PDB or CIF | The loaded model is written through gemmi. PDB chain ids are shortened when needed, and PDB is refused above 62 chains. |

The geometry toolbar outgrew its panel, so the display toggles moved into a *show ▾* menu.

### 3.4 Second round of user feedback: binding lines and analysis modes

| Request | What was done |
|---|---|
| Dotted lines between residues and ligands or ions | `features.binding_contacts`: metal coordination, hydrogen bonds, water bridges and covalent links, drawn as dotted lines with a legend key. It finds the haem Fe–His bond in all four haemoglobin chains and the Zn(His)₃(H₂O) site in carbonic anhydrase. |
| Tidy the analysis buttons; "Pockets/cavities" | An even 2 × 3 grid with single-line labels. |
| Choose before loading: all-atom or Cα, each with coordinating waters, ligands and ions | An *All-atom / Cα* switch. Ligand, ion and coordinated-water atoms join the point set at full detail in both modes. Since a bound ligand would then fill its own pocket, **pockets are computed on the polymer alone**; a test checks they match the polymer-only job exactly. |
| Ligands, ions and waters in the same colours in both views, different from everything else | A shared palette: gold ligand carbons, ions in Jmol element colours, red waters. Ligands are ball-and-stick in both views. |
| Highlight the nodes and edges of the tetrahedra, not the simplices | The filled faces were replaced by the *star* of the selection: every tetrahedron with a selected corner, drawn as its nodes and edges. |
| Note that not every structure has pockets | Added to the Pockets tab and the Analyses panel. |
| The selection is too glaring | The selection colour is now mixed more lightly. Unselected parts lose their hue (fade to grey) instead of keeping a darkened scheme colour, which could be mistaken for the selection. |
| What is the second highlight in the geometry view? | A bug: the selection-geometry overlay used raw coordinates while every other layer is centred on the structure, so it was drawn one centroid-offset away. Fixed and covered by a test. |

During this round the end-to-end console check caught a regression: a refactor had deleted the overlay renderer. It was restored before release.

### 3.5 Clean-up and documentation

- **Removed:** `__pycache__` folders, two reference copies of compGeoD3 JavaScript used while porting, the empty `uploaded_files/`, and 87 MB of layer caches from superseded schemas (they rebuild on demand).
- **Moved:** the user's design references to `docs/reference/` (removed in 3.19).
- **Rewritten:** `README.md`, `.gitignore` and `requirements.txt`.
- **Added:** these notes.

### 3.6 Third round of user feedback: rendering, reset, tetrahedra table

| Request | What was done |
|---|---|
| The Plotly charts looked sharper. Can bioGeometry match them? | The reference dashboard draws fixed-size screen markers (4 px) and lines 2–6 px wide, with no fog. bioGeometry now matches it on four points: (1) supersampling to at least 1.5× the screen's pixel density; (2) anti-aliased sphere edges (alpha-to-coverage under MSAA); (3) no depth fog; (4) thick lines using three.js's `LineSegments2` / `LineMaterial` add-ons, vendored at r170, since WebGL lines are always 1 px. |
| Smaller points in the geometry chart | Points are now fixed-size **screen-space dots** (about 4 px all-atom, 6.5 px Cα), so they no longer swell when you zoom in. *show ▾* offers small / medium / large. |
| Selection like the Plotly dashboard: thicker lines, a faint highlight, no glowing gold | The glow shader was removed. The selection is now a faint cartoon tint, a 5 px backbone line in the selection colour, and thin sticks. In the geometry view, the tetrahedra around the selection are drawn as 4 px lines (between selected nodes) and 2 px lines (to neighbours). The default colour is cyan, as in the reference; the gold the user saw was the amber preset they had picked. |
| Pockets cannot be unselected | A second click now clears the pocket's residues too. The tab gained *pockets shown / hidden* and *✕ clear*, *show ▾* gained an overlay switch, and Esc drops the pocket highlight. |
| After *Measure selection*, or after any analysis, the Residues and Tetrahedra tabs were blank | The dock decided whether a table was still mounted by checking the first child of the dock body, and an analysis tab's charts satisfied that check. Tables now keep a reference to their own root element. |
| The reset button should give a clean slate | `reset_workbench` (Reflex) and `resetWorkbench` (workbench, also `Shift+X`) clear the selection, measurements, overlays and table filters; restore the default layers, parameters, display and toggles; and reload the structure from cache. |
| Tetrahedra with the selected residues should come first | The table partitions its rows (stably) into "involves the selection" and the rest, re-orders on every selection change, and reports the count. |
| Clearer CSV file names | `bioGeometry_<id>-<chains>_<atom set>_<content>.csv`, for example `bioGeometry_4HHB-ABCD_all-atom_delaunay-tetrahedra.csv`; the dock's view exports follow the same pattern. |
| A button for the coordinating waters and residues of a HETATM | **⌬** in the HETATM list selects a ligand or ion's coordination sphere from the binding contacts: the group, the residues bound to it, its coordinating waters, and the residues those waters bridge to. It then flies to the site and measures it. For the haem in 4HHB this is HIS45, HIS58, LYS61, HIS87 and two waters. |
| The ideal tetrahedra output (59 named columns plus two text columns) | A new table in the user's column order and names, over the compact Delaunay tetrahedra of the point set, with ligand, ion and water vertices. New calculations: EC numbers (mmCIF `_entity.pdbx_ec`, else UniProt) and their class; per-vertex residue classes with the user's name sets; Voronoi cell areas; per-vertex SASA; residue depth classes; tetrahedrality; Tropsha class; and interaction types. |
| Residue-name sets for prt / dna / rna / hoh | Adopted as given (`biochem._AMINO`, `_DNA`, `_RNA`, `_WATER`; any water written as HOH). Residues outside the sets fall back to the kind recorded in the file, so a modified amino acid in a protein chain stays `prt`. |

**Two definitions were checked against sources rather than assumed:**

- **Tropsha class:** the five classes of Singh, Tropsha & Vaisman (J. Comput. Biol. 3:213, 1996) are the sizes of the sequence-contiguous runs among the four residues ({4}, {3,1}, {2,2}, {2,1,1}, {1,1,1,1}). This is finer than *simplex-type*, the longest run, which cannot tell {2,2} from {2,1,1}.
- **Tetrahedrality:** Σ(lᵢ − lⱼ)² / 15 l̄², which is 0 for a regular tetrahedron.

**Residue depth was redesigned after a test failed.** The first definition (distance to the nearest exposed atom of another residue) left 4HHB with no buried residues: in a packed protein almost every residue touches some exposed atom. The final definition:

- **surface:** RSA ≥ 20%
- **depth:** the mean distance of the residue's atoms to the nearest exposed atom (in the spirit of Chakravarty & Varadarajan's residue depth)
- **just under the surface:** depth ≤ 2 Å; **buried:** depth > 2 Å

That threshold gives 43–48% surface, 27–37% just under the surface and 20–25% buried on haemoglobin and carbonic anhydrase.

### 3.7 Fourth round of user feedback: depth, Voronoi, downloads

| Request | What was done |
|---|---|
| Residue depth as S / U / B | The depth classes are now `S` (surface), `U` (just under the surface) and `B` (buried), in both the residue and tetrahedra tables. |
| The view after selecting a domain looks ghastly | Two causes. In the structure view, the thick Cα line cut across the ribbons in zig-zags, and 31 residues of side-chain sticks added clutter. In the geometry view, thick star edges and the measurement's hull wireframe overlapped. Now selected residues get a light tint plus an outline around their own cartoon (an enlarged back-face shell); sticks appear only for up to 16 residues; star edges are thinner, and finer still for large selections; and the measurement envelope is faint, with no wireframe, shown only while the Selection tab is open. |
| Voronoi volumes for all residues, using py_della's bounding box | `voronoi_cells(..., box_margin)` mirrors the points across the faces of the expanded bounding box: first those near each face, then all of them if any cell is still open. Every cell closes. The cells tile the box exactly (tested), and interior cells are identical to the unclipped diagram. A *Voronoi interior* flag keeps the Richards test; the packing ratio uses interior residues only; the panel gained a *box margin* slider and an *all cells* option. |
| Downloads flicker (a new tab opens briefly) | Downloads are now fetched by the page and saved through a temporary link, from both the panel and the dock's *files* menu. Because the API is on another port, its responses now expose `Content-Disposition` (`Access-Control-Expose-Headers`), so the browser keeps the proper file name. |
| Is the pockets *✕ clear* button needed? | No. It did the same as clicking the selected pocket again or pressing Esc, so it was removed. *pockets shown / hidden* stays. |
| Ligands, waters and ions are hidden by Delaunay edges and nodes | In the geometry view they are now drawn as in the structure view (ball-and-stick, spheres) in a second pass with a fresh depth buffer: always visible, and still depth-correct among themselves. Their points in the tessellation are not drawn twice. |

**A regression introduced and caught in this round.** Replacing the geometry view's hetero-drawing method cut up to a comment further down the file, silently deleting `setCartoon`, `buildLayer` and `updateSelection`. The first end-to-end run failed on the missing Delaunay layer, and the three methods were restored from this session's record of them.

### 3.8 Fifth round of user feedback: multi-select and the selection's tetrahedra

| Request | What was done |
|---|---|
| Multi-residue selection from the sequence viewer and every other tool; **clear** must clear | Ctrl- and shift-clicks already worked, but nothing showed they existed. A **＋ multi** toggle (`M`) now makes plain clicks toggle and drags / brushes add, in every view: sequence squares and drags, domain bars, chains, 3D picks, table rows, Ramachandran and matrix brushes, pocket and interface rows. All of these now go through `store.toggling / additive / addOrSet / pick`. **clear** and Esc empty the selection, its measurement and the pocket / interface highlights. |
| A selection should select every tetrahedron the residues take part in (edges and nodes), and *measure* should use the same set | Previously the edges and nodes reaching neighbouring residues were drawn muted and read as unselected. Now every node and edge of those tetrahedra is in the selection colour. The measurement adds `participating_tetrahedra` (computed on the full tessellation, with the same max-edge cutoff the view uses): count, wholly-inside count, nodes, edges, volume, mean edge, neighbouring residues. A test shows it equals the Tetrahedra-table rows that involve the selection. |
| A `depth` column joining the four residue depths (SSUU, SUSU, BUBU …) | Added after `residue-4-depth`, in the table and the CSV. |
| The whole workbench changed shade under a selection | Every view faded or desaturated what was not selected. Now nothing outside the selection changes colour. Dimming became the optional *dim* toggle, off by default and applied to all views together. |
| Delaunay nodes 2 px larger | 6.2 px all-atom, 8.5 px Cα (S / M / L still scale them). |

### 3.9 Sixth round of user feedback: selected-edge colour, no gold

| Request | What was done |
|---|---|
| Colour the selected tetrahedra's edges (#e9a3c9, #8a66d2, #87627a tried in turn), 1 px thicker, nodes 1 px larger; "or try a few colours" | Edges are 2.6 px (3.6 px between selected residues, finer above 300 tetrahedra) and nodes are 1 px larger. Ten candidates were rendered on the same 4HHB haem selection. The user then tried **#9a45d1** and then **#922472** (plum), which is now `STAR_EDGE` in `view3d.js`. The comparison shots were saved in `docs/screenshots/edge-colours/` (removed in 3.10): one per colour (11), `edge-colours-sheet.png`, and `workbench-selection-<hex>.png` for the two built-in colours. Candidates are recoloured in linear colour space (as three.js stores them), so each shot matches what the app draws: each recoloured shot differs from the built-in rendering by at most 11/255 per pixel. |
| **No golden colour anywhere** | Removed every gold, amber and yellow colour. The `butter` token became `teal` (#5eead4) in the CSS, JS and Reflex palettes and in the logo / favicon vertex. Ligand carbons became emerald (#34d399), later yellow ochre (see 3.10); hydrogen-bond beads are pink. Sulfur and phosphorus use lime / salmon. The pLDDT 50–70 band is violet. The `amber` selection preset became `violet`. Chain colours past ten skip the 20°–80° hue band. The d3 schemes that pass through yellow (viridis, plasma, cividis, turbo, YlOrRd, RdYlBu, Tableau10, rainbow) were replaced by `RAMP` / `CAT10` in `store.js`, and the contact-map ramp ends in pink. A scan of every colour literal in `assets/`, `biogeometry_app/` and `biogeometry/` finds no hue between 25° and 70°, apart from the ochre ligand carbons the user asked for in 3.10. |

### 3.10 Seventh round of user feedback: sky-blue mesh, ochre ligands

| Request | What was done |
|---|---|
| "The Delaunay mesh does not appear sky-blue" | Each edge was coloured by length, blending sky (#7dd3fc) into lavender (#c4b5fd). Most Cα edges are long, so the mesh looked lavender-white. The mesh is now one sky blue, `MESH_EDGE = #38bdf8`, with the legend swatch to match. |
| Selected edges a shade lighter or darker than the mesh | `STAR_EDGE = #bae6fd` (sky-200, lighter) is built in. The darker alternative #0369a1 (sky-700) is rendered for comparison. The lighter shade separates the selection more clearly on the dark background, and the darker one gives the ligand the most contrast. |
| Ligands: a slightly dark yellow ochre instead of emerald | `LIGAND_CARBON` and the *ligand* class colour are now #b5872f (hue 39°, 59 % saturation, 45 % lightness). This is the user's one explicit exception to "no gold". It is muted and earthy, not bright gold. |
| A few screenshots with this combination; remove the older ones | `docs/screenshots/edge-colours/` was deleted. `docs/screenshots/sky-mesh/` held 11 shots (removed in 3.11) of 4HHB chain A with the haem pocket selected: Cα and heavy-atom workbench views, lighter vs darker selected edges (`00-lighter-vs-darker.png`), and haem close-ups in both views. A `~/.cache` cleanup had removed the Playwright browsers. They were reinstalled (newer Chromium build), and the e2e scripts now pick the newest `chromium-*` build instead of a fixed version. |

### 3.11 Eighth round of user feedback: seaborn palette, identical ligands

| Request | What was done |
|---|---|
| Delaunay edges #42cbee, nodes #5eabd9 | `MESH_EDGE = #42cbee`. Selected edges are a shade lighter, `STAR_EDGE = #aae8f7`. The user chose to make #5eabd9 the first chain colour, so chain A's nodes, cartoon and chips share it. Other chains and other colour schemes still colour the nodes. |
| HEM and other complex ligands: #c9b26e, #bd8b64, #c15b2a | The user chose colouring by role: carbons `LIGAND_CARBON = #c9b26e`, all ligand bonds `LIGAND_BOND = #bd8b64` (one colour instead of per-atom halves), and metal atoms inside a ligand (haem Fe, Mg, Co …) `LIGAND_METAL = #c15b2a`. N, O, S, P and halogens keep element colours. The *ligand* class colour is #c9b26e. |
| Ions and waters: "refer to the url and optimize" | Colours from seaborn's categorical palettes (deep / muted / pastel / colorblind / Set2), grouped by chemistry: group 12 Zn pastel violet #d0bbff (Cd/Hg pastel grey); alkaline earth greens (Mg #6acc64, Ca #8de5a1, Sr/Ba #55a868); alkali purples (Na #956cb4, K #8172b3); Fe = haem iron #c15b2a; Cu #ffb482, Mn #dc7ec0, Co #fab0e4, Ni #029e73; halides Cl #a6d854, Br #8c8c8c, I #e78ac3, F #66c2a5. Waters are seaborn muted red #d65f5f. Each was checked against the scene colours (mesh, nodes, selected edges, selection cyan, ligand roles, N, O): CIELAB ΔE > 20, and contrast ≥ 4 : 1 against the background. The geometry legend lists each ion element present, plus *bonds* and *metal centre*. |
| Ligands must look the same in the structure and geometry viewers | The colour arrays were already equal, but three things made the ligand look different. (1) The geometry view tinted a selected ligand 20 % towards the selection cyan, which turned ochre olive; the structure view only enlarged it. (2) The geometry view's overlay scene (ligands drawn on top) had different lights: ambient 0.55 and key 1.4 with no rim light, against 0.5, 1.5 and rim 0.35. (3) Two near-copies of the builder could drift apart. Now `buildHetBallStick` / `paintHetBallStick` serve both views (selected ligands grow, colours stay) and `addLights` gives both scenes the same lights. The e2e check compares the painted colours, radii, bond colours and lights of a selected haem in both views (all identical). The screenshot script compares rendered pixels at every hetero atom centre with the same camera: the difference is 0/255 at 45 + 39 + 47 atoms, none of them on background. |
| (found while checking) Spheres rendered darker than their colour | The sphere impostor shader wrote three.js's linear colour values straight to the screen without the sRGB step, so every sphere (geometry nodes, ligand / ion / water atoms, sphere style) was darker and more saturated than its hex. For example, ochre #b5872f looked olive-brown and #c9b26e came out near #957228. Adding `#include <colorspace_fragment>` fixed it. A Zn sphere centre now renders #c0adec = 0.92 × #d0bbff, which is the shader's lighting at the centre. Fat lines and Phong meshes already converted correctly. |
| New screenshots; clear the old ones | `docs/screenshots/sky-mesh/` was deleted. `docs/screenshots/seaborn-palette/` has 10 shots: 4HHB chain A (Cα and heavy atoms, haem pocket selected, haem close-up structure vs geometry) and thermolysin 1TLP (Zn, 4 Ca, coordinated waters, inhibitor RDF; zinc-site and calcium-site close-ups, structure vs geometry, captioned with the pixel comparison). |

### 3.12 Ninth round of user feedback: selected edges recoloured

| Request | What was done |
|---|---|
| Start the application and open it in a browser | The app was already running (3100 / 8100), so it was opened in Firefox, the default browser. |
| On selection, do not thicken the Delaunay edges; only change their colour | Selected edges used to be separate, thicker lines drawn over the mesh (2.6 / 3.6 px). Now `updateSelection` recolours the mesh's own edges in place to `STAR_EDGE` (#aae8f7), so width and opacity cannot differ (1.6 px / 0.85 Cα, 1.1 px / 0.55 all-atom, from `meshLine()`). Edges of those tetrahedra that are longer than the mesh's cutoff, if any, are drawn with the same width and opacity. The e2e check confirms that every selected edge is a recoloured mesh edge (6,743 of 6,743, none extra) and that only one width / opacity is in use. The nodes keep their +1 px. The screenshots in `docs/screenshots/seaborn-palette/` were re-taken. |

### 3.13 Tenth round of user feedback: selected edges 0.75 px wider

| Request | What was done |
|---|---|
| Selected Delaunay edges 0.75 px wider | A line set has one width, so the selected edges are drawn separately, `STAR_WIDEN = 0.75` px wider than `meshLine()` (Cα 1.6 → 2.35 px, all-atom 1.1 → 1.85 px), at the same opacity, in `STAR_EDGE`. To avoid drawing them twice, which would double their opacity, `updateSelection` rewrites the mesh's segment buffer with the other edges first and sets `instanceCount` to their number. The selected edges are drawn only by the overlay. Clearing the selection restores the full mesh. The e2e checks confirm: the moved edges plus any beyond the cutoff equal the selection's edges, each is drawn once, the overlay is exactly 0.75 px wider at the same opacity, and after *clear* the mesh draws all its edges again (28,848 of 28,848) with no overlay. The screenshots were re-taken. |

### 3.14 Eleventh round of user feedback: selection list in the sidebar

| Request | What was done |
|---|---|
| The selection panel in the sequence viewer's corner is congested: move it below *Load structure* and above *Geometry parameters*, with the code and logic unchanged | The sidebar (`biogeometry_app.py`) has a **SELECTION** section (`selection_panel()`: a heading and an empty `#bg-selection-panel` in a 26vh scroll box) between the source / status block and the accordion. The app page and the workbench iframe are same-origin, so `sequence.js` renders the same list into `window.parent.document.getElementById('bg-selection-panel')` instead of `#seq-het`, with the same `renderHetero` / `markHetero` code, checkbox, ◎ and ⌬ handlers, and store. The sequence viewer then drops its third column (`.seq-body.no-het`), and the selection accent (`--sel`) is copied to the host on `selstyle`. The panel's CSS is mirrored in `assets/biogeometry.css` under `#bg-selection-panel`. With no parent page, or from another origin, the panel stays in the workbench, so the standalone server (8770) and its e2e tests are unchanged. `e2e_reflex` step 1b checks that the sidebar lists the groups, that the iframe's own panel is hidden, the position (Load structure < SELECTION < Geometry parameters), that ticking HEM in the sidebar selects it in the workbench, and that clearing in the workbench unticks it. |

### 3.15 Twelfth round of user feedback: Ramachandran axes

| Request | What was done |
|---|---|
| Ramachandran axes and axis labels hard to read: use the colour of the *select / filters / cols* buttons | Those chips use `--label` (#9a9ab8). The panel's axis lines (previously `--edge-hi` #3a3a58), tick labels and the φ / ψ titles (previously `--dim` #62627e, now class `axis-title`) use it too. The rule is scoped to `#rama-svg .axis`, so the contact-map and PAE axes in the same panel follow. Computed colours were checked: chip text, axis line, tick label and title are all rgb(154, 154, 184). |
| The same for the 0 / 0 grid lines, so the quadrants are clearly visible | The φ = 0 and ψ = 0 grid lines get class `zero`, in `--label`: first solid, then dotted at the user's request (`stroke-dasharray: 1 2.5`). The ±60° and ±120° lines stay faint and dashed (`--edge`). |

### 3.16 Thirteenth round of user feedback: reference dashboard, tagline, DNA

| Request | What was done |
|---|---|
| Start `protein-geometry-dashboard` (project B) and open it in another tab | Its `dashboard/requirements.txt` pins dash-mantine-components 0.12.1, while `biopyc` has 2.8, whose API is incompatible. The project's own `venv/` does not run on this machine. Rather than downgrade `biopyc`, a uv virtual environment was rebuilt from the requirements in the session's scratch folder; it must be rebuilt after `/tmp` is cleaned. The dashboard runs on :3000 under a 3 GB cap and was opened in a new Firefox tab. |
| Tagline "computational geometry workbench for biomolecules" | Changed in the sidebar brand (`biogeometry_app.py`); it was the only copy. |
| Load a DNA structure to examine | 1BNA (B-DNA dodecamer, 2 chains, 24 nt, 486 atoms) loads cleanly in the standalone workbench: 3,015 Delaunay tetrahedra, every residue classed purine / pyrimidine, no console errors. Its Ramachandran plot is empty because nucleotides have no φ / ψ. It was opened for the user at `localhost:8770/workbench/index.html?api=http://localhost:8770&pdb=1BNA`. |
| How to load local DNA structures, given the standalone workbench | Both ways already existed. In the app: **Upload** or **Local** (a path, **List**, pick a file). In the standalone workbench: the `&path=` URL parameter (`~` is expanded, spaces are written `%20`). A local copy of 1BNA loaded this way. The user was told that DNA is not limited to the standalone workbench. |

### 3.17 Fourteenth round of user feedback: Backbone mode label

| Request | What was done |
|---|---|
| Rename the "Cα" toggle to "backbone" (Cα is protein-only; DNA has a sugar-phosphate backbone) and its description to "Backbone residue(s) + ligands, ions & coordinated waters · applies on Load" | `mode_picker()`: the toggle reads **All-atom / Backbone** and the description is exactly as requested. For consistency, the sidebar status line ("… backbone + hetero …") and the 1AON / 7K00 example tooltips also say "backbone". The internal mode value is still `ca` (level `ca`: Cα for protein, P for nucleic acids), so nothing else changes. `e2e_reflex` step 6 clicks "Backbone" and passes. |

### 3.18 Fifteenth round of user feedback: rename, README, user guide

| Request | What was done |
|---|---|
| Name the project / folder "bioGeometry-workbench"; the interface stays "bioGeometry" | The folder `protein-geometry-dashboard-2` was renamed to `bioGeometry-workbench` with the servers stopped, then restarted from there. Nothing depended on the path: no editable install, and no absolute paths in the code or the `.web` build. |
| Screenshots: the DNA structure 8UCW and several proteins / enzymes | `docs/screenshots/guide/` (12 images, taken through the real app with Playwright / Firefox). 8UCW is the DNA polymerase α catalytic subunit with its DNA template and RNA–DNA primer (plus dGTP, GTP and Mg²⁺). It contains the protein, so its Ramachandran plot is not empty; it was shot in the standalone workbench whole, and with only its nucleic-acid chains (`&chains=C,D`, an empty Ramachandran panel). The others: 4HHB overview, sidebar, and haem coordination sphere measured; 1CA2 zinc site with the Tetrahedra tab; 1TLP pockets; 1HSG chain interfaces, the four geometry layers, and α-filtration; 1AON in Backbone mode; AlphaFold p53 with pLDDT and PAE. |
| New README (old one removed): description, algorithms, installation, basic usage; the DNA / standalone note | `README.md` rewritten. The DNA note is worded accurately: the standalone workbench has no sidebar; only a nucleic-acid-only structure has an empty Ramachandran panel (8UCW contains the polymerase); DNA also loads in the full app. The old README's feature reference moved, condensed, into the user guide. |
| A note for nucleic-acid structures after the 8UCW screenshot (sent while the README was being written) | Added after the 8UCW screenshots in the README, and in the user guide's DNA section. The wording: the workbench loads DNA and RNA and computes their geometry as for proteins, but many residue features and analysis labels (e.g. charge, secondary structure, Ramachandran basins) are still defined for proteins, so for nucleic acids they are not yet accurate; nucleic-acid-specific analysis is work in progress. |
| A user guide with screenshots | `docs/USER_GUIDE.md`: screen layout, loading, selecting, the structure and geometry views, Ramachandran / contacts / PAE, the dock and its filters, measuring, analyses, large assemblies and AlphaFold, DNA and the standalone workbench (with the nucleic-acid caveat), downloads, keys. Every image link resolves; the 79 tests pass from the renamed folder. |
| Record these changes here, and put this file in order | The rounds above, plus: the header covers both days and points to the user guide; the quoted original request is restored, with the rename noted after it; `SCHEMA` reads 10; the code layout lists `serve_standalone.py` and `docs/`; a rounds-at-a-glance table opens section 3 and every heading names its topics; deleted screenshot folders are marked; sections 4–8 are brought up to date. |

### 3.19 Sixteenth round of user feedback: design references removed

| Request | What was done |
|---|---|
| `docs/reference/` and its files are not needed: remove them | Deleted the folder (1.7 MB): `seq-viewer.jpeg`, `seq-viewer-instructions.md`, `bio_dashboard_image_sample.png` and `bio_dashboard_image_sample_02.png`, the design references from 1.3. The mentions in the README's project layout and in these notes (1.3, the code layout, 3.5) were updated. Nothing else used the folder. |

### 3.20 Seventeenth round of user feedback: GitHub and licence

| Request | What was done |
|---|---|
| Commit to GitHub as a public repository `bioGeometry-workbench` (user `y-v-j`) | `git init` on `main`; the existing `.gitignore` keeps out `data/`, `.web/`, `.states/`, `reflex.lock`, `uploaded_files/` and Python caches, leaving 72 files. Before publishing, a scan found no tokens, keys, passwords or email addresses. These notes do contain local paths, which the user was told about. The commit author is the user's GitHub no-reply address (set for this repository only), so no personal email is published. `gh repo create --public --push` created https://github.com/y-v-j/bioGeometry-workbench. |
| Add an appropriate licence | **MIT** (`LICENSE`, © 2026 Yogesh (y-v-j)). It is permissive like the dependencies (Reflex Apache-2.0, gemmi MPL-2.0, numpy / scipy / Starlette / uvicorn BSD) and compatible with the bundled code. The bundled three.js r170 with its add-ons (MIT) and D3 v7.9.0 (ISC) require their notices to travel with copies, so their full texts are in `assets/workbench/vendor/THIRD-PARTY-LICENSES.md`. The README gained a License section. |

---

## 4. Defects found, and how

| Defect | Found by | Fix |
|---|---|---|
| gemmi `Structure` has no `__dict__` | first run | module-level name map |
| Chain description shows "protein" | screenshot review | read `_entity.pdbx_description` |
| Stale layers after format changes | testing | `SCHEMA` bumps (3 → 10 over the project) |
| "layer not found" after switching structure | screenshot review | job pinned through requests |
| `_chain_info` crash on hetero-only chains | new unit tests | walk subchains |
| Unvalidated path in `read_meta` | code review | shared `_safe()` guard |
| e2e tetrahedra check passing on stale rows | review of the test itself | wait for the tetrahedra header |
| Reflex serving stale workbench JS | served-file check | stop excluding `assets/` from hot reload |
| Overlay scrollbars covering table rows | user (Firefox) | gutter scrollbars |
| *show tetrahedra* never drawn; selected ligand invisible in geometry view | user question → inspection | render them |
| Selection overlay displaced by the centroid | user screenshot | centre the layer; regression test |
| Overlay renderer deleted by a refactor | e2e console-error check | restored |
| Test-order artefacts (focus switched chain, off-screen virtual row) | e2e failures | tests fixed; the app was correct |
| Residues / Tetrahedra blank after an analysis tab | user report → reproduction script | tables track their own root element |
| A pocket could not be unselected | user report → reproduction | second click clears the residues; show / hide / clear controls |
| The table info text covered the dock tabs | the new end-to-end test (a click intercepted) | compact text with the details in a tooltip, shorter toolbar labels |
| No buried residues under the first depth definition | new unit test | RSA-based surface + mean residue depth |
| The water group label showed the structure-wide count | screenshot review | label counts the selected waters |
| Three geometry-view methods deleted by an over-wide edit | end-to-end test (Delaunay layer missing) | methods restored |
| Downloads saved under fallback names in the Reflex app | end-to-end test (cross-port header hidden) | expose `Content-Disposition` |
| Colour-comparison candidates rendered lighter than the app draws them | pixel comparison with the built-in colour | recolour in linear colour space, as three.js stores colours (screenshot script) |
| e2e scripts pinned to one Chromium build; the browsers deleted by a `~/.cache` cleanup | e2e launch failure | reinstall; pick the newest `chromium-*` build |
| Delaunay mesh looked lavender-white, not sky blue | user report | one edge colour instead of a length ramp |
| A selected ligand tinted towards the selection colour in the geometry view only (ochre turned olive) | user report → inspection | one shared builder and painter for both views; no tint |
| The geometry view's ligand overlay lit differently from the structure view | inspection | the same lights in both scenes (`addLights`) |
| Spheres darker and more saturated than their colour | rendered-pixel check | the sphere shader gained the sRGB step (`colorspace_fragment`) |

---

## 5. Verification

- **pytest:** 79 tests in about 21 s.
  - Exact checks on synthetic structures built from internal coordinates: ideal helix φ/ψ = −57/−47, DSSP helix, a zinc site.
  - Analytic checks: sphere and two-sphere SASA; regular and right-corner tetrahedra.
  - References: SASA within 3% of Biopython (−0.77% and r = 0.999 when first compared); haem and zinc sites; every layer kind; the API; exports; the memory budget; both analysis modes; polymer-only pockets; overlay centring.
  - Bounding-box Voronoi: cells close and tile the box exactly, interior cells unchanged, a volume for every residue.
  - Tetrahedra table: the reference column order and CSV header, tetrahedrality, Tropsha class and simplex type re-derived independently for every row, EC and interaction types at the zinc site, the residue-class name sets, and depth classes.
- **End-to-end (Playwright, software WebGL):**
  - `e2e_selection.py`: coupled selection across all views.
  - `e2e_features.py`, in Chromium **and Firefox**: linked cameras, gutter scrollbars, domain + ligand → ⌬ measure, binding lines, tetrahedra highlight, selection colour, downloads. Also: a selected haem painted identically in both views (colours, sizes, bond colours, lights); the selected edges drawn once, 0.75 px wider than the mesh at its opacity; the full mesh restored after *clear*.
  - `e2e_reflex.py`: through the Reflex sidebar, including both analysis modes (the *Backbone* toggle), the downloads, the measure button, and the sidebar's SELECTION list (its position; ticking selects in the workbench; clearing unticks).
  - The workbench suites fail on any console error; the Reflex suite fails on page errors and error toasts.
- **Rendered-pixel checks** (screenshot scripts): with the same camera, every ligand, ion and water atom centre renders identically in the two 3D views (0/255 difference over 131 atoms, none on background); a Zn sphere centre renders 0.92 × its colour, the shader's lighting at the centre.
- **Screenshot review** at the user's iframe size (about 1590 × 1000) after each round. Several of the fixes above (toolbar overflow, overlay colour clash, wrapping tables, legend overlap, hue confusion) came from this review rather than from tests.

---

## 6. Resources and stability

An exploratory benchmark (a KD-tree query with k = n on 4HHB, about 2×10⁹ neighbour entries) exhausted memory. The kernel killed it, the session ended and the background servers stopped. Measures taken since:

- **Monitoring during development:** a watch on `/proc/meminfo`, memory pressure (PSI) and the server ports, reporting only state changes (warn below 3 GB, critical below 1.5 GB). While it ran, available memory never fell below 6.0 GB (out of 15.6 GB), and no warning fired. In the later rounds (three servers and a headless browser), free memory, checked before each heavy run, stayed above about 4.5 GB; every browser run was capped at 2–3 GB.
- **Hard caps:** heavy runs and the servers run in cgroup scopes, e.g. `systemd-run --user --scope -p MemoryMax=… -p MemorySwapMax=0`, so the kernel stops only the offending process.
- **Measured costs:** `scripts/memprobe.py` runs one layer per process on 1AON and 7K00. Results: Voronoi and features about 8 KB per heavy point, Delaunay about 1.4 KB per point, parsing about 1.5 KB per atom; all linear.
- **An engine memory budget** (`resources.py`): each estimate is checked against the memory still available (`MemAvailable`, capped by any cgroup limit) minus a 768 MB reserve. A request that does not fit is refused with advice instead of crashing the server, and computations above 300 MB run one at a time.

---

## 7. Development environment

These notes are specific to the machine bioGeometry was built on.

| Item | Detail |
|---|---|
| Folder | `github-new/gemini/bioGeometry-workbench/` (renamed from `protein-geometry-dashboard-2/`) |
| Conda environment | `biopyc`: Python 3.14, reflex 0.9.5.post2, gemmi 0.7.5, numpy 2.5, scipy 1.17, Biopython 1.87 |
| Ports | 3100 Reflex frontend, 8100 Reflex backend + API, 8770 standalone server, 3000 the reference Dash dashboard (project B) |
| Playwright | from `compGeoD3/.venv`; Chromium uses `executable_path` to the newest `chromium-*` build in the cache with `--use-gl=angle --use-angle=swiftshader --enable-unsafe-swiftshader`; Firefox build installed into `~/.cache/ms-playwright`. A `~/.cache` cleanup can delete the browsers: reinstall with `playwright install firefox chromium` |
| Reflex hot reload | `rxconfig.py` excludes `data/`, `scripts/`, `tests/` and `docs/`; Python or asset edits recompile the app |
| Reference dashboard (project B) | `github-new/gemini/protein-geometry-dashboard/`, run from a uv virtual environment built from `dashboard/requirements.txt` in the session scratch folder (dash-mantine-components 0.12.1 is incompatible with `biopyc`'s 2.8); rebuild it after `/tmp` is cleaned |
| Screenshots | `docs/screenshots/guide/` (README and user guide) and `docs/screenshots/seaborn-palette/` (the palette decisions, 3.11–3.13), taken with Playwright's Firefox (`webgl.force-enabled`) |

---

## 8. Known limitations and next steps

- **SASA cost:** SASA is most of the feature-table time on large structures (7K00: about 7 s of 11 s). An adaptive sphere sampling for very large assemblies would be about 30% faster, at about 1% per-residue noise.
- **All-atom on huge assemblies:** this works within the memory budget but is slow (7K00 features about 44 s). *Backbone* (Cα / P) is the default for such examples.
- **Nucleic acids:** their geometry is computed as for proteins, but many residue features and labels are still protein-oriented (e.g. charge, DSSP secondary structure, Ramachandran basins), so for DNA and RNA they are not yet accurate. Nucleic-acid-specific analysis is work in progress (README and user guide carry the same note).
- **Mode label:** the sidebar says *Backbone*, but the workbench's geometry legend and status bar still show the internal level name `ca` (e.g. "ca points").
- **Standalone workbench:** configured only through the URL; it has no sources, parameters, analyses or downloads panel.
- **Binding rules are geometric:** distances only, with no hydrogens or angles. The metal cutoff is 3.0 Å for the dotted lines and 2.8 Å in selection geometry.
- **Selection geometry limits:** pairwise metrics cover the first 16 groups. A domain counts as a group when at least half of it is selected.
- **PDB export:** limited to 62 chains by the format; larger structures are exported as mmCIF only.
- **Domain annotations** depend on PDBe and UniProt being reachable on first load; the results are cached in `data/annotations/`.
- **Tests:** they need 4HHB once (downloaded into `data/structures/`), and the browser checks run on software WebGL.
- **Box-clipped Voronoi cells:** surface cells depend on the box margin and can be very large where the axis-aligned box has empty corners (up to about 3×10⁶ Å³ on the ribosome); use the *Voronoi interior* flag when comparing volumes.
- **Tetrahedra table size:** capped at 400,000 rows (sampled beyond that, with the count shown). Depth classes and binding contacts are geometric heuristics.
- **Possible additions:**
  - selection-geometry history to compare measurements
  - exporting the measurement as CSV
  - distance labels on binding lines
  - hydrogen placement for proper H-bond angles
