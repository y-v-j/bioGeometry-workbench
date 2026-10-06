"""bioGeometry — interactive computational-geometry workbench for biological macromolecules.

One screen, no page scroll: a fixed-width control panel (sources, examples, geometry parameters,
analyses, display, downloads) and the workbench iframe (three.js + d3) filling the rest.
"""
import reflex as rx

from biogeometry.api import api

from .state import EXAMPLES, State

C = {"bg": "#191926", "card": "#20202f", "edge": "#2b2b40", "fg": "#e8e8f2", "label": "#9a9ab8", "dim": "#62627e",
     "sky": "#7dd3fc", "lavender": "#c4b5fd", "rose": "#f9a8d4", "mint": "#86efac", "teal": "#5eead4", "coral": "#fca5a5"}
FONT = '"FantasqueSansM Nerd Font", "Fantasque Sans Mono", "JetBrainsMono Nerd Font", monospace'

ANALYSIS_BUTTONS = [("pockets", "Pockets/cavities", "circle-dot"), ("contacts", "Chain interfaces", "network"),
                    ("filtration", "α-filtration", "chart-line"), ("fourbody", "Four-body stats", "triangle"),
                    ("scaling", "Scaling benchmark", "gauge")]
STYLES = [("cartoon", "cartoon"), ("trace", "trace"), ("spheres", "spheres (vdW)"), ("sticks", "sticks")]
COLORS = [("chain", "chain"), ("ss", "secondary structure"), ("class", "residue class"), ("conf", "B-factor / pLDDT"),
          ("hydropathy", "hydropathy"), ("rainbow", "N → C"), ("element", "element"), ("rsa", "burial (RSA)"),
          ("vor_volume", "Voronoi volume"), ("hull_depth", "hull depth"), ("del_degree", "Delaunay contacts")]


def label(text: str, **kw) -> rx.Component:
    return rx.text(text, size="1", color=C["label"], **kw)


def heading(text: str, color: str) -> rx.Component:
    return rx.hstack(rx.box(width="3px", height="12px", border_radius="2px", background=color, box_shadow=f"0 0 8px {color}"),
                     rx.text(text, size="1", weight="bold", letter_spacing="0.14em", color=C["label"]),
                     spacing="2", align="center")


def select(items, value, on_change, **kw) -> rx.Component:
    return rx.select.root(rx.select.trigger(width="100%"),
                          rx.select.content(*[rx.select.item(lbl, value=val) for val, lbl in items]),
                          value=value, on_change=on_change, size="1", **kw)


def field(text: str, control: rx.Component) -> rx.Component:
    return rx.grid(label(text), control, columns="92px 1fr", align="center", gap="2", width="100%")


def slider(text: str, name: str, var, lo, hi, step, handler) -> rx.Component:
    return rx.vstack(
        rx.hstack(label(text), rx.spacer(), rx.text(var, size="1", font_family=FONT, color=C["sky"]), width="100%"),
        rx.slider(default_value=[var], min=lo, max=hi, step=step, size="1", key=State.job_id + name,
                  on_value_commit=lambda v: handler(name, v)),
        spacing="1", width="100%")


def toggle(text: str, var, on_change) -> rx.Component:
    return rx.hstack(rx.switch(checked=var, on_change=on_change, size="1"), label(text), spacing="2", align="center")


def number(text: str, name: str, var, step: str = "0.5") -> rx.Component:
    return field(text, rx.input(value=var.to_string(), type="number", step=step, size="1",
                                on_change=lambda v: State.set_value(name, v)))


# ------------------------------------------------------------------ source block
def mode_picker() -> rx.Component:
    """The two kinds of analysis, chosen before loading: the polymer at all-atom or backbone level (Cα for
    protein, P for nucleic acids), with the ligands, ions and coordinated waters as atoms in both."""
    return rx.vstack(
        rx.segmented_control.root(
            rx.segmented_control.item("All-atom", value="allatom"), rx.segmented_control.item("Backbone", value="ca"),
            value=State.mode, on_change=State.set_mode, size="1", width="100%"),
        rx.text(rx.cond(State.mode == "ca", "Backbone residue(s)", "every heavy atom"),
                " + ligands, ions & coordinated waters · applies on Load", size="1", color=C["dim"]),
        spacing="1", width="100%")


def source_block() -> rx.Component:
    sv = State.set_text
    rcsb = rx.grid(
        rx.input(placeholder="PDB id (4HHB)", value=State.pdb_id, on_change=lambda v: sv("pdb_id", v), size="1",
                 on_key_down=lambda k: rx.cond(k == "Enter", State.load_structure, rx.noop())),
        rx.input(placeholder="chains: all", value=State.chains, on_change=lambda v: sv("chains", v), size="1",
                 title="Comma-separated author chain ids; blank = whole entry"),
        columns="1fr 96px", gap="2", width="100%")
    alphafold = rx.vstack(
        rx.input(placeholder="UniProt accession (P69905) or AF-P69905-F1", value=State.uniprot,
                 on_change=lambda v: sv("uniprot", v), size="1",
                 on_key_down=lambda k: rx.cond(k == "Enter", State.load_structure, rx.noop())),
        rx.text("AlphaFold DB model · pLDDT in the B column · PAE in the right panel", size="1", color=C["dim"]),
        spacing="1", width="100%")
    local = rx.vstack(
        rx.hstack(rx.input(placeholder="/path/to/dir or file.cif", value=State.local_path,
                           on_change=lambda v: sv("local_path", v), size="1", flex="1"),
                  rx.button("List", size="1", variant="soft", on_click=State.list_local), width="100%", spacing="2"),
        rx.cond(State.local_files.length() > 0,
                rx.select.root(rx.select.trigger(placeholder="pick a file…", width="100%"),
                               rx.select.content(rx.foreach(State.local_file_names, lambda n: rx.select.item(n, value=n))),
                               on_change=State.choose_local, size="1")),
        rx.input(placeholder="chains: all", value=State.chains, on_change=lambda v: sv("chains", v), size="1"),
        spacing="2", width="100%")
    upload = rx.vstack(
        rx.upload(rx.vstack(rx.icon("upload", size=15), rx.text("Drop .pdb / .cif (.gz) or click", size="1"),
                            align="center", spacing="1"),
                  id="structure_upload", multiple=False, max_files=1,
                  on_drop=State.handle_upload(rx.upload_files(upload_id="structure_upload")),
                  border=f"1px dashed {C['edge']}", border_radius="6px", padding="8px", width="100%", cursor="pointer"),
        rx.cond(State.uploaded_name != "", rx.text("file: ", State.uploaded_name, size="1", font_family=FONT)),
        rx.input(placeholder="chains: all", value=State.chains, on_change=lambda v: sv("chains", v), size="1"),
        spacing="1", width="100%")
    return rx.vstack(
        rx.segmented_control.root(
            rx.segmented_control.item("PDB", value="rcsb"), rx.segmented_control.item("AlphaFold", value="alphafold"),
            rx.segmented_control.item("Upload", value="upload"), rx.segmented_control.item("Local", value="local"),
            value=State.source, on_change=State.set_source, size="1", width="100%"),
        rx.match(State.source, ("rcsb", rcsb), ("alphafold", alphafold), ("local", local), upload),
        mode_picker(),
        rx.button(rx.cond(State.busy, rx.spinner(size="1"), rx.icon("download", size=14)),
                  rx.cond(State.source == "local", "Load selected file", "Load structure"),
                  on_click=State.load_structure, size="1", width="100%", disabled=State.busy),
        rx.flex(*[rx.badge(ident, title=desc, variant="surface", color_scheme="gray", size="1", cursor="pointer",
                           on_click=State.load_example(src, ident, mode),
                           style={"font_family": FONT, "border": f"1px solid {C['edge']}"})
                  for src, ident, desc, mode in EXAMPLES],
                wrap="wrap", gap="1", width="100%"),
        spacing="2", width="100%")


def status() -> rx.Component:
    return rx.vstack(
        rx.cond(State.loaded, rx.vstack(
            rx.text(State.title_line, size="1", weight="medium", color=C["fg"],
                    style={"display": "-webkit-box", "-webkit-line-clamp": "2", "-webkit-box-orient": "vertical",
                           "overflow": "hidden"}),
            rx.text(State.meta_line, size="1", color=C["lavender"]),
            rx.text(State.summary_line, size="1", color=C["dim"], font_family=FONT),
            spacing="0", width="100%")),
        rx.cond(State.busy, rx.hstack(rx.spinner(size="1"), rx.text(State.busy_msg, size="1", color=C["teal"]), spacing="2")),
        rx.cond(State.error != "", rx.callout(State.error, icon="triangle-alert", color_scheme="red", size="1", width="100%")),
        rx.cond((State.info != "") & (State.error == ""), rx.text(State.info, size="1", color=C["mint"])),
        spacing="1", width="100%")


# ------------------------------------------------------------------ accordion sections
def geometry_params() -> rx.Component:
    sp = State.set_param
    return rx.vstack(
        rx.text("Toggle hull / Delaunay / α-shape / Voronoi in the geometry panel; these parameters apply live.",
                size="1", color=C["dim"]),
        heading("DELAUNAY", C["sky"]),
        slider("max edge (Å)", "d_max_edge", State.d_max_edge, 3, 30, 0.5, sp),
        field("edges", select([("all", "all"), ("intra", "intra-chain"), ("inter", "inter-chain (interface)")],
                              State.d_filter, lambda v: sp("d_filter", v))),
        toggle("show tetrahedra (shrunken)", State.d_tets, lambda v: sp("d_tets", v)),
        heading("VORONOI", C["lavender"]),
        slider("max cell reach (Å)", "v_radius", State.v_radius, 3, 25, 0.5, sp),
        slider("box margin (Å)", "v_margin", State.v_margin, 1, 40, 0.5, sp),
        slider("shrink", "v_shrink", State.v_shrink, 0.3, 1.0, 0.05, sp),
        field("cells", select([("bounded", "interior cells"), ("all", "all cells (box-closed)"), ("buried", "buried only"),
                               ("surface", "first shell")],
                              State.v_subset, lambda v: sp("v_subset", v))),
        field("style", select([("cells", "solid cells"), ("skeleton", "ridge skeleton")], State.v_style,
                              lambda v: sp("v_style", v))),
        heading("α-SHAPE", C["mint"]),
        slider("α radius (Å)", "a_alpha", State.a_alpha, 1, 40, 0.25, sp),
        heading("CONVEX HULL", C["teal"]),
        field("scope", select([("global", "whole selection"), ("chains", "global + per chain")], State.h_scope,
                              lambda v: sp("h_scope", v))),
        toggle("convex layers (onion peeling)", State.h_layers, lambda v: sp("h_layers", v)),
        spacing="2", width="100%")


def analysis_panel() -> rx.Component:
    def cell(icon, text, **kw):     # equal-width cells, icon + label left-aligned on one line
        return rx.button(rx.icon(icon, size=13, flex_shrink="0"),
                         rx.text(text, size="1", style={"white_space": "nowrap", "overflow": "hidden", "text_overflow": "ellipsis"}),
                         size="1", width="100%", min_width="0", justify_content="flex-start", gap="6px", **kw)
    buttons = rx.grid(
        *[cell(icon, text, variant=rx.cond(State.analysis_done.contains(kind), "solid", "soft"),
               on_click=State.run_analysis(kind), disabled=~State.loaded | State.busy)
          for kind, text, icon in ANALYSIS_BUTTONS],
        rx.tooltip(cell("box", "Measure selection", variant="soft", color_scheme="pink", on_click=State.measure_selection,
                        disabled=~State.loaded),
                   content="Select domains, ligands or ions in the workbench, then measure them: hulls, α-shape, contacts, buried surface"),
        columns="repeat(2, minmax(0, 1fr))", gap="2", width="100%")
    return rx.vstack(
        buttons,
        rx.text("Not every structure has pockets, cavities or channels: small, flat or fibrous chains often have none "
                "above the minimum volume - an empty result is an answer, not an error.", size="1", color=C["dim"]),
        number("pocket α small", "p_small", State.p_small),
        number("pocket α large", "p_large", State.p_large),
        number("min pocket Å³", "p_minvol", State.p_minvol, "50"),
        number("contact cutoff Å", "c_cutoff", State.c_cutoff),
        number("filtration α min", "f_min", State.f_min),
        number("filtration α max", "f_max", State.f_max),
        rx.text("Results open in the workbench dock; clicking a pocket or interface selects its residues everywhere.",
                size="1", color=C["dim"]),
        spacing="2", width="100%")


def display_panel() -> rx.Component:
    sd = State.set_display
    return rx.vstack(
        field("structure", select(STYLES, State.mol_style, lambda v: sd("mol_style", v))),
        field("colour by", select(COLORS, State.mol_color, lambda v: sd("mol_color", v))),
        slider("geometry opacity", "geo_opacity", State.geo_opacity, 0.05, 1.0, 0.05, sd),
        rx.text("Selection — click: one · ctrl: toggle · shift: range/add · drag in the sequence or brush the plots · "
                "Esc clears · F focuses · L links cameras · R resets · S spin · C slab", size="1", color=C["dim"]),
        spacing="2", width="100%")


def load_options() -> rx.Component:
    sv = State.set_value
    return rx.vstack(
        field("assembly", rx.input(value=State.assembly, placeholder="asym. unit (or 1, 2…)", size="1",
                                   on_change=lambda v: sv("assembly", v))),
        field("waters", select([("coordinated", "coordinated (ion / ligand)"), ("none", "none")], State.waters,
                               lambda v: sv("waters", v))),
        rx.text("Changes apply on the next Load.", size="1", color=C["dim"]),
        spacing="2", width="100%")


def downloads() -> rx.Component:
    def link(r, icon, color):         # fetched in the page and saved: no new tab, no flicker
        return rx.link(rx.hstack(rx.icon(icon, size=13), rx.text(r[0], size="1", weight="medium"), spacing="2", align="center"),
                       href="#", on_click=[rx.prevent_default, State.download_file(r[2])], title=r[1], color=color)
    return rx.cond(
        State.loaded,
        rx.vstack(label("structure file"),
                  rx.foreach(State.structure_links, lambda r: link(r, "box", C["rose"])),
                  label("feature tables (CSV)", margin_top="4px"),
                  rx.foreach(State.csv_links, lambda r: link(r, "file-down", C["sky"])),
                  rx.text("Any table view can also be saved as shown (⤓ view CSV in the dock).", size="1", color=C["dim"]),
                  spacing="1", width="100%"),
        rx.text("Load a structure first.", size="1", color=C["dim"]))


def selection_panel() -> rx.Component:
    """Domains, ligands, ions and coordinated waters (tick to add to the selection). The workbench renders the
    list into #bg-selection-panel itself (sequence.js; same origin), so its logic stays in one place."""
    return rx.vstack(
        heading("SELECTION", C["coral"]),
        rx.box(rx.el.div(id="bg-selection-panel"), id="bg-selection-scroll", width="100%"),
        spacing="2", width="100%", padding_bottom="6px", border_bottom=f"1px solid {C['edge']}")


def accordion() -> rx.Component:
    def item(value, title, color, body):
        return rx.accordion.item(header=heading(title, color), content=body, value=value)
    return rx.accordion.root(
        item("params", "GEOMETRY PARAMETERS", C["mint"], geometry_params()),
        item("analysis", "ANALYSES", C["rose"], analysis_panel()),
        item("display", "DISPLAY", C["lavender"], display_panel()),
        item("options", "LOAD OPTIONS", C["teal"], load_options()),
        item("downloads", "DOWNLOADS", C["sky"], downloads()),
        type="single", collapsible=True, default_value="params", variant="ghost", width="100%")


def logo(size: str = "28") -> rx.Component:
    """A Delaunay tetrahedron: two lit faces, the hidden edge dashed, atoms at the vertices."""
    return rx.el.svg(
        rx.el.polygon(points="12,2.5 3,17 14.5,21.5", fill=C["sky"], fill_opacity="0.2"),
        rx.el.polygon(points="12,2.5 14.5,21.5 21,14.5", fill=C["rose"], fill_opacity="0.2"),
        rx.el.line(x1="3", y1="17", x2="21", y2="14.5", stroke=C["label"], stroke_width="1", stroke_dasharray="1.6 1.4"),
        rx.el.path(d="M12 2.5 L3 17 L14.5 21.5 L21 14.5 Z M12 2.5 L14.5 21.5", fill="none", stroke=C["lavender"],
                   stroke_width="1.5", stroke_linejoin="round"),
        rx.el.circle(cx="12", cy="2.5", r="1.9", fill=C["rose"]), rx.el.circle(cx="3", cy="17", r="1.9", fill=C["sky"]),
        rx.el.circle(cx="14.5", cy="21.5", r="1.9", fill=C["lavender"]), rx.el.circle(cx="21", cy="14.5", r="1.9", fill=C["teal"]),
        view_box="0 0 24 24", width=size, height=size, style={"flex_shrink": "0", "overflow": "visible"})


def brand() -> rx.Component:
    return rx.hstack(
        logo(),
        rx.vstack(rx.hstack(rx.text("bio", size="5", color=C["label"]), rx.text("Geometry", size="5", weight="bold", color=C["rose"]),
                            spacing="0"),
                  rx.text("computational geometry workbench for biomolecules", size="1", color=C["dim"], letter_spacing="0.04em"),
                  spacing="0"),
        rx.spacer(),
        rx.tooltip(rx.icon_button(rx.icon("rotate-ccw", size=14), size="1", variant="ghost", on_click=State.reset_workbench),
                   content="Reset the workbench: clear the selection, analyses and overlays, and show the loaded structure afresh"),
        width="100%", align="center", spacing="2")


def sidebar() -> rx.Component:
    return rx.vstack(
        brand(),
        rx.box(height="2px", width="100%", border_radius="1px",
               background=f"linear-gradient(90deg, {C['sky']}, {C['lavender']}, {C['rose']})"),
        source_block(),
        status(),
        selection_panel(),
        rx.scroll_area(accordion(), type="auto", scrollbars="vertical", flex="1", min_height="0"),
        width="330px", min_width="330px", height="100vh", padding="10px 12px", spacing="2",
        border_right=f"1px solid {C['edge']}", background=f"linear-gradient(180deg, #232336, {C['card']} 140px)",
        overflow="hidden")


def index() -> rx.Component:
    return rx.flex(
        sidebar(),
        rx.box(rx.el.iframe(id="bg-workbench", src="/workbench/index.html", title="bioGeometry workbench",
                            style={"border": "none", "width": "100%", "height": "100%", "display": "block"}),
               flex="1", min_width="0", height="100vh"),
        width="100vw", height="100vh", overflow="hidden",
        style={"background": C["bg"], "font_family": FONT},
    )


app = rx.App(api_transformer=api, stylesheets=["/biogeometry.css"], style={"font_family": FONT, "background": C["bg"]},
             head_components=[rx.el.link(rel="icon", type="image/svg+xml", href="/favicon.svg")])
app.add_page(index, title="bioGeometry · computational geometry workbench")
