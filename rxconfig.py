import os
from pathlib import Path

import reflex as rx

# Runtime data (downloads, caches, layer buffers) lives in data/; keep it out of the dev hot-reload so
# loading a structure never restarts the backend mid-request. assets/ stays watched: Reflex copies it
# into .web/public only when it recompiles, so a workbench edit must trigger a reload to be served.
_root = Path(__file__).resolve().parent
os.environ.setdefault(
    "REFLEX_HOT_RELOAD_EXCLUDE_PATHS",
    ":".join(str(_root / d) for d in ("data", "scripts", "tests", "docs") if (_root / d).exists()),
)

config = rx.Config(
    app_name="biogeometry_app",
    frontend_port=3100,
    backend_port=8100,
    api_url="http://localhost:8100",       # the workbench iframe fetches layers from the backend (/bg/...)
    plugins=[rx.plugins.RadixThemesPlugin(theme=rx.theme(appearance="dark", accent_color="cyan",
                                                         gray_color="slate", radius="medium", scaling="95%"))],
    disable_plugins=[rx.plugins.SitemapPlugin],
    telemetry_enabled=False,
    show_built_with_reflex=False,
)
