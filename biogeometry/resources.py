"""Memory budget for heavy computations.

A request that would exhaust RAM is refused with a readable message instead of letting the kernel
OOM-kill the server, and large computations run one at a time (Starlette serves layer requests from a
thread pool, so two big structures could otherwise be tessellated at once).

Per-unit costs are measured peaks (``scripts/memprobe.py`` on 1AON and 7K00, heavy-atom level) with
about 25 % headroom; they scale linearly with the number of points or atoms.
"""
from __future__ import annotations

import threading
from contextlib import contextmanager
from pathlib import Path

RESERVE_MB = 768            # always left for the OS, the browser and the rest of the app
HEAVY_MB = 300              # estimates above this are serialised
PARSE_KB_PER_ATOM = 2.0     # gemmi model + numpy arrays + python lists during conversion
SASA_KB_PER_ATOM = 4.5      # Shrake-Rupley chunks + neighbour queries + interaction search
KB_PER_POINT = {            # tessellation-based layers
    "delaunay": 1.8, "alpha": 1.8, "contacts": 1.8, "pockets": 1.8, "filtration": 2.0, "hull": 0.5,
    "voronoi": 10.5, "scaling": 11.0, "tetrahedra": 12.0,     # the tetrahedra table builds the Voronoi cells too
}

_heavy = threading.Lock()


class MemoryBudgetError(ValueError):
    """User-facing: the computation would not fit in the memory that is currently free."""


CGROUP_ROOT = Path("/sys/fs/cgroup")


def _meminfo_mb() -> float | None:
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) / 1024
    except OSError:
        return None
    return None


def _cgroup_headroom_mb() -> float | None:
    """Room left under the tightest cgroup v2 memory.max above this process (systemd scope, container)."""
    try:
        rel = next(ln[3:].strip() for ln in Path("/proc/self/cgroup").read_text().splitlines() if ln.startswith("0::"))
    except (OSError, StopIteration):
        return None
    best = None
    d = CGROUP_ROOT / rel.lstrip("/")
    while True:
        try:
            limit = (d / "memory.max").read_text().strip()
            if limit != "max":
                room = (int(limit) - int((d / "memory.current").read_text())) / 2 ** 20
                best = room if best is None else min(best, room)
        except (OSError, ValueError):
            pass
        if d == CGROUP_ROOT or d.parent == d:
            return best
        d = d.parent


def available_mb() -> float | None:
    """Memory this process can still use: MemAvailable, capped by any cgroup limit.
    None where neither can be read (no check is made then)."""
    vals = [v for v in (_meminfo_mb(), _cgroup_headroom_mb()) if v is not None]
    return min(vals) if vals else None


def require(need_mb: float, what: str) -> None:
    avail = available_mb()
    if avail is None or need_mb + RESERVE_MB <= avail:
        return
    raise MemoryBudgetError(
        f"{what} needs about {need_mb / 1024:.1f} GB but only {avail / 1024:.1f} GB of memory is free; "
        f"use a coarser point level (ca / centroid), select fewer chains, or close other applications")


@contextmanager
def budget(need_mb: float, what: str):
    """Check the estimate, waiting for any other large computation to finish first."""
    if need_mb < HEAVY_MB:
        require(need_mb, what)
        yield
        return
    with _heavy:
        require(need_mb, what)
        yield


def layer_estimate_mb(kind: str, n_points: int, n_atoms: int, coarse: bool) -> float:
    """Peak memory (MB) to build one layer from scratch."""
    if kind == "features":
        # SASA runs before the Voronoi cells are built, so the peaks do not add up
        vor = (2.0 if coarse else KB_PER_POINT["voronoi"]) * n_points
        return max(SASA_KB_PER_ATOM * n_atoms, vor) / 1024
    if kind == "structure":
        return 1.0 * n_atoms / 1024          # bonds, cartoon mesh, buffers
    return KB_PER_POINT.get(kind, 0.5) * n_points / 1024


def selection_estimate_mb(n_atoms: int) -> float:
    """Selection geometry: tessellation + SASA of the selected atoms (alone, per group, per pair)."""
    return (KB_PER_POINT["delaunay"] + SASA_KB_PER_ATOM) * n_atoms / 1024


def parse_estimate_mb(n_atoms: int) -> float:
    return PARSE_KB_PER_ATOM * n_atoms / 1024
