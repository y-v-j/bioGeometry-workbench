"""Peak memory and time of one layer (plus its prerequisites) on a fresh job.

    python scripts/memprobe.py 1AON heavy voronoi       # KIND = open | structure | delaunay | ... | features

Used to calibrate biogeometry/resources.py. Peak RSS is per process, so run one layer per process,
ideally under a memory cap:

    systemd-run --user --scope -p MemoryMax=6G python scripts/memprobe.py 7K00 heavy features
"""
import resource
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from biogeometry import engine  # noqa: E402


def main(pdb_id, level, kind):
    engine.CACHE_DIR = Path(tempfile.mkdtemp())          # never reuse cached layers
    t = time.perf_counter()
    job = engine.open_job(engine.JobSpec(source="rcsb", pdb_id=pdb_id, level=level))
    t_open = time.perf_counter() - t
    base = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    t = time.perf_counter()
    if kind != "open":
        engine.build_layer(job, kind)
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    print(f"{pdb_id} {job.points.level:7s} points={len(job.points.coords):7d} atoms={job.structure.n_atoms:7d} "
          f"open={t_open:5.1f}s {kind:11s} {time.perf_counter() - t:5.1f}s "
          f"open_peak={base:6.0f}MB layer_peak={peak:6.0f}MB delta={peak - base:6.0f}MB")


if __name__ == "__main__":
    main(*sys.argv[1:4])
