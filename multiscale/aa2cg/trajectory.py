"""Internal aa2cg trajectory helpers for RDF analysis."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from collections.abc import Iterator

from multiscale.aa2cg.rdf import CGFrame
from multiscale.aa2cg.gromacs import parse_gro


@dataclass(frozen=True)
class TrajectoryInfo:
    topology_source: str
    n_atoms: int
    total_frames: int | None
    selected_frames: int | None


def _import_mdanalysis():
    try:
        import MDAnalysis as mda  # type: ignore[import-not-found]
    except ImportError as exc:
        raise RuntimeError("Please install MDAnalysis to use `multiscale aa2cg rdf`.") from exc
    except Exception as exc:
        raise RuntimeError(
            "Failed to import MDAnalysis. This usually means the MDAnalysis/NumPy/SciPy "
            f"environment is inconsistent: {exc}"
        ) from exc
    return mda


def _open_universe(*, mda, topology_tpr: str | Path, coordinates: str | Path, trajectory: str | Path):
    try:
        return mda.Universe(str(topology_tpr), str(trajectory)), "topology_tpr"
    except Exception:
        return mda.Universe(str(coordinates), str(trajectory)), "coordinates"


def _selected_frame_count(
    total_frames: int | None,
    *,
    start: int = 0,
    stop: int | None = None,
    stride: int = 1,
    max_frames: int | None = None,
) -> int | None:
    if total_frames is None:
        return max_frames
    effective_stop = total_frames if stop is None else min(stop, total_frames)
    if start >= effective_stop:
        selected = 0
    else:
        selected = len(range(start, effective_stop, stride))
    if max_frames is not None:
        selected = min(selected, max_frames)
    return selected


def summarize_mdanalysis_trajectory(
    *,
    topology_tpr: str | Path,
    coordinates: str | Path,
    trajectory: str | Path,
    start: int = 0,
    stop: int | None = None,
    stride: int = 1,
    max_frames: int | None = None,
) -> TrajectoryInfo:
    mda = _import_mdanalysis()
    universe, topology_source = _open_universe(
        mda=mda,
        topology_tpr=topology_tpr,
        coordinates=coordinates,
        trajectory=trajectory,
    )
    try:
        total_frames = len(universe.trajectory)
    except Exception:
        total_frames = None
    selected_frames = _selected_frame_count(
        total_frames,
        start=start,
        stop=stop,
        stride=stride,
        max_frames=max_frames,
    )
    return TrajectoryInfo(
        topology_source=topology_source,
        n_atoms=int(universe.atoms.n_atoms),
        total_frames=total_frames,
        selected_frames=selected_frames,
    )


def iter_mdanalysis_aa_frames(
    *,
    topology_tpr: str | Path,
    coordinates: str | Path,
    trajectory: str | Path,
    start: int = 0,
    stop: int | None = None,
    stride: int = 1,
    max_frames: int | None = None,
) -> Iterator[tuple[list[tuple[float, float, float]], tuple[float, float, float]]]:
    mda = _import_mdanalysis()
    universe, _topology_source = _open_universe(
        mda=mda,
        topology_tpr=topology_tpr,
        coordinates=coordinates,
        trajectory=trajectory,
    )

    gro_box = parse_gro(coordinates).box
    emitted = 0
    for ts in universe.trajectory[start:stop:stride]:
        positions_nm = [(float(x) * 0.1, float(y) * 0.1, float(z) * 0.1) for x, y, z in ts.positions]
        dims = getattr(ts, "dimensions", None)
        if dims is not None and len(dims) >= 3 and all(float(v) > 0 for v in dims[:3]):
            if len(dims) >= 6 and any(abs(float(angle) - 90.0) > 1e-6 for angle in dims[3:6]):
                raise ValueError("multiscale aa2cg rdf currently supports only orthorhombic trajectory boxes")
            box_nm = (float(dims[0]) * 0.1, float(dims[1]) * 0.1, float(dims[2]) * 0.1)
        elif gro_box:
            box_nm = gro_box
        else:
            raise ValueError("RDF requires a box from the trajectory frame or fallback GRO coordinates")
        yield positions_nm, box_nm
        emitted += 1
        if max_frames is not None and emitted >= max_frames:
            break


def iter_cg_frames_from_mdanalysis(*, template, **kwargs) -> Iterator[CGFrame]:
    from multiscale.aa2cg.rdf import map_aa_frame_to_cg

    for positions_nm, box_nm in iter_mdanalysis_aa_frames(**kwargs):
        yield map_aa_frame_to_cg(positions_nm, box_nm, template)
