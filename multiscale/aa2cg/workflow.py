"""Execution of the unified AA-to-CG JSON workflow."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Callable

from multiscale.core.cg_system import CGSystem, load_cgsystem_json

from .exporters.galamost import export_galamost_auxiliary, export_galamost_topology
from .exporters.gromacs import export_gromacs_scaffold
from .exporters.lammps import export_lammps_scaffold
from .gromacs_dump import run_gmx_dump
from .json_config import AA2CGConfig, load_aa2cg_config
from .mapper import StructureMapReport, map_structure_to_cgsystem, write_mapped_structure
from .rdf import (
    build_cg_mapping_template,
    compute_rdf_histograms,
    make_rdf_spec_from_topology,
    resolve_rdf_mode_from_spec,
    write_rdf_outputs,
)
from .trajectory import iter_cg_frames_from_mdanalysis, summarize_mdanalysis_trajectory
from .tpr_dump import TprTopologyInput, parse_mdanalysis_tpr, render_tpr_dump


def _write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _load_mapping_topology(
    config: AA2CGConfig,
    *,
    gmx_bin: str,
    log: Callable[[str], None] | None = None,
) -> tuple[TprTopologyInput, dict[str, Any], Path | None]:
    """Return topology metadata and a temporary dump path to clean, if any.

    GROMACS remains the first choice when it is available, preserving the
    established output byte-for-byte.  Missing or failing ``gmx dump`` falls
    back to MDAnalysis' binary TPR parser, which is available on Windows.
    """

    keep_path = (
        config.output_root / "intermediate" / f"{config.topology_tpr.name}.dmp"
        if config.keep_topology_dump
        else None
    )
    try:
        generated_dump = run_gmx_dump(config.topology_tpr, gmx_bin=gmx_bin, keep_path=keep_path)
    except (OSError, RuntimeError) as gmx_error:
        if log:
            log(f"[aa2cg] gmx dump unavailable; reading TPR with MDAnalysis ({type(gmx_error).__name__})")
        try:
            topology = parse_mdanalysis_tpr(config.topology_tpr)
        except Exception as mda_error:
            raise RuntimeError(
                "Unable to read AA topology with either gmx dump or MDAnalysis.\n"
                f"  gmx error: {gmx_error}\n"
                f"  MDAnalysis error: {mda_error}"
            ) from mda_error

        normalized_dump: Path | None = None
        if keep_path is not None:
            keep_path.parent.mkdir(parents=True, exist_ok=True)
            keep_path.write_text(render_tpr_dump(topology), encoding="utf-8")
            normalized_dump = keep_path
        topology_source = {
            "mode": "tpr",
            "tpr_path": str(config.topology_tpr),
            "reader": "mdanalysis",
            "gmx_bin": gmx_bin,
            "gmx_attempted": True,
            "gmx_error_type": type(gmx_error).__name__,
            "gmx_error": str(gmx_error),
            "dump_generated": False,
            "dump_kept": normalized_dump is not None,
            "dump_path": str(normalized_dump) if normalized_dump is not None else None,
            "normalized_dump": normalized_dump is not None,
        }
        return topology, topology_source, None

    topology_source = {
        "mode": "tpr",
        "tpr_path": str(config.topology_tpr),
        "reader": "gmx_dump",
        "gmx_bin": gmx_bin,
        "gmx_attempted": True,
        "dump_generated": True,
        "dump_kept": bool(config.keep_topology_dump),
        "dump_path": str(generated_dump) if config.keep_topology_dump else None,
        "normalized_dump": False,
    }
    temporary_dump = generated_dump if keep_path is None else None
    return generated_dump, topology_source, temporary_dump


def _export_galamost(config: AA2CGConfig, system: CGSystem) -> dict[str, Any]:
    if not config.galamost.enabled:
        return {"enabled": False, "status": "skipped"}
    topology = export_galamost_topology(system, config.galamost.directory, eps=config.galamost.epsilon)
    auxiliary = export_galamost_auxiliary(system, config.galamost.directory)
    topology_report = topology.to_dict()
    _write_json(config.galamost.directory / "pygamd_export_report.json", topology_report)
    _write_json(config.galamost.directory / "pygamd_aux_report.json", auxiliary.to_dict())
    inbox_report = dict(topology_report.pop("topology_inbox"))
    topology_report["ok"] = topology.topology_ok
    return {
        "enabled": True,
        "ok": bool(topology.ok and auxiliary.ok),
        "directory": str(config.galamost.directory),
        "topology": topology_report,
        "topology_inbox": inbox_report,
        "auxiliary": auxiliary.to_dict(),
    }


def _export_gromacs(config: AA2CGConfig, system: CGSystem) -> dict[str, Any]:
    if not config.gromacs.enabled:
        return {"enabled": False, "status": "skipped"}
    report = export_gromacs_scaffold(system, config.gromacs.directory)
    _write_json(config.gromacs.directory / "gromacs_export_report.json", report.to_dict())
    return {
        "enabled": True,
        "ok": report.ok,
        "directory": str(config.gromacs.directory),
        "report": report.to_dict(),
    }


def _export_lammps(config: AA2CGConfig, system: CGSystem) -> dict[str, Any]:
    if not config.lammps.enabled:
        return {"enabled": False, "status": "skipped"}
    report = export_lammps_scaffold(system, config.lammps.directory)
    _write_json(config.lammps.directory / "lammps_export_report.json", report.to_dict())
    return {
        "enabled": True,
        "ok": report.ok,
        "directory": str(config.lammps.directory),
        "report": report.to_dict(),
    }


def _mapping_payload(
    report: StructureMapReport,
    *,
    kept: bool,
    directory: Path | None,
    topology_source: dict[str, Any],
    config_path: Path,
) -> dict[str, Any]:
    """Return an audit payload without exposing temporary file paths."""

    if kept:
        payload = report.to_dict()
        payload["artifacts"] = {
            "kept": True,
            "directory": str(directory),
            "files": dict(report.files),
        }
        return payload

    return {
        "artifacts": {
            "kept": False,
            "directory": None,
            "files": {},
        },
        "counts": report.counts,
        "checks": report.checks,
        "warnings": report.warnings,
        "topology_source": topology_source,
        "entrypoint": "run",
        "config": str(config_path),
        "backend": "internal",
        "ok": report.ok,
    }


def _run_rdf(
    config: AA2CGConfig,
    *,
    topology_dump: TprTopologyInput,
    log: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    rdf = config.rdf
    if not rdf.enabled:
        return {
            "enabled": False,
            "status": "skipped",
            "reason": "rdf is false",
            "trajectory_required": False,
        }
    if rdf.trajectory is None:
        raise ValueError("RDF is enabled but no trajectory was configured")

    spec = make_rdf_spec_from_topology(
        config.topology,
        pairs=rdf.pairs,
        rmax_nm=rdf.rmax_nm,
        dr_nm=rdf.bin_width_nm,
        exclude_bonded_pairs=rdf.exclude_bonded_pairs,
        exclude_angle_endpoints=rdf.exclude_angle_endpoints,
        source=str(config.path),
    )
    if rdf.mode == "auto":
        effective_mode, mode_source, warnings = resolve_rdf_mode_from_spec(spec)
    else:
        effective_mode, mode_source, warnings = rdf.mode, "json_config", []

    if log:
        log(f"[aa2cg] RDF trajectory: {rdf.trajectory}")
        log(
            f"[aa2cg] RDF pairs={len(spec.pairs)}, rmax={spec.rmax_nm:.6f} nm, "
            f"bin_width={spec.dr_nm:.6f} nm, mode={effective_mode}"
        )

    template = build_cg_mapping_template(
        topology_dump=topology_dump,
        coordinates=config.coordinates_gro,
        mapping=config.mapping,
        cg_topology=config.topology,
        rdf_spec=spec,
    )
    trajectory_info = summarize_mdanalysis_trajectory(
        topology_tpr=config.topology_tpr,
        coordinates=config.coordinates_gro,
        trajectory=rdf.trajectory,
        start=rdf.start,
        stop=rdf.stop,
        stride=rdf.stride,
        max_frames=rdf.max_frames,
    )
    frames = iter_cg_frames_from_mdanalysis(
        template=template,
        topology_tpr=config.topology_tpr,
        coordinates=config.coordinates_gro,
        trajectory=rdf.trajectory,
        start=rdf.start,
        stop=rdf.stop,
        stride=rdf.stride,
        max_frames=rdf.max_frames,
    )
    rdf_data = compute_rdf_histograms(
        frames,
        bead_types=template.bead_types,
        molecule_ids=template.bead_molecule_ids,
        bead_names=template.bead_names,
        spec=spec,
        mode=effective_mode,
        rdf_engine=rdf.histogram_engine,
        progress_every=rdf.progress_every,
        total_frames=trajectory_info.selected_frames,
        progress_pairs=rdf.progress_pairs,
        log=log,
    )
    inputs = {
        "config": str(config.path),
        "topology_tpr": str(config.topology_tpr),
        "coordinates_gro": str(config.coordinates_gro),
        "trajectory": str(rdf.trajectory),
        "gmx_bin": config.gmx_bin,
        "frame_selection": {
            "start": rdf.start,
            "stop": rdf.stop,
            "stride": rdf.stride,
            "max_frames": rdf.max_frames,
        },
        "rdf_mode_requested": rdf.mode,
        "rdf_mode_effective": effective_mode,
        "rdf_mode_source": mode_source,
        "rdf_engine": rdf.histogram_engine,
        "warnings": warnings,
        "trajectory_summary": {
            "topology_source": trajectory_info.topology_source,
            "n_atoms": trajectory_info.n_atoms,
            "total_frames": trajectory_info.total_frames,
            "selected_frames": trajectory_info.selected_frames,
        },
    }
    result = write_rdf_outputs(
        outdir=rdf.directory,
        spec=spec,
        rdf_data=rdf_data,
        mode=effective_mode,
        source_trajectory=rdf.trajectory,
        inputs=inputs,
    )
    payload = dict(result.report)
    payload.update({"enabled": True, "ok": True, "directory": str(rdf.directory), "backend": "internal"})
    return payload


def run_aa2cg_config(
    path: str | Path,
    *,
    gmx_bin: str | None = None,
    log: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Run mapping and selected exporters from one ``multiscale.aa2cg.v1`` JSON file.

    Mapping diagnostics are intermediate artifacts and are not retained by
    default. Set ``output.mapping_artifacts`` to ``true`` (or an object with a
    custom directory) when ``cg_system_mapped.json``, ``cg.xmol``, the mapping
    CSV, and the VMD helper script are needed for debugging or auditing.
    """

    config = load_aa2cg_config(path, require_files=True)
    selected_gmx = gmx_bin or config.gmx_bin
    config.output_root.mkdir(parents=True, exist_ok=True)
    if config.keep_mapping_artifacts:
        config.mapped_directory.mkdir(parents=True, exist_ok=True)

    resolved_config_path = config.output_root / "aa2cg_resolved.json"
    _write_json(resolved_config_path, config.to_resolved_dict(gmx_bin=selected_gmx))
    if log:
        log(f"[aa2cg] config: {config.path}")
        log(f"[aa2cg] output: {config.output_root}")
        log("[aa2cg] mapping AA coordinates to CGSystem")

    topology_input, topology_source, temporary_dump = _load_mapping_topology(
        config,
        gmx_bin=selected_gmx,
        log=log,
    )

    cgsystem_path: Path | None = None
    try:
        if config.keep_mapping_artifacts:
            map_report = write_mapped_structure(
                topology_dump=topology_input,
                coordinates=config.coordinates_gro,
                mapping=config.mapping,
                cg_topology=config.topology,
                outdir=config.mapped_directory,
                topology_source=topology_source,
                entrypoint="run",
            )
            cgsystem_path = Path(map_report.files["cgsystem"])
            system = load_cgsystem_json(cgsystem_path)
            if log:
                log(f"[aa2cg] mapping artifacts retained: {config.mapped_directory}")
        else:
            system, map_report = map_structure_to_cgsystem(
                topology_dump=topology_input,
                coordinates=config.coordinates_gro,
                mapping=config.mapping,
                cg_topology=config.topology,
            )
            if log:
                log("[aa2cg] mapping artifacts are intermediate and were not retained")

        exports = {
            "galamost": _export_galamost(config, system),
            "gromacs": _export_gromacs(config, system),
            "lammps": _export_lammps(config, system),
        }
        rdf_report = _run_rdf(config, topology_dump=topology_input, log=log)
    finally:
        if temporary_dump is not None:
            try:
                os.unlink(temporary_dump)
            except FileNotFoundError:
                pass

    step_ok = [map_report.ok]
    for item in exports.values():
        if item.get("enabled"):
            step_ok.append(bool(item.get("ok")))
    if rdf_report.get("enabled"):
        step_ok.append(bool(rdf_report.get("ok")))

    report_path = config.output_root / "aa2cg_run_report.json"
    mapping_payload = _mapping_payload(
        map_report,
        kept=config.keep_mapping_artifacts,
        directory=config.mapped_directory if config.keep_mapping_artifacts else None,
        topology_source=topology_source,
        config_path=config.path,
    )
    report = {
        "schema": "multiscale.aa2cg.run.v1",
        "ok": all(step_ok),
        "name": config.name,
        "config": str(config.path),
        "output_directory": str(config.output_root),
        "input_contract": {
            "coordinates_used_as_is": True,
            "periodic_coordinate_preprocessing_performed": False,
            "user_responsible_for_coordinate_suitability": True,
        },
        "mapping": mapping_payload,
        "exports": exports,
        "rdf": rdf_report,
        "files": {
            "run_report": str(report_path),
            "resolved_config": str(resolved_config_path),
            "cgsystem": str(cgsystem_path) if cgsystem_path is not None else None,
            "galamost_topology": (
                str(config.galamost.directory / "topology.xml") if config.galamost.enabled else None
            ),
            "galamost_topology_inbox": (
                str(config.galamost.directory / "topology_inbox.xml") if config.galamost.enabled else None
            ),
        },
    }
    _write_json(report_path, report)
    return report
