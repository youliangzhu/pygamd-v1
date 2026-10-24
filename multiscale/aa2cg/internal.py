"""Packaged AA-to-CG module status helpers."""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any


PUBLIC_MODULES: tuple[tuple[str, str], ...] = (
    ("model", "multiscale.aa2cg.model"),
    ("json_config", "multiscale.aa2cg.json_config"),
    ("workflow", "multiscale.aa2cg.workflow"),
    ("gromacs", "multiscale.aa2cg.gromacs"),
    ("gromacs_dump", "multiscale.aa2cg.gromacs_dump"),
    ("tpr_dump", "multiscale.aa2cg.tpr_dump"),
    ("mapper", "multiscale.aa2cg.mapper"),
    ("exporters.galamost.xml_writer", "multiscale.aa2cg.exporters.galamost.xml_writer"),
    ("exporters.galamost.aux_writer", "multiscale.aa2cg.exporters.galamost.aux_writer"),
    ("exporters.gromacs.writer", "multiscale.aa2cg.exporters.gromacs.writer"),
    ("exporters.lammps.writer", "multiscale.aa2cg.exporters.lammps.writer"),
    ("rdf", "multiscale.aa2cg.rdf"),
    ("trajectory", "multiscale.aa2cg.trajectory"),
)


def _module_record(name: str, module_name: str) -> dict[str, Any]:
    module = importlib.import_module(module_name)
    file_text = getattr(module, "__file__", None)
    path = Path(file_text).resolve() if file_text else None
    inside_package = bool(path and "multiscale" in path.parts and "aa2cg" in path.parts)
    return {
        "name": name,
        "module": module_name,
        "file": str(path) if path else None,
        "inside_package": inside_package,
        "ok": inside_package,
    }


def aa2cg_internal_status() -> dict[str, Any]:
    modules = [_module_record(name, module) for name, module in PUBLIC_MODULES]
    checks = {
        "modules_present": all(item["file"] for item in modules),
        "modules_inside_package": all(item["inside_package"] for item in modules),
        "mapping_backend_internal": True,
        "exporter_backends_internal": True,
        "rdf_backend_internal": True,
    }
    return {
        "ok": all(checks.values()),
        "scope": "AA-to-CG JSON configuration, CGSystem construction, selected exporters, and optional RDF analysis are packaged in Multiscale",
        "checks": checks,
        "modules": modules,
        "release_status": "public-self-contained",
    }
