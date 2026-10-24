"""Packaged CG-to-AA module status helpers."""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

PUBLIC_MODULES: tuple[str, ...] = (
    "cg_xml",
    "chem",
    "geometry",
    "mapping",
    "molecules",
    "pdb_template",
    "pipeline",
    "pipeline_config",
    "pipeline_forcefield",
    "reactions",
    "refinement",
    "runner",
    "templates",
    "writer",
)


def _module_record(name: str) -> dict[str, Any]:
    module = import_module(f"multiscale.cg2aa.{name}")
    file_text = getattr(module, "__file__", None)
    path = Path(file_text).resolve() if file_text else None
    inside_package = bool(path and "multiscale" in path.parts and "cg2aa" in path.parts)
    return {
        "name": name,
        "module": f"multiscale.cg2aa.{name}",
        "file": str(path) if path else None,
        "inside_package": inside_package,
        "ok": inside_package,
    }


def internal_module_status() -> dict[str, Any]:
    records = [_module_record(name) for name in PUBLIC_MODULES]
    try:
        from multiscale.forcefield.assign import internal_module_status as forcefield_status

        forcefield = forcefield_status()
    except Exception as exc:  # pragma: no cover - defensive status path
        forcefield = {"ok": False, "error_type": type(exc).__name__, "error": str(exc)}
    checks = {
        "modules_present": all(item["file"] for item in records),
        "modules_inside_package": all(item["inside_package"] for item in records),
        "forcefield_backend_internal": bool(forcefield.get("ok")) and forcefield.get("source") == "internal-multiscale",
    }
    return {
        "ok": all(checks.values()),
        "scope": "CG-to-AA template, reaction, output, and bundled force-field modules are packaged in Multiscale",
        "checks": checks,
        "modules": records,
        "forcefield": forcefield,
        "release_status": "public-self-contained",
    }
