"""Environment diagnostics for the self-contained Multiscale package."""

from __future__ import annotations

import importlib
import json
import os
import platform
import shutil
import subprocess
import sys
from typing import Any

from .core.runtime import discover_runtime_layout

REQUIRES_PYTHON = ">=3.12,<3.13"
_REQUIRED_PACKAGES = ("numpy", "networkx")
_OPTIONAL_PACKAGES = ("MDAnalysis",)


def _package_info(name: str) -> dict[str, Any]:
    try:
        mod = importlib.import_module(name)
    except Exception as exc:  # pragma: no cover - diagnostic path
        return {"ok": False, "error_type": type(exc).__name__, "error": str(exc)}
    return {
        "ok": True,
        "version": getattr(mod, "__version__", None),
        "file": getattr(mod, "__file__", None),
    }


def _gmx_info(gmx_bin: str | None) -> dict[str, Any]:
    if not gmx_bin:
        return {"ok": False, "required": False, "error": "gmx not found in PATH"}
    try:
        proc = subprocess.run([gmx_bin, "--version"], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=20)
    except Exception as exc:  # pragma: no cover - diagnostic path
        return {"ok": False, "required": False, "path": gmx_bin, "error_type": type(exc).__name__, "error": str(exc)}
    version_line = None
    for line in proc.stdout.splitlines():
        if "GROMACS version:" in line:
            version_line = line.strip()
            break
    return {
        "ok": proc.returncode == 0,
        "required": False,
        "path": shutil.which(gmx_bin) or gmx_bin,
        "returncode": proc.returncode,
        "version_line": version_line,
    }


def collect_doctor() -> dict[str, Any]:
    layout = discover_runtime_layout()
    python_supported = sys.version_info[:2] == (3, 12)
    required_packages = {name: _package_info(name) for name in _REQUIRED_PACKAGES}
    optional_packages = {name: _package_info(name) for name in _OPTIONAL_PACKAGES}
    gromacs = _gmx_info(layout.gmx_bin)
    required_runtime = {
        "multiscale_importable": _package_info("multiscale").get("ok", False),
        "self_contained_layout": layout.self_contained,
        "python_supported": python_supported,
        "required_python_packages": all(info.get("ok") for info in required_packages.values()),
    }
    optional_tools = {
        "gromacs_available": bool(gromacs.get("ok")),
    }
    return {
        "ok": all(required_runtime.values()),
        "multiscale_version": _package_info("multiscale").get("version"),
        "python": {
            "executable": sys.executable,
            "version": sys.version,
            "requires_python": REQUIRES_PYTHON,
            "supported": python_supported,
            "platform": platform.platform(),
            "conda_default_env": os.environ.get("CONDA_DEFAULT_ENV"),
            "conda_prefix": os.environ.get("CONDA_PREFIX"),
        },
        "layout": layout.to_dict(),
        "self_contained": {
            "ok": True,
            "external_projects_required": False,
            "default_backends": {
                "aa2cg_map": "internal",
                "aa2cg_exporters": "internal",
                "aa2cg_rdf": "internal",
                "cg2aa_run": "internal",
                "forcefield": "internal-multiscale",
            },
        },
        "packages": {**required_packages, **optional_packages},
        "required_packages": required_packages,
        "optional_packages": optional_packages,
        "gromacs": gromacs,
        "optional_tools": optional_tools,
        "checks": required_runtime,
    }


def print_doctor(json_output: bool = False) -> int:
    data = collect_doctor()
    if json_output:
        print(json.dumps(data, indent=2, ensure_ascii=False))
    else:
        print("Multiscale doctor")
        print(f"ok: {data['ok']}")
        print(f"python: {data['python']['executable']} ({data['python']['requires_python']}, supported={data['python']['supported']})")
        print(f"self_contained: {data['self_contained']['ok']}")
        gromacs_label = data["gromacs"].get("path") or data["gromacs"].get("error")
        print(f"gmx optional: {gromacs_label}")
        print("required packages:")
        for name, info in data["required_packages"].items():
            status = "OK" if info.get("ok") else "MISSING/ERROR"
            version = info.get("version") or info.get("error", "")
            print(f"  {name}: {status} {version}")
        print("optional packages/tools:")
        for name, info in data["optional_packages"].items():
            status = "OK" if info.get("ok") else "MISSING/ERROR"
            version = info.get("version") or info.get("error", "")
            print(f"  {name}: {status} {version}")
        print(f"  gromacs: {'OK' if data['optional_tools']['gromacs_available'] else 'MISSING/ERROR'}")
    return 0 if data["ok"] else 1
