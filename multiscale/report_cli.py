"""CLI for Multiscale report schema/wrapping/validation."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from .core.report import emit_json, load_json, make_report, report_schema, validate_report, write_json


def _collect_run_dir(run_dir: str | Path) -> dict[str, Any]:
    root = Path(run_dir).expanduser().resolve()
    files: dict[str, str] = {}
    known = {
        "aa2cg_run_report": ["aa2cg_run_report.json"],
        "aa2cg_map_cgsystem": ["mapped/cg_system_mapped.json", "aa2cg_map/cg_system_mapped.json"],
        "aa2cg_map_report": ["mapped/map_structure_report.json", "aa2cg_map/map_structure_report.json"],
        "galamost_topology": ["galamost/topology.xml", "aa2cg_export_galamost/topology.xml"],
        "galamost_topology_inbox": ["galamost/topology_inbox.xml", "aa2cg_export_galamost/topology_inbox.xml"],
        "gromacs_gro": ["gromacs/cg.gro", "aa2cg_export_gromacs/cg.gro"],
        "gromacs_itp": ["gromacs/cg.itp", "aa2cg_export_gromacs/cg.itp"],
        "gromacs_top": ["gromacs/topol.top", "aa2cg_export_gromacs/topol.top"],
        "lammps_data": ["lammps/system.data", "aa2cg_export_lammps/system.data"],
        "lammps_input": ["lammps/in.cg", "aa2cg_export_lammps/in.cg"],
        "rdf_report": ["rdf/rdf_report.json", "aa2cg_rdf/rdf_report.json"],
    }
    for key, candidates in known.items():
        for rel in candidates:
            candidate = root / rel
            if candidate.exists():
                files[key] = str(candidate)
                break
    if "aa2cg_run_report" not in files:
        matches = sorted(root.rglob("aa2cg_run_report.json")) if root.is_dir() else []
        if matches:
            files["aa2cg_run_report"] = str(matches[0])
    payload = {
        "ok": root.is_dir(),
        "run_dir": str(root),
        "files": files,
        "counts": {"known_files_found": len(files)},
        "checks": {"run_dir_exists": root.is_dir()},
    }
    return make_report(module="report", command="collect", payload=payload, ok=payload["ok"])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="multiscale report", description="Multiscale report protocol tools")
    sub = parser.add_subparsers(dest="subcommand")

    p_schema = sub.add_parser("schema", help="print the Multiscale report schema")
    p_schema.add_argument("--json", action="store_true")

    p_wrap = sub.add_parser("wrap", help="wrap an existing JSON payload in a Multiscale report envelope")
    p_wrap.add_argument("payload_json")
    p_wrap.add_argument("--module", required=True)
    p_wrap.add_argument("--command", dest="wrapped_command", required=True)
    p_wrap.add_argument("--out")
    p_wrap.add_argument("--json", action="store_true")

    p_validate = sub.add_parser("validate", help="validate a Multiscale report JSON file")
    p_validate.add_argument("report_json")
    p_validate.add_argument("--json", action="store_true")

    p_collect = sub.add_parser("collect", help="collect common Multiscale output files from a run directory")
    p_collect.add_argument("run_dir")
    p_collect.add_argument("--out")
    p_collect.add_argument("--json", action="store_true")

    args = parser.parse_args(argv)
    if args.subcommand in (None, "schema"):
        data = report_schema()
        emit_json(data)
        return 0

    if args.subcommand == "wrap":
        try:
            payload = load_json(args.payload_json)
            data = make_report(module=args.module, command=args.wrapped_command, payload=payload)
            if args.out:
                write_json(args.out, data)
        except Exception as exc:
            data = make_report(
                module=args.module,
                command=args.wrapped_command,
                ok=False,
                payload={},
                errors=[{"type": type(exc).__name__, "message": str(exc)}],
            )
        emit_json(data)
        return 0 if data.get("ok") else 1

    if args.subcommand == "validate":
        try:
            data = validate_report(load_json(args.report_json))
        except Exception as exc:
            data = {"ok": False, "error_type": type(exc).__name__, "error": str(exc)}
        emit_json(data)
        return 0 if data.get("ok") else 1

    if args.subcommand == "collect":
        data = _collect_run_dir(args.run_dir)
        if args.out:
            write_json(args.out, data)
        emit_json(data)
        return 0 if data.get("ok") else 1

    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
