"""Multiscale AA-to-CG command-line interface."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from multiscale.core.cg_system import load_cgsystem_json
from multiscale.core.report import make_report

from .exporters.galamost import export_galamost_auxiliary, export_galamost_topology
from .exporters.gromacs import export_gromacs_scaffold
from .exporters.lammps import export_lammps_scaffold
from .internal import aa2cg_internal_status
from .json_config import validate_aa2cg_config
from .workflow import run_aa2cg_config

_COMMANDS = {
    "run",
    "validate-config",
    "export-pygamd",
    "export-gromacs",
    "export-lammps",
    "internal-status",
    "-h",
    "--help",
}


def _print(data: dict[str, Any], *, json_output: bool, multiscale_report: bool = False, command: str) -> None:
    if multiscale_report:
        data = make_report(module="aa2cg", command=command, payload=data)
        json_output = True
    if json_output:
        print(json.dumps(data, indent=2, ensure_ascii=False))
    else:
        for key, value in data.items():
            print(f"{key}: {value}")


def _help() -> None:
    print(
        """usage: multiscale aa2cg [-h] {run,validate-config,export-pygamd,export-gromacs,export-lammps,internal-status} ...

AA-to-CG mapping, selected format export, and optional RDF analysis from one
formal multiscale.aa2cg.v1 JSON configuration.

commands:
  run              execute mapping, selected exporters, and optional RDF
  validate-config  validate a formal multiscale.aa2cg.v1 JSON configuration
  export-pygamd    export PyGAMD topology and support files from a retained CGSystem JSON file
  export-gromacs   export GROMACS files from a retained CGSystem JSON file
  export-lammps    export a LAMMPS topology scaffold from a retained CGSystem JSON file
  internal-status  report availability of packaged AA-to-CG modules

A configuration path may be supplied directly as shorthand for `run`.
Add --multiscale-report to wrap command output in multiscale.report.v1.
"""
    )


def _load_system(source: str | Path):
    source_path = Path(source).expanduser().resolve()
    if not source_path.is_file():
        raise ValueError(f"CGSystem JSON input not found: {source_path}")
    return load_cgsystem_json(source_path)


def _export_pygamd(source: str, outdir: str, epsilon: float) -> dict[str, Any]:
    system = _load_system(source)
    topology = export_galamost_topology(system, outdir, eps=epsilon)
    auxiliary = export_galamost_auxiliary(system, outdir)
    output = Path(outdir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "pygamd_export_report.json").write_text(
        json.dumps(topology.to_dict(), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (output / "pygamd_aux_report.json").write_text(
        json.dumps(auxiliary.to_dict(), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return {
        **topology.to_dict(),
        "auxiliary": auxiliary.to_dict(),
        "backend": "internal",
        "ok": bool(topology.ok and auxiliary.ok),
    }


def _export_gromacs(source: str, outdir: str) -> dict[str, Any]:
    report = export_gromacs_scaffold(_load_system(source), outdir)
    output = Path(outdir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "gromacs_export_report.json").write_text(
        json.dumps(report.to_dict(), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return {**report.to_dict(), "backend": "internal"}


def _export_lammps(source: str, outdir: str) -> dict[str, Any]:
    report = export_lammps_scaffold(_load_system(source), outdir)
    output = Path(outdir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "lammps_export_report.json").write_text(
        json.dumps(report.to_dict(), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return {**report.to_dict(), "backend": "internal"}


def main(argv: list[str] | None = None) -> int:
    args_list = list([] if argv is None else argv)
    if not args_list or args_list[0] in {"-h", "--help"}:
        _help()
        return 0
    if args_list[0] not in _COMMANDS:
        return main(["run", *args_list])

    command, rest = args_list[0], args_list[1:]

    if command == "run":
        parser = argparse.ArgumentParser(prog="multiscale aa2cg run")
        parser.add_argument("config")
        parser.add_argument("--gmx-bin", default=None, help="override the configured GROMACS executable")
        parser.add_argument("--json", action="store_true")
        parser.add_argument("--multiscale-report", action="store_true")
        args = parser.parse_args(rest)
        log = None if args.json or args.multiscale_report else lambda message: print(message, file=sys.stderr, flush=True)
        try:
            data = run_aa2cg_config(args.config, gmx_bin=args.gmx_bin, log=log)
        except (RuntimeError, ValueError) as exc:
            print(f"aa2cg run failed: {exc}", file=sys.stderr)
            return 2
        except Exception as exc:  # pragma: no cover - defensive CLI boundary
            print(f"aa2cg run failed: {exc}", file=sys.stderr)
            return 1
        _print(data, json_output=args.json, multiscale_report=args.multiscale_report, command="run")
        return 0 if data.get("ok") else 1

    if command == "validate-config":
        parser = argparse.ArgumentParser(prog="multiscale aa2cg validate-config")
        parser.add_argument("config")
        parser.add_argument("--no-file-check", action="store_true")
        parser.add_argument("--json", action="store_true")
        parser.add_argument("--multiscale-report", action="store_true")
        args = parser.parse_args(rest)
        data = validate_aa2cg_config(args.config, require_files=not args.no_file_check)
        _print(data, json_output=args.json, multiscale_report=args.multiscale_report, command="validate-config")
        return 0 if data.get("ok") else 1

    if command == "export-pygamd":
        parser = argparse.ArgumentParser(prog="multiscale aa2cg export-pygamd")
        parser.add_argument("source", help="retained CGSystem JSON input")
        parser.add_argument("--outdir", required=True)
        parser.add_argument("--eps", type=float, default=1.0)
        parser.add_argument("--json", action="store_true")
        parser.add_argument("--multiscale-report", action="store_true")
        args = parser.parse_args(rest)
        try:
            data = _export_pygamd(args.source, args.outdir, args.eps)
        except Exception as exc:
            print(f"export-pygamd failed: {exc}", file=sys.stderr)
            return 1
        _print(data, json_output=args.json, multiscale_report=args.multiscale_report, command="export-pygamd")
        return 0 if data.get("ok") else 1

    if command == "export-gromacs":
        parser = argparse.ArgumentParser(prog="multiscale aa2cg export-gromacs")
        parser.add_argument("source", help="retained CGSystem JSON input")
        parser.add_argument("--outdir", required=True)
        parser.add_argument("--json", action="store_true")
        parser.add_argument("--multiscale-report", action="store_true")
        args = parser.parse_args(rest)
        try:
            data = _export_gromacs(args.source, args.outdir)
        except Exception as exc:
            print(f"export-gromacs failed: {exc}", file=sys.stderr)
            return 1
        _print(data, json_output=args.json, multiscale_report=args.multiscale_report, command="export-gromacs")
        return 0 if data.get("ok") else 1

    if command == "export-lammps":
        parser = argparse.ArgumentParser(prog="multiscale aa2cg export-lammps")
        parser.add_argument("source", help="retained CGSystem JSON input")
        parser.add_argument("--outdir", required=True)
        parser.add_argument("--json", action="store_true")
        parser.add_argument("--multiscale-report", action="store_true")
        args = parser.parse_args(rest)
        try:
            data = _export_lammps(args.source, args.outdir)
        except Exception as exc:
            print(f"export-lammps failed: {exc}", file=sys.stderr)
            return 1
        _print(data, json_output=args.json, multiscale_report=args.multiscale_report, command="export-lammps")
        return 0 if data.get("ok") else 1

    if command == "internal-status":
        parser = argparse.ArgumentParser(prog="multiscale aa2cg internal-status")
        parser.add_argument("--json", action="store_true")
        parser.add_argument("--multiscale-report", action="store_true")
        args = parser.parse_args(rest)
        data = aa2cg_internal_status()
        _print(data, json_output=args.json, multiscale_report=args.multiscale_report, command="internal-status")
        return 0 if data.get("ok") else 1

    _help()
    return 2
