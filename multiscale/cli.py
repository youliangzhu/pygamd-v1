"""Top-level Multiscale CLI."""

from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .aa2cg.cli import main as aa2cg_main
from .audit.cli import main as audit_main
from .clean import main as clean_main
from .doctor import print_doctor
from .forcefield.cli import main as forcefield_main
from .cg2aa.cli import main as cg2aa_main
from .report_cli import main as report_main
from .core.runtime import discover_runtime_layout
from .core.cg_system import load_cgsystem_json
from .core.galamost_xml import compare_galamost_counts, inspect_galamost_xml
from .core.report import make_report


def _print_main_help() -> None:
    print("""Multiscale CLI

Usage:
  multiscale doctor [--json]
  multiscale clean-cache [PATH] [--dry-run] [--include-build] [--json]
  multiscale core info [--json] [--multiscale-report]
  multiscale core cgsystem-summary <cg_system.json> [--json] [--multiscale-report]
  multiscale core inspect-pygamd <topology.xml> [--json] [--multiscale-report]
  multiscale core compare-pygamd <reference.xml> <generated.xml> [--json] [--multiscale-report]
  multiscale report schema [--json]
  multiscale report wrap <payload.json> --module MODULE --command COMMAND [--out report.json] [--json]
  multiscale report validate <report.json> [--json]
  multiscale report collect <run_dir> [--out report.json] [--json]
  multiscale audit status [--json] [--multiscale-report]
  multiscale audit source-files [--json] [--out-json source_files.json]
  multiscale audit references [--json] [--out-json refs.json]
  multiscale audit matrix [--json] [--out-json matrix.json] [--out-md component_matrix.md] [--multiscale-report]
  multiscale aa2cg <aa2cg.json>
  multiscale aa2cg validate-config <aa2cg.json> [--json] [--multiscale-report]
  multiscale aa2cg run <aa2cg.json> [--gmx-bin PATH] [--json] [--multiscale-report]
  multiscale aa2cg export-pygamd <cg_system.json> --outdir <outdir> [--json] [--multiscale-report]
  multiscale aa2cg export-gromacs <cg_system.json> --outdir <outdir> [--json] [--multiscale-report]
  multiscale aa2cg export-lammps <cg_system.json> --outdir <outdir> [--json] [--multiscale-report]
  multiscale aa2cg internal-status [--json] [--multiscale-report]
  multiscale cg2aa <cg2aa.json>
  multiscale cg2aa run <cg2aa.json> [--multiscale-report]
  multiscale cg2aa summary <cg2aa.json> [--json] [--multiscale-report]
  multiscale cg2aa validate <cg2aa.json> [--json] [--multiscale-report]
  multiscale cg2aa internal-status [--json] [--multiscale-report]
  multiscale forcefield list [--json] [--multiscale-report]
  multiscale forcefield inspect opls [--json] [--multiscale-report]
  multiscale forcefield validate opls [--json] [--multiscale-report]
  multiscale forcefield internal-status opls [--json] [--multiscale-report]
  multiscale forcefield validate-config <cg2aa.json> [--json] [--multiscale-report]
  multiscale version [--json]

Multiscale is a self-contained AA→CG / CG→AA workflow package. AA→CG and CG→AA projects use the formal multiscale.aa2cg.v1 and multiscale.cg2aa.v1 JSON schemas.
""")


def _print_data(data: dict, *, json_output: bool, multiscale_report: bool = False, module: str = "core", command: str = "info") -> None:
    if multiscale_report:
        data = make_report(module=module, command=command, payload=data)
        json_output = True
    if json_output:
        print(json.dumps(data, indent=2, ensure_ascii=False))
    else:
        for key, value in data.items():
            print(f"{key}: {value}")


def _core(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="multiscale core")
    sub = parser.add_subparsers(dest="command")

    p_info = sub.add_parser("info", help="show runtime layout")
    p_info.add_argument("--json", action="store_true")
    p_info.add_argument("--multiscale-report", action="store_true")

    p_cg = sub.add_parser("cgsystem-summary", help="summarize a CGSystem JSON file")
    p_cg.add_argument("path")
    p_cg.add_argument("--json", action="store_true")
    p_cg.add_argument("--multiscale-report", action="store_true")

    p_inspect = sub.add_parser("inspect-pygamd", help="inspect a topology XML file used by PyGAMD")
    p_inspect.add_argument("path")
    p_inspect.add_argument("--json", action="store_true")
    p_inspect.add_argument("--multiscale-report", action="store_true")

    p_compare = sub.add_parser("compare-pygamd", help="compare two topology XML files used by PyGAMD")
    p_compare.add_argument("reference")
    p_compare.add_argument("generated")
    p_compare.add_argument("--json", action="store_true")
    p_compare.add_argument("--multiscale-report", action="store_true")

    args = parser.parse_args(argv)
    if args.command in (None, "info"):
        data = {"ok": True, "layout": discover_runtime_layout().to_dict()}
        _print_data(data, json_output=getattr(args, "json", False), multiscale_report=getattr(args, "multiscale_report", False), module="core", command="info")
        return 0

    if args.command == "cgsystem-summary":
        try:
            system = load_cgsystem_json(args.path)
            data = {"ok": True, "path": args.path, "summary": system.summary()}
        except Exception as exc:
            data = {"ok": False, "path": args.path, "error_type": type(exc).__name__, "error": str(exc)}
        _print_data(data, json_output=args.json, multiscale_report=args.multiscale_report, module="core", command="cgsystem-summary")
        return 0 if data.get("ok") else 1

    if args.command == "inspect-pygamd":
        try:
            data = inspect_galamost_xml(args.path)
        except Exception as exc:
            data = {"ok": False, "path": args.path, "error_type": type(exc).__name__, "error": str(exc)}
        _print_data(data, json_output=args.json, multiscale_report=args.multiscale_report, module="core", command="inspect-pygamd")
        return 0 if data.get("ok") else 1

    if args.command == "compare-pygamd":
        try:
            data = compare_galamost_counts(args.reference, args.generated)
        except Exception as exc:
            data = {"ok": False, "reference": args.reference, "generated": args.generated, "error_type": type(exc).__name__, "error": str(exc)}
        _print_data(data, json_output=args.json, multiscale_report=args.multiscale_report, module="core", command="compare-pygamd")
        return 0 if data.get("ok") else 1

    return 0


def _version(argv: list[str]) -> int:
    json_output = "--json" in argv
    data = {"ok": True, "version": __version__, "release": "public"}
    if json_output:
        print(json.dumps(data, indent=2, ensure_ascii=False))
    else:
        print(__version__)
    return 0


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in {"-h", "--help"}:
        _print_main_help()
        return 0

    command, rest = argv[0], argv[1:]
    if command == "doctor":
        return print_doctor(json_output="--json" in rest)
    if command == "clean-cache":
        return clean_main(rest)
    if command == "core":
        return _core(rest)
    if command == "report":
        return report_main(rest)
    if command == "audit":
        return audit_main(rest)
    if command == "aa2cg":
        return aa2cg_main(rest)
    if command == "cg2aa":
        return cg2aa_main(rest)
    if command == "forcefield":
        return forcefield_main(rest)
    if command in {"version", "--version"}:
        return _version(rest)

    print(f"ERROR: unknown multiscale command: {command}", file=sys.stderr)
    _print_main_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
