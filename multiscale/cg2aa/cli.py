"""Multiscale CG-to-AA command-line interface."""

from __future__ import annotations

import argparse
import json
from typing import Any

from multiscale.core.report import make_report

from .config import summarize_config, validate_config
from .internal import internal_module_status
from .runner import run_cg2aa

_COMMANDS = {"run", "summary", "validate", "internal-status", "-h", "--help"}


def _print(data: dict[str, Any], json_output: bool, *, multiscale_report: bool = False, command: str) -> None:
    if multiscale_report:
        data = make_report(module="cg2aa", command=command, payload=data)
        json_output = True
    if json_output:
        print(json.dumps(data, indent=2, ensure_ascii=False))
    else:
        for key, value in data.items():
            print(f"{key}: {value}")


def main(argv: list[str] | None = None) -> int:
    args_list = list([] if argv is None else argv)
    if args_list and args_list[0] not in _COMMANDS:
        return main(["run", *args_list])

    parser = argparse.ArgumentParser(
        prog="multiscale cg2aa",
        description="CG-to-AA template placement, reaction processing, and OPLS-style force-field export from one formal multiscale.cg2aa.v1 JSON configuration.",
    )
    sub = parser.add_subparsers(dest="command")

    p_run = sub.add_parser("run", help="run the CG-to-AA pipeline")
    p_run.add_argument("config")
    p_run.add_argument("--multiscale-report", action="store_true")

    p_summary = sub.add_parser("summary", help="summarize a formal CG-to-AA configuration")
    p_summary.add_argument("config")
    p_summary.add_argument("--json", action="store_true")
    p_summary.add_argument("--multiscale-report", action="store_true")

    p_validate = sub.add_parser("validate", help="validate configuration paths and force-field settings")
    p_validate.add_argument("config")
    p_validate.add_argument("--json", action="store_true")
    p_validate.add_argument("--multiscale-report", action="store_true")
    p_validate.add_argument("--no-forcefield", action="store_true")

    p_status = sub.add_parser("internal-status", help="report availability of packaged CG-to-AA modules")
    p_status.add_argument("--json", action="store_true")
    p_status.add_argument("--multiscale-report", action="store_true")

    args = parser.parse_args(args_list)
    if args.command is None:
        parser.print_help()
        return 0

    if args.command == "run":
        try:
            data = run_cg2aa(args.config)
        except Exception as exc:
            data = {"ok": False, "config": args.config, "error_type": type(exc).__name__, "error": str(exc)}
            _print(data, True, multiscale_report=args.multiscale_report, command="run")
            return 1
        data.setdefault("ok", True)
        _print(data, True, multiscale_report=args.multiscale_report, command="run")
        return 0

    if args.command == "summary":
        try:
            data = {"ok": True, "summary": summarize_config(args.config).to_dict()}
        except Exception as exc:
            data = {"ok": False, "config": args.config, "error_type": type(exc).__name__, "error": str(exc)}
        _print(data, args.json, multiscale_report=args.multiscale_report, command="summary")
        return 0 if data.get("ok") else 1

    if args.command == "validate":
        try:
            data = validate_config(args.config, check_forcefield=not args.no_forcefield)
        except Exception as exc:
            data = {"ok": False, "config": args.config, "error_type": type(exc).__name__, "error": str(exc)}
        _print(data, args.json, multiscale_report=args.multiscale_report, command="validate")
        return 0 if data.get("ok") else 1

    data = internal_module_status()
    _print(data, args.json, multiscale_report=args.multiscale_report, command="internal-status")
    return 0 if data.get("ok") else 1


if __name__ == "__main__":
    raise SystemExit(main())
