"""Multiscale forcefield CLI namespace."""

from __future__ import annotations

import argparse
import json
from typing import Any

from multiscale.core.report import make_report
from .registry import list_forcefields
from .opls import atom_type_lookup, summarize_opls, validate_config_forcefield
from .assign import internal_module_status


def _print(data: dict[str, Any], json_output: bool, *, multiscale_report: bool = False, command: str = "list") -> None:
    if multiscale_report:
        data = make_report(module="forcefield", command=command, payload=data)
        json_output = True
    if json_output:
        print(json.dumps(data, indent=2, ensure_ascii=False))
    else:
        for key, value in data.items():
            print(f"{key}: {value}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="multiscale forcefield", description="Multiscale self-contained force-field namespace with bundled OPLS assets")
    sub = parser.add_subparsers(dest="command")

    p_list = sub.add_parser("list", help="list discoverable force fields")
    p_list.add_argument("--json", action="store_true")
    p_list.add_argument("--multiscale-report", action="store_true")

    p_inspect = sub.add_parser("inspect", help="inspect one force field")
    p_inspect.add_argument("name")
    p_inspect.add_argument("--json", action="store_true")
    p_inspect.add_argument("--multiscale-report", action="store_true")
    p_inspect.add_argument("--atom-type")

    p_validate = sub.add_parser("validate", help="validate force-field assets")
    p_validate.add_argument("name")
    p_validate.add_argument("--json", action="store_true")
    p_validate.add_argument("--multiscale-report", action="store_true")

    p_internal = sub.add_parser("internal-status", help="check bundled internal force-field backend")
    p_internal.add_argument("name")
    p_internal.add_argument("--json", action="store_true")
    p_internal.add_argument("--multiscale-report", action="store_true")

    p_config = sub.add_parser("validate-config", help="validate a cg2aa force_field config block")
    p_config.add_argument("config")
    p_config.add_argument("--json", action="store_true")
    p_config.add_argument("--multiscale-report", action="store_true")

    args = parser.parse_args(argv)
    if args.command in (None, "list"):
        data = {"ok": True, "forcefields": list_forcefields()}
        _print(data, getattr(args, "json", False), multiscale_report=getattr(args, "multiscale_report", False), command="list")
        return 0

    if args.command == "inspect":
        if args.name != "opls":
            data = {"ok": False, "forcefield": None, "error": f"force field not found: {args.name}"}
            _print(data, args.json, multiscale_report=args.multiscale_report, command="inspect")
            return 1
        summary = summarize_opls().to_dict()
        data: dict[str, Any] = {"ok": summary.get("ok", False), "forcefield": "opls", "summary": summary}
        if args.atom_type:
            data["atom_type"] = atom_type_lookup().get(args.atom_type)
            if data["atom_type"] is None:
                data["ok"] = False
                data["error"] = f"atom type not found: {args.atom_type}"
        _print(data, args.json, multiscale_report=args.multiscale_report, command="inspect")
        return 0 if data.get("ok") else 1

    if args.command == "validate":
        if args.name != "opls":
            data = {"ok": False, "forcefield": args.name, "error": f"unsupported force field: {args.name}"}
        else:
            summary = summarize_opls().to_dict()
            data = {"ok": summary.get("ok", False), "forcefield": "opls", "summary": summary}
        _print(data, args.json, multiscale_report=args.multiscale_report, command="validate")
        return 0 if data.get("ok") else 1

    if args.command == "internal-status":
        data = internal_module_status() if args.name == "opls" else {"ok": False, "forcefield": args.name, "error": f"unsupported force field: {args.name}"}
        _print(data, args.json, multiscale_report=args.multiscale_report, command="internal-status")
        return 0 if data.get("ok") else 1

    if args.command == "validate-config":
        try:
            data = validate_config_forcefield(args.config)
        except Exception as exc:
            data = {"ok": False, "config": args.config, "error_type": type(exc).__name__, "error": str(exc)}
        _print(data, args.json, multiscale_report=args.multiscale_report, command="validate-config")
        return 0 if data.get("ok") else 1

    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
