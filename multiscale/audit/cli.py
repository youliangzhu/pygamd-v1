"""CLI for Multiscale self-contained audit tools."""

from __future__ import annotations

import argparse
from typing import Any

from multiscale.core.report import emit_json, make_report, write_json
from .component import (
    build_component_matrix,
    matrix_to_markdown,
    scan_external_dependency_references,
    self_contained_status,
    source_file_inventory,
    source_project_summary,
    write_markdown_matrix,
)


def _emit(data: dict[str, Any], *, json_output: bool, multiscale_report: bool = False, command: str = "status") -> None:
    if multiscale_report:
        data = make_report(module="audit", command=command, payload=data)
        json_output = True
    if json_output:
        emit_json(data)
    else:
        for key, value in data.items():
            print(f"{key}: {value}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="multiscale audit", description="Multiscale self-contained source/runtime audit tools")
    sub = parser.add_subparsers(dest="command")

    p_status = sub.add_parser("status", help="show self-contained package status")
    p_status.add_argument("--json", action="store_true")
    p_status.add_argument("--multiscale-report", action="store_true")

    p_source = sub.add_parser("source-files", help="inventory source files by public category")
    p_source.add_argument("--json", action="store_true")
    p_source.add_argument("--multiscale-report", action="store_true")
    p_source.add_argument("--full", action="store_true", help="include per-file inventory instead of summary only")
    p_source.add_argument("--source-root")
    p_source.add_argument("--out-json")

    p_refs = sub.add_parser("references", help="scan current multiscale source for external-runtime dependency references")
    p_refs.add_argument("--json", action="store_true")
    p_refs.add_argument("--multiscale-report", action="store_true")
    p_refs.add_argument("--source-root")
    p_refs.add_argument("--out-json")

    p_matrix = sub.add_parser("matrix", help="build the public Multiscale component matrix")
    p_matrix.add_argument("--json", action="store_true")
    p_matrix.add_argument("--multiscale-report", action="store_true")
    p_matrix.add_argument("--out-json")
    p_matrix.add_argument("--out-md")

    args = parser.parse_args(argv)

    if args.command in (None, "status"):
        data = self_contained_status()
        _emit(data, json_output=getattr(args, "json", False), multiscale_report=getattr(args, "multiscale_report", False), command="status")
        return 0

    if args.command == "source-files":
        data = source_file_inventory(args.source_root) if args.full else source_project_summary(args.source_root)
        if args.out_json:
            write_json(args.out_json, data)
        _emit(data, json_output=args.json, multiscale_report=args.multiscale_report, command="source-files")
        return 0 if data.get("ok", True) else 1

    if args.command == "references":
        data = scan_external_dependency_references(args.source_root) if args.source_root else scan_external_dependency_references()
        if args.out_json:
            write_json(args.out_json, data)
        _emit(data, json_output=args.json, multiscale_report=args.multiscale_report, command="references")
        return 0 if data.get("ok", True) else 1

    if args.command == "matrix":
        data = build_component_matrix()
        if args.out_json:
            write_json(args.out_json, data)
        if args.out_md:
            path = write_markdown_matrix(args.out_md, data)
            data = dict(data)
            data.setdefault("files", {})
            data["files"] = {"component_matrix": str(path), **data.get("files", {})}
        if args.json or args.multiscale_report:
            _emit(data, json_output=args.json, multiscale_report=args.multiscale_report, command="matrix")
        else:
            print(matrix_to_markdown(data))
        return 0 if data.get("ok", True) else 1

    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
