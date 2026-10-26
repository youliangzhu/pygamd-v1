"""Unified Multiscale report protocol helpers.

Multiscale standard commands may emit the ``multiscale.report.v1`` envelope with
``--multiscale-report``. Arbitrary JSON output can also be wrapped with
``multiscale report wrap``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

from multiscale import __version__

REPORT_SCHEMA = "multiscale.report.v1"


def _now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


@dataclass
class MultiscaleReport:
    ok: bool
    module: str
    command: str
    payload: dict[str, Any] = field(default_factory=dict)
    files: dict[str, Any] = field(default_factory=dict)
    counts: dict[str, Any] = field(default_factory=dict)
    checks: dict[str, Any] = field(default_factory=dict)
    warnings: list[Any] = field(default_factory=list)
    errors: list[Any] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    schema: str = REPORT_SCHEMA
    version: str = __version__
    timestamp: str = field(default_factory=_now_iso)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "ok": self.ok,
            "module": self.module,
            "command": self.command,
            "version": self.version,
            "timestamp": self.timestamp,
            "files": self.files,
            "counts": self.counts,
            "checks": self.checks,
            "warnings": self.warnings,
            "errors": self.errors,
            "metadata": self.metadata,
            "payload": self.payload,
        }


def extract_common(payload: dict[str, Any]) -> dict[str, Any]:
    """Extract common Agent-friendly fields from an existing command payload."""
    return {
        "files": _as_dict(payload.get("files") or payload.get("outputs")),
        "counts": _as_dict(payload.get("counts")),
        "checks": _as_dict(payload.get("checks")),
        "warnings": _as_list(payload.get("warnings")),
    }


def make_report(
    *,
    module: str,
    command: str,
    payload: dict[str, Any] | None = None,
    ok: bool | None = None,
    files: dict[str, Any] | None = None,
    counts: dict[str, Any] | None = None,
    checks: dict[str, Any] | None = None,
    warnings: list[Any] | None = None,
    errors: list[Any] | None = None,
    metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    payload = payload or {}
    common = extract_common(payload)
    inferred_ok = bool(payload.get("ok", False)) if "ok" in payload else not bool(errors)
    report = MultiscaleReport(
        ok=inferred_ok if ok is None else bool(ok),
        module=module,
        command=command,
        payload=payload,
        files=common["files"] if files is None else files,
        counts=common["counts"] if counts is None else counts,
        checks=common["checks"] if checks is None else checks,
        warnings=common["warnings"] if warnings is None else warnings,
        errors=[] if errors is None else errors,
        metadata={} if metadata is None else metadata,
    )
    return report.to_dict()


def ok_report(module: str, command: str, **extra: Any) -> dict[str, Any]:
    result: dict[str, Any] = {"ok": True, "module": module, "command": command}
    result.update(extra)
    return result


def fail_report(module: str, command: str, message: str, **extra: Any) -> dict[str, Any]:
    result: dict[str, Any] = {"ok": False, "module": module, "command": command, "error": message}
    result.update(extra)
    return result


def validate_report(data: dict[str, Any]) -> dict[str, Any]:
    checks = {
        "has_schema": data.get("schema") == REPORT_SCHEMA,
        "has_ok": isinstance(data.get("ok"), bool),
        "has_module": isinstance(data.get("module"), str) and bool(data.get("module")),
        "has_command": isinstance(data.get("command"), str) and bool(data.get("command")),
        "has_payload": isinstance(data.get("payload"), dict),
        "has_files": isinstance(data.get("files"), dict),
        "has_counts": isinstance(data.get("counts"), dict),
        "has_checks": isinstance(data.get("checks"), dict),
        "has_warnings": isinstance(data.get("warnings"), list),
        "has_errors": isinstance(data.get("errors"), list),
    }
    return {"ok": all(checks.values()), "schema": REPORT_SCHEMA, "checks": checks}


def report_schema() -> dict[str, Any]:
    return {
        "ok": True,
        "schema": REPORT_SCHEMA,
        "version": __version__,
        "required_fields": [
            "schema",
            "ok",
            "module",
            "command",
            "version",
            "timestamp",
            "files",
            "counts",
            "checks",
            "warnings",
            "errors",
            "metadata",
            "payload",
        ],
        "field_notes": {
            "payload": "Original command JSON payload, preserved without mutation.",
            "files": "Agent-friendly file/output map, inferred from payload.files or payload.outputs when possible.",
            "counts": "Agent-friendly numeric/count summary, inferred from payload.counts when possible.",
            "checks": "Agent-friendly validation booleans, inferred from payload.checks when possible.",
            "warnings": "List of non-fatal issues.",
            "errors": "List of fatal issues or structured error records.",
        },
    }


def load_json(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).expanduser().resolve().read_text())


def write_json(path: str | Path, data: dict[str, Any]) -> Path:
    p = Path(path).expanduser().resolve()
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    return p


def emit_json(data: dict[str, Any]) -> None:
    print(json.dumps(data, indent=2, ensure_ascii=False, sort_keys=False))
