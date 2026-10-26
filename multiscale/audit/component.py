"""Public self-contained audit helpers for Multiscale."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from collections.abc import Iterable
from typing import Any
import json
import re

from multiscale import __version__

COMPONENTS: tuple[dict[str, Any], ...] = (
    {
        "component": "aa2cg",
        "title": "AA-to-CG mapping, exporters, and RDF",
        "packages": ["multiscale/aa2cg"],
        "commands": ["multiscale aa2cg run", "multiscale aa2cg validate-config", "multiscale aa2cg export-pygamd", "multiscale aa2cg export-gromacs", "multiscale aa2cg export-lammps"],
        "status": "implemented",
    },
    {
        "component": "cg2aa",
        "title": "CG-to-AA backmapping",
        "packages": ["multiscale/cg2aa"],
        "commands": ["multiscale cg2aa run", "multiscale cg2aa summary", "multiscale cg2aa validate"],
        "status": "implemented",
    },
    {
        "component": "forcefield",
        "title": "Bundled OPLS force-field helpers",
        "packages": ["multiscale/forcefield"],
        "commands": ["multiscale forcefield list", "multiscale forcefield validate"],
        "status": "implemented",
    },
    {
        "component": "core",
        "title": "Shared CGSystem, GALAMOST XML, runtime, and report helpers",
        "packages": ["multiscale/core"],
        "commands": ["multiscale core info", "multiscale core inspect-pygamd", "multiscale report schema"],
        "status": "implemented",
    },
    {
        "component": "cleanup",
        "title": "Source/runtime cleanup",
        "packages": ["multiscale/clean.py"],
        "commands": ["multiscale clean-cache"],
        "status": "implemented",
    },
)

# Operational dependency markers that should not appear in the public
# self-contained runtime.
DEPENDENCY_PATTERNS: tuple[str, ...] = (
    "delegate_aa2cg",
    "delegate_cg2aa",
    "from multiscale._delegate",
    "run_cg2aa.py",
    "cg2aa.forcefield",
)


@dataclass(frozen=True)
class SourceFileRecord:
    path: str
    relpath: str
    bytes: int
    suffix: str
    category: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _safe_size(path: Path) -> int:
    try:
        return path.stat().st_size
    except OSError:
        return 0


def _iter_files(root: Path) -> Iterable[Path]:
    skip_dirs = {".git", "__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache", "build", "dist", "htmlcov"}
    for path in root.rglob("*"):
        if any(part in skip_dirs for part in path.parts):
            continue
        if path.is_file():
            yield path


def classify_source_file(relpath: str) -> dict[str, str]:
    rel = relpath.replace("\\", "/")
    suffix = Path(rel).suffix.lower()
    if rel.startswith("multiscale/") and suffix == ".py":
        return {"category": "runtime_python"}
    if rel.startswith("tests/") and suffix == ".py":
        return {"category": "test_python"}
    if rel.startswith("docs/") or rel in {"README.md", "NOTICE.md", "LICENSE"}:
        return {"category": "documentation_or_license"}
    if rel.startswith("scripts/"):
        return {"category": "script"}
    if rel == "pyproject.toml" or rel == "MANIFEST.in":
        return {"category": "packaging_metadata"}
    return {"category": "data_or_unknown"}


def source_file_inventory(source_root: str | Path | None = None) -> dict[str, Any]:
    """Return a source-tree inventory. Public audit API."""
    root = Path(source_root).expanduser().resolve() if source_root else Path(__file__).resolve().parents[2]
    files: list[SourceFileRecord] = []
    if root.exists():
        for path in _iter_files(root):
            rel = path.relative_to(root).as_posix()
            files.append(SourceFileRecord(
                path=str(path),
                relpath=rel,
                bytes=_safe_size(path),
                suffix=path.suffix,
                category=classify_source_file(rel)["category"],
            ))
    return {
        "ok": root.exists(),
        "source_root": str(root),
        "files_total": len(files),
        "bytes_total": sum(item.bytes for item in files),
        "files": [item.to_dict() for item in files],
    }


def source_project_summary(source_root: str | Path | None = None) -> dict[str, Any]:
    """Return a compact source-tree summary. Public audit API."""
    inv = source_file_inventory(source_root)
    categories: dict[str, int] = {}
    for item in inv["files"]:
        cat = item["category"]
        categories[cat] = categories.get(cat, 0) + 1
    return {
        "ok": inv["ok"],
        "source_root": inv["source_root"],
        "files_total": inv["files_total"],
        "bytes_total": inv["bytes_total"],
        "categories": categories,
    }


def _line_is_effective_dependency(line: str, pattern: str) -> bool:
    text = line.strip()
    if not text or text.startswith(("#", "\'\'\'", "\"\"\"", "\"", "\'")):
        return False
    if text.startswith(("if pattern in", "return pattern in")):
        return False
    if "DEPENDENCY_PATTERNS" in text:
        return False
    if "not_from_external" in text or "self_contained" in text:
        return False
    if "schema" in text and pattern in text:
        return False
    return pattern in text
def _scan_text_file(path: Path) -> list[dict[str, Any]]:
    findings: list[dict[str, Any]] = []
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return findings
    except OSError:
        return findings
    for lineno, line in enumerate(text.splitlines(), start=1):
        for pattern in DEPENDENCY_PATTERNS:
            if pattern in line and _line_is_effective_dependency(line, pattern):
                findings.append({
                    "path": str(path),
                    "line": lineno,
                    "pattern": pattern,
                    "text": line.strip()[:240],
                })
    if re.search(r"^\s*(from|import)\s+cg2aa(\.|\s|$)", text, flags=re.MULTILINE):
        findings.append({"path": str(path), "line": 0, "pattern": "top_level_cg2aa_import", "text": "top-level cg2aa import"})
    return findings


def scan_external_dependency_references(source_root: str | Path | None = None) -> dict[str, Any]:
    """Scan the public Multiscale source tree for operational external-runtime references."""
    root = Path(source_root).expanduser().resolve() if source_root else Path(__file__).resolve().parents[1]
    paths = sorted(root.rglob("*.py")) if root.exists() else []
    findings: list[dict[str, Any]] = []
    for path in paths:
        findings.extend(_scan_text_file(path))
    files = sorted({item["path"] for item in findings})
    return {
        "ok": True,
        "source_root": str(root),
        "dependency_reference_count": len(findings),
        "dependency_files": files,
        "findings": findings,
    }


def self_contained_status() -> dict[str, Any]:
    refs = scan_external_dependency_references()
    complete = refs["dependency_reference_count"] == 0
    blockers = []
    if refs["dependency_reference_count"]:
        blockers.append("Multiscale source still contains operational external-runtime dependency references")
    return {
        "ok": True,
        "self_contained_complete": complete,
        "self_contained_default": complete,
        "dependency_reference_count": refs["dependency_reference_count"],
        "dependency_files": refs["dependency_files"],
        "blockers": blockers,
        "release_status": "public-self-contained" if complete else "needs-cleanup",
        "version": __version__,
    }



def build_component_matrix() -> dict[str, Any]:
    refs = scan_external_dependency_references()
    status = self_contained_status()
    return {
        "ok": True,
        "schema": "multiscale.component_matrix.v1",
        "version": __version__,
        "status": status,
        "components": list(COMPONENTS),
        "dependency_references": refs,
        "notes": [
            "Multiscale is expected to use internal package backends by default.",
            "The audit matrix is a public component inventory, not an internal development plan.",
        ],
    }


def _md_table(rows: list[list[str]]) -> str:
    if not rows:
        return ""
    lines = ["| " + " | ".join(rows[0]) + " |", "| " + " | ".join("---" for _ in rows[0]) + " |"]
    for row in rows[1:]:
        lines.append("| " + " | ".join(str(cell).replace("|", "\\|") for cell in row) + " |")
    return "\n".join(lines)


def matrix_to_markdown(matrix: dict[str, Any]) -> str:
    status = matrix["status"]
    lines: list[str] = ["# Multiscale component matrix", ""]
    lines.append("## Self-contained status")
    lines.append("")
    lines.append(_md_table([
        ["Check", "Value"],
        ["self_contained_complete", str(status["self_contained_complete"])],
        ["dependency_reference_count", str(status["dependency_reference_count"])],
        ["release_status", str(status["release_status"])],
        ["version", str(status["version"])],
    ]))
    lines.append("")
    lines.append("## Components")
    lines.append("")
    rows = [["Component", "Title", "Packages", "Commands", "Status"]]
    for item in matrix["components"]:
        rows.append([
            item["component"],
            item["title"],
            ", ".join(item["packages"]),
            ", ".join(item["commands"]),
            item["status"],
        ])
    lines.append(_md_table(rows))
    lines.append("")
    refs = matrix["dependency_references"]["findings"]
    lines.append("## External-runtime dependency references")
    lines.append("")
    if not refs:
        lines.append("No operational external-runtime dependency references found.")
    else:
        for item in refs[:120]:
            lines.append(f"- `{item['path']}:{item['line']}` `{item['pattern']}` — `{item['text']}`")
        if len(refs) > 120:
            lines.append(f"- ... {len(refs) - 120} more references")
    lines.append("")
    return "\n".join(lines)


def write_markdown_matrix(path: str | Path, matrix: dict[str, Any] | None = None) -> Path:
    if matrix is None:
        matrix = build_component_matrix()
    p = Path(path).expanduser().resolve()
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(matrix_to_markdown(matrix) + "\n", encoding="utf-8")
    return p


def write_json(path: str | Path, data: dict[str, Any]) -> Path:
    p = Path(path).expanduser().resolve()
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return p
