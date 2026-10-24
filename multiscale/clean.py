"""Cache and generated-file cleanup helpers for Multiscale source/runtime trees."""

from __future__ import annotations

import argparse
import fnmatch
import json
import os
import shutil
from dataclasses import asdict, dataclass, field
from pathlib import Path

CACHE_DIR_NAMES: frozenset[str] = frozenset({
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    ".ruff_cache",
})

BUILD_DIR_NAMES: frozenset[str] = frozenset({
    "build",
    "dist",
    "htmlcov",
})

CACHE_FILE_PATTERNS: tuple[str, ...] = (
    "*.pyc",
    "*.pyo",
    "*.pyd",
    ".coverage",
    "coverage.xml",
    "run.log",
)

CACHE_DIR_PATTERNS: tuple[str, ...] = (
    "*.egg-info",
    "multiscale_validation_run_*",
)


@dataclass(frozen=True)
class CleanResult:
    """Structured result returned by ``multiscale clean-cache``."""

    ok: bool
    root: str
    dry_run: bool
    include_build: bool
    candidates: list[str] = field(default_factory=list)
    removed: list[str] = field(default_factory=list)
    errors: list[dict[str, str]] = field(default_factory=list)

    def to_dict(self) -> dict[str, object]:
        data = asdict(self)
        data["counts"] = {
            "candidates": len(self.candidates),
            "removed": len(self.removed),
            "errors": len(self.errors),
        }
        return data


def _matches_any(name: str, patterns: tuple[str, ...]) -> bool:
    return any(fnmatch.fnmatch(name, pattern) for pattern in patterns)


def _is_cleanable_dir(name: str, *, include_build: bool) -> bool:
    if name in CACHE_DIR_NAMES or _matches_any(name, CACHE_DIR_PATTERNS):
        return True
    return include_build and name in BUILD_DIR_NAMES


def _is_cleanable_file(name: str) -> bool:
    return _matches_any(name, CACHE_FILE_PATTERNS)


def collect_cleanable_paths(root: str | Path, *, include_build: bool = False) -> list[Path]:
    """Return cache/generated paths that can be removed safely.

    The default set targets Python/runtime noise only: bytecode caches,
    pytest caches, editable-install metadata, run logs, and Multiscale validation-run
    directories.  Build directories are only included when ``include_build`` is
    true so that a normal cleanup does not remove a user's intentionally built
    distribution artifacts.
    """

    base = Path(root).expanduser().resolve()
    if not base.exists():
        return []

    found: list[Path] = []
    if base.is_file():
        return [base] if _is_cleanable_file(base.name) else []

    for dirpath_str, dirnames, filenames in os.walk(base, topdown=True):
        dirpath = Path(dirpath_str)
        for dirname in tuple(dirnames):
            if _is_cleanable_dir(dirname, include_build=include_build):
                path = dirpath / dirname
                found.append(path)
                dirnames.remove(dirname)
        for filename in filenames:
            if _is_cleanable_file(filename):
                found.append(dirpath / filename)

    return sorted(found, key=lambda path: (len(path.parts), str(path)))


def clean_cache(root: str | Path = ".", *, dry_run: bool = False, include_build: bool = False) -> CleanResult:
    """Remove Multiscale cache/generated files below ``root`` and report the result."""

    base = Path(root).expanduser().resolve()
    if not base.exists():
        return CleanResult(
            ok=False,
            root=str(base),
            dry_run=dry_run,
            include_build=include_build,
            errors=[{"path": str(base), "error": "path does not exist"}],
        )

    candidates = collect_cleanable_paths(base, include_build=include_build)
    removed: list[str] = []
    errors: list[dict[str, str]] = []

    if not dry_run:
        for path in candidates:
            try:
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink(missing_ok=True)
                removed.append(str(path))
            except OSError as exc:
                errors.append({"path": str(path), "error": str(exc)})

    return CleanResult(
        ok=not errors,
        root=str(base),
        dry_run=dry_run,
        include_build=include_build,
        candidates=[str(path) for path in candidates],
        removed=removed,
        errors=errors,
    )


def _print_result(result: CleanResult, *, json_output: bool) -> None:
    data = result.to_dict()
    if json_output:
        print(json.dumps(data, indent=2, ensure_ascii=False))
        return

    counts = data["counts"]
    print(f"ok: {data['ok']}")
    print(f"root: {data['root']}")
    print(f"dry_run: {data['dry_run']}")
    print(f"candidates: {counts['candidates']}")
    print(f"removed: {counts['removed']}")
    if result.errors:
        print("errors:")
        for item in result.errors:
            print(f"  {item['path']}: {item['error']}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="multiscale clean-cache",
        description="Remove Python/runtime cache files and Multiscale validation-run artifacts from a source tree.",
    )
    parser.add_argument("path", nargs="?", default=".", help="source/runtime tree to clean; defaults to current directory")
    parser.add_argument("--dry-run", action="store_true", help="list cleanable paths without deleting them")
    parser.add_argument("--include-build", action="store_true", help="also remove build/, dist/, and htmlcov/ directories")
    parser.add_argument("--json", action="store_true", help="emit a structured JSON report")
    args = parser.parse_args(argv)

    result = clean_cache(args.path, dry_run=args.dry_run, include_build=args.include_build)
    _print_result(result, json_output=args.json)
    return 0 if result.ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
