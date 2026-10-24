"""Runtime discovery for the self-contained Multiscale package."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import os
import shutil
import sys
from collections.abc import Iterable
from typing import Any


@dataclass(frozen=True)
class RuntimeLayout:
    package_root: str | None
    project_root: str | None
    gmx_bin: str | None
    python_executable: str
    self_contained: bool = True

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _candidate_roots(start: Path | None = None) -> Iterable[Path]:
    seen: set[Path] = set()

    env_root = os.environ.get("MULTISCALE_ROOT")
    if env_root:
        p = Path(env_root).expanduser()
        if p not in seen:
            seen.add(p)
            yield p

    package_root = Path(__file__).resolve().parents[2]
    for p in [package_root, package_root.parent]:
        if p not in seen:
            seen.add(p)
            yield p

    cwd = Path.cwd().resolve()
    for p in [cwd, *cwd.parents]:
        if p not in seen:
            seen.add(p)
            yield p

    if start is not None:
        s = start.resolve()
        for p in [s, *s.parents]:
            if p not in seen:
                seen.add(p)
                yield p


def _looks_like_multiscale_project(root: Path) -> bool:
    return (root / "pyproject.toml").is_file() and (root / "multiscale" / "cli.py").is_file()


def _find_project_root(package_root: Path) -> str | None:
    for root in _candidate_roots(package_root):
        if _looks_like_multiscale_project(root) or root == package_root:
            try:
                return str(root.resolve())
            except OSError:
                return str(root)
    return str(package_root.parent.resolve())


def discover_runtime_layout() -> RuntimeLayout:
    package_root = Path(__file__).resolve().parents[2]
    gmx_bin = os.environ.get("GMX_BIN") or shutil.which("gmx")
    return RuntimeLayout(
        package_root=str(package_root.resolve()),
        project_root=_find_project_root(package_root),
        gmx_bin=gmx_bin,
        python_executable=os.environ.get("MULTISCALE_PYTHON") or sys.executable,
        self_contained=True,
    )


def prepend_pythonpath(env: dict[str, str], path: str | None) -> dict[str, str]:
    """Return ``env`` with ``path`` prepended to PYTHONPATH."""
    new_env = env.copy()
    if path:
        current = new_env.get("PYTHONPATH", "")
        new_env["PYTHONPATH"] = path if not current else f"{path}:{current}"
    return new_env
