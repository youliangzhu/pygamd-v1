"""Path resolution helpers."""
from __future__ import annotations

from pathlib import Path
from typing import Any


def resolve_path(base_file_or_dir: str | Path, value: str | Path) -> Path:
    """Resolve ``value`` relative to a config file or directory.

    If ``base_file_or_dir`` is a file, its parent is used as the base.
    """
    p = Path(value).expanduser()
    if p.is_absolute():
        return p.resolve()
    base = Path(base_file_or_dir).expanduser()
    if base.suffix and not base.is_dir():
        base = base.parent
    return (base / p).resolve()


def absolutize_config_paths(obj: Any, base_file_or_dir: str | Path, keys: set[str]) -> Any:
    """Recursively convert selected path-valued dict keys to absolute strings."""
    if isinstance(obj, dict):
        out: dict[str, Any] = {}
        for key, value in obj.items():
            if key in keys and isinstance(value, str) and value:
                out[key] = str(resolve_path(base_file_or_dir, value))
            else:
                out[key] = absolutize_config_paths(value, base_file_or_dir, keys)
        return out
    if isinstance(obj, list):
        return [absolutize_config_paths(item, base_file_or_dir, keys) for item in obj]
    return obj
