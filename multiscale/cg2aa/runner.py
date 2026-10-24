"""Internal cg2aa runner for self-contained Multiscale."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from .pipeline import run_pipeline


def _with_backend(summary: dict[str, Any], backend: str) -> dict[str, Any]:
    data = dict(summary)
    data["backend"] = backend
    return data


def run_cg2aa(config: str | Path) -> dict[str, Any]:
    """Run the packaged pipeline and annotate the summary."""
    return _with_backend(run_pipeline(config), "internal")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="multiscale cg2aa run",
        description="Run the self-contained Multiscale cg2aa pipeline.",
    )
    parser.add_argument("config", help="Path to a JSON config.")
    parser.add_argument("--backend", choices=["internal"], default="internal")
    args = parser.parse_args(argv)
    summary = run_cg2aa(args.config)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
