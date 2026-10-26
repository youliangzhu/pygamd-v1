"""Multiscale cg2aa package."""

from __future__ import annotations

from typing import Any

__all__ = ["run_cg2aa"]


def __getattr__(name: str) -> Any:
    if name == "run_cg2aa":
        from .runner import run_cg2aa
        return run_cg2aa
    raise AttributeError(name)
