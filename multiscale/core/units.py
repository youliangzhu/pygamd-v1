"""Unit conversion helpers shared by Multiscale modules."""
from __future__ import annotations

import math
from collections.abc import Iterable

COULOMB_GALAMOST_CHARGE_FACTOR = math.sqrt(138.935456)


def normalize_length_unit(unit: str | None, *, default: str = "angstrom") -> str:
    value = (unit or default).strip().lower()
    aliases = {
        "a": "angstrom",
        "å": "angstrom",
        "angstroms": "angstrom",
        "ang": "angstrom",
        "nanometer": "nm",
        "nanometers": "nm",
        "nanometre": "nm",
        "nanometres": "nm",
    }
    return aliases.get(value, value)


def length_to_nm(value: float, unit: str | None) -> float:
    normalized = normalize_length_unit(unit)
    if normalized == "angstrom":
        return float(value) * 0.1
    if normalized == "nm":
        return float(value)
    raise ValueError(f"Unsupported length unit: {unit!r}")


def length_from_nm(value: float, unit: str | None) -> float:
    normalized = normalize_length_unit(unit)
    if normalized == "angstrom":
        return float(value) * 10.0
    if normalized == "nm":
        return float(value)
    raise ValueError(f"Unsupported length unit: {unit!r}")


def vector_to_nm(values: Iterable[float], unit: str | None) -> tuple[float, ...]:
    return tuple(length_to_nm(float(v), unit) for v in values)
