"""PyGAMD auxiliary-file writer for Multiscale CGSystem.

This writer intentionally focuses on the small Python include files that
connect an existing topology XML to PyGAMD exclusion rules and table-potential file
references.  This module does not generate or modify potential table data.
"""
from __future__ import annotations

from dataclasses import dataclass
import ast
import re
from pathlib import Path
from collections.abc import Iterable
from typing import Any

from multiscale.core.cg_system import CGSystem


@dataclass(frozen=True)
class TablePotentialRef:
    family: str  # pair, bond, angle
    types: tuple[str, ...]
    filename: str
    rmin: float | int = 0
    rmax: float | int = 1

    def to_dict(self) -> dict[str, Any]:
        return {
            "family": self.family,
            "types": list(self.types),
            "filename": self.filename,
            "rmin": self.rmin,
            "rmax": self.rmax,
        }


@dataclass(frozen=True)
class GalamostAuxExportReport:
    outdir: str
    files: dict[str, str]
    counts: dict[str, int]
    reference_counts: dict[str, int] | None
    checks: dict[str, bool | None]

    @property
    def ok(self) -> bool:
        return all(value is not False for value in self.checks.values())

    def to_dict(self) -> dict[str, Any]:
        return {
            "outdir": self.outdir,
            "files": self.files,
            "counts": self.counts,
            "reference_counts": self.reference_counts,
            "checks": self.checks,
            "ok": self.ok,
        }


_FAMILY_BY_VAR = {"pairs": "pair", "bonds": "bond", "angles": "angle"}
_VAR_BY_FAMILY = {"pair": "pairs", "bond": "bonds", "angle": "angles"}


def _read_text(path: str | Path) -> str:
    return Path(path).read_text(encoding="utf-8", errors="replace")


def _normalize_text(text: str) -> str:
    return "\n".join(line.rstrip() for line in text.strip().splitlines())


def _format_number(value: float | int) -> str:
    try:
        value_float = float(value)
    except (TypeError, ValueError):
        return str(value)
    if value_float.is_integer():
        return str(int(value_float))
    return f"{value_float:g}"


def _table_filename(types: Iterable[str]) -> str:
    return "table_" + "-".join(str(t).replace(":", "-") for t in types) + ".dat"


def _type_sort_key(value: str) -> tuple[str, int]:
    prefix, sep, suffix = value.rpartition(":")
    if sep:
        try:
            return (prefix, int(suffix))
        except ValueError:
            return (prefix, 0)
    return (value, 0)


def parse_tables_inc(path: str | Path) -> list[TablePotentialRef]:
    """Parse GALAMOST setPotential references from tables.inc.py."""
    refs: list[TablePotentialRef] = []
    pattern = re.compile(r"^(pairs|bonds|angles)\.setPotential\((.*)\)\s*$")
    for raw in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        match = pattern.match(line)
        if not match:
            continue
        var_name, args_text = match.groups()
        family = _FAMILY_BY_VAR[var_name]
        args = ast.literal_eval(f"({args_text},)")
        if family == "pair":
            if len(args) < 5:
                raise ValueError(f"Invalid pair setPotential line: {line!r}")
            refs.append(TablePotentialRef(family, (str(args[0]), str(args[1])), str(args[2]), args[3], args[4]))
        else:
            if len(args) < 4:
                raise ValueError(f"Invalid {family} setPotential line: {line!r}")
            refs.append(TablePotentialRef(family, (str(args[0]),), str(args[1]), args[2], args[3]))
    return refs


def summarize_table_refs(refs: Iterable[TablePotentialRef]) -> dict[str, int]:
    counts = {"pair_potentials": 0, "bond_potentials": 0, "angle_potentials": 0}
    for ref in refs:
        if ref.family == "pair":
            counts["pair_potentials"] += 1
        elif ref.family == "bond":
            counts["bond_potentials"] += 1
        elif ref.family == "angle":
            counts["angle_potentials"] += 1
    return counts


def parse_galamost_aux_counts(directory: str | Path) -> dict[str, int]:
    directory = Path(directory)
    refs = parse_tables_inc(directory / "tables.inc.py") if (directory / "tables.inc.py").exists() else []
    counts = summarize_table_refs(refs)
    counts["has_exclusionsSR"] = int((directory / "exclusionsSR.inc.py").exists())
    counts["has_exclusionsEL"] = int((directory / "exclusionsEL.inc.py").exists())
    counts["has_tables"] = int((directory / "tables.inc.py").exists())
    return counts


def _infer_pair_refs(system: CGSystem) -> list[TablePotentialRef]:
    # Preserve the bead-type order declared by the JSON model so generated
    # pair-table references remain deterministic.
    type_groups = system.metadata.get("bead_type_groups") or {}
    if isinstance(type_groups, dict) and type_groups:
        order = [str(key) for key in type_groups.keys()]
    else:
        order = system.summary().get("bead_types") or []
    refs: list[TablePotentialRef] = []
    for i, left in enumerate(order):
        for right in order[i:]:
            refs.append(TablePotentialRef("pair", (str(left), str(right)), _table_filename((str(left), str(right)))))
    return refs


def _infer_bond_refs(system: CGSystem) -> list[TablePotentialRef]:
    return [
        TablePotentialRef("bond", (bond_type,), _table_filename((bond_type,)))
        for bond_type in sorted(system.summary().get("bond_types") or [], key=_type_sort_key)
    ]


def _infer_angle_refs(system: CGSystem) -> list[TablePotentialRef]:
    return [
        TablePotentialRef("angle", (angle_type,), _table_filename((angle_type,)))
        for angle_type in sorted(system.summary().get("angle_types") or [], key=_type_sort_key)
    ]


def infer_table_refs(system: CGSystem, reference_tables: str | Path | None = None) -> list[TablePotentialRef]:
    """Return deterministic pair, bond, and angle table references."""
    if reference_tables and Path(reference_tables).exists():
        refs = parse_tables_inc(reference_tables)
        if refs:
            return refs
    return _infer_pair_refs(system) + _infer_bond_refs(system) + _infer_angle_refs(system)


def render_exclusions_sr() -> str:
    return (
        "neighbor_listSR = gala.NeighborList(all_info, rcut, rbuffer)\n"
        "neighbor_listSR.addExclusionsFromBonds()\n"
        "neighbor_listSR.addExclusionsFromAngles()\n"
    )


def render_exclusions_el() -> str:
    return "neighbor_listEL = neighbor_listSR\n"


def render_tables_inc(refs: Iterable[TablePotentialRef]) -> str:
    grouped = {"pair": [], "bond": [], "angle": []}
    for ref in refs:
        grouped[ref.family].append(ref)

    lines: list[str] = []
    if grouped["pair"]:
        lines.append("pairs = gala.PairForceTable(all_info, neighbor_listSR, npoints_in_table)")
        for ref in grouped["pair"]:
            lines.append(
                "pairs.setPotential("
                f"'{ref.types[0]}','{ref.types[1]}', '{ref.filename}', {_format_number(ref.rmin)} ,{_format_number(ref.rmax)})"
            )
        lines.append("app.add(pairs)")
    if grouped["bond"]:
        lines.append("bonds = gala.BondForceTable(all_info, npoints_in_table)")
        for ref in grouped["bond"]:
            lines.append(
                "bonds.setPotential("
                f"'{ref.types[0]}', '{ref.filename}', {_format_number(ref.rmin)} ,{_format_number(ref.rmax)})"
            )
        lines.append("app.add(bonds)")
    if grouped["angle"]:
        lines.append("angles = gala.AngleForceTable(all_info, npoints_in_table-1)")
        for ref in grouped["angle"]:
            lines.append(
                "angles.setPotential("
                f"'{ref.types[0]}', '{ref.filename}', {_format_number(ref.rmin)} ,{_format_number(ref.rmax)})"
            )
        lines.append("app.add(angles)")
    return "\n".join(lines) + "\n"


def _reference_dir_from_system(system: CGSystem) -> Path | None:
    source_root = system.metadata.get("source_root")
    if not source_root:
        return None
    root = Path(source_root)
    candidates = [root / "coor", root]
    for candidate in candidates:
        if (candidate / "tables.inc.py").exists():
            return candidate
    return None


def export_galamost_auxiliary(
    system: CGSystem,
    outdir: str | Path,
    *,
    reference_dir: str | Path | None = None,
) -> GalamostAuxExportReport:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    reference = Path(reference_dir).resolve() if reference_dir else _reference_dir_from_system(system)
    reference_tables = reference / "tables.inc.py" if reference else None

    refs = infer_table_refs(system, reference_tables)
    files = {
        "exclusionsSR": outdir / "exclusionsSR.inc.py",
        "exclusionsEL": outdir / "exclusionsEL.inc.py",
        "tables": outdir / "tables.inc.py",
    }
    files["exclusionsSR"].write_text(render_exclusions_sr(), encoding="utf-8")
    files["exclusionsEL"].write_text(render_exclusions_el(), encoding="utf-8")
    files["tables"].write_text(render_tables_inc(refs), encoding="utf-8")

    counts = parse_galamost_aux_counts(outdir)
    reference_counts = parse_galamost_aux_counts(reference) if reference and reference.exists() else None
    checks: dict[str, bool | None] = {
        "has_exclusionsSR": counts["has_exclusionsSR"] == 1,
        "has_exclusionsEL": counts["has_exclusionsEL"] == 1,
        "has_tables": counts["has_tables"] == 1,
        "has_pair_potentials": counts["pair_potentials"] > 0,
        "has_bond_potentials": counts["bond_potentials"] > 0,
        "has_angle_potentials": counts["angle_potentials"] > 0,
    }
    if reference_counts:
        checks.update(
            {
                "pair_potentials_match_reference": counts["pair_potentials"] == reference_counts["pair_potentials"],
                "bond_potentials_match_reference": counts["bond_potentials"] == reference_counts["bond_potentials"],
                "angle_potentials_match_reference": counts["angle_potentials"] == reference_counts["angle_potentials"],
            }
        )
    return GalamostAuxExportReport(
        outdir=str(outdir),
        files={key: str(path) for key, path in files.items()},
        counts=counts,
        reference_counts=reference_counts,
        checks=checks,
    )


def compare_galamost_auxiliary(reference_dir: str | Path, generated_dir: str | Path) -> dict[str, Any]:
    reference_dir = Path(reference_dir)
    generated_dir = Path(generated_dir)
    reference_counts = parse_galamost_aux_counts(reference_dir)
    generated_counts = parse_galamost_aux_counts(generated_dir)
    checks = {
        "pair_potentials_match": generated_counts["pair_potentials"] == reference_counts["pair_potentials"],
        "bond_potentials_match": generated_counts["bond_potentials"] == reference_counts["bond_potentials"],
        "angle_potentials_match": generated_counts["angle_potentials"] == reference_counts["angle_potentials"],
        "has_exclusionsSR": generated_counts["has_exclusionsSR"] == 1,
        "has_exclusionsEL": generated_counts["has_exclusionsEL"] == 1,
        "has_tables": generated_counts["has_tables"] == 1,
    }
    exact_checks: dict[str, bool | None] = {}
    for filename in ("exclusionsSR.inc.py", "exclusionsEL.inc.py", "tables.inc.py"):
        ref_file = reference_dir / filename
        gen_file = generated_dir / filename
        exact_checks[f"{filename}_normalized_text_match"] = (
            _normalize_text(_read_text(ref_file)) == _normalize_text(_read_text(gen_file))
            if ref_file.exists() and gen_file.exists()
            else None
        )
    checks.update(exact_checks)
    return {
        "reference_counts": reference_counts,
        "generated_counts": generated_counts,
        "checks": checks,
        "ok": all(value is not False for value in checks.values()),
    }
