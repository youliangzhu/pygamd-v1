"""GROMACS exporter scaffold for Multiscale CGSystem.

This exporter intentionally focuses on a coordinate file (`cg.gro`) and a
minimal topology scaffold (`cg.itp`, `topol.top`) from an existing CGSystem.
It does not generate nonbonded parameters, bonded parameters, IBI/IMC data, or
potential tables.  Parameter assignment is left to a future potential module or
user-provided templates.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from collections.abc import Iterable
from typing import Any

from multiscale.core.cg_system import CGBead, CGSystem


@dataclass(frozen=True)
class GromacsExportReport:
    outdir: str
    files: dict[str, str]
    counts: dict[str, Any]
    checks: dict[str, bool | None]

    @property
    def ok(self) -> bool:
        return all(value is not False for value in self.checks.values())

    def to_dict(self) -> dict[str, Any]:
        return {
            "outdir": self.outdir,
            "files": self.files,
            "counts": self.counts,
            "checks": self.checks,
            "ok": self.ok,
        }


def _length_to_nm(value: float, unit: str | None) -> float:
    normalized = (unit or "angstrom").lower()
    if normalized in {"angstrom", "angstroms", "a", "å"}:
        return float(value) * 0.1
    if normalized in {"nm", "nanometer", "nanometers"}:
        return float(value)
    raise ValueError(f"Unsupported length unit for GROMACS export: {unit!r}")


def _safe_name(value: str | None, *, fallback: str, max_len: int | None = None) -> str:
    text = str(value or fallback).strip() or fallback
    text = text.replace(" ", "_")
    if max_len is not None and len(text) > max_len:
        return text[:max_len]
    return text


def _molecule_name(system: CGSystem) -> str:
    return _safe_name(system.name, fallback="CGMOL")


def _residue_name(system: CGSystem) -> str:
    mapping = system.metadata.get("mapping") or {}
    parent = mapping.get("bead_mapping", {}).get("parent_type") if isinstance(mapping, dict) else None
    return _safe_name(parent or system.name, fallback="CG", max_len=5)


def _box_to_nm(system: CGSystem) -> tuple[float, float, float]:
    if not system.box:
        return (0.0, 0.0, 0.0)
    return tuple(_length_to_nm(float(v), system.box_unit) for v in system.box)  # type: ignore[return-value]


def _position_to_nm(system: CGSystem, bead: CGBead) -> tuple[float, float, float]:
    if bead.position is None:
        return (0.0, 0.0, 0.0)
    return tuple(_length_to_nm(float(v), system.coordinate_unit) for v in bead.position)  # type: ignore[return-value]


def _template_beads(system: CGSystem) -> list[CGBead]:
    if system.beads_per_molecule:
        first_mol = [b for b in system.beads if b.molecule_index == 1]
        if first_mol:
            return sorted(first_mol, key=lambda b: b.molecule_bead_index or b.index)
        return sorted(system.beads, key=lambda b: b.index)[: system.beads_per_molecule]
    return sorted(system.beads, key=lambda b: b.index)


def _local_index(global_index: int, molecule_index: int | None, beads_per_molecule: int | None) -> int:
    if molecule_index and beads_per_molecule:
        return global_index - (molecule_index - 1) * beads_per_molecule
    return global_index


def _template_bonds(system: CGSystem) -> list[tuple[int, int, str | None]]:
    selected = [b for b in system.bonds if b.molecule_index == 1]
    if not selected and system.beads_per_molecule:
        selected = [b for b in system.bonds if b.i <= system.beads_per_molecule and b.j <= system.beads_per_molecule]
    if not selected:
        selected = list(system.bonds)
    return [
        (
            _local_index(b.i, b.molecule_index, system.beads_per_molecule),
            _local_index(b.j, b.molecule_index, system.beads_per_molecule),
            b.bond_type,
        )
        for b in selected
    ]


def _template_angles(system: CGSystem) -> list[tuple[int, int, int, str | None]]:
    selected = [a for a in system.angles if a.molecule_index == 1]
    if not selected and system.beads_per_molecule:
        selected = [
            a
            for a in system.angles
            if a.i <= system.beads_per_molecule and a.j <= system.beads_per_molecule and a.k <= system.beads_per_molecule
        ]
    if not selected:
        selected = list(system.angles)
    return [
        (
            _local_index(a.i, a.molecule_index, system.beads_per_molecule),
            _local_index(a.j, a.molecule_index, system.beads_per_molecule),
            _local_index(a.k, a.molecule_index, system.beads_per_molecule),
            a.angle_type,
        )
        for a in selected
    ]


def _gro_atom_line(system: CGSystem, bead: CGBead) -> str:
    resnr = int(bead.molecule_index or 1) % 100000
    atomnr = int(bead.index) % 100000
    resname = _residue_name(system)
    atomname = _safe_name(bead.name, fallback=f"B{bead.index}", max_len=5)
    x, y, z = _position_to_nm(system, bead)
    return f"{resnr:5d}{resname:<5}{atomname:>5}{atomnr:5d}{x:8.3f}{y:8.3f}{z:8.3f}"


def render_gro(system: CGSystem) -> str:
    lx, ly, lz = _box_to_nm(system)
    stage = str((system.metadata or {}).get("stage") or "")
    if "step7" in stage.lower():
        title = f"Generated by Multiscale AA-to-CG structure mapping: {system.name}"
    else:
        title = f"Generated by Multiscale GROMACS scaffold: {system.name}"
    lines = [
        title,
        f"{len(system.beads):5d}",
    ]
    lines.extend(_gro_atom_line(system, bead) for bead in sorted(system.beads, key=lambda b: b.index))
    lines.append(f"{lx:10.5f}{ly:10.5f}{lz:10.5f}")
    return "\n".join(lines) + "\n"


def render_itp(system: CGSystem) -> str:
    molname = _molecule_name(system)
    resname = _residue_name(system)
    beads = _template_beads(system)
    bonds = _template_bonds(system)
    angles = _template_angles(system)

    lines: list[str] = [
        "; Generated by Multiscale GROMACS scaffold.",
        "; This file intentionally does not contain force-field parameters.",
        "; Fill atomtypes, nonbonded, bond and angle parameters externally before production MD.",
        "",
        "[ moleculetype ]",
        "; name  nrexcl",
        f"{molname:<16} 3",
        "",
        "[ atoms ]",
        "; nr  type  resnr residue atom  cgnr  charge       mass",
    ]
    for local_idx, bead in enumerate(beads, start=1):
        atom_type = _safe_name(bead.bead_type or bead.name, fallback=f"B{local_idx}")
        atom_name = _safe_name(bead.name, fallback=f"B{local_idx}", max_len=8)
        charge = float(bead.charge or 0.0)
        mass = float(bead.mass or 0.0)
        lines.append(f"{local_idx:5d} {atom_type:<8} {1:5d} {resname:<6} {atom_name:<8} {local_idx:5d} {charge:12.6f} {mass:12.6f}")

    lines.extend(["", "[ bonds ]", "; ai    aj  funct  ; original_type; parameters are not generated"])
    for i, j, bond_type in bonds:
        lines.append(f"{i:5d} {j:5d} {1:5d}  ; {bond_type or 'bond'}")

    lines.extend(["", "[ angles ]", "; ai    aj    ak  funct  ; original_type; parameters are not generated"])
    for i, j, k, angle_type in angles:
        lines.append(f"{i:5d} {j:5d} {k:5d} {1:5d}  ; {angle_type or 'angle'}")

    return "\n".join(lines) + "\n"


def render_top(system: CGSystem, *, itp_name: str = "cg.itp") -> str:
    molname = _molecule_name(system)
    nmol = int(system.nmol or 1)
    return (
        "; Generated by Multiscale GROMACS scaffold.\n"
        "; This topology is a structural scaffold only; it is not parameter-complete.\n\n"
        "[ defaults ]\n"
        "; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ\n"
        "1         1          no         1.0      1.0\n\n"
        f'#include "{itp_name}"\n\n'
        "[ system ]\n"
        f"{system.name} coarse-grained system\n\n"
        "[ molecules ]\n"
        "; name            number\n"
        f"{molname:<16} {nmol}\n"
    )


def write_gromacs_files(system: CGSystem, outdir: str | Path) -> dict[str, Path]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = {
        "gro": outdir / "cg.gro",
        "itp": outdir / "cg.itp",
        "top": outdir / "topol.top",
    }
    paths["gro"].write_text(render_gro(system), encoding="utf-8")
    paths["itp"].write_text(render_itp(system), encoding="utf-8")
    paths["top"].write_text(render_top(system), encoding="utf-8")
    return paths


def _count_section_lines(lines: Iterable[str], section_name: str) -> int:
    active = False
    count = 0
    wanted = section_name.strip().lower()
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith("[") and line.endswith("]"):
            active = line.strip("[] ").lower() == wanted
            continue
        if active:
            count += 1
    return count


def parse_gromacs_export_counts(outdir: str | Path) -> dict[str, Any]:
    outdir = Path(outdir)
    counts: dict[str, Any] = {
        "gro_natoms": None,
        "gro_has_box": False,
        "itp_atoms": 0,
        "itp_bonds": 0,
        "itp_angles": 0,
        "top_molecules": 0,
        "has_gro": int((outdir / "cg.gro").exists()),
        "has_itp": int((outdir / "cg.itp").exists()),
        "has_top": int((outdir / "topol.top").exists()),
    }

    gro = outdir / "cg.gro"
    if gro.exists():
        gro_lines = gro.read_text(encoding="utf-8", errors="replace").splitlines()
        if len(gro_lines) >= 3:
            try:
                counts["gro_natoms"] = int(gro_lines[1].strip())
            except ValueError:
                counts["gro_natoms"] = None
            counts["gro_has_box"] = bool(gro_lines[-1].split())

    itp = outdir / "cg.itp"
    if itp.exists():
        itp_lines = itp.read_text(encoding="utf-8", errors="replace").splitlines()
        counts["itp_atoms"] = _count_section_lines(itp_lines, "atoms")
        counts["itp_bonds"] = _count_section_lines(itp_lines, "bonds")
        counts["itp_angles"] = _count_section_lines(itp_lines, "angles")

    top = outdir / "topol.top"
    if top.exists():
        top_lines = top.read_text(encoding="utf-8", errors="replace").splitlines()
        counts["top_molecules"] = _count_section_lines(top_lines, "molecules")
    return counts


def export_gromacs_scaffold(system: CGSystem, outdir: str | Path) -> GromacsExportReport:
    paths = write_gromacs_files(system, outdir)
    counts = parse_gromacs_export_counts(outdir)
    expected_template_atoms = int(system.beads_per_molecule or len(system.beads))
    expected_template_bonds = len(_template_bonds(system))
    expected_template_angles = len(_template_angles(system))
    checks: dict[str, bool | None] = {
        "has_gro": counts["has_gro"] == 1,
        "has_itp": counts["has_itp"] == 1,
        "has_top": counts["has_top"] == 1,
        "gro_natoms_match_cgsystem": counts["gro_natoms"] == len(system.beads),
        "gro_has_box": bool(counts["gro_has_box"]),
        "itp_atoms_match_template": counts["itp_atoms"] == expected_template_atoms,
        "itp_bonds_match_template": counts["itp_bonds"] == expected_template_bonds,
        "itp_angles_match_template": counts["itp_angles"] == expected_template_angles,
        "top_has_one_molecule_entry": counts["top_molecules"] == 1,
        "no_forcefield_parameters_generated": True,
    }
    return GromacsExportReport(
        outdir=str(Path(outdir)),
        files={key: str(path) for key, path in paths.items()},
        counts=counts,
        checks=checks,
    )
