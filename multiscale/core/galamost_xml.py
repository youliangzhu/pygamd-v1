"""GALAMOST XML reader/writer helpers for Multiscale core.

The functions are deliberately conservative and preserve Multiscale CGSystem behavior:

* GALAMOST topology XML indices are 0-based in ``bond``/``angle`` sections.
* CGSystem indices remain 1-based internally.
* Coordinates in GALAMOST XML are centered and expressed in nm.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET

from .cg_system import CGSystem
from .units import COULOMB_GALAMOST_CHARGE_FACTOR, length_to_nm


@dataclass
class GalamostXml:
    path: str | None
    timestep: int
    dimensions: int
    natoms: int
    box: tuple[float, float, float, float, float, float]
    sections: dict[str, list[list[str]]] = field(default_factory=dict)

    def rows(self, section: str) -> list[list[str]]:
        return self.sections.get(section, [])

    @property
    def positions(self) -> list[tuple[float, float, float]]:
        return [tuple(float(v) for v in row[:3]) for row in self.rows("position") if len(row) >= 3]

    @property
    def images(self) -> list[tuple[int, int, int]]:
        return [tuple(int(v) for v in row[:3]) for row in self.rows("image") if len(row) >= 3]

    @property
    def masses(self) -> list[float]:
        return [float(row[0]) for row in self.rows("mass") if row]

    @property
    def types(self) -> list[str]:
        return [str(row[0]) for row in self.rows("type") if row]

    @property
    def charges(self) -> list[float]:
        return [float(row[0]) for row in self.rows("charge") if row]

    @property
    def bonds(self) -> list[tuple[str, int, int]]:
        out: list[tuple[str, int, int]] = []
        for row in self.rows("bond"):
            if len(row) >= 3:
                out.append((str(row[0]), int(row[1]), int(row[2])))
        return out

    @property
    def angles(self) -> list[tuple[str, int, int, int]]:
        out: list[tuple[str, int, int, int]] = []
        for row in self.rows("angle"):
            if len(row) >= 4:
                out.append((str(row[0]), int(row[1]), int(row[2]), int(row[3])))
        return out

    def counts(self) -> dict[str, Any]:
        return {
            "natoms": self.natoms,
            "box": [self.box[0], self.box[1], self.box[2]],
            "positions": len(self.rows("position")),
            "images": len(self.rows("image")),
            "masses": len(self.rows("mass")),
            "types": len(self.rows("type")),
            "charges": len(self.rows("charge")),
            "bonds": len(self.rows("bond")),
            "angles": len(self.rows("angle")),
            "unique_types": sorted(set(self.types)),
        }

    def consistency_checks(self) -> dict[str, bool]:
        counts = self.counts()
        return {
            "positions_match_natoms": counts["positions"] == self.natoms,
            "images_match_natoms_or_absent": counts["images"] in (0, self.natoms),
            "masses_match_natoms": counts["masses"] in (0, self.natoms),
            "types_match_natoms": counts["types"] == self.natoms,
            "charges_match_natoms": counts["charges"] in (0, self.natoms),
            "has_box": bool(self.box[:3]),
        }

    def to_dict(self) -> dict[str, Any]:
        checks = self.consistency_checks()
        return {
            "path": self.path,
            "timestep": self.timestep,
            "dimensions": self.dimensions,
            "natoms": self.natoms,
            "box": list(self.box),
            "counts": self.counts(),
            "checks": checks,
            "ok": all(checks.values()),
        }


def _rows(node: ET.Element | None) -> list[list[str]]:
    if node is None or node.text is None:
        return []
    rows: list[list[str]] = []
    for line in node.text.splitlines():
        fields = line.split()
        if fields:
            rows.append(fields)
    return rows


def read_galamost_xml(path: str | Path) -> GalamostXml:
    path = Path(path)
    root = ET.parse(path).getroot()
    conf = root.find("configuration")
    if conf is None:
        raise ValueError(f"No <configuration> found in {path}")

    natoms = int(conf.attrib["natoms"])
    timestep = int(conf.attrib.get("time_step", 0))
    dimensions = int(conf.attrib.get("dimensions", 3))

    box_node = conf.find("box")
    if box_node is None:
        raise ValueError(f"No <box> found in {path}")
    box = (
        float(box_node.attrib.get("lx", 0.0)),
        float(box_node.attrib.get("ly", 0.0)),
        float(box_node.attrib.get("lz", 0.0)),
        float(box_node.attrib.get("xy", 0.0)),
        float(box_node.attrib.get("xz", 0.0)),
        float(box_node.attrib.get("yz", 0.0)),
    )
    sections = {child.tag: _rows(child) for child in conf}
    return GalamostXml(str(path.resolve()), timestep, dimensions, natoms, box, sections)


def inspect_galamost_xml(path: str | Path) -> dict[str, Any]:
    xml = read_galamost_xml(path)
    return {"ok": True, "galamost": xml.to_dict()}


def compare_galamost_counts(reference_xml: str | Path, generated_xml: str | Path) -> dict[str, Any]:
    reference = read_galamost_xml(reference_xml)
    generated = read_galamost_xml(generated_xml)
    ref_counts = reference.counts()
    gen_counts = generated.counts()
    checks = {
        "natoms_match": gen_counts.get("natoms") == ref_counts.get("natoms"),
        "positions_match": gen_counts.get("positions") == ref_counts.get("positions"),
        "images_match": gen_counts.get("images") == ref_counts.get("images"),
        "masses_match": gen_counts.get("masses") == ref_counts.get("masses"),
        "types_match": gen_counts.get("types") == ref_counts.get("types"),
        "charges_match": gen_counts.get("charges") == ref_counts.get("charges"),
        "bonds_match": gen_counts.get("bonds") == ref_counts.get("bonds"),
        "angles_match": gen_counts.get("angles") == ref_counts.get("angles"),
        "box_match": gen_counts.get("box") == ref_counts.get("box"),
        "unique_types_match": gen_counts.get("unique_types") == ref_counts.get("unique_types"),
    }
    return {
        "ok": all(checks.values()),
        "reference": {"path": str(Path(reference_xml).resolve()), "counts": ref_counts},
        "generated": {"path": str(Path(generated_xml).resolve()), "counts": gen_counts},
        "checks": checks,
    }


def _section(lines: list[str], tag: str, values: list[str]) -> None:
    lines.append(f'<{tag} num="{len(values)}">')
    lines.extend(values)
    lines.append(f'</{tag}>')


def _type_sort_key(value: str | None) -> tuple[str, int]:
    if not value:
        return ("", 0)
    prefix, sep, suffix = value.rpartition(":")
    if sep:
        try:
            return (prefix, int(suffix))
        except ValueError:
            return (prefix, 0)
    return (value, 0)


def _position_to_galamost_nm(system: CGSystem, position: tuple[float, float, float] | None) -> tuple[float, float, float]:
    if position is None:
        return (0.0, 0.0, 0.0)
    if system.box:
        centered = tuple(float(p) - float(b) / 2.0 for p, b in zip(position, system.box))
    else:
        centered = tuple(float(p) for p in position)
    return tuple(length_to_nm(v, system.coordinate_unit) for v in centered)  # type: ignore[return-value]


def _box_to_galamost_nm(system: CGSystem) -> tuple[float, float, float]:
    if not system.box:
        return (0.0, 0.0, 0.0)
    return tuple(length_to_nm(float(v), system.box_unit) for v in system.box)  # type: ignore[return-value]


def render_galamost_xml(system: CGSystem, *, eps: float = 1.0) -> str:
    """Render an Multiscale CGSystem as GALAMOST topology XML."""
    if eps <= 0:
        raise ValueError("eps must be positive")
    lx, ly, lz = _box_to_galamost_nm(system)
    charge_factor = COULOMB_GALAMOST_CHARGE_FACTOR / (eps ** 0.5)

    lines: list[str] = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<galamost_xml version="1.3">',
        f'<configuration time_step="0" dimensions="3" natoms="{len(system.beads)}" >',
        f'<box lx="{lx:.5f}" ly="{ly:.5f}" lz="{lz:.5f}"/>',
    ]

    position_values: list[str] = []
    mass_values: list[str] = []
    type_values: list[str] = []
    charge_values: list[str] = []
    for bead in system.beads:
        x, y, z = _position_to_galamost_nm(system, bead.position)
        position_values.append(f'{x:11.5f} {y:11.5f} {z:11.5f}')
        mass_values.append(f'{float(bead.mass or 0.0):11.5f}')
        type_values.append(str(bead.bead_type or bead.name))
        charge_values.append(f'{float(bead.charge or 0.0) * charge_factor:11.5f}')

    _section(lines, "position", position_values)
    _section(lines, "mass", mass_values)
    _section(lines, "type", type_values)
    _section(lines, "charge", charge_values)

    lines.append("<bond>")
    for bond in sorted(system.bonds, key=lambda b: (_type_sort_key(b.bond_type), b.i, b.j)):
        lines.append(f'{bond.bond_type or "bond"} {bond.i - system.index_base} {bond.j - system.index_base}')
    lines.append("</bond>")

    lines.append("<angle>")
    for angle in sorted(system.angles, key=lambda a: (_type_sort_key(a.angle_type), a.i, a.j, a.k)):
        lines.append(f'{angle.angle_type or "angle"} {angle.i - system.index_base} {angle.j - system.index_base} {angle.k - system.index_base}')
    lines.append("</angle>")
    lines.append("</configuration>")
    lines.append("</galamost_xml>")
    return "\n".join(lines) + "\n"


def write_galamost_xml(system: CGSystem, path: str | Path, *, eps: float = 1.0) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(render_galamost_xml(system, eps=eps), encoding="utf-8")
    return out
