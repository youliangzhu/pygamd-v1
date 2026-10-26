"""GALAMOST topology XML writers for Multiscale CGSystem.

The standard ``topology.xml`` preserves the mapped coordinates exactly as before.
The additional ``topology_inbox.xml`` contains coordinates wrapped into the
centered periodic box together with integer image vectors that reconstruct the
standard coordinates.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET

from multiscale.core.cg_system import CGAngle, CGBond, CGSystem

COULOMB_GALAMOST_CHARGE_FACTOR = math.sqrt(138.935456)
_XML_DECIMALS = 5
_RECONSTRUCTION_TOLERANCE = 2.0 * 10 ** (-_XML_DECIMALS)


def _read_lines(path: str | Path) -> list[str]:
    return Path(path).read_text(encoding="utf-8", errors="replace").splitlines()


def parse_galamost_xml_counts(path: str | Path) -> dict[str, Any]:
    """Return GALAMOST XML section counts without loading large arrays."""
    import re

    counts = {
        "natoms": None,
        "box": None,
        "positions": 0,
        "images": 0,
        "masses": 0,
        "types": 0,
        "charges": 0,
        "bonds": 0,
        "angles": 0,
    }
    section: str | None = None
    tag_to_count = {
        "position": "positions",
        "image": "images",
        "mass": "masses",
        "type": "types",
        "charge": "charges",
        "bond": "bonds",
        "angle": "angles",
    }
    for raw in _read_lines(path):
        line = raw.strip()
        if line.startswith("<configuration"):
            match = re.search(r'natoms="(\d+)"', line)
            if match:
                counts["natoms"] = int(match.group(1))
        elif line.startswith("<box"):
            vals = []
            for key in ("lx", "ly", "lz"):
                match = re.search(rf'{key}="([-+0-9.eE]+)"', line)
                vals.append(float(match.group(1)) if match else None)
            counts["box"] = vals
        else:
            opened = next((tag for tag in tag_to_count if line.startswith(f"<{tag}") and not line.startswith(f"</{tag}")), None)
            closed = next((tag for tag in tag_to_count if line.startswith(f"</{tag}")), None)
            if opened is not None:
                section = tag_to_count[opened]
            elif closed is not None:
                section = None
            elif section in tag_to_count.values() and line and not line.startswith("<"):
                counts[section] += 1
    return counts


@dataclass(frozen=True)
class GalamostExportReport:
    output: str
    counts: dict[str, Any]
    reference_counts: dict[str, Any] | None
    checks: dict[str, bool | None]
    topology_inbox: dict[str, Any]

    @property
    def topology_ok(self) -> bool:
        return all(value is not False for value in self.checks.values())

    @property
    def ok(self) -> bool:
        return bool(self.topology_ok and self.topology_inbox.get("ok"))

    def to_dict(self) -> dict[str, Any]:
        return {
            "output": self.output,
            "counts": self.counts,
            "reference_counts": self.reference_counts,
            "checks": self.checks,
            "topology_inbox": self.topology_inbox,
            "ok": self.ok,
        }


def load_cgsystem_json(path: str | Path) -> CGSystem:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return CGSystem.from_dict(data)


def _length_to_nm(value: float, unit: str | None) -> float:
    normalized = (unit or "angstrom").lower()
    if normalized in {"angstrom", "angstroms", "a", "å"}:
        return value * 0.1
    if normalized in {"nm", "nanometer", "nanometers"}:
        return value
    raise ValueError(f"Unsupported length unit for GALAMOST export: {unit!r}")


def _position_to_galamost_nm(
    system: CGSystem,
    position: tuple[float, float, float] | None,
) -> tuple[float, float, float]:
    if position is None:
        return (0.0, 0.0, 0.0)
    if system.box:
        centered = tuple(float(p) - float(b) / 2.0 for p, b in zip(position, system.box))
    else:
        centered = tuple(float(p) for p in position)
    return tuple(_length_to_nm(value, system.coordinate_unit) for value in centered)  # type: ignore[return-value]


def _box_to_galamost_nm(system: CGSystem) -> tuple[float, float, float]:
    if not system.box:
        return (0.0, 0.0, 0.0)
    return tuple(_length_to_nm(float(value), system.box_unit) for value in system.box)  # type: ignore[return-value]


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


def _bond_sort_key(bond: CGBond) -> tuple[tuple[str, int], int, int]:
    return (_type_sort_key(bond.bond_type), bond.i, bond.j)


def _angle_sort_key(angle: CGAngle) -> tuple[tuple[str, int], int, int, int]:
    return (_type_sort_key(angle.angle_type), angle.i, angle.j, angle.k)


def _xml_float(value: float) -> float:
    """Return the numeric value represented by the XML five-decimal format."""
    return float(f"{value:.{_XML_DECIMALS}f}")


def _wrap_centered_coordinate(value: float, box_length: float) -> tuple[float, int]:
    """Wrap one centered coordinate into ``[-L/2, L/2)``.

    The calculation uses the same five-decimal values that are written to XML,
    so ``wrapped + image * box`` reconstructs the coordinate in
    ``topology.xml`` within the XML precision.
    """
    value_xml = _xml_float(value)
    box_xml = _xml_float(box_length)
    if box_xml <= 0.0:
        raise ValueError("topology_inbox.xml requires positive box lengths")

    image = math.floor(value_xml / box_xml + 0.5)
    wrapped = value_xml - image * box_xml
    half = box_xml / 2.0

    if wrapped >= half:
        wrapped -= box_xml
        image += 1
    elif wrapped < -half:
        wrapped += box_xml
        image -= 1

    wrapped = _xml_float(wrapped)
    if wrapped >= half:
        wrapped = _xml_float(wrapped - box_xml)
        image += 1
    elif wrapped < -half:
        wrapped = _xml_float(wrapped + box_xml)
        image -= 1
    if abs(wrapped) < 0.5 * 10 ** (-_XML_DECIMALS):
        wrapped = 0.0
    return wrapped, int(image)


def _append_topology_sections(
    lines: list[str],
    system: CGSystem,
    *,
    charge_factor: float,
) -> None:
    mass_values = [f'{float(bead.mass or 0.0):11.5f}' for bead in system.beads]
    type_values = [str(bead.bead_type or bead.name) for bead in system.beads]
    charge_values = [f'{float(bead.charge or 0.0) * charge_factor:11.5f}' for bead in system.beads]

    _section(lines, "mass", mass_values)
    _section(lines, "type", type_values)
    _section(lines, "charge", charge_values)

    lines.append("<bond>")
    for bond in sorted(system.bonds, key=_bond_sort_key):
        lines.append(f'{bond.bond_type or "bond"} {bond.i - system.index_base} {bond.j - system.index_base}')
    lines.append("</bond>")

    lines.append("<angle>")
    for angle in sorted(system.angles, key=_angle_sort_key):
        lines.append(
            f'{angle.angle_type or "angle"} '
            f'{angle.i - system.index_base} {angle.j - system.index_base} {angle.k - system.index_base}'
        )
    lines.append("</angle>")
    lines.append("</configuration>")
    lines.append("</galamost_xml>")


def render_galamost_xml(system: CGSystem, *, eps: float = 1.0) -> str:
    """Render the unchanged mapped-coordinate ``topology.xml`` payload."""
    if eps <= 0:
        raise ValueError("eps must be positive")
    lx, ly, lz = _box_to_galamost_nm(system)
    charge_factor = COULOMB_GALAMOST_CHARGE_FACTOR / math.sqrt(eps)

    lines: list[str] = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<galamost_xml version="1.3">',
        f'<configuration time_step="0" dimensions="3" natoms="{len(system.beads)}" >',
        f'<box lx="{lx:.5f}" ly="{ly:.5f}" lz="{lz:.5f}"/>',
    ]

    position_values = []
    for bead in system.beads:
        x, y, z = _position_to_galamost_nm(system, bead.position)
        position_values.append(f'{x:11.5f} {y:11.5f} {z:11.5f}')
    _section(lines, "position", position_values)
    _append_topology_sections(lines, system, charge_factor=charge_factor)
    return "\n".join(lines) + "\n"


def render_galamost_inbox_xml(system: CGSystem, *, eps: float = 1.0) -> str:
    """Render centered box coordinates and integer image vectors."""
    if eps <= 0:
        raise ValueError("eps must be positive")
    lx, ly, lz = _box_to_galamost_nm(system)
    box = (lx, ly, lz)
    if any(length <= 0.0 for length in box):
        raise ValueError("topology_inbox.xml requires a positive three-dimensional box")
    charge_factor = COULOMB_GALAMOST_CHARGE_FACTOR / math.sqrt(eps)

    lines: list[str] = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<galamost_xml version="1.3">',
        f'<configuration time_step="0" dimensions="3" natoms="{len(system.beads)}" >',
        f'<box lx="{lx:.5f}" ly="{ly:.5f}" lz="{lz:.5f}"/>',
    ]

    position_values: list[str] = []
    image_values: list[str] = []
    for bead in system.beads:
        original = _position_to_galamost_nm(system, bead.position)
        wrapped_and_images = [
            _wrap_centered_coordinate(value, length)
            for value, length in zip(original, box)
        ]
        wrapped = [item[0] for item in wrapped_and_images]
        images = [item[1] for item in wrapped_and_images]
        position_values.append(f'{wrapped[0]:11.5f} {wrapped[1]:11.5f} {wrapped[2]:11.5f}')
        image_values.append(f'{images[0]} {images[1]} {images[2]}')

    _section(lines, "position", position_values)
    _section(lines, "image", image_values)
    _append_topology_sections(lines, system, charge_factor=charge_factor)
    return "\n".join(lines) + "\n"


def write_galamost_topology_xml(
    system: CGSystem,
    outpath: str | Path,
    *,
    eps: float = 1.0,
) -> Path:
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    outpath.write_text(render_galamost_xml(system, eps=eps), encoding="utf-8")
    return outpath


def write_galamost_inbox_xml(
    system: CGSystem,
    outpath: str | Path,
    *,
    eps: float = 1.0,
) -> Path:
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    outpath.write_text(render_galamost_inbox_xml(system, eps=eps), encoding="utf-8")
    return outpath


def _xml_rows(path: Path, section: str) -> list[list[str]]:
    root = ET.parse(path).getroot()
    configuration = root.find("configuration")
    if configuration is None:
        raise ValueError(f"No <configuration> found in {path}")
    node = configuration.find(section)
    if node is None or node.text is None:
        return []
    return [line.split() for line in node.text.splitlines() if line.split()]


def _xml_box(path: Path) -> tuple[float, float, float]:
    root = ET.parse(path).getroot()
    configuration = root.find("configuration")
    if configuration is None:
        raise ValueError(f"No <configuration> found in {path}")
    box = configuration.find("box")
    if box is None:
        raise ValueError(f"No <box> found in {path}")
    return tuple(float(box.attrib[key]) for key in ("lx", "ly", "lz"))  # type: ignore[return-value]


def _validate_inbox_xml(original_path: Path, inbox_path: Path) -> dict[str, Any]:
    original_counts = parse_galamost_xml_counts(original_path)
    inbox_counts = parse_galamost_xml_counts(inbox_path)
    box = _xml_box(inbox_path)
    original_positions = [[float(value) for value in row[:3]] for row in _xml_rows(original_path, "position")]
    inbox_positions = [[float(value) for value in row[:3]] for row in _xml_rows(inbox_path, "position")]
    images = [[int(value) for value in row[:3]] for row in _xml_rows(inbox_path, "image")]

    all_positions_inside = len(inbox_positions) == int(inbox_counts.get("natoms") or 0)
    if all_positions_inside:
        for position in inbox_positions:
            for value, length in zip(position, box):
                if not (-length / 2.0 <= value < length / 2.0):
                    all_positions_inside = False
                    break
            if not all_positions_inside:
                break

    reconstructs = (
        len(original_positions) == len(inbox_positions) == len(images)
    )
    if reconstructs:
        for original, wrapped, image in zip(original_positions, inbox_positions, images):
            for original_value, wrapped_value, image_value, length in zip(original, wrapped, image, box):
                reconstructed = wrapped_value + image_value * length
                if not math.isclose(
                    reconstructed,
                    original_value,
                    rel_tol=0.0,
                    abs_tol=_RECONSTRUCTION_TOLERANCE,
                ):
                    reconstructs = False
                    break
            if not reconstructs:
                break

    topology_sections = ("mass", "type", "charge", "bond", "angle")
    topology_matches = all(
        _xml_rows(original_path, section) == _xml_rows(inbox_path, section)
        for section in topology_sections
    )
    checks = {
        "has_positions": inbox_counts.get("positions") == inbox_counts.get("natoms"),
        "has_images": inbox_counts.get("images") == inbox_counts.get("natoms"),
        "image_count_matches_natoms": inbox_counts.get("images") == inbox_counts.get("natoms"),
        "all_positions_inside_box": all_positions_inside,
        "unwrapped_positions_reconstruct": reconstructs,
        "topology_matches_original": topology_matches,
        "box_matches_original": inbox_counts.get("box") == original_counts.get("box"),
        "natoms_matches_original": inbox_counts.get("natoms") == original_counts.get("natoms"),
    }
    return {
        "output": str(inbox_path),
        "counts": inbox_counts,
        "checks": checks,
        "coordinate_interval": "[-L/2, L/2)",
        "reconstruction": "topology = topology_inbox + image * box",
        "ok": all(checks.values()),
    }


def compare_galamost_xml_counts(
    reference_xml: str | Path,
    generated_xml: str | Path,
) -> dict[str, Any]:
    reference = parse_galamost_xml_counts(reference_xml)
    generated = parse_galamost_xml_counts(generated_xml)
    checks = {
        "natoms_match": generated.get("natoms") == reference.get("natoms"),
        "positions_match": generated.get("positions") == reference.get("positions"),
        "masses_match": generated.get("masses") == reference.get("masses"),
        "types_match": generated.get("types") == reference.get("types"),
        "charges_match": generated.get("charges") == reference.get("charges"),
        "bonds_match": generated.get("bonds") == reference.get("bonds"),
        "angles_match": generated.get("angles") == reference.get("angles"),
        "box_match": generated.get("box") == reference.get("box"),
    }
    return {
        "reference_counts": reference,
        "generated_counts": generated,
        "checks": checks,
        "ok": all(checks.values()),
    }


def export_galamost_topology(
    system: CGSystem,
    outdir: str | Path,
    *,
    eps: float = 1.0,
    reference_xml: str | Path | None = None,
) -> GalamostExportReport:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    xml_path = write_galamost_topology_xml(system, outdir / "topology.xml", eps=eps)
    inbox_path = write_galamost_inbox_xml(system, outdir / "topology_inbox.xml", eps=eps)

    counts = parse_galamost_xml_counts(xml_path)
    reference_counts = parse_galamost_xml_counts(reference_xml) if reference_xml else None
    checks: dict[str, bool | None] = {
        "has_positions": counts.get("positions") == counts.get("natoms"),
        "has_masses": counts.get("masses") == counts.get("natoms"),
        "has_types": counts.get("types") == counts.get("natoms"),
        "has_charges": counts.get("charges") == counts.get("natoms"),
        "has_bonds": int(counts.get("bonds") or 0) > 0,
        "has_angles": int(counts.get("angles") or 0) > 0,
    }
    if reference_counts:
        checks.update(
            {
                "natoms_match_reference": counts.get("natoms") == reference_counts.get("natoms"),
                "bonds_match_reference": counts.get("bonds") == reference_counts.get("bonds"),
                "angles_match_reference": counts.get("angles") == reference_counts.get("angles"),
                "box_match_reference": counts.get("box") == reference_counts.get("box"),
            }
        )
    inbox_report = _validate_inbox_xml(xml_path, inbox_path)
    return GalamostExportReport(
        str(xml_path),
        counts,
        reference_counts,
        checks,
        inbox_report,
    )
