"""Neutral coarse-grained system data model for Multiscale.

The model is shared by ``aa2cg``/``cg2aa``/``forcefield`` modules. Bead, bond,
and angle indices are 1-based unless an IO backend converts explicitly.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Any


@dataclass
class CGBead:
    index: int
    name: str
    bead_type: str | None = None
    molecule_index: int | None = None
    molecule_bead_index: int | None = None
    mass: float | None = None
    charge: float | None = None
    position: tuple[float, float, float] | None = None

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "CGBead":
        position = data.get("position")
        if position is not None:
            position = tuple(float(x) for x in position)
        return cls(
            index=int(data["index"]),
            name=str(data.get("name") or data.get("bead_type") or data["index"]),
            bead_type=data.get("bead_type"),
            molecule_index=data.get("molecule_index"),
            molecule_bead_index=data.get("molecule_bead_index"),
            mass=data.get("mass"),
            charge=data.get("charge"),
            position=position,  # type: ignore[arg-type]
        )


@dataclass
class CGBond:
    i: int
    j: int
    bond_type: str | None = None
    molecule_index: int | None = None

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "CGBond":
        return cls(
            i=int(data["i"]),
            j=int(data["j"]),
            bond_type=data.get("bond_type"),
            molecule_index=data.get("molecule_index"),
        )


@dataclass
class CGAngle:
    i: int
    j: int
    k: int
    angle_type: str | None = None
    molecule_index: int | None = None

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "CGAngle":
        return cls(
            i=int(data["i"]),
            j=int(data["j"]),
            k=int(data["k"]),
            angle_type=data.get("angle_type"),
            molecule_index=data.get("molecule_index"),
        )


@dataclass
class CGSystem:
    name: str
    box: tuple[float, float, float] | None = None
    box_unit: str | None = None
    coordinate_unit: str | None = None
    index_base: int = 1
    nmol: int | None = None
    beads_per_molecule: int | None = None
    beads: list[CGBead] = field(default_factory=list)
    bonds: list[CGBond] = field(default_factory=list)
    angles: list[CGAngle] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def summary(self) -> dict[str, Any]:
        bead_types = sorted({b.bead_type for b in self.beads if b.bead_type})
        bond_types = sorted({b.bond_type for b in self.bonds if b.bond_type})
        angle_types = sorted({a.angle_type for a in self.angles if a.angle_type})
        return {
            "name": self.name,
            "nmol": self.nmol,
            "beads_per_molecule": self.beads_per_molecule,
            "natoms": len(self.beads),
            "nbonds": len(self.bonds),
            "nangles": len(self.angles),
            "bead_types": bead_types,
            "bond_types": bond_types,
            "angle_types": angle_types,
            "box": list(self.box) if self.box else None,
            "box_unit": self.box_unit,
            "coordinate_unit": self.coordinate_unit,
            "index_base": self.index_base,
        }

    def to_dict(self, *, include_records: bool = True) -> dict[str, Any]:
        data: dict[str, Any] = {
            "schema": "multiscale.cgsystem.v0",
            "summary": self.summary(),
            "metadata": self.metadata,
        }
        if include_records:
            data["beads"] = [asdict(b) for b in self.beads]
            data["bonds"] = [asdict(b) for b in self.bonds]
            data["angles"] = [asdict(a) for a in self.angles]
        return data

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "CGSystem":
        summary = data.get("summary") or {}
        box = summary.get("box")
        if box is not None:
            box = tuple(float(x) for x in box)
        return cls(
            name=str(summary.get("name") or "cg_system"),
            box=box,  # type: ignore[arg-type]
            box_unit=summary.get("box_unit"),
            coordinate_unit=summary.get("coordinate_unit"),
            index_base=int(summary.get("index_base") or 1),
            nmol=summary.get("nmol"),
            beads_per_molecule=summary.get("beads_per_molecule"),
            beads=[CGBead.from_dict(item) for item in data.get("beads", [])],
            bonds=[CGBond.from_dict(item) for item in data.get("bonds", [])],
            angles=[CGAngle.from_dict(item) for item in data.get("angles", [])],
            metadata=data.get("metadata") or {},
        )


def load_cgsystem_json(path: str | Path) -> CGSystem:
    """Load an Multiscale CGSystem JSON file."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return CGSystem.from_dict(data)


def write_cgsystem_json(system: CGSystem, path: str | Path, *, include_records: bool = True) -> Path:
    """Write a CGSystem JSON file."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(system.to_dict(include_records=include_records), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return out


def normalized_cgsystem_dict(data: dict[str, Any]) -> dict[str, Any]:
    """Normalize non-semantic fields for regression comparisons.

    The aa2cg mapper records a temporary GROMACS dump path under
    ``metadata.source_files.topology_dump``.  That path is expected to change
    between runs and should not affect semantic comparisons.
    """
    clone = json.loads(json.dumps(data))
    try:
        clone["metadata"]["source_files"].pop("topology_dump", None)
    except Exception:
        pass
    return clone
