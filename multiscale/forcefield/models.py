"""Data models for the Multiscale forcefield namespace."""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

@dataclass(frozen=True)
class ForceFieldAssets:
    name: str
    source: str
    asset_dir: Path | None
    template_path: Path | None
    bonded_path: Path | None
    missing: tuple[str, ...] = ()
    @property
    def ok(self) -> bool: return not self.missing
    def to_dict(self) -> dict[str, Any]:
        return {"name": self.name, "source": self.source,
                "asset_dir": str(self.asset_dir) if self.asset_dir else None,
                "template_path": str(self.template_path) if self.template_path else None,
                "bonded_path": str(self.bonded_path) if self.bonded_path else None,
                "missing": list(self.missing), "ok": self.ok}

@dataclass(frozen=True)
class OPLSAtomType:
    element: str
    atom_name: str
    gromacs_type: str
    smarts: str
    moltemplate_type: str
    charge: float
    description: str
    def to_dict(self) -> dict[str, Any]:
        return {"element": self.element, "atom_name": self.atom_name,
                "gromacs_type": self.gromacs_type, "smarts": self.smarts,
                "moltemplate_type": self.moltemplate_type, "charge": self.charge,
                "description": self.description}

@dataclass(frozen=True)
class OPLSSummary:
    assets: ForceFieldAssets
    atom_type_count: int
    unique_moltemplate_types: int
    elements: tuple[str, ...]
    bonded_section_counts: dict[str, int] = field(default_factory=dict)
    required_defaults_present: dict[str, bool] = field(default_factory=dict)
    @property
    def ok(self) -> bool:
        return self.assets.ok and self.atom_type_count > 0 and all(self.required_defaults_present.values())
    def to_dict(self) -> dict[str, Any]:
        return {"ok": self.ok, "assets": self.assets.to_dict(),
                "atom_type_count": self.atom_type_count,
                "unique_moltemplate_types": self.unique_moltemplate_types,
                "elements": list(self.elements),
                "bonded_section_counts": dict(self.bonded_section_counts),
                "required_defaults_present": dict(self.required_defaults_present)}

@dataclass
class ForceFieldAssignment:
    molecule: Any
    bonds: list[str]
    angles: list[str]
    dihedrals: list[str]
    atom_count: int
    atom_type_counts: dict[str, int]
    opls_type_counts: dict[str, int]
    total_charge: float
    warnings: list[str] = field(default_factory=list)
    @property
    def ok(self) -> bool: return self.atom_count >= 0
    def to_dict(self) -> dict[str, Any]:
        return {"ok": self.ok, "atom_count": self.atom_count,
                "bond_count": len(self.bonds), "angle_count": len(self.angles),
                "dihedral_count": len(self.dihedrals),
                "atom_type_counts": dict(self.atom_type_counts),
                "opls_type_counts": dict(self.opls_type_counts),
                "total_charge": self.total_charge, "warnings": list(self.warnings)}
