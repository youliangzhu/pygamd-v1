"""In-memory AA-to-CG mapping and topology models."""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class BeadMapping:
    """One coarse-grained bead and its 1-based atom membership."""

    name: str
    atom_indices_1based: list[int]

    @property
    def atom_count(self) -> int:
        return len(self.atom_indices_1based)

    @property
    def is_count_consistent(self) -> bool:
        return self.atom_count > 0


@dataclass(frozen=True)
class MappingSpec:
    """Mapping definition for one repeated molecular type."""

    name: str
    parent_type: str
    beads: list[BeadMapping]
    source: str | None = None

    def to_dict(self) -> dict:
        return {
            "source": self.source,
            "name": self.name,
            "parent_type": self.parent_type,
            "beads": [
                asdict(bead) | {"atom_count": bead.atom_count}
                for bead in self.beads
            ],
        }


@dataclass(frozen=True)
class BondGroup:
    molecular_type: str
    type_id: str
    pairs: list[tuple[str, str]]

    def to_dict(self) -> dict:
        return {
            "molecular_type": self.molecular_type,
            "type_id": self.type_id,
            "pairs": [list(pair) for pair in self.pairs],
        }


@dataclass(frozen=True)
class AngleGroup:
    molecular_type: str
    type_id: str
    triples: list[tuple[str, str, str]]

    def to_dict(self) -> dict:
        return {
            "molecular_type": self.molecular_type,
            "type_id": self.type_id,
            "triples": [list(triple) for triple in self.triples],
        }


@dataclass(frozen=True)
class TopologySpec:
    """Bead-type groups plus per-molecule bond and angle templates."""

    cgtypes: dict[str, list[str]]
    bond_groups: list[BondGroup]
    angle_groups: list[AngleGroup]
    source: str | None = None

    @property
    def has_bond_angle_topology(self) -> bool:
        return bool(self.bond_groups or self.angle_groups)

    def to_dict(self) -> dict:
        return {
            "source": self.source,
            "cgtypes": self.cgtypes,
            "bond_groups": [group.to_dict() for group in self.bond_groups],
            "angle_groups": [group.to_dict() for group in self.angle_groups],
        }
