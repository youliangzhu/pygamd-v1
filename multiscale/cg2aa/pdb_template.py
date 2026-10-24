from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class AtomRecord:
    serial: int
    name: str
    residue: str
    residue_id: int
    coord: np.ndarray
    element: str


@dataclass
class PdbTemplate:
    atoms: list[AtomRecord]
    bonds: list[tuple[int, int]]

    @property
    def coords(self) -> np.ndarray:
        return np.asarray([atom.coord for atom in self.atoms], dtype=float)


def _guess_element(atom_name: str, raw_element: str) -> str:
    element = raw_element.strip()
    if element:
        return element
    stripped = "".join(ch for ch in atom_name if ch.isalpha())
    if not stripped:
        return "X"
    if len(stripped) >= 2 and stripped[:2].title() in {"Cl", "Br", "Si"}:
        return stripped[:2].title()
    return stripped[0].upper()


def read_pdb(path: str | Path) -> PdbTemplate:
    atoms: list[AtomRecord] = []
    serial_to_index: dict[int, int] = {}
    bonds: set[tuple[int, int]] = set()
    path = Path(path)

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            if record in {"ATOM", "HETATM"}:
                serial = int(line[6:11])
                name = line[12:16].strip() or f"A{serial}"
                residue = line[17:20].strip() or "MOL"
                residue_raw = line[22:26].strip()
                residue_id = int(residue_raw) if residue_raw else 1
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = _guess_element(name, line[76:78] if len(line) >= 78 else "")
                serial_to_index[serial] = len(atoms)
                atoms.append(
                    AtomRecord(
                        serial=serial,
                        name=name,
                        residue=residue,
                        residue_id=residue_id,
                        coord=np.asarray([x, y, z], dtype=float),
                        element=element,
                    )
                )
            elif record == "CONECT":
                fields = line.split()
                if len(fields) < 3:
                    continue
                src_serial = int(fields[1])
                for value in fields[2:]:
                    dst_serial = int(value)
                    if src_serial != dst_serial:
                        bonds.add(tuple(sorted((src_serial, dst_serial))))

    if not atoms:
        raise ValueError(f"No ATOM/HETATM records found in {path}")

    indexed_bonds = []
    for a_serial, b_serial in sorted(bonds):
        if a_serial in serial_to_index and b_serial in serial_to_index:
            indexed_bonds.append((serial_to_index[a_serial], serial_to_index[b_serial]))
    return PdbTemplate(atoms=atoms, bonds=indexed_bonds)
