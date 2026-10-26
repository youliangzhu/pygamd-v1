"""Small GROMACS text-file readers used by native multiscale workflows.

The scope is intentionally limited to text artifacts that are useful during
Single-frame AA structure mapping from GROMACS ``.gro`` coordinates.  Binary ``.tpr``
and ``.xtc`` files remain outside this module.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class GroAtom:
    index: int
    residue_number: int
    residue_name: str
    atom_name: str
    atom_number: int
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class GroFrame:
    title: str
    natoms: int
    atoms: list[GroAtom]
    box: tuple[float, float, float] | None
    box_raw: tuple[float, ...] | None = None


def _parse_gro_atom_line(line: str, index: int) -> GroAtom:
    """Parse one GRO atom line.

    GROMACS ``.gro`` atom records are fixed-width, but permissive text inputs may
    occasionally be whitespace-separated.  The fixed-width parser is tried first
    and a conservative split fallback is used only if that fails.
    """
    try:
        return GroAtom(
            index=index,
            residue_number=int(line[0:5]),
            residue_name=line[5:10].strip(),
            atom_name=line[10:15].strip(),
            atom_number=int(line[15:20]),
            x=float(line[20:28]),
            y=float(line[28:36]),
            z=float(line[36:44]),
        )
    except (ValueError, IndexError):
        parts = line.split()
        if len(parts) < 6:
            raise ValueError(f"Invalid GRO atom line {index}: {line!r}") from None
        # Common split form: resnr resname atomname atomnr x y z
        if len(parts) >= 7 and parts[0].lstrip("+-").isdigit():
            return GroAtom(
                index=index,
                residue_number=int(parts[0]),
                residue_name=parts[1],
                atom_name=parts[2],
                atom_number=int(parts[3]),
                x=float(parts[4]),
                y=float(parts[5]),
                z=float(parts[6]),
            )
        # Fallback for fixed-width merged residue field split by whitespace, e.g.
        # Some text inputs omit fixed-width spacing. Keep this parser deliberately narrow.
        res_field = parts[0]
        digits = ""
        letters = ""
        for char in res_field:
            if char.isdigit() and not letters:
                digits += char
            else:
                letters += char
        if not digits or not letters or len(parts) < 6:
            raise ValueError(f"Invalid GRO atom line {index}: {line!r}") from None
        return GroAtom(
            index=index,
            residue_number=int(digits),
            residue_name=letters,
            atom_name=parts[1],
            atom_number=int(parts[2]),
            x=float(parts[3]),
            y=float(parts[4]),
            z=float(parts[5]),
        )


def parse_gro(path: str | Path) -> GroFrame:
    path = Path(path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 3:
        raise ValueError(f"GRO file is too short: {path}")
    title = lines[0].rstrip("\n")
    try:
        natoms = int(lines[1].strip())
    except ValueError as exc:
        raise ValueError(f"Invalid GRO atom count line in {path}: {lines[1]!r}") from exc
    if len(lines) < natoms + 3:
        raise ValueError(f"GRO atom count mismatch: header={natoms}, available={max(0, len(lines) - 3)}")
    atoms = [_parse_gro_atom_line(line, idx) for idx, line in enumerate(lines[2 : 2 + natoms], start=1)]
    box_parts = [float(x) for x in lines[2 + natoms].split()]
    box_raw = tuple(box_parts) if box_parts else None
    box = tuple(box_parts[:3]) if len(box_parts) >= 3 else None
    return GroFrame(title=title, natoms=natoms, atoms=atoms, box=box, box_raw=box_raw)
