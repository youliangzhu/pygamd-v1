"""Utilities for reading GROMACS .gro/.top/.itp files.

This module keeps the light interface used by the newer ``gro2xml.py`` while
covering the sections that the old ``gro_to_xml.py`` handled:

- [atomtypes], [bondtypes], [constrainttypes], [angletypes], [dihedraltypes]
- [nonbond_params] and [pairtypes]
- per-molecule [atoms], [bonds], [pairs], [constraints], [virtual_sites3],
  [angles], [dihedrals], and [exclusions]
- simple ``#define`` aliases used by bonded/non-bonded parameters
"""

from __future__ import annotations

import os
import re
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# Basic text helpers
# ---------------------------------------------------------------------------

def _strip_comment(line: str) -> str:
    """Remove GROMACS ';' comments and strip whitespace."""
    idx = line.find(";")
    if idx != -1:
        line = line[:idx]
    return line.strip()


def _parse_section_header(line: str) -> Optional[str]:
    """Return a normalized section name for a GROMACS section header."""
    match = re.match(r"^\[\s*([^\]]+?)\s*\]", line)
    if match:
        return match.group(1).strip().lower()
    return None


def _parse_include(line: str) -> Optional[str]:
    """Return the path inside a #include directive, if present."""
    match = re.match(r'^\s*#include\s+["<]([^">]+)[">]', line)
    if match:
        return match.group(1).strip()
    return None


def _parse_define(line: str) -> Optional[Tuple[str, List[str]]]:
    """Return (name, values) for a simple '#define name value...' line."""
    parts = _strip_comment(line).split()
    if len(parts) >= 3 and parts[0] == "#define":
        return parts[1], parts[2:]
    return None


def read1(line: str, list1: List[List[str]]) -> None:
    """Compatibility helper: append a split non-empty data line."""
    line = _strip_comment(line)
    parts = line.split()
    if parts and parts[0] != ";":
        list1.append(parts)


# ---------------------------------------------------------------------------
# .gro reader
# ---------------------------------------------------------------------------

def read_gro_simple(filename: str):
    """Read a GROMACS .gro file.

    Returns
    -------
    tuple
        ``(num, box, pos, vel, gro_types)`` where ``num`` is an int, ``box`` is
        ``[lx, ly, lz]``, ``pos`` and ``vel`` are coordinate lists, and
        ``gro_types`` uses the atom name field from the .gro file.
    """
    pos: List[List[float]] = []
    vel: List[List[float]] = []
    gro_types: List[str] = []

    with open(filename, encoding="utf-8", errors="replace") as rf:
        _title = rf.readline()
        count_line = rf.readline()
        if not count_line:
            raise RuntimeError(f"Invalid .gro file: missing atom count: {filename}")
        num = int(count_line.strip())

        for _ in range(num):
            line = rf.readline().rstrip("\n")
            if not line:
                raise RuntimeError(f"Invalid .gro file: missing atom line in {filename}")

            # Prefer fixed-width .gro fields; fall back to split for loose files.
            atom_name = line[10:15].strip() if len(line) >= 15 else ""
            parts = line.split()
            if len(parts) >= 6:
                try:
                    x, y, z = float(parts[-6]), float(parts[-5]), float(parts[-4])
                    vx, vy, vz = float(parts[-3]), float(parts[-2]), float(parts[-1])
                    pos.append([x, y, z])
                    vel.append([vx, vy, vz])
                except ValueError:
                    x, y, z = float(parts[-3]), float(parts[-2]), float(parts[-1])
                    pos.append([x, y, z])
            else:
                x = float(line[20:28])
                y = float(line[28:36])
                z = float(line[36:44])
                pos.append([x, y, z])
                if len(line) >= 68:
                    vel.append([float(line[44:52]), float(line[52:60]), float(line[60:68])])

            if not atom_name:
                atom_name = parts[1] if len(parts) > 1 else parts[0]
            gro_types.append(atom_name)

        box_line = rf.readline()
        if not box_line:
            raise RuntimeError(f"Invalid .gro file: missing box line: {filename}")
        box_parts = box_line.split()
        if len(box_parts) < 3:
            raise RuntimeError(f"Only orthorhombic .gro boxes are supported: {box_line!r}")
        box = [float(box_parts[0]), float(box_parts[1]), float(box_parts[2])]

    if len(pos) != num:
        raise RuntimeError(f".gro atom count mismatch: header={num}, read={len(pos)}")
    if len(vel) != num:
        vel = []

    return num, box, pos, vel, gro_types


# ---------------------------------------------------------------------------
# .top reader
# ---------------------------------------------------------------------------

def _resolve_include_path(include_name: str, base_dir: str) -> str:
    """Resolve a local include relative to base_dir."""
    if os.path.isabs(include_name):
        return os.path.normpath(include_name)
    return os.path.normpath(os.path.join(base_dir, include_name))


def _collect_includes(filename: str, visited: Optional[set] = None) -> List[str]:
    """Recursively collect local #include files, preserving order."""
    if visited is None:
        visited = set()

    abs_path = os.path.abspath(filename)
    if abs_path in visited:
        return []
    visited.add(abs_path)

    base_dir = os.path.dirname(abs_path)
    result: List[str] = []
    try:
        with open(abs_path, encoding="utf-8", errors="replace") as rf:
            for raw_line in rf:
                line = _strip_comment(raw_line)
                inc = _parse_include(line)
                if not inc:
                    continue
                inc_path = _resolve_include_path(inc, base_dir)
                result.append(inc_path)
                if os.path.isfile(inc_path):
                    result.extend(_collect_includes(inc_path, visited))
    except FileNotFoundError:
        print(f"[warning] include/top file not found while collecting includes: {abs_path}")

    # Deduplicate while preserving order.
    out: List[str] = []
    seen = set()
    for path in result:
        if path not in seen:
            seen.add(path)
            out.append(path)
    return out


def read_top(filename: str):
    """Read a .top file.

    Returns
    -------
    tuple
        ``(files_l, defaults_l, molecules_l)``.
    """
    files_l = _collect_includes(filename)
    defaults_l: List[List[str]] = []
    molecules_l: List[List[str]] = []

    cur_section: Optional[str] = None
    try:
        with open(filename, encoding="utf-8", errors="replace") as rf:
            for raw_line in rf:
                line = _strip_comment(raw_line)
                if not line:
                    continue
                if line.startswith("#"):
                    continue

                section = _parse_section_header(line)
                if section is not None:
                    cur_section = section
                    continue

                parts = line.split()
                if not parts:
                    continue
                if cur_section == "defaults":
                    defaults_l.append(parts)
                elif cur_section == "molecules":
                    if len(parts) >= 2:
                        molecules_l.append([parts[0], parts[1]])
    except FileNotFoundError:
        raise RuntimeError(f"top file not found: {filename}")

    return files_l, defaults_l, molecules_l


# ---------------------------------------------------------------------------
# .itp reader
# ---------------------------------------------------------------------------

def _empty_forcefield_dict() -> Dict[str, object]:
    return {
        "defines": {},
        "defaults": [],
        "atomtypes": [],
        "bondtypes": [],
        "constrainttypes": [],
        "angletypes": [],
        "dihedraltypes": [],
        "nonbond_params": [],
        "pairtypes": [],
    }


def read_itp_full(filename: str):
    """Read molecule and force-field sections from an .itp/.top-like file.

    Returns
    -------
    tuple
        ``(mol_d, forcefield_d)``.
    """
    mol_d: Dict[str, Dict[str, List[List[str]]]] = {}
    ff_d = _empty_forcefield_dict()
    cur_section: Optional[str] = None
    cur_mol: Optional[str] = None

    molecule_sections = {
        "atoms": "atoms",
        "bonds": "bonds",
        "pairs": "pairs",
        "pairs_nb": "pairs",
        "constraints": "constraints",
        "virtual_sites3": "virtual_sites3",
        "angles": "angles",
        "dihedrals": "dihedrals",
        "exclusions": "exclusions",
    }
    global_sections = {
        "defaults": "defaults",
        "atomtypes": "atomtypes",
        "bondtypes": "bondtypes",
        "constrainttypes": "constrainttypes",
        "angletypes": "angletypes",
        "dihedraltypes": "dihedraltypes",
        "nonbond_params": "nonbond_params",
        "pairtypes": "pairtypes",
    }

    try:
        with open(filename, encoding="utf-8", errors="replace") as rf:
            for raw_line in rf:
                define = _parse_define(raw_line)
                if define is not None:
                    name, values = define
                    ff_d["defines"][name] = values
                    continue

                line = _strip_comment(raw_line)
                if not line:
                    continue
                if line.startswith("#"):
                    # Skip #ifdef/#ifndef/#else/#endif and non-local includes here.
                    continue

                section = _parse_section_header(line)
                if section is not None:
                    cur_section = section
                    continue

                parts = line.split()
                if not parts:
                    continue

                if cur_section == "moleculetype":
                    cur_mol = parts[0]
                    mol_d.setdefault(cur_mol, {})
                    continue

                if cur_section in global_sections:
                    ff_d[global_sections[cur_section]].append(parts)
                    continue

                if cur_section in molecule_sections and cur_mol:
                    key = molecule_sections[cur_section]
                    mol_d[cur_mol].setdefault(key, []).append(parts)
                    continue

                # Unknown or currently unsupported sections are intentionally ignored.
    except FileNotFoundError:
        raise RuntimeError(f"itp file not found: {filename}")

    return mol_d, ff_d


def merge_forcefields(forcefields: Iterable[Dict[str, object]]) -> Dict[str, object]:
    """Merge force-field dictionaries from multiple topology files."""
    merged = _empty_forcefield_dict()
    for ff in forcefields:
        merged["defines"].update(ff.get("defines", {}))
        for key in (
            "defaults", "atomtypes", "bondtypes", "constrainttypes",
            "angletypes", "dihedraltypes", "nonbond_params", "pairtypes",
        ):
            merged[key].extend(ff.get(key, []))
    return merged


def read_atomtypes(filename: str):
    """Compatibility helper returning only [atomtypes] rows."""
    _mol_d, ff_d = read_itp_full(filename)
    return ff_d["atomtypes"]


def read_itp_all(filename: str):
    """Compatibility interface used by older gro2xml code.

    Returns ``(mol_d, atomtypes_l)``.
    """
    mol_d, ff_d = read_itp_full(filename)
    return mol_d, ff_d["atomtypes"]


def read_line(filename: str):
    """Compatibility helper for simple single-file molecule sections."""
    mol_d, _ff_d = read_itp_full(filename)
    if not mol_d:
        return [], [], [], [], []
    first = next(iter(mol_d.values()))
    return (
        first.get("atoms", []),
        first.get("bonds", []),
        first.get("angles", []),
        first.get("dihedrals", []),
        first.get("pairs", []),
    )
