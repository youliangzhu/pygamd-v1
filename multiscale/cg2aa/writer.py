from __future__ import annotations

from pathlib import Path

import numpy as np

from .geometry import cell_lengths_and_angles


XML_TEMPLATE = """<?xml version ="1.0" encoding ="UTF-8" ?>
<{program}_xml version="{version}">
<configuration time_step="0" dimensions="3" natoms="{n_atoms:d}" >
<box lx="{lx:.8f}" ly="{ly:.8f}" lz="{lz:.8f}" xy="{xy:.8f}" xz="{xz:.8f}" yz="{yz:.8f}"/>
<position num="{n_atoms:d}">
{positions}</position>
{image_block}<type num="{n_atoms:d}">
{types}</type>
<opls_type num="{n_atoms:d}">
{opls_type}</opls_type>
<monomer_id num="{n_atoms:d}">
{monomer_id}</monomer_id>
<charge num="{n_atoms:d}">
{charge}</charge>
<mass num="{n_atoms:d}">
{mass}</mass>
<bond num="{n_bonds:d}">
{bond}
</bond>
<angle num="{n_angles:d}">
{angle}
</angle>
<dihedral num="{n_dihedrals:d}">
{dihedral}
</dihedral>
</configuration>
</{program}_xml>"""


def _atom_prop(atom, key: str, default: str) -> str:
    try:
        return atom.GetProp(key)
    except KeyError:
        return default


def _atom_int_prop(atom, key: str, default: int) -> int:
    try:
        return atom.GetIntProp(key)
    except KeyError:
        return default


def write_forcefield_xml(
    path: str | Path,
    molecule,
    coords: list[np.ndarray],
    box,
    bonds: list[str],
    angles: list[str],
    dihedrals: list[str],
    position_scale: float = 1.0,
    program: str = "galamost",
    version: str = "1.3",
    images: np.ndarray | list | None = None,
) -> str:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    n_atoms = molecule.GetNumAtoms()
    positions = ""
    types = ""
    opls_type = ""
    monomer_id = ""
    charge = ""
    mass = ""
    for atom in molecule.GetAtoms():
        pos = np.asarray(coords[atom.GetIdx()], dtype=float) * position_scale
        positions += "%.6f %.6f %.6f\n" % (pos[0], pos[1], pos[2])
        types += "%s\n" % _atom_prop(atom, "OplsType", _atom_prop(atom, "AtomType", ""))
        opls_type += "%s\n" % _atom_prop(atom, "ElementType", atom.GetSymbol())
        monomer_id += "%d\n" % _atom_int_prop(atom, "molecule_id", 0)
        charge += "%.6f\n" % float(_atom_prop(atom, "Charge", "0.0"))
        mass += "%.6f\n" % atom.GetMass()

    image_block = ""
    if images is not None:
        image_array = np.asarray(images, dtype=int)
        if image_array.shape != (n_atoms, 3):
            raise ValueError(f"images must have shape ({n_atoms}, 3), got {image_array.shape}")
        image_text = "".join(f"{row[0]} {row[1]} {row[2]}\n" for row in image_array)
        image_block = f'<image num="{n_atoms:d}">\n{image_text}</image>\n'

    values = [float(value) for value in box]
    values.extend([0.0] * (6 - len(values)))
    lx, ly, lz, xy, xz, yz = values[:6]
    text = XML_TEMPLATE.format(
        program=program,
        version=version,
        n_atoms=n_atoms,
        lx=lx,
        ly=ly,
        lz=lz,
        xy=xy,
        xz=xz,
        yz=yz,
        positions=positions,
        image_block=image_block,
        types=types,
        opls_type=opls_type,
        monomer_id=monomer_id,
        charge=charge,
        mass=mass,
        n_bonds=len(bonds),
        bond="\n".join(bonds),
        n_angles=len(angles),
        angle="\n".join(angles),
        n_dihedrals=len(dihedrals),
        dihedral="\n".join(dihedrals),
    )
    path.write_text(text, encoding="utf-8")
    return str(path)


def write_pdb(
    path: str | Path,
    molecule,
    coords: list[np.ndarray],
    *,
    box=None,
    position_scale: float = 1.0,
) -> str:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    adjacency: dict[int, list[int]] = {}
    for bond in molecule.GetBonds():
        i = bond.GetBeginAtomIdx() + 1
        j = bond.GetEndAtomIdx() + 1
        adjacency.setdefault(i, []).append(j)
        adjacency.setdefault(j, []).append(i)

    with path.open("w", encoding="utf-8") as handle:
        handle.write("REMARK Multiscale cg2aa generated PDB\n")
        if box is not None:
            scaled_box = [float(value) * float(position_scale) for value in box]
            a, b, c, alpha, beta, gamma = cell_lengths_and_angles(scaled_box)
            handle.write(
                f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1           1\n"
            )
        for idx, atom in enumerate(molecule.GetAtoms(), start=1):
            pos = np.asarray(coords[idx - 1], dtype=float) * float(position_scale)
            name = _atom_prop(atom, "AtomName", atom.GetSymbol())[:4]
            residue = _atom_prop(atom, "ResidueName", "MOL")[:3]
            molecule_id = _atom_int_prop(atom, "molecule_id", 0) + 1
            handle.write(
                f"ATOM  {idx:5d} {name:>4s} {residue:>3s} A{molecule_id:4d}    "
                f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00  0.00          {atom.GetSymbol()[:2]:>2s}\n"
            )
        for atom_id in sorted(adjacency):
            neighbors = "".join(f"{neighbor:5d}" for neighbor in sorted(set(adjacency[atom_id])))
            handle.write(f"CONECT{atom_id:5d}{neighbors}\n")
        handle.write("END\n")
    return str(path)
