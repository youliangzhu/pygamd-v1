"""Read the small subset of TPR topology needed by the AA-to-CG mapper.

Two inputs are supported:

* text produced by ``gmx dump`` / ``gmxdump``; and
* a binary ``.tpr`` read directly with MDAnalysis.

The direct reader is important on platforms where a GROMACS executable is not
installed (notably a normal Windows Python environment).  Both readers return
the same :class:`TprDumpTopology` value, so mapping and RDF code do not need to
care which backend supplied atom names, masses, and charges.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import warnings
from typing import TypeAlias


_FLOAT = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?|[-+]?\d+(?:[eE][-+]?\d+)?"


@dataclass(frozen=True)
class TprDumpAtom:
    index: int
    name: str | None
    mass: float
    charge: float
    atom_type: int | None = None


@dataclass(frozen=True)
class TprDumpTopology:
    path: str
    name: str | None
    total_atoms: int | None
    molecules: int | None
    atoms_per_molecule: int | None
    atoms: list[TprDumpAtom]
    box: tuple[float, float, float] | None
    reader: str = "gmx_dump"


TprTopologyInput: TypeAlias = str | Path | TprDumpTopology


def _import_mdanalysis():
    try:
        import MDAnalysis as mda  # type: ignore[import-not-found]
    except ImportError as exc:
        raise RuntimeError(
            "MDAnalysis is required to read a binary TPR without GROMACS; "
            "install pygamd[aa2cg] (or pygamd[rdf]), or install MDAnalysis>=2.9."
        ) from exc
    except Exception as exc:
        raise RuntimeError(
            "Failed to import MDAnalysis while preparing the binary TPR fallback. "
            f"The MDAnalysis/NumPy/SciPy environment may be inconsistent: {exc}"
        ) from exc
    return mda


def _atom_attribute(atoms, name: str, *, path: Path) -> list[object]:
    """Return a required MDAnalysis atom attribute with an actionable error."""

    try:
        values = getattr(atoms, name)
    except Exception as exc:
        raise ValueError(f"TPR topology {path} does not provide required atom attribute {name!r}") from exc
    result = list(values)
    if len(result) != int(atoms.n_atoms):
        raise ValueError(
            f"TPR topology {path} returned {len(result)} {name} values for {atoms.n_atoms} atoms"
        )
    return result


def parse_mdanalysis_tpr(path: str | Path) -> TprDumpTopology:
    """Read a homogeneous molecular TPR directly through MDAnalysis.

    The current mapper intentionally models a repeated molecule template.  A
    mixed TPR (for example solute plus explicit solvent) cannot be represented
    by that contract and is rejected instead of silently assigning the first
    molecule's masses to unrelated atoms.
    """

    tpr_path = Path(path).expanduser().resolve()
    if not tpr_path.is_file():
        raise FileNotFoundError(f"TPR topology not found: {tpr_path}")
    mda = _import_mdanalysis()
    try:
        # A TPR contains topology but MDAnalysis does not expose it as a
        # coordinate trajectory.  Silence that expected, harmless warning.
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=r"No coordinate reader found for .*\.tpr\. Skipping this file\.",
                category=UserWarning,
            )
            universe = mda.Universe(str(tpr_path))
    except Exception as exc:
        raise RuntimeError(f"MDAnalysis could not read binary TPR topology {tpr_path}: {exc}") from exc

    aa_atoms = universe.atoms
    total_atoms = int(aa_atoms.n_atoms)
    if total_atoms <= 0:
        raise ValueError(f"TPR topology contains no atoms: {tpr_path}")

    names = _atom_attribute(aa_atoms, "names", path=tpr_path)
    masses = _atom_attribute(aa_atoms, "masses", path=tpr_path)
    charges = _atom_attribute(aa_atoms, "charges", path=tpr_path)
    molnums = _atom_attribute(aa_atoms, "molnums", path=tpr_path)
    try:
        moltypes = _atom_attribute(aa_atoms, "moltypes", path=tpr_path)
    except ValueError:
        moltypes = [None] * total_atoms
    try:
        atom_types = _atom_attribute(aa_atoms, "types", path=tpr_path)
    except ValueError:
        atom_types = [None] * total_atoms

    # Preserve the TPR atom order.  Offset-based mapping requires every
    # molecule to occupy one contiguous, equally-sized atom block.
    ordered_molnums: list[object] = []
    groups: list[list[int]] = []
    group_by_molnum: dict[object, list[int]] = {}
    for atom_index, raw_molnum in enumerate(molnums):
        molnum = raw_molnum.item() if hasattr(raw_molnum, "item") else raw_molnum
        try:
            group = group_by_molnum.get(molnum)
        except TypeError as exc:
            raise ValueError(f"TPR topology has an unusable molecule id at atom {atom_index + 1}") from exc
        if group is None:
            ordered_molnums.append(molnum)
            group = []
            group_by_molnum[molnum] = group
            groups.append(group)
        group.append(atom_index)

    if not groups:
        raise ValueError(f"TPR topology has no molecule assignments: {tpr_path}")
    atoms_per_molecule = len(groups[0])
    if atoms_per_molecule <= 0 or any(len(group) != atoms_per_molecule for group in groups):
        counts = ", ".join(str(len(group)) for group in groups[:8])
        raise ValueError(
            "AA-to-CG mapping requires one repeated molecule template, but the TPR "
            f"contains different molecule sizes (first counts: {counts})"
        )
    for molecule_index, group in enumerate(groups):
        expected_start = molecule_index * atoms_per_molecule
        if group != list(range(expected_start, expected_start + atoms_per_molecule)):
            raise ValueError(
                "AA-to-CG mapping requires contiguous molecule atom blocks, but the "
                f"TPR atom order is interleaved near molecule {ordered_molnums[molecule_index]!r}"
            )

    nonempty_moltypes = {str(value) for value in moltypes if value is not None and str(value)}
    if len(nonempty_moltypes) > 1:
        listed = ", ".join(sorted(nonempty_moltypes)[:8])
        raise ValueError(
            "AA-to-CG mapping currently supports a homogeneous repeated molecule TPR; "
            f"found molecule types: {listed}"
        )

    template_indices = groups[0]
    # Detect a surprising heterogeneous template even when moltype metadata is
    # absent or unreliable.  Float tolerances cover normal TPR precision.
    for group in groups[1:]:
        for local_index, (reference, candidate) in enumerate(zip(template_indices, group), start=1):
            if str(names[reference]) != str(names[candidate]):
                raise ValueError(
                    "TPR molecules do not share one atom-name template; "
                    f"mismatch at local atom {local_index}"
                )
            if abs(float(masses[reference]) - float(masses[candidate])) > 1.0e-5:
                raise ValueError(
                    "TPR molecules do not share one mass template; "
                    f"mismatch at local atom {local_index}"
                )
            if abs(float(charges[reference]) - float(charges[candidate])) > 1.0e-5:
                raise ValueError(
                    "TPR molecules do not share one charge template; "
                    f"mismatch at local atom {local_index}"
                )

    type_ids: dict[str, int] = {}
    template: list[TprDumpAtom] = []
    for local_index, global_index in enumerate(template_indices, start=1):
        raw_type = atom_types[global_index]
        atom_type = None
        if raw_type is not None:
            type_key = str(raw_type)
            if type_key not in type_ids:
                type_ids[type_key] = len(type_ids)
            atom_type = type_ids[type_key]
        template.append(
            TprDumpAtom(
                index=local_index,
                name=str(names[global_index]) if names[global_index] is not None else None,
                mass=float(masses[global_index]),
                charge=float(charges[global_index]),
                atom_type=atom_type,
            )
        )

    molecule_name = next(iter(nonempty_moltypes), None)
    if molecule_name is None:
        try:
            first_resname = aa_atoms[template_indices[0]].resname
        except Exception:
            first_resname = None
        molecule_name = str(first_resname) if first_resname else None

    return TprDumpTopology(
        path=str(tpr_path),
        name=molecule_name,
        total_atoms=total_atoms,
        molecules=len(groups),
        atoms_per_molecule=atoms_per_molecule,
        atoms=template,
        box=None,
        reader="mdanalysis",
    )


def parse_tpr_dump(path: str | Path) -> TprDumpTopology:
    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    total_atoms = None
    total_match = re.search(r"^\s*#atoms\s*=\s*(\d+)\s*$", text, re.MULTILINE)
    if not total_match:
        total_match = re.search(r"^\s*natoms\s*=\s*(\d+)\s*$", text, re.MULTILINE)
    if total_match:
        total_atoms = int(total_match.group(1))

    molecules = None
    mol_match = re.search(r"^\s*#molecules\s*=\s*(\d+)\s*$", text, re.MULTILINE)
    if mol_match:
        molecules = int(mol_match.group(1))

    name = None
    name_match = re.search(r'^\s*name="([^"]+)"\s*$', text, re.MULTILINE)
    if name_match:
        name = name_match.group(1)
    else:
        moltype_match = re.search(r'^\s*moltype\s*=\s*\d+\s+"([^"]+)"\s*$', text, re.MULTILINE)
        if moltype_match:
            name = moltype_match.group(1)

    atom_data: dict[int, dict[str, float | int | str | None]] = {}
    mass_re = re.compile(
        rf"atom\[\s*(\d+)\s*\]=\{{.*?type=\s*(\d+).*?\bm=\s*({_FLOAT}),\s*q=\s*({_FLOAT})",
    )
    name_re = re.compile(r'atom\[\s*(\d+)\s*\]=\{name="([^"]+)"\}')
    for line in lines:
        mass_match = mass_re.search(line)
        if mass_match:
            idx0 = int(mass_match.group(1))
            atom_data.setdefault(idx0 + 1, {})
            atom_data[idx0 + 1].update(
                {
                    "atom_type": int(mass_match.group(2)),
                    "mass": float(mass_match.group(3)),
                    "charge": float(mass_match.group(4)),
                }
            )
            continue
        name_match = name_re.search(line)
        if name_match:
            idx0 = int(name_match.group(1))
            atom_data.setdefault(idx0 + 1, {})["name"] = name_match.group(2)

    atoms: list[TprDumpAtom] = []
    for idx in sorted(atom_data):
        item = atom_data[idx]
        if "mass" not in item or "charge" not in item:
            continue
        atoms.append(
            TprDumpAtom(
                index=idx,
                name=item.get("name") if isinstance(item.get("name"), str) else None,
                mass=float(item["mass"]),
                charge=float(item["charge"]),
                atom_type=int(item["atom_type"]) if item.get("atom_type") is not None else None,
            )
        )

    atoms_per_molecule = len(atoms) or None
    if atoms_per_molecule is None and total_atoms is not None and molecules:
        atoms_per_molecule = total_atoms // molecules

    box_rows: dict[int, tuple[float, float, float]] = {}
    box_re = re.compile(rf"^\s*box\[\s*(\d+)\s*\]=\{{\s*({_FLOAT}),\s*({_FLOAT}),\s*({_FLOAT})\s*\}}")
    for line in lines:
        match = box_re.match(line)
        if match:
            row = int(match.group(1))
            box_rows[row] = (float(match.group(2)), float(match.group(3)), float(match.group(4)))
    box = None
    if all(i in box_rows for i in (0, 1, 2)):
        # Orthorhombic target; keep diagonal lengths.
        box = (box_rows[0][0], box_rows[1][1], box_rows[2][2])

    return TprDumpTopology(
        path=str(path),
        name=name,
        total_atoms=total_atoms,
        molecules=molecules,
        atoms_per_molecule=atoms_per_molecule,
        atoms=atoms,
        box=box,
        reader="gmx_dump",
    )


def load_tpr_topology(source: TprTopologyInput) -> TprDumpTopology:
    """Load an in-memory topology, a binary TPR, or a text GROMACS dump."""

    if isinstance(source, TprDumpTopology):
        return source
    path = Path(source)
    if path.suffix.lower() == ".tpr":
        return parse_mdanalysis_tpr(path)
    return parse_tpr_dump(path)


def render_tpr_dump(topology: TprDumpTopology) -> str:
    """Render a normalized, parseable text dump for optional diagnostics."""

    name = topology.name or "MOL"
    lines = [
        f"natoms = {topology.total_atoms or 0}",
        "topology:",
        f'   name="{name}"',
        f"   #atoms = {topology.total_atoms or 0}",
        "   molblock (0):",
        f'      moltype = 0 "{name}"',
        f"      #molecules = {topology.molecules or 0}",
        "   moltype (0):",
        f'      name="{name}"',
        "      atoms:",
        f"         atom ({len(topology.atoms)}):",
    ]
    for atom in topology.atoms:
        idx0 = atom.index - 1
        atom_type = atom.atom_type if atom.atom_type is not None else 0
        lines.append(
            f"            atom[{idx0:6d}]={{type={atom_type:3d}, typeB={atom_type:3d}, "
            f"ptype=Atom, m={atom.mass:.12e}, q={atom.charge:.12e}}}"
        )
    lines.append(f"         atom ({len(topology.atoms)}):")
    for atom in topology.atoms:
        idx0 = atom.index - 1
        atom_name = (atom.name or f"A{atom.index}").replace('"', "'")
        lines.append(f'            atom[{idx0}]={{name="{atom_name}"}}')
    if topology.box is not None:
        x, y, z = topology.box
        lines.extend(
            [
                "box (3x3):",
                f"   box[0]={{ {x:.12e}, 0.000000000000e+00, 0.000000000000e+00}}",
                f"   box[1]={{ 0.000000000000e+00, {y:.12e}, 0.000000000000e+00}}",
                f"   box[2]={{ 0.000000000000e+00, 0.000000000000e+00, {z:.12e}}}",
            ]
        )
    return "\n".join(lines) + "\n"
