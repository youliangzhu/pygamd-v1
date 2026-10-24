"""Convert GROMACS .gro/.top/.itp files to GALAMOST/pygamd XML.

Typical CLI use:
    python gro2xml.py --top topol.top --gro eq_NPT.gro --er 1.0
    python gro2xml.py --top topol.top --gro eq_NPT.gro --no-force-field --no-exclude
    
    gro2xml(gro="eq_NPT.gro", top="topol.top", er=1.0)
"""

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .read_itp_sel import merge_forcefields, read_gro_simple, read_itp_full, read_top
from shape_model.xml_io import View_, reimage

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _is_number(text: object) -> bool:
    try:
        float(text)
        return True
    except (TypeError, ValueError):
        return False



def _dedup_rows(rows: Iterable[Sequence[str]]) -> List[List[str]]:
    out: List[List[str]] = []
    seen = set()
    for row in rows:
        key = tuple(row)
        if key not in seen:
            seen.add(key)
            out.append(list(row))
    return out


def _default_combination_rule(defaults_l: Sequence[Sequence[str]]) -> str:
    if defaults_l and len(defaults_l[0]) >= 2:
        return str(defaults_l[0][1])
    return "2"


def _macro_values(macros: Dict[str, List[str]], name: str, *, context: str) -> List[str]:
    if name not in macros:
        raise RuntimeError(f"Parameter alias {name!r} not found while parsing {context}.")
    return list(macros[name])


def _atomtype_bond_name(atom_type: str, name_to_bondtype: Dict[str, str]) -> str:
    return name_to_bondtype.get(atom_type, atom_type)


def _format_row(row: Sequence[object]) -> str:
    return " ".join(str(x) for x in row)


# ---------------------------------------------------------------------------
# Atomtype and parameter lookup helpers
# ---------------------------------------------------------------------------

@dataclass
class AtomTypeInfo:
    name: str
    bond_type: str
    mass: Optional[float]
    sigma: Optional[float]
    epsilon: Optional[float]
    raw: List[str] = field(default_factory=list)


def _parse_atomtype(row: Sequence[str]) -> AtomTypeInfo:
    """Parse common GROMACS atomtypes layouts."""
    row = list(row)
    if not row:
        raise RuntimeError("empty atomtypes row")

    name = row[0]
    bond_type = name
    mass = sigma = epsilon = None

    # Common layouts:
    # name  mass charge ptype sigma epsilon                         (6 cols)
    # name  at.num mass charge ptype sigma epsilon                  (7 cols)
    # name  bond_type at.num mass charge ptype sigma epsilon        (8 cols)
    if len(row) >= 8 and not _is_number(row[1]):
        bond_type = row[1]
        if _is_number(row[3]):
            mass = float(row[3])
    elif len(row) >= 7 and _is_number(row[2]):
        if _is_number(row[2]):
            mass = float(row[2])
    elif len(row) >= 6 and _is_number(row[1]):
        mass = float(row[1])

    if len(row) >= 2 and _is_number(row[-2]):
        sigma = float(row[-2])
    if len(row) >= 1 and _is_number(row[-1]):
        epsilon = float(row[-1])

    return AtomTypeInfo(name=name, bond_type=bond_type, mass=mass, sigma=sigma, epsilon=epsilon, raw=row)


def _build_atomtype_maps(atomtypes_l: Sequence[Sequence[str]]):
    atomtypes: Dict[str, AtomTypeInfo] = {}
    name_to_bondtype: Dict[str, str] = {}
    for row in _dedup_rows(atomtypes_l):
        info = _parse_atomtype(row)
        atomtypes[info.name] = info
        name_to_bondtype[info.name] = info.bond_type
    return atomtypes, name_to_bondtype


def _lookup_mass(atom_type: str, atomtypes: Dict[str, AtomTypeInfo]) -> float:
    info = atomtypes.get(atom_type)
    if info and info.mass is not None:
        return info.mass
    raise RuntimeError(f"Particle mass for atom type {atom_type!r} was not found.")


def _lookup_bondtype(rows: Sequence[Sequence[str]], a: str, b: str, func: str):
    for row in rows:
        if len(row) >= 4 and row[2] == func:
            if (row[0] == a and row[1] == b) or (row[0] == b and row[1] == a):
                return list(row[2:])
    return None


def _lookup_angletype(rows: Sequence[Sequence[str]], a: str, b: str, c: str, func: str):
    for row in rows:
        if len(row) >= 6 and row[3] == func:
            forward = row[0] == a and row[1] == b and row[2] == c
            reverse = row[0] == c and row[1] == b and row[2] == a
            if forward or reverse:
                return list(row[3:])
    return None


def _match_dihedral_atom(pattern: str, atom: str) -> bool:
    return pattern == atom or pattern == "X"


def _lookup_dihedraltype(rows: Sequence[Sequence[str]], a: str, b: str, c: str, d: str, func: str):
    for row in rows:
        if len(row) >= 8 and row[4] == func:
            forward = (
                _match_dihedral_atom(row[0], a)
                and _match_dihedral_atom(row[1], b)
                and _match_dihedral_atom(row[2], c)
                and _match_dihedral_atom(row[3], d)
            )
            reverse = (
                _match_dihedral_atom(row[3], a)
                and _match_dihedral_atom(row[2], b)
                and _match_dihedral_atom(row[1], c)
                and _match_dihedral_atom(row[0], d)
            )
            if forward or reverse:
                return list(row[4:])
    return None


# ---------------------------------------------------------------------------
# Bonded type registries
# ---------------------------------------------------------------------------

class BondedRegistry:
    """Register bonded parameter rows and assign unique type names."""

    def __init__(self):
        self.rows: List[List[object]] = []
        self._by_base: Dict[str, List[int]] = {}
        self._counter: Dict[str, int] = {}

    def register(self, base_name: str, params: Sequence[object]) -> str:
        params_l = [str(x) for x in params]
        for idx in self._by_base.get(base_name, []):
            if [str(x) for x in self.rows[idx][1:]] == params_l:
                return str(self.rows[idx][0])

        if base_name not in self._by_base:
            name = base_name
            self._by_base[base_name] = []
            self._counter[base_name] = 0
        else:
            self._counter[base_name] += 1
            name = f"{base_name}_{self._counter[base_name]}"

        self.rows.append([name] + list(params))
        self._by_base[base_name].append(len(self.rows) - 1)
        return name


class DihedralRegistry:
    """Register dihedral parameters, merging function-9/function-4 terms."""

    def __init__(self):
        self.rows: List[List[object]] = []
        self._base_to_index: Dict[str, int] = {}
        self._counter: Dict[str, int] = {}

    def register(self, base_name: str, params: Sequence[str]) -> str:
        params = [str(x) for x in params]
        func = params[0] if params else ""

        if func == "9":
            if len(params) < 4:
                raise RuntimeError(f"Dihedral function 9 needs phi/k/multiplicity: {base_name} {params}")
            row = [base_name, "9", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0"]
            multiplicity = params[3]
            if multiplicity == "1":
                row[2], row[6] = params[2], params[1]
            elif multiplicity == "2":
                row[3], row[7] = params[2], params[1]
            elif multiplicity == "3":
                row[4], row[8] = params[2], params[1]
            elif multiplicity == "4":
                row[5], row[9] = params[2], params[1]
            else:
                # Keep unsupported multiplicities as a normal unique type.
                return self._register_unique(base_name, params)

            if base_name in self._base_to_index:
                idx = self._base_to_index[base_name]
                if str(self.rows[idx][1]) != "9":
                    raise RuntimeError(f"Dihedral type conflict for {base_name}")
                for j in range(2, 10):
                    if row[j] != "0.0":
                        self.rows[idx][j] = row[j]
                return str(self.rows[idx][0])

            self.rows.append(row)
            self._base_to_index[base_name] = len(self.rows) - 1
            return base_name

        if func == "4" and len(params) >= 4:
            row = [base_name, "4", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0"]
            if params[3] == "2":
                row[2], row[3] = params[1], params[2]

            if base_name in self._base_to_index:
                idx = self._base_to_index[base_name]
                if str(self.rows[idx][1]) != "4":
                    raise RuntimeError(f"Dihedral type conflict for {base_name}")
                if params[3] == "2":
                    self.rows[idx][2], self.rows[idx][3] = params[1], params[2]
                return str(self.rows[idx][0])

            self.rows.append(row)
            self._base_to_index[base_name] = len(self.rows) - 1
            return base_name

        return self._register_unique(base_name, params)

    def _register_unique(self, base_name: str, params: Sequence[object]) -> str:
        params_l = [str(x) for x in params]
        for idx, row in enumerate(self.rows):
            if str(row[0]).startswith(base_name) and [str(x) for x in row[1:]] == params_l:
                return str(row[0])

        if base_name not in self._counter:
            name = base_name
            self._counter[base_name] = 0
        else:
            self._counter[base_name] += 1
            name = f"{base_name}_{self._counter[base_name]}"
        self.rows.append([name] + list(params_l))
        return name


# ---------------------------------------------------------------------------
# Molecule template building
# ---------------------------------------------------------------------------

@dataclass
class MoleculeTemplate:
    name: str
    atom_types: List[str]
    charges: List[float]
    masses: List[float]
    bonds: List[List[object]] = field(default_factory=list)
    pairs: List[List[object]] = field(default_factory=list)
    constraints: List[List[object]] = field(default_factory=list)
    vsites: List[List[object]] = field(default_factory=list)
    angles: List[List[object]] = field(default_factory=list)
    dihedrals: List[List[object]] = field(default_factory=list)
    exclusions: List[List[int]] = field(default_factory=list)


def _atom_type_list(atom_rows: Sequence[Sequence[str]]) -> List[str]:
    return [row[1] for row in atom_rows]


def _type_name(types: Sequence[str], indices: Sequence[int]) -> str:
    return "-".join(types[i - 1] for i in indices)



def _normalize_bond_line(row, atom_types, ff, name_to_bondtype):
    if len(row) < 3:
        raise RuntimeError(f"Invalid bond row: {row}")
    i, j, func = row[0], row[1], row[2]
    base = _type_name(atom_types, [int(i), int(j)])

    if len(row) == 3:
        params = _lookup_bondtype(ff["bondtypes"], atom_types[int(i) - 1], atom_types[int(j) - 1], func)
        if params is None:
            mapped_i = _atomtype_bond_name(atom_types[int(i) - 1], name_to_bondtype)
            mapped_j = _atomtype_bond_name(atom_types[int(j) - 1], name_to_bondtype)
            params = _lookup_bondtype(ff["bondtypes"], mapped_i, mapped_j, func)
            if params is not None:
                base = f"{mapped_i}-{mapped_j}"
        if params is None:
            raise RuntimeError(f"Bond parameters not found for row: {row}")
        return base, [i, j] + params

    if len(row) == 4:
        values = _macro_values(ff["defines"], row[3], context="bond")
        return base, [i, j, func] + values

    return base, list(row)


def _normalize_constraint_line(row, atom_types, ff, name_to_bondtype):
    if len(row) < 3:
        raise RuntimeError(f"Invalid constraint row: {row}")
    i, j, func = row[0], row[1], row[2]
    base = _type_name(atom_types, [int(i), int(j)])
    if len(row) == 3:
        params = _lookup_bondtype(ff["constrainttypes"], atom_types[int(i) - 1], atom_types[int(j) - 1], func)
        if params is None:
            mapped_i = _atomtype_bond_name(atom_types[int(i) - 1], name_to_bondtype)
            mapped_j = _atomtype_bond_name(atom_types[int(j) - 1], name_to_bondtype)
            params = _lookup_bondtype(ff["constrainttypes"], mapped_i, mapped_j, func)
            if params is not None:
                base = f"{mapped_i}-{mapped_j}"
        if params is None:
            raise RuntimeError(f"Constraint parameters not found for row: {row}")
        return base, [i, j] + params
    if len(row) == 4 and row[3] in ff["defines"]:
        return base, [i, j, func] + _macro_values(ff["defines"], row[3], context="constraint")
    return base, list(row)


def _normalize_angle_line(row, atom_types, ff, name_to_bondtype):
    if len(row) < 4:
        raise RuntimeError(f"Invalid angle row: {row}")
    i, j, k, func = row[0], row[1], row[2], row[3]
    base = _type_name(atom_types, [int(i), int(j), int(k)])
    if len(row) == 4:
        params = _lookup_angletype(ff["angletypes"], atom_types[int(i) - 1], atom_types[int(j) - 1], atom_types[int(k) - 1], func)
        if params is None:
            mapped = [_atomtype_bond_name(atom_types[int(x) - 1], name_to_bondtype) for x in (i, j, k)]
            params = _lookup_angletype(ff["angletypes"], mapped[0], mapped[1], mapped[2], func)
            if params is not None:
                base = "-".join(mapped)
        if params is None:
            raise RuntimeError(f"Angle parameters not found for row: {row}")
        return base, [i, j, k] + params
    if len(row) == 5:
        values = _macro_values(ff["defines"], row[4], context="angle")
        return base, [i, j, k, func] + values
    return base, list(row)


def _normalize_dihedral_line(row, atom_types, ff, name_to_bondtype):
    if len(row) < 5:
        raise RuntimeError(f"Invalid dihedral row: {row}")
    i, j, k, l, func = row[0], row[1], row[2], row[3], row[4]
    base = _type_name(atom_types, [int(i), int(j), int(k), int(l)]) + f"_F{func}"
    if len(row) == 5:
        params = _lookup_dihedraltype(
            ff["dihedraltypes"],
            atom_types[int(i) - 1], atom_types[int(j) - 1], atom_types[int(k) - 1], atom_types[int(l) - 1], func,
        )
        if params is None:
            mapped = [_atomtype_bond_name(atom_types[int(x) - 1], name_to_bondtype) for x in (i, j, k, l)]
            params = _lookup_dihedraltype(ff["dihedraltypes"], mapped[0], mapped[1], mapped[2], mapped[3], func)
            if params is not None:
                base = "-".join(mapped) + f"_F{func}"
        if params is None:
            raise RuntimeError(f"Dihedral parameters not found for row: {row}")
        return base, [i, j, k, l] + params
    return base, list(row)


def build_templates(mols_d, molecules_l, ff, er: float):
    atomtypes, name_to_bondtype = _build_atomtype_maps(ff["atomtypes"])
    charge_scale = math.sqrt(138.935 / er)

    bond_reg = BondedRegistry()
    constraint_reg = BondedRegistry()
    vsite_reg = BondedRegistry()
    angle_reg = BondedRegistry()
    dihedral_reg = DihedralRegistry()

    templates: Dict[str, MoleculeTemplate] = {}

    for mol_name, _mol_count in molecules_l:
        if mol_name in templates:
            continue
        if mol_name not in mols_d:
            raise RuntimeError(f"Molecule type {mol_name!r} not found in included topology files.")

        mol_data = mols_d[mol_name]
        atom_rows = mol_data.get("atoms", [])
        if not atom_rows:
            raise RuntimeError(f"Molecule type {mol_name!r} has no [atoms] section.")

        atom_types = _atom_type_list(atom_rows)
        charges: List[float] = []
        masses: List[float] = []
        for row in atom_rows:
            if len(row) < 7:
                raise RuntimeError(f"Invalid [atoms] row in {mol_name}: {row}")
            charges.append(float(row[6]) * charge_scale)
            if len(row) >= 8 and _is_number(row[7]):
                masses.append(float(row[7]))
            else:
                masses.append(_lookup_mass(row[1], atomtypes))

        tmpl = MoleculeTemplate(name=mol_name, atom_types=atom_types, charges=charges, masses=masses)

        for row in mol_data.get("bonds", []):
            base, norm = _normalize_bond_line(row, atom_types, ff, name_to_bondtype)
            type_name = bond_reg.register(base, norm[2:])
            tmpl.bonds.append([type_name, int(norm[0]) - 1, int(norm[1]) - 1])

        for row in mol_data.get("pairs", []):
            if len(row) >= 3 and row[2] == "1":
                tmpl.pairs.append(["p1", int(row[0]) - 1, int(row[1]) - 1])

        for row in mol_data.get("constraints", []):
            base, norm = _normalize_constraint_line(row, atom_types, ff, name_to_bondtype)
            type_name = constraint_reg.register(base, norm[2:])
            tmpl.constraints.append([type_name, int(norm[0]) - 1, int(norm[1]) - 1])

        for row in mol_data.get("virtual_sites3", []):
            if len(row) < 7:
                raise RuntimeError(f"Invalid [virtual_sites3] row in {mol_name}: {row}")
            row = list(row)
            if len(row) == 7:
                row.append("0.0")
            base = _type_name(atom_types, [int(row[0]), int(row[1]), int(row[2]), int(row[3])])
            type_name = vsite_reg.register(base, row[4:])
            tmpl.vsites.append([type_name, int(row[0]) - 1, int(row[1]) - 1, int(row[2]) - 1, int(row[3]) - 1])

        for row in mol_data.get("angles", []):
            base, norm = _normalize_angle_line(row, atom_types, ff, name_to_bondtype)
            type_name = angle_reg.register(base, norm[3:])
            tmpl.angles.append([type_name, int(norm[0]) - 1, int(norm[1]) - 1, int(norm[2]) - 1])

        for row in mol_data.get("dihedrals", []):
            base, norm = _normalize_dihedral_line(row, atom_types, ff, name_to_bondtype)
            type_name = dihedral_reg.register(base, norm[4:])
            tmpl.dihedrals.append([type_name, int(norm[0]) - 1, int(norm[1]) - 1, int(norm[2]) - 1, int(norm[3]) - 1])

        for row in mol_data.get("exclusions", []):
            if len(row) < 2:
                continue
            i0 = int(row[0]) - 1
            for item in row[1:]:
                tmpl.exclusions.append([i0, int(item) - 1])

        templates[mol_name] = tmpl

    registries = {
        "bond": bond_reg,
        "constraint": constraint_reg,
        "vsite": vsite_reg,
        "angle": angle_reg,
        "dihedral": dihedral_reg,
        "atomtypes": atomtypes,
    }
    return templates, registries


# ---------------------------------------------------------------------------
# Force-field output
# ---------------------------------------------------------------------------

def _pair_params_from_forcefield(ff, atom_type_order, atomtypes, defaults_l):
    comb = _default_combination_rule(defaults_l)
    macros = ff["defines"]

    explicit_rows = list(ff.get("nonbond_params", [])) + list(ff.get("pairtypes", []))
    explicit_pairs: List[Tuple[str, str]] = []
    explicit_values: List[List[object]] = []

    for row in explicit_rows:
        if len(row) >= 5 and row[0] in atom_type_order and row[1] in atom_type_order:
            p0, p1 = row[0], row[1]
            c6, c12 = float(row[-2]), float(row[-1])
        elif len(row) == 4 and row[0] in atom_type_order and row[1] in atom_type_order:
            values = _macro_values(macros, row[3], context="nonbond_params/pairtypes")
            if len(values) < 2:
                raise RuntimeError(f"nonbond macro needs two values: {row}")
            p0, p1 = row[0], row[1]
            c6, c12 = float(values[0]), float(values[1])
        else:
            continue

        if c6 != 0.0:
            sigma6 = c12 / c6
            sigma = sigma6 ** (1.0 / 6.0)
            epsilon = c6 / sigma6 / 4.0
            if comb == "2":
                value = [c12, c6]
            elif comb == "1":
                value = [epsilon, sigma]
            else:
                value = [epsilon, sigma]
        else:
            value = [0, 0.47]
        explicit_pairs.append((p0, p1))
        explicit_values.append(value)

    if explicit_pairs:
        rows = [[a, b, value[0], value[1], 1.0] for (a, b), value in zip(explicit_pairs, explicit_values)]
        covered = {tuple(sorted(pair)) for pair in explicit_pairs}
        for i, a in enumerate(atom_type_order):
            for b in atom_type_order[i:]:
                if tuple(sorted((a, b))) not in covered:
                    print(f"[warning] no explicit nonbond/pairtype parameters for {a}-{b}")
        return rows

    rows = []
    for i, a in enumerate(atom_type_order):
        for b in atom_type_order[i:]:
            info_a = atomtypes.get(a)
            info_b = atomtypes.get(b)
            if info_a is None or info_b is None or info_a.sigma is None or info_b.sigma is None or info_a.epsilon is None or info_b.epsilon is None:
                print(f"[warning] missing atomtype sigma/epsilon for pair {a}-{b}; skipped")
                continue
            if comb == "3":
                sigma = (info_a.sigma * info_b.sigma) ** 0.5
            else:
                sigma = (info_a.sigma + info_b.sigma) / 2.0
            epsilon = (info_a.epsilon * info_b.epsilon) ** 0.5
            rows.append([a, b, epsilon, sigma, 1.0])
    return rows


def write_force_field(filename: str, ff, templates, registries, defaults_l):
    atom_type_order: List[str] = []
    for tmpl in templates.values():
        for atom_type in tmpl.atom_types:
            if atom_type not in atom_type_order:
                atom_type_order.append(atom_type)

    pair_params = _pair_params_from_forcefield(ff, atom_type_order, registries["atomtypes"], defaults_l)

    def write_block(out, label, rows):
        if not rows:
            return
        out.write(f"<{label}>\n")
        for row in rows:
            out.write(_format_row(row) + "\n")
        out.write(f"</{label}>\n\n")

    # Match old gro_to_xml force-field ordering and row reshuffling.
    constraint_rows = [[r[0], r[2], r[1]] for r in registries["constraint"].rows if len(r) >= 3]
    vsite_rows = [[r[0], r[2], r[3], r[4], r[1]] for r in registries["vsite"].rows if len(r) >= 5]
    bond_rows = [[r[0], r[3], r[2], r[1]] for r in registries["bond"].rows if len(r) >= 4]

    angle_rows = []
    for r in registries["angle"].rows:
        if len(r) == 4:
            angle_rows.append([r[0], r[3], r[2], r[1]])
        elif len(r) == 6:
            angle_rows.append([r[0], r[3], r[2], r[5], r[4], r[1]])
        elif len(r) > 1:
            angle_rows.append(r)

    dihedral_rows = []
    for r in registries["dihedral"].rows:
        if len(r) >= 2:
            dihedral_rows.append([r[0]] + list(r[2:]) + [r[1]])

    with open(filename, "w", encoding="utf-8") as out:
        write_block(out, "pair_params", pair_params)
        write_block(out, "constraint_params", constraint_rows)
        write_block(out, "vsite_params", vsite_rows)
        write_block(out, "bond_params", bond_rows)
        write_block(out, "angle_params", angle_rows)
        write_block(out, "dihedral_params", dihedral_rows)


# ---------------------------------------------------------------------------
# Expansion and conversion
# ---------------------------------------------------------------------------

def _shift_positions_to_galamost(pos: List[List[float]], box: Sequence[float]):
    out = [[p[0] - box[0] / 2.0, p[1] - box[1] / 2.0, p[2] - box[2] / 2.0] for p in pos]
    images = [reimage(p, box) for p in out]
    return out, images


def _offset_rows(rows: Sequence[Sequence[object]], offset: int) -> List[List[object]]:
    out = []
    for row in rows:
        out.append([row[0]] + [int(x) + offset for x in row[1:]])
    return out


def expand_system(templates: Dict[str, MoleculeTemplate], molecules_l):
    types: List[str] = []
    charges: List[float] = []
    masses: List[float] = []
    bonds: List[List[object]] = []
    pairs: List[List[object]] = []
    constraints: List[List[object]] = []
    vsites: List[List[object]] = []
    angles: List[List[object]] = []
    dihedrals: List[List[object]] = []
    exclusions: List[List[int]] = []
    first_offsets: Dict[str, int] = {}

    offset = 0
    for mol_name, mol_count_s in molecules_l:
        mol_count = int(mol_count_s)
        tmpl = templates[mol_name]
        n = len(tmpl.atom_types)
        if mol_name not in first_offsets:
            first_offsets[mol_name] = offset
        for _ in range(mol_count):
            types.extend(tmpl.atom_types)
            charges.extend(tmpl.charges)
            masses.extend(tmpl.masses)
            bonds.extend(_offset_rows(tmpl.bonds, offset))
            pairs.extend(_offset_rows(tmpl.pairs, offset))
            constraints.extend(_offset_rows(tmpl.constraints, offset))
            vsites.extend(_offset_rows(tmpl.vsites, offset))
            angles.extend(_offset_rows(tmpl.angles, offset))
            dihedrals.extend(_offset_rows(tmpl.dihedrals, offset))
            exclusions.extend([[i + offset, j + offset] for i, j in tmpl.exclusions])
            offset += n

    return {
        "type": types,
        "charge": charges,
        "mass": masses,
        "bond": bonds,
        "pair": pairs,
        "constraint": constraints,
        "vsite": vsites,
        "angle": angles,
        "dihedral": dihedrals,
        "exclusions": exclusions,
        "first_offsets": first_offsets,
    }


def _load_topology(top_path: str):
    files_l, top_defaults, molecules_l = read_top(top_path)
    all_files: List[str] = []
    for path in files_l + [top_path]:
        abs_path = os.path.abspath(path)
        if abs_path not in all_files:
            all_files.append(abs_path)

    mols_d = {}
    ff_parts = []
    for path in all_files:
        if not os.path.isfile(path):
            print(f"[warning] topology include not found, skipped: {path}")
            continue
        mol_d, ff_d = read_itp_full(path)
        for key, value in mol_d.items():
            mols_d[key] = value
        ff_parts.append(ff_d)

    ff = merge_forcefields(ff_parts)
    if top_defaults:
        ff["defaults"] = top_defaults + ff.get("defaults", [])
    ff["atomtypes"] = _dedup_rows(ff.get("atomtypes", []))
    ff["bondtypes"] = _dedup_rows(ff.get("bondtypes", []))
    ff["constrainttypes"] = _dedup_rows(ff.get("constrainttypes", []))
    ff["angletypes"] = _dedup_rows(ff.get("angletypes", []))
    ff["dihedraltypes"] = _dedup_rows(ff.get("dihedraltypes", []))
    ff["nonbond_params"] = _dedup_rows(ff.get("nonbond_params", []))
    ff["pairtypes"] = _dedup_rows(ff.get("pairtypes", []))

    if not molecules_l:
        raise RuntimeError(f"No [molecules] section found in top file: {top_path}")
    return mols_d, ff, molecules_l


def gro2xml(
    *,
    gro: str,
    top: str,
    er: float = 1.0,
    output_dir: Optional[str] = None,
    xml: bool = True,
    force_field: bool = True,
    exclude: bool = True,
    molecule_xml: bool = True,
):
    if er <= 0:
        raise RuntimeError("Dielectric constant er must be positive.")

    top_path = os.path.abspath(top)
    gro_path = os.path.abspath(gro)
    if output_dir is None:
        output_dir = os.path.dirname(gro_path)
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    mols_d, ff, molecules_l = _load_topology(top_path)
    num, box, pos_raw, vel, _ = read_gro_simple(gro_path)
    pos, images = _shift_positions_to_galamost(pos_raw, box)

    templates, registries = build_templates(mols_d, molecules_l, ff, er=er)
    expanded = expand_system(templates, molecules_l)

    expected_natoms = len(expanded["type"])
    if expected_natoms != num:
        raise RuntimeError(f"Atom number mismatch: .gro has {num}, topology expands to {expected_natoms}.")

    out_stem = os.path.splitext(os.path.basename(gro_path))[0]
    outputs: Dict[str, str] = {}

    base_info = {
        "position": pos,
        "velocity": vel,
        "image": images,
        "box": box,
        "type": expanded["type"],
        "charge": expanded["charge"],
        "mass": expanded["mass"],
        "bond": expanded["bond"],
        "pair": expanded["pair"],
        "constraint": expanded["constraint"],
        "vsite": expanded["vsite"],
        "angle": expanded["angle"],
        "dihedral": expanded["dihedral"],
        "time_step": 0,
        "dimensions": 3,
    }

    if xml:
        out_xml = os.path.join(output_dir, out_stem + ".xml")
        View_(out_xml, base_info)
        outputs["xml"] = out_xml
        print(os.path.basename(out_xml), "success!")

    if molecule_xml:
        written = set()
        for mol_name, _mol_count in molecules_l:
            if mol_name in written:
                continue
            written.add(mol_name)
            tmpl = templates[mol_name]
            start = expanded["first_offsets"][mol_name]
            end = start + len(tmpl.atom_types)
            mol_info = {
                "position": pos[start:end],
                "velocity": vel[start:end] if len(vel) == len(pos) else [],
                "image": images[start:end],
                "box": box,
                "type": tmpl.atom_types,
                "charge": tmpl.charges,
                "mass": tmpl.masses,
                "bond": tmpl.bonds,
                "pair": tmpl.pairs,
                "constraint": tmpl.constraints,
                "vsite": tmpl.vsites,
                "angle": tmpl.angles,
                "dihedral": tmpl.dihedrals,
                "time_step": 0,
                "dimensions": 3,
            }
            mol_xml = os.path.join(output_dir, mol_name + ".xml")
            View_(mol_xml, mol_info)
            outputs.setdefault("molecule_xml", []).append(mol_xml)
            print(os.path.basename(mol_xml), "success!")

    if force_field:
        ff_path = os.path.join(output_dir, out_stem + ".force_field")
        defaults_l = ff.get("defaults", [])
        write_force_field(ff_path, ff, templates, registries, defaults_l)
        outputs["force_field"] = ff_path
        print(os.path.basename(ff_path), "success!")

    if exclude and expanded["exclusions"]:
        excl_path = os.path.join(output_dir, "exclusion_list.dat")
        with open(excl_path, "w", encoding="utf-8") as out:
            for i, j in expanded["exclusions"]:
                out.write(f"{i} {j}\n")
        outputs["exclusion_list"] = excl_path
        print("exclusion_list.dat success!")

    print("success!")
    return outputs


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

__all__ = ["gro2xml", "main"]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Convert GROMACS gro/top/itp files to GALAMOST/pygamd XML.")
    parser.add_argument("--top", default=None, help="GROMACS .top file")
    parser.add_argument("--gro", default=None, help="GROMACS .gro file")
    parser.add_argument("--er", type=float, default=1.0, help="dielectric constant; default 1.0")
    parser.add_argument("--out-dir", default=None, help="output directory; default is the .gro directory")

    parser.set_defaults(
        write_xml=True,
        write_forcefield=True,
        write_exclusions=True,
        write_molecule_xml=True,
    )

    parser.add_argument("--xml", dest="write_xml", action="store_true", help="write main XML")
    parser.add_argument("--no-xml", dest="write_xml", action="store_false", help="do not write main XML")

    parser.add_argument("--force-field", dest="write_forcefield", action="store_true", help="write .force_field")
    parser.add_argument("--no-force-field", dest="write_forcefield", action="store_false", help="do not write .force_field")

    parser.add_argument("--exclude", dest="write_exclusions", action="store_true", help="write exclusion_list.dat")
    parser.add_argument("--no-exclude", dest="write_exclusions", action="store_false", help="do not write exclusion_list.dat")

    parser.add_argument("--molecule-xml", "--mol-xml", dest="write_molecule_xml", action="store_true",
                        help="write one XML file for each molecule type")
    parser.add_argument("--no-molecule-xml", "--no-mol-xml", dest="write_molecule_xml", action="store_false",
                        help="do not write per-molecule XML files")
    return parser.parse_args(argv)


def main(argv=None):
    args = _parse_args(argv)
    top_path = args.top
    gro_path = args.gro
    er = args.er

    if not top_path or not os.path.isfile(top_path):
        raise RuntimeError("Please set a valid top file with --top=....")
    if not gro_path or not os.path.isfile(gro_path):
        raise RuntimeError("Please set a valid gro file with --gro=....")
    if er <= 0:
        raise RuntimeError("Dielectric constant --er must be positive.")

    print("dielectric constant =", er)
    return gro2xml(
        gro=gro_path,
        top=top_path,
        er=er,
        output_dir=args.out_dir,
        xml=args.write_xml,
        force_field=args.write_forcefield,
        exclude=args.write_exclusions,
        molecule_xml=args.write_molecule_xml,
    )


if __name__ == "__main__":
    main()
