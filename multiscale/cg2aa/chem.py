"""Small chemistry graph toolkit used as an RDKit replacement.

This module implements the subset of the RDKit API that this project uses:
SMILES/SMARTS parsing, simple substructure matching, reaction SMARTS mapping,
hydrogen addition/removal, connected-fragment extraction, feature matching and
initial 3D coordinates for downstream placement and refinement.

It is intentionally lightweight. It supports the organic SMILES/SMARTS syntax
used by the bundled back-mapping scripts and OPLS templates, but it is not a
drop-in replacement for full RDKit cheminformatics.
"""

from __future__ import annotations

from dataclasses import dataclass
import itertools
import math
import re
from collections.abc import Iterable, Sequence


class BondType:
    SINGLE = "SINGLE"
    DOUBLE = "DOUBLE"
    TRIPLE = "TRIPLE"
    AROMATIC = "AROMATIC"
    ANY = "ANY"


BOND_FROM_SYMBOL = {
    "-": BondType.SINGLE,
    "=": BondType.DOUBLE,
    "#": BondType.TRIPLE,
    ":": BondType.AROMATIC,
    "~": BondType.ANY,
}

BOND_TO_SYMBOL = {
    BondType.SINGLE: "",
    BondType.DOUBLE: "=",
    BondType.TRIPLE: "#",
    BondType.AROMATIC: ":",
    BondType.ANY: "~",
}

ATOMIC_NUMBERS = {
    "H": 1,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Br": 35,
    "I": 53,
}

ELEMENTS = {v: k for k, v in ATOMIC_NUMBERS.items()}

ATOMIC_MASSES = {
    "H": 1.008,
    "B": 10.81,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "P": 30.974,
    "S": 32.06,
    "Cl": 35.45,
    "Br": 79.904,
    "I": 126.904,
}

ORGANIC_TWO_LETTER = {"Cl", "Br"}
ORGANIC_ONE_LETTER = set("BCNOPSFHI")
AROMATIC_TO_ELEMENT = {
    "b": "B",
    "c": "C",
    "n": "N",
    "o": "O",
    "p": "P",
    "s": "S",
}


def _bond_order(bond_type: str, aromatic_as_one: bool = True) -> float:
    if bond_type == BondType.DOUBLE:
        return 2.0
    if bond_type == BondType.TRIPLE:
        return 3.0
    if bond_type == BondType.AROMATIC:
        return 1.0 if aromatic_as_one else 1.5
    return 1.0


def _normal_symbol(symbol: str) -> tuple[str, bool]:
    if not symbol:
        return "*", False
    if symbol in AROMATIC_TO_ELEMENT:
        return AROMATIC_TO_ELEMENT[symbol], True
    if len(symbol) == 1:
        return symbol.upper(), False
    return symbol[0].upper() + symbol[1:].lower(), False


class Atom:
    def __init__(
        self,
        symbol: str = "C",
        aromatic: bool = False,
        atom_map_num: int = 0,
        query: AtomQuery | None = None,
    ):
        symbol, inferred_aromatic = _normal_symbol(symbol)
        self.symbol = symbol
        self.aromatic = aromatic or inferred_aromatic
        self.atom_map_num = atom_map_num
        self.query = query
        self.props: dict[str, object] = {}
        self._idx = -1
        self._mol: RWMol | None = None

    def copy(self) -> "Atom":
        atom = Atom(self.symbol, self.aromatic, self.atom_map_num, self.query)
        atom.props = dict(self.props)
        return atom

    def GetIdx(self) -> int:
        return self._idx

    def GetAtomicNum(self) -> int:
        if self.symbol == "*":
            return 0
        return ATOMIC_NUMBERS.get(self.symbol, 0)

    def GetSymbol(self) -> str:
        return self.symbol

    def GetMass(self) -> float:
        return ATOMIC_MASSES.get(self.symbol, 0.0)

    def GetAtomMapNum(self) -> int:
        return self.atom_map_num

    def SetAtomMapNum(self, value: int) -> None:
        self.atom_map_num = int(value)

    def GetIsAromatic(self) -> bool:
        return self.aromatic

    def SetIsAromatic(self, value) -> None:
        self.aromatic = bool(value)

    def SetIntProp(self, key: str, value: int) -> None:
        self.props[key] = int(value)

    def GetIntProp(self, key: str) -> int:
        return int(self.props[key])

    def SetProp(self, key: str, value: object) -> None:
        self.props[key] = str(value)

    def GetProp(self, key: str) -> str:
        if key not in self.props:
            raise KeyError(key)
        return str(self.props[key])

    def GetPropsAsDict(self) -> dict[str, object]:
        return dict(self.props)

    def GetNeighbors(self) -> list[Atom]:
        if self._mol is None:
            return []
        neighbors = []
        for bond in self._mol.GetBonds():
            if bond.begin == self._idx:
                neighbors.append(self._mol.GetAtomWithIdx(bond.end))
            elif bond.end == self._idx:
                neighbors.append(self._mol.GetAtomWithIdx(bond.begin))
        return neighbors

    def GetBonds(self) -> list[Bond]:
        if self._mol is None:
            return []
        return [b for b in self._mol.GetBonds() if b.begin == self._idx or b.end == self._idx]


class Bond:
    def __init__(self, begin: int, end: int, bond_type: str = BondType.SINGLE):
        self.begin = int(begin)
        self.end = int(end)
        self.bond_type = bond_type
        self._mol: RWMol | None = None

    def copy(self) -> "Bond":
        return Bond(self.begin, self.end, self.bond_type)

    def GetBeginAtomIdx(self) -> int:
        return self.begin

    def GetEndAtomIdx(self) -> int:
        return self.end

    def GetBondType(self) -> str:
        return self.bond_type

    def SetBondType(self, value: str) -> None:
        self.bond_type = value


class RWMol:
    def __init__(self, mol: RWMol | None = None):
        self.atoms: list[Atom] = []
        self.bonds: list[Bond] = []
        self._conformers: dict[int, Conformer] = {}
        if mol is not None:
            mapping = {}
            for atom in mol.GetAtoms():
                mapping[atom.GetIdx()] = self.AddAtom(atom)
            for bond in mol.GetBonds():
                self.AddBond(mapping[bond.GetBeginAtomIdx()], mapping[bond.GetEndAtomIdx()], bond.GetBondType())

    def _refresh(self) -> None:
        for idx, atom in enumerate(self.atoms):
            atom._idx = idx
            atom._mol = self
        for bond in self.bonds:
            bond._mol = self

    def AddAtom(self, atom: Atom) -> int:
        new_atom = atom.copy() if isinstance(atom, Atom) else Atom(str(atom))
        new_atom._idx = len(self.atoms)
        new_atom._mol = self
        self.atoms.append(new_atom)
        return new_atom._idx

    def AddBond(self, begin: int, end: int, bond_type: str = BondType.SINGLE) -> None:
        if begin == end:
            return
        if self.GetBondBetweenAtoms(begin, end) is not None:
            self.GetBondBetweenAtoms(begin, end).SetBondType(bond_type)
            return
        bond = Bond(begin, end, bond_type)
        bond._mol = self
        self.bonds.append(bond)

    def RemoveBond(self, begin: int, end: int) -> None:
        self.bonds = [
            b for b in self.bonds
            if {b.begin, b.end} != {int(begin), int(end)}
        ]
        self._refresh()

    def RemoveAtom(self, idx: int) -> None:
        idx = int(idx)
        self.atoms.pop(idx)
        new_bonds = []
        for bond in self.bonds:
            if bond.begin == idx or bond.end == idx:
                continue
            begin = bond.begin - 1 if bond.begin > idx else bond.begin
            end = bond.end - 1 if bond.end > idx else bond.end
            new_bonds.append(Bond(begin, end, bond.bond_type))
        self.bonds = new_bonds
        self._refresh()

    def GetAtomWithIdx(self, idx: int) -> Atom:
        return self.atoms[int(idx)]

    def GetAtoms(self) -> list[Atom]:
        return list(self.atoms)

    def GetBonds(self) -> list[Bond]:
        return list(self.bonds)

    def GetNumAtoms(self) -> int:
        return len(self.atoms)

    def GetNumBonds(self) -> int:
        return len(self.bonds)

    def GetBondBetweenAtoms(self, begin: int, end: int) -> Bond | None:
        begin = int(begin)
        end = int(end)
        for bond in self.bonds:
            if (bond.begin == begin and bond.end == end) or (bond.begin == end and bond.end == begin):
                return bond
        return None

    def GetConformer(self, conf_id: int = 0) -> "Conformer":
        return self._conformers[conf_id]

    def AddConformer(self, conf: "Conformer", assignId: bool = True) -> int:
        conf_id = 0 if not self._conformers else max(self._conformers) + 1
        self._conformers[conf_id] = conf
        return conf_id


class QueryMol(RWMol):
    pass


@dataclass
class Position:
    x: float
    y: float
    z: float


class Conformer:
    def __init__(self, positions: dict[int, Position] | None = None):
        self.positions = positions or {}

    def GetAtomPosition(self, idx: int) -> Position:
        return self.positions[int(idx)]

    def SetAtomPosition(self, idx: int, pos: Position) -> None:
        self.positions[int(idx)] = pos


class Feature:
    def __init__(self, feature_type: str, atom_ids: Sequence[int]):
        self.feature_type = feature_type
        self.atom_ids = tuple(atom_ids)

    def GetType(self) -> str:
        return self.feature_type

    def GetAtomIds(self) -> tuple[int, ...]:
        return self.atom_ids


class FeatureFactory:
    def __init__(self, rules: Sequence[tuple[str, str]]):
        self.rules = []
        for feature_type, smarts in rules:
            try:
                self.rules.append((feature_type, MolFromSmarts(smarts)))
            except ValueError:
                continue

    def GetFeaturesForMol(self, molecule: RWMol) -> list[Feature]:
        features = []
        for feature_type, query in self.rules:
            for match in find_substructure_matches(query, molecule):
                if match:
                    ordered = tuple(match[i] for i in range(query.GetNumAtoms()))
                    features.append(Feature(feature_type, ordered))
        return features


class ChemicalFeatures:
    @staticmethod
    def BuildFeatureFactoryFromString(fdefn: str) -> FeatureFactory:
        rules = []
        current_type = None
        current_pattern = None
        for raw_line in fdefn.splitlines():
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("DefineFeature"):
                parts = line.split(None, 2)
                if len(parts) == 3:
                    current_type, current_pattern = parts[1], parts[2]
            elif line.startswith("EndFeature") and current_type and current_pattern:
                rules.append((current_type, current_pattern))
                current_type = None
                current_pattern = None
        return FeatureFactory(rules)


class AtomQuery:
    def __init__(self, condition: Condition | None = None):
        self.condition = condition or AnyCondition()

    def match(self, atom: Atom) -> bool:
        return self.condition.match(atom)


class Condition:
    def match(self, atom: Atom) -> bool:
        return True


class AnyCondition(Condition):
    pass


class AndCondition(Condition):
    def __init__(self, conditions: Iterable[Condition]):
        self.conditions = [c for c in conditions if c is not None]

    def match(self, atom: Atom) -> bool:
        return all(c.match(atom) for c in self.conditions)


class OrCondition(Condition):
    def __init__(self, conditions: Iterable[Condition]):
        self.conditions = [c for c in conditions if c is not None]

    def match(self, atom: Atom) -> bool:
        return any(c.match(atom) for c in self.conditions)


class NotCondition(Condition):
    def __init__(self, condition: Condition):
        self.condition = condition

    def match(self, atom: Atom) -> bool:
        return not self.condition.match(atom)


class ElementCondition(Condition):
    def __init__(self, symbol: str, aromatic: bool | None = None):
        self.symbol, inferred_aromatic = _normal_symbol(symbol)
        self.aromatic = inferred_aromatic if aromatic is None and inferred_aromatic else aromatic

    def match(self, atom: Atom) -> bool:
        if self.symbol != "*" and atom.GetSymbol() != self.symbol:
            return False
        if self.aromatic is not None and atom.GetIsAromatic() != self.aromatic:
            return False
        return True


class AtomicNumberCondition(Condition):
    def __init__(self, atomic_num: int):
        self.atomic_num = int(atomic_num)

    def match(self, atom: Atom) -> bool:
        return atom.GetAtomicNum() == self.atomic_num


class AromaticAnyCondition(Condition):
    def match(self, atom: Atom) -> bool:
        return atom.GetIsAromatic()


class AliphaticAnyCondition(Condition):
    def match(self, atom: Atom) -> bool:
        return not atom.GetIsAromatic()


class HCountCondition(Condition):
    def __init__(self, count: int):
        self.count = int(count)

    def match(self, atom: Atom) -> bool:
        return _total_h_count(atom) == self.count


class DegreeCondition(Condition):
    def __init__(self, degree: int):
        self.degree = int(degree)

    def match(self, atom: Atom) -> bool:
        return len(atom.GetNeighbors()) == self.degree


class RingCondition(Condition):
    def __init__(self, size: int | None = None):
        self.size = size

    def match(self, atom: Atom) -> bool:
        return _atom_in_ring(atom, self.size)


class ChargeCondition(Condition):
    def __init__(self, charge: int):
        self.charge = int(charge)

    def match(self, atom: Atom) -> bool:
        return int(atom.props.get("formal_charge", 0)) == self.charge


class RecursiveCondition(Condition):
    def __init__(self, smarts: str):
        self.smarts = smarts
        self.query = MolFromSmarts(smarts)

    def match(self, atom: Atom) -> bool:
        if atom._mol is None:
            return False
        return bool(find_substructure_matches(self.query, atom._mol, rooted=(0, atom.GetIdx())))


def _split_top_level(text: str, delimiter: str) -> list[str]:
    parts = []
    start = 0
    paren = 0
    bracket = 0
    i = 0
    while i < len(text):
        ch = text[i]
        if ch == "[":
            bracket += 1
        elif ch == "]":
            bracket -= 1
        elif ch == "(":
            paren += 1
        elif ch == ")":
            paren -= 1
        elif ch == delimiter and paren == 0 and bracket == 0:
            parts.append(text[start:i])
            start = i + 1
        i += 1
    parts.append(text[start:])
    return [p for p in parts if p != ""]


def _strip_outer_recursive(text: str) -> str | None:
    if not text.startswith("$(") or not text.endswith(")"):
        return None
    depth = 0
    for i, ch in enumerate(text):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0 and i != len(text) - 1:
                return None
    return text[2:-1]


def _read_number(text: str, start: int) -> tuple[int | None, int]:
    i = start
    while i < len(text) and text[i].isdigit():
        i += 1
    if i == start:
        return None, start
    return int(text[start:i]), i


def _parse_atom_part(part: str) -> Condition:
    part = part.strip()
    if not part:
        return AnyCondition()
    if part.startswith("!"):
        return NotCondition(_parse_atom_part(part[1:]))
    recursive = _strip_outer_recursive(part)
    if recursive is not None:
        return RecursiveCondition(recursive)

    conditions: list[Condition] = []
    i = 0
    element_seen = False

    if part == "*":
        return AnyCondition()
    if part == "a":
        return AromaticAnyCondition()
    if part == "A":
        return AliphaticAnyCondition()

    while i < len(part):
        ch = part[i]
        if ch == "*":
            conditions.append(AnyCondition())
            i += 1
        elif ch == "#":
            number, i2 = _read_number(part, i + 1)
            if number is not None:
                conditions.append(AtomicNumberCondition(number))
                element_seen = True
                i = i2
            else:
                i += 1
        elif part[i:i + 2] in ORGANIC_TWO_LETTER:
            conditions.append(ElementCondition(part[i:i + 2], aromatic=False))
            element_seen = True
            i += 2
        elif ch in AROMATIC_TO_ELEMENT:
            # Lower-case r means ring; all other lower-case organic symbols are aromatic atoms.
            if ch == "r":
                number, i2 = _read_number(part, i + 1)
                conditions.append(RingCondition(number))
                i = i2 if i2 != i + 1 else i + 1
            else:
                conditions.append(ElementCondition(ch, aromatic=True))
                element_seen = True
                i += 1
        elif ch in ORGANIC_ONE_LETTER:
            if ch == "H" and element_seen:
                number, i2 = _read_number(part, i + 1)
                conditions.append(HCountCondition(1 if number is None else number))
                i = i2 if number is not None else i + 1
            elif ch == "D":
                number, i2 = _read_number(part, i + 1)
                if number is not None:
                    conditions.append(DegreeCondition(number))
                    i = i2
                else:
                    i += 1
            elif ch == "R":
                number, i2 = _read_number(part, i + 1)
                conditions.append(RingCondition(None))
                i = i2 if number is not None else i + 1
            else:
                conditions.append(ElementCondition(ch, aromatic=False))
                element_seen = True
                i += 1
        elif ch == "H":
            number, i2 = _read_number(part, i + 1)
            conditions.append(HCountCondition(1 if number is None else number))
            i = i2 if number is not None else i + 1
        elif ch == "D":
            number, i2 = _read_number(part, i + 1)
            if number is not None:
                conditions.append(DegreeCondition(number))
                i = i2
            else:
                i += 1
        elif ch in "+-":
            sign = 1 if ch == "+" else -1
            number, i2 = _read_number(part, i + 1)
            conditions.append(ChargeCondition(sign * (1 if number is None else number)))
            i = i2 if number is not None else i + 1
        else:
            # Ignore unsupported SMARTS decorators such as @ stereochemistry.
            i += 1

    return AndCondition(conditions) if conditions else AnyCondition()


def _parse_atom_query(spec: str) -> AtomQuery:
    and_parts = _split_top_level(spec, ";")
    and_conditions = []
    for and_part in and_parts:
        or_parts = _split_top_level(and_part, ",")
        if len(or_parts) == 1:
            and_conditions.append(_parse_atom_part(or_parts[0]))
        else:
            and_conditions.append(OrCondition(_parse_atom_part(p) for p in or_parts))
    return AtomQuery(AndCondition(and_conditions))


def _extract_atom_map(spec: str) -> tuple[str, int]:
    depth = 0
    for i in range(len(spec) - 1, -1, -1):
        ch = spec[i]
        if ch == ")":
            depth += 1
        elif ch == "(":
            depth -= 1
        elif ch == ":" and depth == 0:
            map_text = spec[i + 1:]
            if map_text.isdigit():
                return spec[:i], int(map_text)
    return spec, 0


def _first_symbol_from_spec(spec: str) -> tuple[str, bool]:
    spec = spec.strip()
    if not spec:
        return "*", False
    if spec.startswith("#"):
        number, _ = _read_number(spec, 1)
        return ELEMENTS.get(number or 0, "*"), False
    if spec.startswith("*"):
        return "*", False
    if spec.startswith("a"):
        return "*", True
    if spec[:2] in ORGANIC_TWO_LETTER:
        return spec[:2], False
    if spec and spec[0] in AROMATIC_TO_ELEMENT:
        return AROMATIC_TO_ELEMENT[spec[0]], True
    if spec and spec[0] in ORGANIC_ONE_LETTER:
        return spec[0], False
    return "*", False


class SmilesParser:
    def __init__(self, text: str, query: bool = False):
        self.text = text.strip()
        self.query = query
        self.i = 0
        self.mol = QueryMol() if query else RWMol()
        self.current: int | None = None
        self.branch_stack: list[int | None] = []
        self.rings: dict[str, tuple[int, str | None]] = {}
        self.pending_bond: str | None = None

    def parse(self) -> RWMol:
        while self.i < len(self.text):
            ch = self.text[self.i]
            if ch in BOND_FROM_SYMBOL:
                self.pending_bond = BOND_FROM_SYMBOL[ch]
                self.i += 1
            elif ch == ".":
                self.current = None
                self.pending_bond = None
                self.i += 1
            elif ch == "(":
                self.branch_stack.append(self.current)
                self.i += 1
            elif ch == ")":
                if not self.branch_stack:
                    raise ValueError(f"Unmatched ')' in {self.text!r}")
                self.current = self.branch_stack.pop()
                self.i += 1
            elif ch.isdigit() or ch == "%":
                self._parse_ring()
            else:
                atom = self._parse_atom()
                new_idx = self.mol.AddAtom(atom)
                if self.current is not None:
                    bond_type = self.pending_bond or self._default_bond(self.current, new_idx)
                    self.mol.AddBond(self.current, new_idx, bond_type)
                self.current = new_idx
                self.pending_bond = None
        if self.branch_stack:
            raise ValueError(f"Unmatched '(' in {self.text!r}")
        if self.rings:
            raise ValueError(f"Unclosed ring in {self.text!r}: {sorted(self.rings)}")
        return self.mol

    def _default_bond(self, begin: int, end: int) -> str:
        a = self.mol.GetAtomWithIdx(begin)
        b = self.mol.GetAtomWithIdx(end)
        if a.GetIsAromatic() and b.GetIsAromatic():
            return BondType.AROMATIC
        return BondType.SINGLE

    def _parse_ring(self) -> None:
        if self.text[self.i] == "%":
            ring = self.text[self.i + 1:self.i + 3]
            self.i += 3
        else:
            ring = self.text[self.i]
            self.i += 1
        if self.current is None:
            raise ValueError(f"Ring closure without atom in {self.text!r}")
        if ring in self.rings:
            other, stored_bond = self.rings.pop(ring)
            bond_type = self.pending_bond or stored_bond or self._default_bond(self.current, other)
            self.mol.AddBond(self.current, other, bond_type)
        else:
            self.rings[ring] = (self.current, self.pending_bond)
        self.pending_bond = None

    def _parse_atom(self) -> Atom:
        ch = self.text[self.i]
        if ch == "[":
            end = self._find_closing_bracket(self.i)
            spec = self.text[self.i + 1:end]
            self.i = end + 1
            spec, atom_map = _extract_atom_map(spec)
            symbol, aromatic = _first_symbol_from_spec(spec)
            query = _parse_atom_query(spec) if self.query else None
            atom = Atom(symbol, aromatic=aromatic, atom_map_num=atom_map, query=query)
            charge_match = re.search(r"([+-])(\d*)", spec)
            if charge_match:
                sign = 1 if charge_match.group(1) == "+" else -1
                atom.props["formal_charge"] = sign * int(charge_match.group(2) or "1")
            return atom

        if self.text[self.i:self.i + 2] in ORGANIC_TWO_LETTER:
            symbol = self.text[self.i:self.i + 2]
            self.i += 2
            return Atom(symbol)
        if self.query and ch == "a":
            self.i += 1
            return Atom("*", aromatic=True, query=AtomQuery(AromaticAnyCondition()))
        if ch in AROMATIC_TO_ELEMENT:
            self.i += 1
            return Atom(ch)
        if ch in ORGANIC_ONE_LETTER:
            self.i += 1
            return Atom(ch)
        if ch == "*" and self.query:
            self.i += 1
            return Atom("*", query=AtomQuery(AnyCondition()))
        raise ValueError(f"Unsupported atom token {ch!r} in {self.text!r}")

    def _find_closing_bracket(self, start: int) -> int:
        i = start + 1
        depth = 0
        while i < len(self.text):
            ch = self.text[i]
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            elif ch == "]" and depth == 0:
                return i
            i += 1
        raise ValueError(f"Unclosed '[' in {self.text!r}")


def MolFromSmiles(smiles: str) -> RWMol:
    return SmilesParser(smiles, query=False).parse()


def MolFromSmarts(smarts: str) -> QueryMol:
    return SmilesParser(smarts, query=True).parse()


def _atom_matches(query_atom: Atom, target_atom: Atom) -> bool:
    if query_atom.query is not None:
        return query_atom.query.match(target_atom)
    if query_atom.GetSymbol() != "*" and query_atom.GetSymbol() != target_atom.GetSymbol():
        return False
    if query_atom.GetIsAromatic() != target_atom.GetIsAromatic():
        return False
    return True


def _bond_matches(query_bond: Bond, target_bond: Bond | None) -> bool:
    if target_bond is None:
        return False
    qtype = query_bond.GetBondType()
    ttype = target_bond.GetBondType()
    if qtype == BondType.ANY:
        return True
    return qtype == ttype


def find_substructure_matches(
    query: RWMol,
    target: RWMol,
    rooted: tuple[int, int] | None = None,
) -> list[dict[int, int]]:
    if query.GetNumAtoms() == 0:
        return []

    q_neighbors = {
        atom.GetIdx(): [n.GetIdx() for n in atom.GetNeighbors()]
        for atom in query.GetAtoms()
    }
    q_order = sorted(range(query.GetNumAtoms()), key=lambda i: len(q_neighbors[i]), reverse=True)
    if rooted is not None:
        q_root, target_root = rooted
        q_order = [q_root] + [i for i in q_order if i != q_root]
    matches: list[dict[int, int]] = []

    def backtrack(mapping: dict[int, int], used_targets: set) -> None:
        if len(mapping) == query.GetNumAtoms():
            matches.append(dict(mapping))
            return

        q_idx = next(i for i in q_order if i not in mapping)
        q_atom = query.GetAtomWithIdx(q_idx)
        if rooted is not None and q_idx == rooted[0]:
            candidates = [rooted[1]]
        else:
            candidates = range(target.GetNumAtoms())

        for t_idx in candidates:
            if t_idx in used_targets:
                continue
            t_atom = target.GetAtomWithIdx(t_idx)
            if not _atom_matches(q_atom, t_atom):
                continue
            ok = True
            for q_other, t_other in mapping.items():
                q_bond = query.GetBondBetweenAtoms(q_idx, q_other)
                if q_bond is None:
                    continue
                t_bond = target.GetBondBetweenAtoms(t_idx, t_other)
                if not _bond_matches(q_bond, t_bond):
                    ok = False
                    break
            if not ok:
                continue
            mapping[q_idx] = t_idx
            used_targets.add(t_idx)
            backtrack(mapping, used_targets)
            used_targets.remove(t_idx)
            del mapping[q_idx]

    backtrack({}, set())
    return matches


def _connected_components(molecule: RWMol) -> list[list[int]]:
    seen = set()
    comps = []
    for atom in molecule.GetAtoms():
        start = atom.GetIdx()
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        comp = []
        while stack:
            idx = stack.pop()
            comp.append(idx)
            for nbr in molecule.GetAtomWithIdx(idx).GetNeighbors():
                nidx = nbr.GetIdx()
                if nidx not in seen:
                    seen.add(nidx)
                    stack.append(nidx)
        comps.append(sorted(comp))
    return comps


def GetMolFrags(molecule: RWMol, asMols: bool = False):
    comps = _connected_components(molecule)
    if not asMols:
        return tuple(tuple(c) for c in comps)
    frags = []
    for comp in comps:
        frag = RWMol()
        mapping = {}
        for old_idx in comp:
            mapping[old_idx] = frag.AddAtom(molecule.GetAtomWithIdx(old_idx))
        for bond in molecule.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a in mapping and b in mapping:
                frag.AddBond(mapping[a], mapping[b], bond.GetBondType())
        frags.append(frag)
    return tuple(frags)


def _default_valence(atom: Atom) -> int:
    sym = atom.GetSymbol()
    if sym == "H":
        return 1
    if sym == "C":
        return 3 if atom.GetIsAromatic() else 4
    if sym == "N":
        return 2 if atom.GetIsAromatic() else 3
    if sym == "O":
        return 2
    if sym in {"F", "Cl", "Br", "I"}:
        return 1
    if sym == "P":
        return 3
    if sym == "S":
        return 2
    return 0


def _bond_valence(atom: Atom) -> float:
    return sum(_bond_order(bond.GetBondType(), aromatic_as_one=True) for bond in atom.GetBonds())


def _total_h_count(atom: Atom) -> int:
    explicit = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 1)
    if explicit:
        return explicit
    if atom.GetAtomicNum() == 1:
        return 0
    target = _default_valence(atom)
    if target <= 0:
        return 0
    return max(0, int(round(target - _bond_valence(atom))))


def AddHs(molecule: RWMol) -> RWMol:
    result = RWMol(molecule)
    heavy_atoms = [atom.GetIdx() for atom in result.GetAtoms() if atom.GetAtomicNum() != 1]
    for atom_idx in heavy_atoms:
        atom = result.GetAtomWithIdx(atom_idx)
        current_h = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 1)
        target_h = _total_h_count(atom)
        for _ in range(max(0, target_h - current_h)):
            h = Atom("H")
            for prop in ("molecule_id",):
                if prop in atom.props:
                    h.props[prop] = atom.props[prop]
            h_idx = result.AddAtom(h)
            result.AddBond(atom_idx, h_idx, BondType.SINGLE)
    return result


def RemoveAllHs(molecule: RWMol) -> RWMol:
    result = RWMol()
    mapping = {}
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        mapping[atom.GetIdx()] = result.AddAtom(atom)
    for bond in molecule.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in mapping and b in mapping:
            result.AddBond(mapping[a], mapping[b], bond.GetBondType())
    return result


def SanitizeMol(molecule: RWMol) -> int:
    for atom in molecule.GetAtoms():
        if atom.GetSymbol() in {"C", "N", "O", "S", "P"}:
            if any(b.GetBondType() == BondType.AROMATIC for b in atom.GetBonds()):
                atom.SetIsAromatic(True)
    return 0


def MolToSmiles(molecule: RWMol) -> str:
    if molecule.GetNumAtoms() == 0:
        return ""
    comps = _connected_components(molecule)
    pieces = []
    for comp in comps:
        seen = set()

        def atom_text(atom: Atom) -> str:
            if atom.GetAtomicNum() == 1:
                return "[H]"
            sym = atom.GetSymbol().lower() if atom.GetIsAromatic() and atom.GetSymbol() != "*" else atom.GetSymbol()
            return sym

        def dfs(idx: int, parent: int | None = None) -> str:
            seen.add(idx)
            atom = molecule.GetAtomWithIdx(idx)
            text = atom_text(atom)
            branches = []
            for nbr in atom.GetNeighbors():
                nidx = nbr.GetIdx()
                if nidx == parent or nidx in seen:
                    continue
                bond = molecule.GetBondBetweenAtoms(idx, nidx)
                branches.append(BOND_TO_SYMBOL.get(bond.GetBondType(), "") + dfs(nidx, idx))
            if branches:
                text += branches[0]
                for branch in branches[1:]:
                    text += f"({branch})"
            return text

        pieces.append(dfs(comp[0]))
    return ".".join(pieces)


def _atom_in_ring(atom: Atom, size: int | None = None) -> bool:
    mol = atom._mol
    if mol is None:
        return False
    start = atom.GetIdx()

    def dfs(current: int, parent: int, depth: int, seen: set) -> bool:
        for nbr in mol.GetAtomWithIdx(current).GetNeighbors():
            nidx = nbr.GetIdx()
            if nidx == parent:
                continue
            if nidx == start and depth >= 2:
                return size is None or depth + 1 == size
            if size is not None and depth + 1 >= size:
                continue
            if nidx in seen:
                continue
            if dfs(nidx, current, depth + 1, seen | {nidx}):
                return True
        return False

    return dfs(start, -1, 0, {start})


class MolOps:
    @staticmethod
    def GetMolFrags(molecule: RWMol, asMols: bool = False):
        return GetMolFrags(molecule, asMols=asMols)


@dataclass
class _Reaction:
    reactant_templates: list[QueryMol]
    product_templates: list[QueryMol]
    smarts: str

    def GetNumReactantTemplates(self) -> int:
        return len(self.reactant_templates)

    def GetReactantTemplate(self, idx: int) -> QueryMol:
        return self.reactant_templates[idx]

    def RunReactants(self, reactants: Sequence[RWMol]):
        if len(reactants) != len(self.reactant_templates):
            raise ValueError("Reaction reactant count does not match template count")

        template_matches = [
            find_substructure_matches(template, reactant)
            for template, reactant in zip(self.reactant_templates, reactants)
        ]
        if any(not matches for matches in template_matches):
            return tuple()

        map_to_reactant_atom: dict[int, tuple[int, int]] = {}
        for ridx, template in enumerate(self.reactant_templates):
            for atom in template.GetAtoms():
                if atom.GetAtomMapNum():
                    map_to_reactant_atom[atom.GetAtomMapNum()] = (ridx, atom.GetIdx())

        products = []
        for match_combo in itertools.product(*template_matches):
            product_mols = []
            for product_template in self.product_templates:
                product = RWMol()
                for p_atom in product_template.GetAtoms():
                    new_atom = Atom(p_atom.GetSymbol(), p_atom.GetIsAromatic(), p_atom.GetAtomMapNum(), p_atom.query)
                    if p_atom.GetAtomMapNum() in map_to_reactant_atom:
                        ridx, reactant_template_atom_idx = map_to_reactant_atom[p_atom.GetAtomMapNum()]
                        reactant_atom_idx = match_combo[ridx][reactant_template_atom_idx]
                        new_atom.props["old_mapno"] = p_atom.GetAtomMapNum()
                        new_atom.props["react_atom_idx"] = reactant_atom_idx
                        new_atom.props["reactant_idx"] = ridx
                    product.AddAtom(new_atom)
                for bond in product_template.GetBonds():
                    product.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType())
                product_mols.append(product)
            products.append(tuple(product_mols))
        return tuple(products)


class rdChemReactions:
    ChemicalReaction = _Reaction

    @staticmethod
    def ReactionFromSmarts(smarts: str) -> _Reaction:
        if ">>" not in smarts:
            raise ValueError(f"Reaction SMARTS must contain '>>': {smarts}")
        reactants_text, products_text = smarts.split(">>", 1)
        reactants = [MolFromSmarts(part) for part in _split_top_level(reactants_text, ".")]
        products = [MolFromSmarts(part) for part in _split_top_level(products_text, ".")]
        return _Reaction(reactants, products, smarts)


class AllChem:
    AddHs = staticmethod(AddHs)

    @staticmethod
    def EmbedMolecule(molecule: RWMol, useRandomCoords: bool = True) -> int:
        positions: dict[int, Position] = {}
        n = max(1, molecule.GetNumAtoms())
        golden = math.pi * (3.0 - math.sqrt(5.0))
        radius = max(1.0, n ** (1.0 / 3.0) * 0.7)
        for idx, atom in enumerate(molecule.GetAtoms()):
            y = 1.0 - (idx / max(1, n - 1)) * 2.0
            r = math.sqrt(max(0.0, 1.0 - y * y))
            theta = golden * idx
            scale = radius + 0.08 * (atom.GetAtomicNum() % 5)
            positions[idx] = Position(
                scale * math.cos(theta) * r,
                scale * y,
                scale * math.sin(theta) * r,
            )
        molecule._conformers[0] = Conformer(positions)
        return 0


class Chem:
    Atom = Atom
    Bond = Bond
    BondType = BondType
    RWMol = RWMol
    Mol = RWMol
    MolFromSmiles = staticmethod(MolFromSmiles)
    MolFromSmarts = staticmethod(MolFromSmarts)
    MolToSmiles = staticmethod(MolToSmiles)
    AddHs = staticmethod(AddHs)
    RemoveAllHs = staticmethod(RemoveAllHs)
    SanitizeMol = staticmethod(SanitizeMol)
    ChemicalFeatures = ChemicalFeatures
    rdmolops = MolOps
