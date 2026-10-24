"""Internal OPLS force-field assignment for Multiscale.

This module reads bundled Multiscale package assets from
``multiscale/forcefield/assets/opls`` and uses the shared ``multiscale.cg2aa.chem`` chemistry
layer. It is the default backend for ``multiscale cg2aa run``. The feature-definition
conversion and atom-typing workflow have historical lineage from
``asteeves/rdlt``; provenance and redistribution notes are recorded in
``NOTICE.md``.
"""

import os
import re
import warnings
from collections import Counter
from typing import Any

from multiscale.cg2aa.chem import Chem
from .models import ForceFieldAssignment


this_dir, this_filename = os.path.split(__file__)
_asset_dir = os.path.join(this_dir, "assets", "opls")
opls_tpl_path = os.path.join(_asset_dir, "STaGE_opls_tomoltemplate_opls.txt")
opls_bonded_itp = os.path.join(_asset_dir, "bonded.itp")


def generate_feature_defn(fpath):
    """Build RDKit feature definitions from a Moltemplate conversion file."""
    info = {}
    fdefn = ''
    with open(fpath, 'r') as infile:
        feat_index = 0
        for line in [line.strip() for line in infile if line.strip()]:
            if line[0] != '*':
                el, atomname, typename, patt, lttype, chg, desc = [el.strip() for el in line.split("|")]
                info[lttype] = (atomname, typename, chg, desc)
                # write lttype, SMARTS to feature definintion file
                # NOTE: using feature family to store the atom names is dangerous
                # because rdkit won't assign mutliple features in same family.
                # So I had to assign an index to make unique names [AHS]
                fdefn += \
                    """
                    DefineFeature {0} {1}
                    Family {2}{3}
                    EndFeature""".format(lttype, patt, feat_index, atomname)
                # add new charge dictionary entry
                feat_index += 1
    return info, fdefn


def _bond_order(bond):
    bond_type = bond.GetBondType()
    if bond_type == Chem.BondType.DOUBLE:
        return 2
    if bond_type == Chem.BondType.TRIPLE:
        return 3
    return 1


def _explicit_h_count(atom):
    return sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1)


def _heavy_neighbors(atom):
    return [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]


def _has_bond_order(atom, order):
    return any(_bond_order(bond) == order for bond in atom.GetBonds())


def _neighbors_by_bond_order(atom, order):
    return [
        nbr for nbr in atom.GetNeighbors()
        if _bond_order_to(atom, nbr) == order
    ]


def _bond_order_to(atom, nbr):
    bond = atom._mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
    return _bond_order(bond) if bond is not None else 1


def _known(info, atom_type):
    return atom_type if atom_type in info else None


def _has_aromatic_neighbor(atom, symbol=None):
    for nbr in atom.GetNeighbors():
        if not nbr.GetIsAromatic():
            continue
        if symbol is None or nbr.GetSymbol() == symbol:
            return True
    return False


def _has_neighbor(atom, symbol=None, bond_order=None, aromatic=None):
    for nbr in atom.GetNeighbors():
        if symbol is not None and nbr.GetSymbol() != symbol:
            continue
        if aromatic is not None and nbr.GetIsAromatic() != aromatic:
            continue
        if bond_order is not None and _bond_order_to(atom, nbr) != bond_order:
            continue
        return True
    return False


def _is_phenol_oxygen(atom):
    return (
        atom.GetSymbol() == "O"
        and _explicit_h_count(atom) > 0
        and _has_aromatic_neighbor(atom, "C")
    )


def _is_triazine_carbon(atom):
    return (
        atom.GetSymbol() == "C"
        and atom.GetIsAromatic()
        and _has_aromatic_neighbor(atom, "N")
    )


def _is_triazine_hydrogen(atom):
    if atom.GetSymbol() != "H" or not atom.GetNeighbors():
        return False
    parent = atom.GetNeighbors()[0]
    return _is_triazine_carbon(parent)


def _is_cyanide_carbon(atom):
    return atom.GetSymbol() == "C" and _has_neighbor(atom, "N", bond_order=3)


def _cyanide_single_neighbor(atom):
    for nbr in atom.GetNeighbors():
        if nbr.GetSymbol() == "N" and _bond_order_to(atom, nbr) == 3:
            continue
        return nbr
    return None


def _neutralize_total_charge(molecule, total_charge):
    atoms = molecule.GetAtoms()
    if not atoms:
        return 0.0
    correction = -total_charge / len(atoms)
    charges = [
        round(float(atom.GetProp("Charge")) + correction, 6)
        for atom in atoms
    ]
    residual = round(-sum(charges), 6)
    if residual:
        idx = min(range(len(charges)), key=lambda i: abs(charges[i]))
        charges[idx] = round(charges[idx] + residual, 6)
    for atom, charge in zip(atoms, charges):
        atom.SetProp("Charge", f"{charge:.6f}")
    return sum(charges)


def _quick_opls_type(atom, info, defaults):
    """Assign common OPLS types for supported polymer environments."""
    symbol = atom.GetSymbol()
    neighbors = _heavy_neighbors(atom)
    h_count = _explicit_h_count(atom)

    if symbol == "H":
        parent = atom.GetNeighbors()[0] if atom.GetNeighbors() else None
        if _is_triazine_hydrogen(atom):
            return _known(info, "@atom:584")
        if parent is not None and parent.GetIsAromatic():
            return _known(info, "@atom:91")
        if parent is not None and _is_phenol_oxygen(parent):
            return _known(info, "@atom:110")
        if parent is not None and parent.GetSymbol() == "O":
            return _known(info, "@atom:97")
        if parent is not None and parent.GetSymbol() == "N":
            return _known(info, "@atom:739")
        if parent is not None and _has_bond_order(parent, 3):
            return _known(info, "@atom:756")
        return _known(info, "@atom:85")

    if symbol == "C":
        if atom.GetIsAromatic():
            if _is_triazine_carbon(atom):
                return _known(info, "@atom:583")
            if any(_is_phenol_oxygen(n) for n in atom.GetNeighbors()):
                return _known(info, "@atom:108")
            if _has_neighbor(atom, "O"):
                return _known(info, "@atom:141")
            if any(_is_cyanide_carbon(n) for n in atom.GetNeighbors()):
                return _known(info, "@atom:202")
            return _known(info, "@atom:90")
        if _is_cyanide_carbon(atom):
            single_nbr = _cyanide_single_neighbor(atom)
            if single_nbr is not None and single_nbr.GetIsAromatic():
                return _known(info, "@atom:203")
            return _known(info, "@atom:695")
        if _has_bond_order(atom, 3):
            if h_count:
                return _known(info, "@atom:755")
            return _known(info, "@atom:759") or _known(info, "@atom:769")
        if any(n.GetSymbol() == "O" and _bond_order_to(atom, n) == 2 for n in neighbors):
            if any(n.GetSymbol() == "N" for n in neighbors):
                return _known(info, "@atom:177")
            if h_count:
                return _known(info, "@atom:219")
            return _known(info, "@atom:222")
        if h_count >= 3:
            if any(n.GetIsAromatic() for n in neighbors):
                return _known(info, "@atom:93")
            return _known(info, "@atom:80")
        if h_count == 2:
            if any(n.GetIsAromatic() for n in neighbors):
                return _known(info, "@atom:94")
            return _known(info, "@atom:81")
        if h_count == 1:
            return _known(info, "@atom:82")
        return _known(info, "@atom:84")

    if symbol == "N":
        if atom.GetIsAromatic():
            return _known(info, "@atom:582")
        if _has_bond_order(atom, 3):
            triple_carbons = _neighbors_by_bond_order(atom, 3)
            if triple_carbons:
                single_nbr = _cyanide_single_neighbor(triple_carbons[0])
                if single_nbr is not None and single_nbr.GetIsAromatic():
                    return _known(info, "@atom:204")
            return _known(info, "@atom:694")
        if h_count >= 2:
            return _known(info, "@atom:730")
        if h_count == 1:
            return _known(info, "@atom:731")
        return _known(info, "@atom:732") or (defaults or {}).get(symbol)

    if symbol == "O":
        if _has_bond_order(atom, 2):
            double_bond_carbons = [
                nbr for nbr in neighbors
                if nbr.GetSymbol() == "C" and _bond_order_to(atom, nbr) == 2
            ]
            amide_like = any(
                any(nb.GetSymbol() == "N" for nb in c.GetNeighbors())
                for c in double_bond_carbons
            )
            if amide_like:
                return _known(info, "@atom:178") or _known(info, "@atom:220")
            return _known(info, "@atom:220") or _known(info, "@atom:178")
        if h_count:
            if _is_phenol_oxygen(atom):
                return _known(info, "@atom:109")
            if any(n.GetIsAromatic() for n in neighbors):
                return _known(info, "@atom:109")
            return _known(info, "@atom:96")
        if _has_aromatic_neighbor(atom, "C") and _has_neighbor(atom, "C"):
            return _known(info, "@atom:121")
        return _known(info, "@atom:122") or (defaults or {}).get(symbol)

    return (defaults or {}).get(symbol)


def ff(molecule, **kwargs):
    info, fdefn = generate_feature_defn(opls_tpl_path)

    defaults = kwargs.get('defaults', {})
    use_smarts = kwargs.get("use_smarts", False)
    neutralize = kwargs.get("neutralize", True)

    for atom in molecule.GetAtoms():
        atom_type = _quick_opls_type(atom, info, defaults)
        if atom_type is not None:
            atom.SetProp("AtomType", atom_type)

    if use_smarts:
        factory = Chem.ChemicalFeatures.BuildFeatureFactoryFromString(fdefn)
        features = factory.GetFeaturesForMol(molecule)
        [molecule.GetAtomWithIdx(f.GetAtomIds()[0]).SetProp('AtomType', f.GetType()) for f in features]

    for atom in molecule.GetAtoms():
        try:
            _ = atom.GetProp("AtomType")
        except KeyError:
            _m = Chem.RWMol()
            i_a = count = 0
            _m.AddAtom(atom)
            count += 1
            for nbr_atom in atom.GetNeighbors():
                j = count
                _m.AddAtom(nbr_atom)
                count += 1
                _m.AddBond(i_a, j, molecule.GetBondBetweenAtoms(atom.GetIdx(), nbr_atom.GetIdx()).GetBondType())
                for _nbr_atom in nbr_atom.GetNeighbors():
                    if _nbr_atom.GetIdx() != atom.GetIdx():
                        k = count
                        _m.AddAtom(_nbr_atom)
                        count += 1
                        _m.AddBond(j, k,
                                   molecule.GetBondBetweenAtoms(_nbr_atom.GetIdx(), nbr_atom.GetIdx()).GetBondType())
            smiles = Chem.MolToSmiles(_m)
            msg = "---------------------------------------------------------\n" \
                  "Chemical env of atom {0}({1}) is\n\t{2}\nType is not found from database, set manually as\n\t{3}"
            s = "*** NO DEFAULT TYPE ***\n"
            if defaults is not None:
                if defaults.get(atom.GetSymbol()) is not None:
                    atom.SetProp("AtomType", defaults.get(atom.GetSymbol()))
                    dt = defaults.get(atom.GetSymbol())
                    s = info[dt][0] + ", " + info[dt][3]
            print(msg.format(atom.GetIdx(), atom.GetSymbol(), smiles, s))


    sum_of_charge = 0
    for atom in molecule.GetAtoms():
        names = info[atom.GetProp("AtomType")]
        atom.SetProp("ElementType", names[0])
        atom.SetProp("OplsType", names[1])
        atom.SetProp("Charge", names[2])
        sum_of_charge += float(names[2])
    if neutralize and abs(sum_of_charge) > 1e-6:
        sum_of_charge = _neutralize_total_charge(molecule, sum_of_charge)
    if abs(sum_of_charge) > 0.01:
        warnings.warn("*** Total Charge is %.6f !" % sum_of_charge)

    bonded_itp = open(opls_bonded_itp, 'r').read()
    bonds = []
    angles = []
    dihedrals = []
    for bond in molecule.GetBonds():  # bonds
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        atom_i = molecule.GetAtomWithIdx(i)
        atom_j = molecule.GetAtomWithIdx(j)

        bonds.append('%s-%s %d %d' % (atom_i.GetProp("ElementType"),
                                      atom_j.GetProp("ElementType"),
                                      i, j))

        for nbr_i in atom_i.GetNeighbors():  # a bond (i, j) and neighbor of i and j define a dihedral
            ni = nbr_i.GetIdx()
            if ni == j:
                continue
            for nbr_j in atom_j.GetNeighbors():
                nj = nbr_j.GetIdx()
                if nj == i or ni == nj:  # avoid 3-ring
                    continue
                dihedral = [nbr_i.GetProp("ElementType"),
                            atom_i.GetProp("ElementType"),
                            atom_j.GetProp("ElementType"),
                            nbr_j.GetProp("ElementType")]
                forward = re.search(r'\s+'.join(dihedral) + r'\s+([1-9])\s', bonded_itp)
                reverse = re.search(r'\s+'.join(dihedral[::-1]) + r'\s+([1-9])\s', bonded_itp)
                match = forward or reverse
                if match:
                    func = match.group(1)
                    dihedrals.append('%s_F%s %d %d %d %d' % ('-'.join(dihedral), func, ni, i, j, nj))

    for atom in molecule.GetAtoms():  # any atom can be a center of an angle
        for nbr_i in atom.GetNeighbors():
            for nbr_j in atom.GetNeighbors():
                if nbr_i.GetIdx() <= nbr_j.GetIdx():
                    continue
                angle = [nbr_i.GetProp("ElementType"), atom.GetProp("ElementType"), nbr_j.GetProp("ElementType")]
                if re.search(r'\s*'.join(angle) + r'\s*[1-9]+', bonded_itp) or re.search(
                        r'\s*'.join(angle[::-1]) + r'\s*[1-9]+', bonded_itp):
                    angles.append('%s %d %d %d' % ('-'.join(angle), nbr_i.GetIdx(), atom.GetIdx(), nbr_j.GetIdx()))

    return molecule, bonds, angles, dihedrals


def _atom_prop(atom: Any, key: str, default: str = "") -> str:
    try:
        return str(atom.GetProp(key))
    except Exception:
        return default


def _summarize_assignment(molecule: Any, bonds: list[str], angles: list[str], dihedrals: list[str], warning_messages: list[str]) -> ForceFieldAssignment:
    atoms = list(molecule.GetAtoms())
    atom_type_counts = Counter(_atom_prop(atom, "AtomType", "") for atom in atoms)
    opls_type_counts = Counter(_atom_prop(atom, "OplsType", "") for atom in atoms)
    total_charge = 0.0
    for atom in atoms:
        try:
            total_charge += float(atom.GetProp("Charge"))
        except Exception:
            pass
    return ForceFieldAssignment(
        molecule,
        list(bonds),
        list(angles),
        list(dihedrals),
        len(atoms),
        dict(atom_type_counts),
        dict(opls_type_counts),
        round(total_charge, 6),
        warning_messages,
    )


def assign_forcefield(molecule: Any, name: str = "opls", **kwargs: Any) -> ForceFieldAssignment:
    """Assign bundled OPLS parameters and return an Multiscale summary object."""
    if name != "opls":
        raise ValueError(f"unsupported force field: {name}")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        typed_mol, bonds, angles, dihedrals = ff(molecule, **kwargs)
    return _summarize_assignment(typed_mol, bonds, angles, dihedrals, [str(w.message) for w in caught])


def internal_module_status() -> dict[str, Any]:
    """Report whether the bundled OPLS implementation/assets are usable."""
    template = os.path.abspath(opls_tpl_path)
    bonded = os.path.abspath(opls_bonded_itp)
    this_file = os.path.abspath(__file__)
    package_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    return {
        "ok": os.path.isfile(template) and os.path.isfile(bonded),
        "module": __name__,
        "file": this_file,
        "source": "internal-multiscale",
        "assets": {
            "template_path": template,
            "bonded_path": bonded,
        },
        "checks": {
            "template_exists": os.path.isfile(template),
            "bonded_exists": os.path.isfile(bonded),
            "module_inside_package": this_file.startswith(package_root),
        },
        "has_ff": callable(ff),
        "has_assign_forcefield": callable(assign_forcefield),
    }
