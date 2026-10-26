from __future__ import annotations

from collections import namedtuple

from .chem import rdChemReactions


atom_info = namedtuple(
    "atom_info",
    ("reactant_id", "reactant_atom_id", "product_id", "product_atom_id", "atomic_number"),
)
bond_info = namedtuple(
    "bond_info",
    ("product_id", "product_atoms_id", "reactants_id", "reactant_atoms_id", "bond_type", "status"),
)


def process_reactants(reactants):
    for r_idx, molecule in enumerate(reactants):
        for atom in molecule.GetAtoms():
            atom.SetIntProp("reactant_idx", r_idx)
    return reactants


def map_reacting_atoms(reaction):
    reacting_map = {}
    for r_idx in range(reaction.GetNumReactantTemplates()):
        template = reaction.GetReactantTemplate(r_idx)
        for atom in template.GetAtoms():
            if atom.GetAtomMapNum():
                reacting_map[atom.GetAtomMapNum()] = r_idx
    return reacting_map


def map_atoms(products, reacting_map):
    atom_map = []
    reacting_atoms = {}
    for product_id, product in enumerate(products):
        for atom in product.GetAtoms():
            old_mapno = atom.GetPropsAsDict().get("old_mapno")
            react_atom_idx = atom.GetPropsAsDict().get("react_atom_idx")
            if old_mapno is None:
                continue
            reactant_id = reacting_map[old_mapno]
            reacting_atoms.setdefault(reactant_id, []).append(react_atom_idx)
            atom_map.append(
                atom_info(
                    reactant_id,
                    react_atom_idx,
                    product_id,
                    atom.GetIdx(),
                    atom.GetAtomicNum(),
                )
            )
    return atom_map, reacting_atoms


def atom_map(products, reaction):
    reacting_map = map_reacting_atoms(reaction)
    return map_atoms(products, reacting_map)


def bond_map(reactants: list, products: list, reaction: rdChemReactions.ChemicalReaction) -> list:
    atom_mapping, _ = atom_map(products, reaction)
    result = []
    for reactant_id, reactant in enumerate(reactants):
        for bond in reactant.GetBonds():
            a_pid = b_pid = None
            a = bond.GetBeginAtomIdx()
            b = bond.GetEndAtomIdx()
            bond_type = bond.GetBondType()
            for atom in atom_mapping:
                if atom.reactant_id != reactant_id:
                    continue
                if atom.reactant_atom_id == a:
                    a_pid = atom.product_id
                    a_paid = atom.product_atom_id
                if atom.reactant_atom_id == b:
                    b_pid = atom.product_id
                    b_paid = atom.product_atom_id
            if a_pid is None or b_pid is None:
                continue
            if a_pid != b_pid:
                result.append(bond_info(a_pid, (a_paid, b_paid), (reactant_id, reactant_id), (a, b), bond_type, "deleted"))
                continue
            product_bond = products[a_pid].GetBondBetweenAtoms(a_paid, b_paid)
            if product_bond is None:
                result.append(bond_info(a_pid, (a_paid, b_paid), (reactant_id, reactant_id), (a, b), bond_type, "deleted"))

    for product_id, product in enumerate(products):
        for bond in product.GetBonds():
            a_rid = b_rid = None
            a = bond.GetBeginAtomIdx()
            b = bond.GetEndAtomIdx()
            bond_type = bond.GetBondType()
            for atom in atom_mapping:
                if atom.product_id != product_id:
                    continue
                if atom.product_atom_id == a:
                    a_rid = atom.reactant_id
                    a_raid = atom.reactant_atom_id
                if atom.product_atom_id == b:
                    b_rid = atom.reactant_id
                    b_raid = atom.reactant_atom_id
            if a_rid is None or b_rid is None:
                continue
            if a_rid != b_rid:
                result.append(bond_info(product_id, (a, b), (a_rid, b_rid), (a_raid, b_raid), bond_type, "new"))
                continue
            reactant_bond = reactants[a_rid].GetBondBetweenAtoms(a_raid, b_raid)
            if reactant_bond is None:
                result.append(bond_info(product_id, (a, b), (a_rid, b_rid), (a_raid, b_raid), bond_type, "new"))
            elif reactant_bond.GetBondType() != bond_type:
                result.append(bond_info(product_id, (a, b), (a_rid, b_rid), (a_raid, b_raid), bond_type, "changed"))
    return result
