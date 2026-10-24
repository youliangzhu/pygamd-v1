from __future__ import annotations

from itertools import combinations

import networkx as nx
import numpy as np

from .chem import Chem, rdChemReactions
from .mapping import atom_map, bond_map, process_reactants


class ReactionChoices:
    def __init__(self, reaction_id: str, reactants_tuple: tuple[str, ...], reaction_meta: dict):
        self.reaction_id = reaction_id
        self.reactants_tuple = reactants_tuple
        self.reaction_meta = reaction_meta
        self.index = list(range(len(reaction_meta["reactions"])))
        probabilities = [entry[1] for entry in reaction_meta["reactions"]]
        self.probabilities = np.asarray(probabilities, dtype=float) / sum(probabilities)
        self.reaction_maps = {idx: [] for idx in self.index}
        self.product_idx = {idx: reaction_meta["reactions"][idx][-1] for idx in self.index}
        self.all_reactions = set(self.index)

    def process(self, reactants_meta: dict[str, dict]) -> None:
        for idx in self.index:
            reaction = self.reaction_meta["reactions"][idx][0]
            molecules = [Chem.MolFromSmiles(reactants_meta[t]["smiles"]) for t in self.reactants_tuple]
            reactants = process_reactants(molecules)
            products = reaction.RunReactants(reactants)
            for product in products:
                amap, reacting_atoms = atom_map(product, reaction)
                bmap = bond_map(reactants, product, reaction)
                self.reaction_maps[idx].append((reacting_atoms, amap, bmap))

    def choose(self, idx: int | None = None):
        if idx is None:
            idx = int(np.random.choice(self.index, p=self.probabilities))
        return idx, self.reaction_maps[idx], self.product_idx[idx]


def _site_allowed(reaction_map, site_constraints: dict[int, set[int]] | None) -> bool:
    if not site_constraints:
        return True
    reacting_atoms = reaction_map[0]
    for reactant_id, allowed_atoms in site_constraints.items():
        atoms = set(reacting_atoms.get(reactant_id, []))
        if atoms and not atoms.issubset(allowed_atoms):
            return False
    return True


def _allowed_reaction(reacted_atoms, reactions: ReactionChoices, site_constraints: dict[int, set[int]] | None = None):
    seen = set()
    while seen != reactions.all_reactions:
        idx, reaction_maps, product_idx = reactions.choose()
        if idx in seen:
            continue
        seen.add(idx)
        for reaction_map in reaction_maps:
            if not _site_allowed(reaction_map, site_constraints):
                continue
            allowed = True
            for reactant_id in reaction_map[0]:
                if set(reaction_map[0][reactant_id]).issubset(reacted_atoms[reactant_id]):
                    allowed = False
                    break
            if allowed:
                return reaction_map, product_idx
    return None, None


def _product_ids(value):
    if value in (None, "all"):
        return None
    return [int(v) for v in value]


def build_reaction_templates(config: dict) -> dict:
    templates = {}
    for reaction in config.get("reactions", []):
        reaction_id = reaction["id"]
        reactants = reaction.get("reactants_tuple", reaction.get("reactants"))
        if not reactants:
            raise ValueError(f"Reaction {reaction_id} must define reactants.")
        entries = reaction.get("templates") or [reaction]
        compiled = []
        for entry in entries:
            compiled.append(
                (
                    rdChemReactions.ReactionFromSmarts(entry["smarts"]),
                    float(entry.get("probability", reaction.get("probability", 1.0))),
                    _product_ids(entry.get("keep_products", reaction.get("keep_products"))),
                )
            )
        templates[reaction_id] = {
            "reactants_tuple": tuple(reactants),
            "reactions": compiled,
        }
    return templates


def build_instance_graph(cg, instances: list[dict], bead_lookup: dict[int, dict]) -> nx.Graph:
    graph = nx.Graph()
    for instance in instances:
        graph.add_node(
            instance["id"],
            type=instance["type"],
            template=instance["template"],
        )
    for bond_type, bead_i, bead_j in cg.bonds:
        if bead_i not in bead_lookup or bead_j not in bead_lookup:
            continue
        left = bead_lookup[bead_i]
        right = bead_lookup[bead_j]
        if left["instance"] == right["instance"]:
            continue
        i = int(left["instance"])
        j = int(right["instance"])
        if graph.has_edge(i, j):
            graph.edges[i, j]["bond_types"].add(str(bond_type))
            graph.edges[i, j]["site_pairs"].append((left["site"], right["site"]))
            graph.edges[i, j]["site_links"].append((i, left["site"], j, right["site"], str(bond_type)))
        else:
            graph.add_edge(
                i,
                j,
                bond_types={str(bond_type)},
                site_pairs=[(left["site"], right["site"])],
                site_links=[(i, left["site"], j, right["site"], str(bond_type))],
            )
    return graph


def _coerce_node(value, bead_lookup: dict[int, dict], index_space: str) -> int:
    value = int(value)
    if index_space == "instance":
        return value
    if value not in bead_lookup:
        raise ValueError(f"CG bead {value} is not assigned to a template instance")
    return int(bead_lookup[value]["instance"])


def _dedupe_events(events):
    seen = set()
    unique = []
    for event in events:
        key = (event[0], tuple(event[1:]))
        if key not in seen:
            seen.add(key)
            unique.append(event)
    return unique


def _events_from_xml_section(rule, cg, bead_lookup):
    section = rule["section"]
    reaction_id = rule["reaction"]
    columns = rule.get("columns")
    if columns is None:
        raise ValueError(f"xml_section rule for {section} must define columns.")
    index_space = rule.get("index_space", "bead")
    events = []
    for row in cg.rows(section):
        nodes = [_coerce_node(row[int(col)], bead_lookup, index_space) for col in columns]
        events.append(tuple([reaction_id] + nodes))
    return events


def _types_match(actual, expected):
    return tuple(actual) == tuple(expected)


def _events_from_bonds(rule, graph):
    reaction_id = rule["reaction"]
    expected = rule.get("reactants")
    allowed_bonds = set(map(str, rule.get("bond_types", [])))
    allow_reverse = rule.get("allow_reverse", True)
    events = []
    for i, j, data in graph.edges(data=True):
        for left_node, left_site, right_node, right_site, bond_type in data.get("site_links", []):
            if allowed_bonds and bond_type not in allowed_bonds:
                continue
            nodes = (left_node, right_node)
            sites = (left_site, right_site)
            if expected is None:
                events.append((reaction_id, nodes[0], nodes[1], sites[0], sites[1]))
                continue
            actual = (graph.nodes[nodes[0]]["type"], graph.nodes[nodes[1]]["type"])
            if _types_match(actual, expected):
                events.append((reaction_id, nodes[0], nodes[1], sites[0], sites[1]))
            elif allow_reverse and _types_match(actual[::-1], expected):
                events.append((reaction_id, nodes[1], nodes[0], sites[1], sites[0]))
    return events


def _role_nodes(center, neighbors, roles):
    nodes = []
    next_neighbor = 0
    for role in roles:
        if role == "center":
            nodes.append(center)
        elif role == "neighbor":
            nodes.append(neighbors[next_neighbor])
            next_neighbor += 1
        else:
            raise ValueError(f"Unsupported node_neighbors role: {role}")
    return nodes


def _events_from_node_neighbors(rule, graph):
    center_type = rule["center_type"]
    neighbor_types = set(rule.get("neighbor_types", []))
    events = []
    for center, attrs in graph.nodes(data=True):
        if attrs["type"] != center_type:
            continue
        neighbors = [
            node for node in sorted(graph.neighbors(center))
            if not neighbor_types or graph.nodes[node]["type"] in neighbor_types
        ]
        for case in rule.get("cases", []):
            count = int(case["neighbor_count"])
            if case.get("exact", True) and len(neighbors) != count:
                continue
            if len(neighbors) < count:
                continue
            for selected in combinations(neighbors, count):
                events.append(tuple([case["reaction"]] + _role_nodes(center, list(selected), case["roles"])))
    return events


def build_reaction_events(config: dict, cg, graph, bead_lookup):
    events = []
    for rule in config.get("reaction_events", []):
        source = rule["source"]
        if source == "xml_section":
            events.extend(_events_from_xml_section(rule, cg, bead_lookup))
        elif source == "bonds":
            events.extend(_events_from_bonds(rule, graph))
        elif source == "node_neighbors":
            events.extend(_events_from_node_neighbors(rule, graph))
        else:
            raise ValueError(f"Unsupported reaction event source: {source}")
    return _dedupe_events(events)


def _make_reaction_choices(reaction_templates, reactants_meta):
    choices = {}
    for reaction_id, meta in reaction_templates.items():
        choices[reaction_id] = ReactionChoices(reaction_id, meta["reactants_tuple"], meta)
        choices[reaction_id].process(reactants_meta)
    return choices


def apply_reactions(molecule, coords: list, instances: list[dict], reaction_templates: dict, events: list[tuple], reactants_meta: dict):
    choices = _make_reaction_choices(reaction_templates, reactants_meta)
    meta = {}
    for instance in instances:
        meta[instance["id"]] = {
            "atom_idx": dict(instance["local_to_global"]),
            "site_chem_atoms": {
                site: set(atoms)
                for site, atoms in instance.get("site_chem_atoms", {}).items()
            },
            "reacting_map": {},
            "rm_atoms": set(),
        }

    for event in events:
        reaction_id = event[0]
        reaction_choices = choices.get(reaction_id)
        if reaction_choices is None:
            raise ValueError(f"Reaction event {event} references undefined reaction {reaction_id}")
        reactant_count = len(reaction_choices.reactants_tuple)
        reactant_ids = tuple(int(v) for v in event[1:1 + reactant_count])
        site_names = tuple(event[1 + reactant_count:])
        if site_names and len(site_names) != reactant_count:
            raise ValueError(f"Reaction event {event} has {len(site_names)} sites for {reactant_count} reactants")

        reactants = [meta[idx] for idx in reactant_ids]
        site_constraints = None
        if site_names:
            site_constraints = {}
            for reactant_index, site in enumerate(site_names):
                allowed_atoms = reactants[reactant_index]["site_chem_atoms"].get(site)
                if allowed_atoms is None:
                    raise ValueError(f"Reaction event {event} references undefined site {site!r}")
                site_constraints[reactant_index] = allowed_atoms
        reacted_atoms = {idx: set() for idx in range(len(reactants))}
        key = tuple(sorted(reactant_ids))
        for reactant_index, reactant in enumerate(reactants):
            for stored_key in reactant["reacting_map"]:
                reacted_atoms[reactant_index].update(reactant["reacting_map"][stored_key])

        reaction_map, product_idx = _allowed_reaction(reacted_atoms, reaction_choices, site_constraints)
        if not reaction_map:
            raise ValueError(f"{reactant_ids} can not react via {reaction_id}")

        reacting_atoms, amap, bmap = reaction_map
        for reactant_index, reactant in enumerate(reactants):
            reactant["reacting_map"].setdefault(key, set())
            for atom_id in reacting_atoms.get(reactant_index, []):
                reactant["reacting_map"][key].add(atom_id)

        for atom in amap:
            if product_idx is None or atom.product_id in product_idx:
                continue
            reactant = reactants[atom.reactant_id]
            if atom.reactant_atom_id in reactant["atom_idx"]:
                reactant["rm_atoms"].add(reactant["atom_idx"][atom.reactant_atom_id])

        for bond in bmap:
            if bond.status == "deleted":
                reactant = reactants[bond.reactants_id[0]]
                if bond.reactant_atoms_id[0] not in reactant["atom_idx"] or bond.reactant_atoms_id[1] not in reactant["atom_idx"]:
                    continue
                i = reactant["atom_idx"][bond.reactant_atoms_id[0]]
                j = reactant["atom_idx"][bond.reactant_atoms_id[1]]
                molecule.RemoveBond(i, j)
            elif bond.status == "changed":
                reactant = reactants[bond.reactants_id[0]]
                if bond.reactant_atoms_id[0] not in reactant["atom_idx"] or bond.reactant_atoms_id[1] not in reactant["atom_idx"]:
                    continue
                i = reactant["atom_idx"][bond.reactant_atoms_id[0]]
                j = reactant["atom_idx"][bond.reactant_atoms_id[1]]
                target_bond = molecule.GetBondBetweenAtoms(i, j)
                if target_bond is not None:
                    target_bond.SetBondType(bond.bond_type)
            elif bond.status == "new":
                reactant0 = reactants[bond.reactants_id[0]]
                reactant1 = reactants[bond.reactants_id[1]]
                if bond.reactant_atoms_id[0] not in reactant0["atom_idx"] or bond.reactant_atoms_id[1] not in reactant1["atom_idx"]:
                    continue
                i = reactant0["atom_idx"][bond.reactant_atoms_id[0]]
                j = reactant1["atom_idx"][bond.reactant_atoms_id[1]]
                molecule.AddBond(i, j, bond.bond_type)

    remove_atoms = sorted({atom for data in meta.values() for atom in data["rm_atoms"]}, reverse=True)
    for atom_idx in remove_atoms:
        molecule.RemoveAtom(atom_idx)
        coords.pop(atom_idx)
    return molecule, coords
