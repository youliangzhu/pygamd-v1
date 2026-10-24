from __future__ import annotations

import math

import numpy as np

from .chem import Chem
from .geometry import minimum_image


def connected_components(molecule) -> list[list[int]]:
    seen = set()
    components = []
    for atom in molecule.GetAtoms():
        start = atom.GetIdx()
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        component = []
        while stack:
            idx = stack.pop()
            component.append(idx)
            for neighbor in molecule.GetAtomWithIdx(idx).GetNeighbors():
                nidx = neighbor.GetIdx()
                if nidx not in seen:
                    seen.add(nidx)
                    stack.append(nidx)
        components.append(sorted(component))
    return components


def split_molecule_with_coords(molecule, coords: list[np.ndarray]) -> list[tuple[object, list[np.ndarray]]]:
    fragments = []
    for component in connected_components(molecule):
        fragment = Chem.RWMol()
        mapping = {}
        fragment_coords = []
        for old_idx in component:
            new_idx = fragment.AddAtom(molecule.GetAtomWithIdx(old_idx))
            mapping[old_idx] = new_idx
            fragment_coords.append(np.asarray(coords[old_idx], dtype=float))
        for bond in molecule.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            if i in mapping and j in mapping:
                fragment.AddBond(mapping[i], mapping[j], bond.GetBondType())
        fragments.append((fragment, fragment_coords))
    return fragments


def _unit(vector: np.ndarray) -> np.ndarray | None:
    norm = float(np.linalg.norm(vector))
    if norm < 1.0e-12:
        return None
    return np.asarray(vector, dtype=float) / norm


def _candidate_directions(parent_idx: int) -> list[np.ndarray]:
    golden = math.pi * (3.0 - math.sqrt(5.0))
    directions = []
    candidate_count = 96
    for index in range(candidate_count):
        z = 1.0 - 2.0 * ((index + 0.5) / candidate_count)
        radius = math.sqrt(max(0.0, 1.0 - z * z))
        theta = golden * (index + parent_idx * 0.38196601125)
        directions.append(np.asarray([math.cos(theta) * radius, math.sin(theta) * radius, z], dtype=float))
    return directions


def _hydrogen_candidates(parent_idx: int, neighbor_vectors: list[np.ndarray]) -> list[np.ndarray]:
    """Return deterministic candidate directions, with local-geometry choices first."""

    occupied = [direction for vector in neighbor_vectors if (direction := _unit(vector)) is not None]
    candidates = _candidate_directions(parent_idx)
    preferred = []
    if occupied:
        opposite = _unit(-np.sum(occupied, axis=0))
        if opposite is not None:
            preferred.append(opposite)
        if len(occupied) >= 2:
            normal = _unit(np.cross(occupied[0], occupied[1]))
            if normal is not None:
                preferred.extend([normal, -normal])

    unique = []
    for candidate in preferred + candidates:
        direction = _unit(candidate)
        if direction is None:
            continue
        if any(float(np.dot(direction, existing)) > 0.999999 for existing in unique):
            continue
        unique.append(direction)
    return unique


def _topology_exclusions(molecule) -> list[set[int]]:
    adjacency = [set() for _ in range(molecule.GetNumAtoms())]
    for bond in molecule.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        adjacency[i].add(j)
        adjacency[j].add(i)

    exclusions = []
    for atom_idx, neighbors in enumerate(adjacency):
        excluded = {atom_idx, *neighbors}
        for neighbor in neighbors:
            excluded.update(adjacency[neighbor])
        exclusions.append(excluded)
    return exclusions


def _place_hydrogens_with_clearance(
    molecule,
    coords: list[np.ndarray],
    hydrogens_by_parent: dict[int, list[int]],
    *,
    box=None,
    angstrom_to_coordinate: float,
    passes: int = 1,
) -> None:
    """Place H atoms at fixed bond lengths while avoiding nonbonded overlaps."""

    exclusions = _topology_exclusions(molecule)
    heavy_indices = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() != 1]
    hydrogen_indices = sorted(
        atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 1
    )
    placed = set(heavy_indices)
    severe_cutoff = 0.80 * float(angstrom_to_coordinate)
    soft_cutoff = 1.20 * float(angstrom_to_coordinate)

    for pass_index in range(max(1, int(passes))):
        if pass_index:
            placed = set(heavy_indices)
        for parent_idx in sorted(hydrogens_by_parent):
            parent = molecule.GetAtomWithIdx(parent_idx)
            hydrogen_group = sorted(hydrogens_by_parent[parent_idx])
            bond_length = _hydrogen_bond_length(parent.GetSymbol()) * float(angstrom_to_coordinate)
            heavy_neighbor_vectors = [
                np.asarray(coords[neighbor.GetIdx()], dtype=float) - np.asarray(coords[parent_idx], dtype=float)
                for neighbor in parent.GetNeighbors()
                if neighbor.GetAtomicNum() != 1
            ]
            selected_directions: list[np.ndarray] = []
            for hydrogen_idx in hydrogen_group:
                local_vectors = heavy_neighbor_vectors + selected_directions
                candidates = _hydrogen_candidates(parent_idx + hydrogen_idx, local_vectors)
                local_directions = [
                    direction
                    for vector in local_vectors
                    if (direction := _unit(vector)) is not None
                ]
                comparison_indices = sorted(placed - exclusions[hydrogen_idx])
                candidate_array = np.asarray(candidates, dtype=float)
                candidate_positions = (
                    np.asarray(coords[parent_idx], dtype=float)[None, :]
                    + bond_length * candidate_array
                )
                if comparison_indices:
                    comparison_coords = np.asarray(
                        [coords[index] for index in comparison_indices],
                        dtype=float,
                    )
                    deltas = comparison_coords[None, :, :] - candidate_positions[:, None, :]
                    if box is not None:
                        deltas = minimum_image(deltas, box)
                    distances = np.linalg.norm(deltas, axis=2)
                    severe_count = np.count_nonzero(distances < severe_cutoff, axis=1)
                    soft_count = np.count_nonzero(distances < soft_cutoff, axis=1)
                    minimum_distance = np.min(distances, axis=1)
                    soft_penalty = np.sum(np.maximum(0.0, soft_cutoff - distances) ** 2, axis=1)
                else:
                    size = len(candidate_array)
                    severe_count = np.zeros(size, dtype=int)
                    soft_count = np.zeros(size, dtype=int)
                    minimum_distance = np.full(size, np.inf)
                    soft_penalty = np.zeros(size, dtype=float)

                if local_directions:
                    dots = candidate_array @ np.asarray(local_directions, dtype=float).T
                    local_overlap = np.count_nonzero(dots > 0.45, axis=1)
                    maximum_dot = np.max(dots, axis=1)
                    angular_penalty = np.sum((dots + 1.0 / 3.0) ** 2, axis=1)
                else:
                    size = len(candidate_array)
                    local_overlap = np.zeros(size, dtype=int)
                    maximum_dot = np.full(size, -1.0)
                    angular_penalty = np.zeros(size, dtype=float)

                order = np.lexsort(
                    (
                        -minimum_distance,
                        angular_penalty,
                        maximum_dot,
                        soft_penalty,
                        soft_count,
                        local_overlap,
                        severe_count,
                    )
                )
                best_index = int(order[0])
                coords[hydrogen_idx] = candidate_positions[best_index]
                selected_directions.append(
                    np.asarray(coords[hydrogen_idx], dtype=float) - np.asarray(coords[parent_idx], dtype=float)
                )
                placed.add(hydrogen_idx)

        # On later passes each H is reconsidered against all H atoms placed earlier in the pass.
        if pass_index + 1 < max(1, int(passes)):
            for hydrogen_idx in hydrogen_indices:
                placed.discard(hydrogen_idx)


def _hydrogen_bond_length(parent_symbol: str) -> float:
    return {
        "C": 1.09,
        "N": 1.01,
        "O": 0.96,
        "S": 1.34,
    }.get(parent_symbol, 1.00)


def rebuild_hydrogens(
    molecule,
    coords: list[np.ndarray],
    *,
    angstrom_to_coordinate: float = 1.0,
    box=None,
) -> tuple[object, list[np.ndarray]]:
    """Rebuild hydrogens with chemically reasonable bond lengths and directions."""

    heavy_coords = []
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() != 1:
            heavy_coords.append(np.asarray(coords[atom.GetIdx()], dtype=float))

    heavy = Chem.RemoveAllHs(molecule)
    Chem.SanitizeMol(heavy)
    rebuilt = Chem.AddHs(heavy)
    new_coords = [np.asarray(coord, dtype=float) for coord in heavy_coords]

    hydrogens_by_parent: dict[int, list[int]] = {}
    for atom in rebuilt.GetAtoms():
        if atom.GetAtomicNum() != 1:
            continue
        parent = atom.GetNeighbors()[0] if atom.GetNeighbors() else None
        if parent is None:
            new_coords.append(np.zeros(3))
            continue
        parent_idx = parent.GetIdx()
        hydrogens_by_parent.setdefault(parent_idx, []).append(atom.GetIdx())
        try:
            atom.SetIntProp("molecule_id", parent.GetIntProp("molecule_id"))
        except KeyError:
            pass
        new_coords.append(np.zeros(3))

    _place_hydrogens_with_clearance(
        rebuilt,
        new_coords,
        hydrogens_by_parent,
        box=box,
        angstrom_to_coordinate=angstrom_to_coordinate,
    )
    return rebuilt, new_coords


def set_missing_h_molecule_ids(molecule) -> object:
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() != 1 or not atom.GetNeighbors():
            continue
        try:
            atom.GetIntProp("molecule_id")
        except KeyError:
            parent = atom.GetNeighbors()[0]
            try:
                atom.SetIntProp("molecule_id", parent.GetIntProp("molecule_id"))
            except KeyError:
                atom.SetIntProp("molecule_id", 0)
    return molecule
