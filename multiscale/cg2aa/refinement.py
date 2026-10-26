from __future__ import annotations

from collections import defaultdict
import math

import numpy as np

from .geometry import box_matrix, fractional_rotation, kabsch_fit, minimum_image


_COVALENT_RADII_ANGSTROM = {
    "H": 0.31,
    "B": 0.85,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Br": 1.20,
    "I": 1.39,
}


def _molecule_id(atom) -> int:
    try:
        return atom.GetIntProp("molecule_id")
    except KeyError:
        return -1


def _reaction_bonds(molecule, molecule_ids: np.ndarray) -> list[tuple[int, int]]:
    result = []
    for bond in molecule.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if molecule_ids[i] != molecule_ids[j]:
            result.append((i, j))
    return result


def _target_bond_length(molecule, i: int, j: int, angstrom_to_coordinate: float) -> float:
    left = molecule.GetAtomWithIdx(i).GetSymbol()
    right = molecule.GetAtomWithIdx(j).GetSymbol()
    value = _COVALENT_RADII_ANGSTROM.get(left, 0.75) + _COVALENT_RADII_ANGSTROM.get(right, 0.75)
    bond = molecule.GetBondBetweenAtoms(i, j)
    if bond is not None:
        bond_type = str(bond.GetBondType())
        if bond_type == "DOUBLE":
            value *= 0.90
        elif bond_type == "TRIPLE":
            value *= 0.84
        elif bond_type == "AROMATIC":
            value *= 0.93
    return value * float(angstrom_to_coordinate)


def _target_center_angle_degrees(molecule, center: int) -> float:
    atom = molecule.GetAtomWithIdx(center)
    heavy_neighbors = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() > 1]
    heavy_degree = len(heavy_neighbors)
    bond_types = {str(bond.GetBondType()) for bond in atom.GetBonds()}
    if atom.GetIsAromatic() or "DOUBLE" in bond_types or "AROMATIC" in bond_types:
        return 120.0
    if heavy_degree >= 4:
        return 109.47122063449069
    if heavy_degree == 3:
        return 120.0
    if heavy_degree == 2 and atom.GetSymbol() in {"C", "N", "O", "S"}:
        return 120.0
    return 109.47122063449069


def _constraint_bond_length(
    molecule,
    i: int,
    j: int,
    reference_coords: np.ndarray,
    box,
    molecule_ids: np.ndarray,
    angstrom_to_coordinate: float,
) -> float:
    if molecule_ids[i] != molecule_ids[j]:
        return _target_bond_length(molecule, i, j, angstrom_to_coordinate)
    return float(np.linalg.norm(minimum_image(reference_coords[j] - reference_coords[i], box)))


def _topology_exclusions(molecule) -> set[tuple[int, int]]:
    adjacency = [set() for _ in range(molecule.GetNumAtoms())]
    for bond in molecule.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        adjacency[i].add(j)
        adjacency[j].add(i)

    exclusions: set[tuple[int, int]] = set()
    for i, neighbors in enumerate(adjacency):
        for j in neighbors:
            exclusions.add(tuple(sorted((i, j))))
            for k in adjacency[j]:
                if k != i:
                    exclusions.add(tuple(sorted((i, k))))
    return exclusions


def _instance_components(instance_ids: list[int], reaction_bonds: list[tuple[int, int]], molecule_ids: np.ndarray) -> list[list[int]]:
    graph = {instance_id: set() for instance_id in instance_ids}
    for i, j in reaction_bonds:
        left = int(molecule_ids[i])
        right = int(molecule_ids[j])
        graph[left].add(right)
        graph[right].add(left)

    components = []
    seen = set()
    for instance_id in instance_ids:
        if instance_id in seen:
            continue
        stack = [instance_id]
        seen.add(instance_id)
        component = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in graph[current]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        components.append(component)
    return components


def _cell_clashes(
    molecule,
    coords: np.ndarray,
    box,
    molecule_ids: np.ndarray,
    exclusions: set[tuple[int, int]],
    cutoff: float,
    *,
    atom_scope: str = "heavy",
    interinstance_only: bool = True,
) -> list[tuple[int, int, np.ndarray, float]]:
    if cutoff <= 0.0:
        return []
    if atom_scope == "heavy":
        selected_indices = np.asarray(
            [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1],
            dtype=int,
        )
    elif atom_scope == "all":
        selected_indices = np.arange(molecule.GetNumAtoms(), dtype=int)
    else:
        raise ValueError("atom_scope must be 'heavy' or 'all'")
    if len(selected_indices) < 2:
        return []

    matrix = box_matrix(box)
    inverse = np.linalg.inv(matrix)
    lengths = np.asarray([np.linalg.norm(vector) for vector in matrix], dtype=float)
    cell_counts = np.maximum(1, np.floor(lengths / cutoff).astype(int))
    fractional = (coords[selected_indices] @ inverse) % 1.0

    cells: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for atom_idx, point in zip(selected_indices, fractional):
        key = tuple((np.floor(point * cell_counts).astype(int) % cell_counts).tolist())
        cells[key].append(int(atom_idx))

    result = []
    seen_pairs: set[tuple[int, int]] = set()
    for key, atoms in cells.items():
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neighbor_key = (
                        (key[0] + dx) % int(cell_counts[0]),
                        (key[1] + dy) % int(cell_counts[1]),
                        (key[2] + dz) % int(cell_counts[2]),
                    )
                    for i in atoms:
                        for j in cells.get(neighbor_key, []):
                            if i >= j:
                                continue
                            if interinstance_only and molecule_ids[i] == molecule_ids[j]:
                                continue
                            pair = (i, j)
                            if pair in seen_pairs or pair in exclusions:
                                continue
                            seen_pairs.add(pair)
                            delta = minimum_image(coords[j] - coords[i], box)
                            distance = float(np.linalg.norm(delta))
                            if distance < cutoff:
                                result.append((i, j, delta, distance))
    return result


def _reaction_angle_triplets(molecule, molecule_ids: np.ndarray) -> list[tuple[int, int, int]]:
    """Return heavy-atom angles that include at least one cross-instance bond."""

    heavy_adjacency = [set() for _ in range(molecule.GetNumAtoms())]
    cross_edges: set[tuple[int, int]] = set()
    for bond in molecule.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if molecule.GetAtomWithIdx(i).GetAtomicNum() <= 1 or molecule.GetAtomWithIdx(j).GetAtomicNum() <= 1:
            continue
        heavy_adjacency[i].add(j)
        heavy_adjacency[j].add(i)
        if molecule_ids[i] != molecule_ids[j]:
            cross_edges.add(tuple(sorted((i, j))))

    triplets = []
    for center, neighbors in enumerate(heavy_adjacency):
        ordered = sorted(neighbors)
        for left_index, left in enumerate(ordered):
            for right in ordered[left_index + 1 :]:
                if (
                    tuple(sorted((left, center))) not in cross_edges
                    and tuple(sorted((center, right))) not in cross_edges
                ):
                    continue
                triplets.append((left, center, right))
    return triplets


def _angle_degrees(left: np.ndarray, center: np.ndarray, right: np.ndarray, box) -> float:
    first = minimum_image(left - center, box)
    second = minimum_image(right - center, box)
    first_norm = float(np.linalg.norm(first))
    second_norm = float(np.linalg.norm(second))
    if first_norm < 1.0e-12 or second_norm < 1.0e-12:
        return 0.0
    cosine = float(np.dot(first, second) / (first_norm * second_norm))
    return math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0))))


def _angle_summary(values: list[float], minimum: float, maximum: float) -> dict:
    if not values:
        return {
            "count": 0,
            "min_deg": None,
            "median_deg": None,
            "max_deg": None,
            "below_minimum": 0,
            "above_maximum": 0,
            "minimum_allowed_deg": float(minimum),
            "maximum_allowed_deg": float(maximum),
        }
    array = np.asarray(values, dtype=float)
    return {
        "count": int(len(array)),
        "min_deg": float(np.min(array)),
        "median_deg": float(np.median(array)),
        "max_deg": float(np.max(array)),
        "below_minimum": int(np.count_nonzero(array < float(minimum))),
        "above_maximum": int(np.count_nonzero(array > float(maximum))),
        "minimum_allowed_deg": float(minimum),
        "maximum_allowed_deg": float(maximum),
    }


def _length_summary(values: list[float], coordinate_to_nm: float) -> dict:
    if not values:
        return {"count": 0, "min_nm": None, "median_nm": None, "max_nm": None}
    array = np.asarray(values, dtype=float) * float(coordinate_to_nm)
    return {
        "count": int(len(array)),
        "min_nm": float(np.min(array)),
        "median_nm": float(np.median(array)),
        "max_nm": float(np.max(array)),
    }


def geometry_audit(
    molecule,
    coords,
    box,
    *,
    coordinate_to_nm: float,
    severe_clash_angstrom: float = 1.00,
    severe_all_atom_clash_angstrom: float = 0.80,
    minimum_reaction_angle_deg: float = 70.0,
    maximum_reaction_angle_deg: float = 170.0,
    angstrom_to_coordinate: float = 1.0,
    check_hydrogens: bool = True,
) -> dict:
    coordinates = np.asarray(coords, dtype=float)
    molecule_ids = np.asarray([_molecule_id(atom) for atom in molecule.GetAtoms()], dtype=int)
    reaction_bonds = _reaction_bonds(molecule, molecule_ids)
    reaction_lengths = [
        float(np.linalg.norm(minimum_image(coordinates[j] - coordinates[i], box)))
        for i, j in reaction_bonds
    ]
    reaction_angle_triplets = _reaction_angle_triplets(molecule, molecule_ids)
    reaction_angles = [
        _angle_degrees(coordinates[left], coordinates[center], coordinates[right], box)
        for left, center, right in reaction_angle_triplets
    ]

    hydrogen_lengths = []
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() != 1 or not atom.GetNeighbors():
            continue
        parent = atom.GetNeighbors()[0]
        hydrogen_lengths.append(float(np.linalg.norm(coordinates[atom.GetIdx()] - coordinates[parent.GetIdx()])))

    exclusions = _topology_exclusions(molecule)
    severe_cutoff = float(severe_clash_angstrom) * float(angstrom_to_coordinate)
    heavy_clashes = _cell_clashes(
        molecule,
        coordinates,
        box,
        molecule_ids,
        exclusions,
        severe_cutoff,
        atom_scope="heavy",
        interinstance_only=True,
    )
    all_atom_cutoff = float(severe_all_atom_clash_angstrom) * float(angstrom_to_coordinate)
    all_atom_clashes = _cell_clashes(
        molecule,
        coordinates,
        box,
        molecule_ids,
        exclusions,
        all_atom_cutoff,
        atom_scope="all",
        interinstance_only=False,
    )

    matrix = box_matrix(box)
    fractional = coordinates @ np.linalg.inv(matrix)
    outside = np.any((fractional < -0.5) | (fractional >= 0.5), axis=1)

    reaction_summary = _length_summary(reaction_lengths, coordinate_to_nm)
    reaction_angle_summary = _angle_summary(
        reaction_angles,
        minimum_reaction_angle_deg,
        maximum_reaction_angle_deg,
    )
    hydrogen_summary = _length_summary(hydrogen_lengths, coordinate_to_nm)
    reaction_ok = (
        reaction_summary["count"] == 0
        or (
            reaction_summary["min_nm"] is not None
            and reaction_summary["min_nm"] >= 0.11
            and reaction_summary["max_nm"] <= 0.20
        )
    )
    hydrogen_ok = (
        hydrogen_summary["count"] == 0
        or (
            hydrogen_summary["min_nm"] is not None
            and hydrogen_summary["min_nm"] >= 0.085
            and hydrogen_summary["max_nm"] <= 0.14
        )
    )
    reaction_angles_ok = (
        reaction_angle_summary["count"] == 0
        or (
            reaction_angle_summary["below_minimum"] == 0
            and reaction_angle_summary["above_maximum"] == 0
        )
    )
    checks = {
        "reaction_bond_lengths_reasonable": bool(reaction_ok),
        "reaction_center_angles_reasonable": bool(reaction_angles_ok),
        "no_severe_interinstance_heavy_atom_clashes": not heavy_clashes,
    }
    if check_hydrogens:
        checks.update(
            {
                "hydrogen_bond_lengths_reasonable": bool(hydrogen_ok),
                "no_severe_nonbonded_all_atom_clashes": not all_atom_clashes,
            }
        )
    return {
        "audit_scope": "full_all_atom" if check_hydrogens else "reactive_heavy_atom_stage",
        "reaction_bonds": reaction_summary,
        "reaction_center_angles": reaction_angle_summary,
        "hydrogen_bonds": hydrogen_summary,
        "severe_interinstance_heavy_atom_clashes": {
            "threshold_nm": severe_cutoff * float(coordinate_to_nm),
            "count": len(heavy_clashes),
            "minimum_nm": (
                min(item[3] for item in heavy_clashes) * float(coordinate_to_nm)
                if heavy_clashes
                else None
            ),
        },
        "severe_nonbonded_all_atom_clashes": {
            "threshold_nm": all_atom_cutoff * float(coordinate_to_nm),
            "count": len(all_atom_clashes),
            "minimum_nm": (
                min(item[3] for item in all_atom_clashes) * float(coordinate_to_nm)
                if all_atom_clashes
                else None
            ),
        },
        "positions_outside_primary_box": int(np.count_nonzero(outside)),
        "checks": checks,
        "ok": all(checks.values()),
    }


def _distance_geometry_relaxation(
    molecule,
    coords: np.ndarray,
    reference_coords: np.ndarray,
    box,
    molecule_ids: np.ndarray,
    *,
    angstrom_to_coordinate: float,
    iterations: int = 120,
) -> np.ndarray:
    heavy = {atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1}
    adjacency = [set() for _ in range(molecule.GetNumAtoms())]
    constraints: dict[tuple[int, int], tuple[float, float]] = {}

    for bond in molecule.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        adjacency[i].add(j)
        adjacency[j].add(i)
        if i not in heavy or j not in heavy:
            continue
        pair = tuple(sorted((i, j)))
        if molecule_ids[i] != molecule_ids[j]:
            target = _target_bond_length(molecule, i, j, angstrom_to_coordinate)
        else:
            target = float(np.linalg.norm(minimum_image(reference_coords[j] - reference_coords[i], box)))
        constraints[pair] = (target, 1.60)

    for i in heavy:
        for neighbor in adjacency[i]:
            for k in adjacency[neighbor]:
                if k == i or k not in heavy or molecule_ids[i] != molecule_ids[k]:
                    continue
                pair = tuple(sorted((i, k)))
                if pair in constraints:
                    continue
                target = float(np.linalg.norm(minimum_image(reference_coords[pair[1]] - reference_coords[pair[0]], box)))
                constraints[pair] = (target, 0.55)

    # Cross-instance 1-3 distances enforce chemically plausible reaction-center
    # angles without requiring a force-field-specific angle table.
    for center in heavy:
        neighbors = sorted(index for index in adjacency[center] if index in heavy)
        if len(neighbors) < 2:
            continue
        target_angle = math.radians(_target_center_angle_degrees(molecule, center))
        for left_index, left in enumerate(neighbors):
            for right in neighbors[left_index + 1 :]:
                if (
                    molecule_ids[center] == molecule_ids[left]
                    and molecule_ids[center] == molecule_ids[right]
                ):
                    continue
                pair = tuple(sorted((left, right)))
                if pair in constraints:
                    continue
                left_length = _constraint_bond_length(
                    molecule,
                    center,
                    left,
                    reference_coords,
                    box,
                    molecule_ids,
                    angstrom_to_coordinate,
                )
                right_length = _constraint_bond_length(
                    molecule,
                    center,
                    right,
                    reference_coords,
                    box,
                    molecule_ids,
                    angstrom_to_coordinate,
                )
                target = math.sqrt(
                    max(
                        0.0,
                        left_length * left_length
                        + right_length * right_length
                        - 2.0 * left_length * right_length * math.cos(target_angle),
                    )
                )
                constraints[pair] = (target, 0.80)

    if constraints:
        constraint_i = np.asarray([pair[0] for pair in constraints], dtype=int)
        constraint_j = np.asarray([pair[1] for pair in constraints], dtype=int)
        targets = np.asarray([value[0] for value in constraints.values()], dtype=float)
        weights = np.asarray([value[1] for value in constraints.values()], dtype=float)
    else:
        constraint_i = np.empty(0, dtype=int)
        constraint_j = np.empty(0, dtype=int)
        targets = np.empty(0, dtype=float)
        weights = np.empty(0, dtype=float)

    exclusions = _topology_exclusions(molecule)
    clash_cutoff = 1.10 * float(angstrom_to_coordinate)
    maximum_step = 0.04 * float(angstrom_to_coordinate)
    coordinates = np.asarray(coords, dtype=float).copy()

    for _ in range(int(iterations)):
        displacement = np.zeros_like(coordinates)
        counts = np.zeros(len(coordinates), dtype=float)

        if len(constraint_i):
            delta = minimum_image(coordinates[constraint_j] - coordinates[constraint_i], box)
            distances = np.linalg.norm(delta, axis=1)
            valid = distances > 1.0e-12
            correction = np.zeros_like(delta)
            correction[valid] = (
                0.5
                * weights[valid, None]
                * (distances[valid] - targets[valid])[:, None]
                * delta[valid]
                / distances[valid, None]
            )
            np.add.at(displacement, constraint_i, correction)
            np.add.at(displacement, constraint_j, -correction)
            np.add.at(counts, constraint_i, weights)
            np.add.at(counts, constraint_j, weights)

        clashes = _cell_clashes(
            molecule,
            coordinates,
            box,
            molecule_ids,
            exclusions,
            clash_cutoff,
        )
        for i, j, delta, distance in clashes:
            if distance < 1.0e-12:
                direction = np.asarray([1.0, 0.0, 0.0], dtype=float)
            else:
                direction = delta / distance
            weight = 1.5
            correction = 0.5 * weight * (distance - clash_cutoff) * direction
            displacement[i] += correction
            displacement[j] -= correction
            counts[i] += weight
            counts[j] += weight

        active = counts > 0.0
        step = np.zeros_like(coordinates)
        step[active] = displacement[active] / counts[active, None]
        magnitudes = np.linalg.norm(step, axis=1)
        large = magnitudes > maximum_step
        step[large] *= (maximum_step / magnitudes[large])[:, None]
        coordinates += 0.65 * step

    return coordinates


def refine_reactive_geometry(
    molecule,
    coords,
    box,
    *,
    angstrom_to_coordinate: float = 1.0,
    coordinate_to_nm: float = 0.1,
    iterations: int = 120,
) -> tuple[list[np.ndarray], dict]:
    """Relax reaction geometry while preserving template bonds and local angles."""

    coordinates = np.asarray(coords, dtype=float).copy()
    reference_coordinates = coordinates.copy()
    molecule_ids = np.asarray([_molecule_id(atom) for atom in molecule.GetAtoms()], dtype=int)
    instance_ids = sorted(set(int(value) for value in molecule_ids))
    groups = {instance_id: np.where(molecule_ids == instance_id)[0] for instance_id in instance_ids}
    reaction_bonds = _reaction_bonds(molecule, molecule_ids)
    if not reaction_bonds:
        report = geometry_audit(
            molecule,
            coordinates,
            box,
            coordinate_to_nm=coordinate_to_nm,
            angstrom_to_coordinate=angstrom_to_coordinate,
            check_hydrogens=False,
        )
        return [np.asarray(point, dtype=float) for point in coordinates], {
            "applied": False,
            "rigid_body_iterations": 0,
            "distance_geometry_iterations": 0,
            "before": report,
            "after": report,
        }

    before = geometry_audit(
        molecule,
        coordinates,
        box,
        coordinate_to_nm=coordinate_to_nm,
        angstrom_to_coordinate=angstrom_to_coordinate,
        check_hydrogens=False,
    )
    exclusions = _topology_exclusions(molecule)
    by_instance: dict[int, list[tuple[int, int]]] = {instance_id: [] for instance_id in instance_ids}
    for i, j in reaction_bonds:
        by_instance[int(molecule_ids[i])].append((i, j))
        by_instance[int(molecule_ids[j])].append((j, i))

    components = _instance_components(instance_ids, reaction_bonds, molecule_ids)
    component_indices = {
        tuple(component): np.concatenate([groups[instance_id] for instance_id in component])
        for component in components
    }
    component_centers = {
        component: coordinates[indices].mean(axis=0)
        for component, indices in component_indices.items()
    }

    fit_fraction = 0.25
    clash_fraction = 0.30
    clash_cutoff = 1.20 * float(angstrom_to_coordinate)
    maximum_shift = 0.15 * float(angstrom_to_coordinate)

    for _ in range(int(iterations)):
        updated = coordinates.copy()
        for instance_id in instance_ids:
            links = by_instance[instance_id]
            if not links:
                continue
            source = []
            target = []
            for own, other in links:
                delta = minimum_image(coordinates[own] - coordinates[other], box)
                distance = float(np.linalg.norm(delta))
                if distance < 1.0e-12:
                    axis = np.asarray([1.0, 0.0, 0.0], dtype=float)
                    if (own + other) % 3 == 1:
                        axis = np.asarray([0.0, 1.0, 0.0], dtype=float)
                    elif (own + other) % 3 == 2:
                        axis = np.asarray([0.0, 0.0, 1.0], dtype=float)
                else:
                    axis = delta / distance
                nearby_other = coordinates[own] - delta
                source.append(coordinates[own])
                target.append(
                    nearby_other
                    + axis * _target_bond_length(molecule, own, other, angstrom_to_coordinate)
                )

            source_array = np.asarray(source, dtype=float)
            target_array = np.asarray(target, dtype=float)
            rotation, _ = kabsch_fit(source_array, target_array)
            partial_rotation = fractional_rotation(rotation, fit_fraction)
            source_center = source_array.mean(axis=0)
            target_center = target_array.mean(axis=0)
            indices = groups[instance_id]
            updated[indices] = (
                (coordinates[indices] - source_center) @ partial_rotation
                + source_center
                + fit_fraction * (target_center - source_center)
            )
        coordinates = updated

        clashes = _cell_clashes(
            molecule,
            coordinates,
            box,
            molecule_ids,
            exclusions,
            clash_cutoff,
        )
        shifts = {instance_id: np.zeros(3, dtype=float) for instance_id in instance_ids}
        counts = {instance_id: 0 for instance_id in instance_ids}
        for i, j, delta, distance in clashes:
            if distance < 1.0e-12:
                direction = np.asarray([1.0, 0.0, 0.0], dtype=float)
            else:
                direction = delta / distance
            correction = clash_fraction * (clash_cutoff - distance) * direction
            left = int(molecule_ids[i])
            right = int(molecule_ids[j])
            shifts[left] -= 0.5 * correction
            shifts[right] += 0.5 * correction
            counts[left] += 1
            counts[right] += 1

        for instance_id in instance_ids:
            if counts[instance_id] == 0:
                continue
            shift = shifts[instance_id] / max(1.0, math.sqrt(counts[instance_id]))
            magnitude = float(np.linalg.norm(shift))
            if magnitude > maximum_shift:
                shift *= maximum_shift / magnitude
            coordinates[groups[instance_id]] += shift

        for component, indices in component_indices.items():
            coordinates[indices] += component_centers[component] - coordinates[indices].mean(axis=0)

    distance_iterations = 160
    coordinates = _distance_geometry_relaxation(
        molecule,
        coordinates,
        reference_coordinates,
        box,
        molecule_ids,
        angstrom_to_coordinate=angstrom_to_coordinate,
        iterations=distance_iterations,
    )

    after = geometry_audit(
        molecule,
        coordinates,
        box,
        coordinate_to_nm=coordinate_to_nm,
        angstrom_to_coordinate=angstrom_to_coordinate,
        check_hydrogens=False,
    )
    return [np.asarray(point, dtype=float) for point in coordinates], {
        "applied": True,
        "rigid_body_iterations": int(iterations),
        "distance_geometry_iterations": distance_iterations,
        "before": before,
        "after": after,
    }
