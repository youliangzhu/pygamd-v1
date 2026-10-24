from __future__ import annotations

import numpy as np

from multiscale.cg2aa.chem import Chem
from multiscale.cg2aa.geometry import apply_transform, kabsch_fit, wrap_positions
from multiscale.cg2aa.molecules import rebuild_hydrogens
from multiscale.cg2aa.refinement import geometry_audit, refine_reactive_geometry


def _add_atom(molecule, symbol: str, molecule_id: int) -> int:
    index = molecule.AddAtom(Chem.Atom(symbol))
    molecule.GetAtomWithIdx(index).SetIntProp("molecule_id", molecule_id)
    return index


def test_kabsch_fit_recovers_row_vector_rigid_transform():
    source = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.2, 1.1, 0.0],
            [0.1, 0.3, 1.4],
        ]
    )
    angle = np.deg2rad(63.0)
    rotation = np.asarray(
        [
            [np.cos(angle), np.sin(angle), 0.0],
            [-np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    translation = np.asarray([2.5, -1.2, 0.7])
    target = source @ rotation + translation

    fitted_rotation, fitted_translation = kabsch_fit(source, target)
    fitted = apply_transform(source, fitted_rotation, fitted_translation)
    assert np.linalg.det(fitted_rotation) > 0.0
    assert np.max(np.abs(fitted - target)) < 1.0e-10


def test_wrap_positions_reconstructs_triclinic_coordinates():
    box = (4.0, 5.0, 6.0, 0.4, -0.2, 0.3)
    positions = np.asarray(
        [
            [0.1, 0.2, 0.3],
            [5.0, -3.0, 7.5],
            [-4.7, 6.2, -8.1],
        ]
    )
    wrapped, images = wrap_positions(positions, box)
    matrix = np.asarray(
        [
            [box[0], 0.0, 0.0],
            [box[3], box[1], 0.0],
            [box[4], box[5], box[2]],
        ]
    )
    fractional = wrapped @ np.linalg.inv(matrix)
    assert np.all(fractional >= -0.5)
    assert np.all(fractional < 0.5)
    assert np.max(np.abs(positions - (wrapped + images @ matrix))) < 1.0e-12


def test_reactive_geometry_repairs_a_collapsed_cross_instance_angle():
    molecule = Chem.RWMol()
    left = _add_atom(molecule, "C", 0)
    center = _add_atom(molecule, "N", 0)
    right = _add_atom(molecule, "C", 1)
    molecule.AddBond(left, center, Chem.BondType.SINGLE)
    molecule.AddBond(center, right, Chem.BondType.SINGLE)

    angle = np.deg2rad(50.0)
    coords = [
        np.asarray([1.40, 0.0, 0.0]),
        np.asarray([0.0, 0.0, 0.0]),
        np.asarray([1.47 * np.cos(angle), 1.47 * np.sin(angle), 0.0]),
    ]
    box = (50.0, 50.0, 50.0)

    before = geometry_audit(molecule, coords, box, coordinate_to_nm=0.1)
    assert before["checks"]["reaction_bond_lengths_reasonable"] is True
    assert before["checks"]["reaction_center_angles_reasonable"] is False

    refined, report = refine_reactive_geometry(
        molecule,
        coords,
        box,
        coordinate_to_nm=0.1,
    )
    after = geometry_audit(molecule, refined, box, coordinate_to_nm=0.1)
    assert report["applied"] is True
    assert after["checks"]["reaction_center_angles_reasonable"] is True
    assert 100.0 <= after["reaction_center_angles"]["min_deg"] <= 140.0
    assert after["ok"] is True


def test_geometry_audit_detects_nonbonded_hydrogen_overlap():
    molecule = Chem.RWMol()
    carbon_left = _add_atom(molecule, "C", 0)
    hydrogen_left = _add_atom(molecule, "H", 0)
    carbon_right = _add_atom(molecule, "C", 1)
    hydrogen_right = _add_atom(molecule, "H", 1)
    molecule.AddBond(carbon_left, hydrogen_left, Chem.BondType.SINGLE)
    molecule.AddBond(carbon_right, hydrogen_right, Chem.BondType.SINGLE)

    coords = [
        np.asarray([0.0, 0.0, 0.0]),
        np.asarray([1.09, 0.0, 0.0]),
        np.asarray([2.18, 0.0, 0.0]),
        np.asarray([1.09, 0.0, 0.0]),
    ]
    report = geometry_audit(
        molecule,
        coords,
        (50.0, 50.0, 50.0),
        coordinate_to_nm=0.1,
    )
    assert report["checks"]["hydrogen_bond_lengths_reasonable"] is True
    assert report["checks"]["no_severe_interinstance_heavy_atom_clashes"] is True
    assert report["checks"]["no_severe_nonbonded_all_atom_clashes"] is False
    assert report["severe_nonbonded_all_atom_clashes"]["count"] >= 1
    assert report["ok"] is False


def test_hydrogen_rebuild_avoids_severe_periodic_nonbonded_clashes():
    molecule = Chem.RWMol()
    _add_atom(molecule, "C", 0)
    _add_atom(molecule, "C", 1)
    rebuilt, coords = rebuild_hydrogens(
        molecule,
        [np.asarray([0.0, 0.0, 0.0]), np.asarray([2.18, 0.0, 0.0])],
        box=(50.0, 50.0, 50.0),
    )
    report = geometry_audit(
        rebuilt,
        coords,
        (50.0, 50.0, 50.0),
        coordinate_to_nm=0.1,
    )
    assert rebuilt.GetNumAtoms() == 10
    assert report["hydrogen_bonds"]["count"] == 8
    assert report["checks"]["hydrogen_bond_lengths_reasonable"] is True
    assert report["checks"]["no_severe_nonbonded_all_atom_clashes"] is True
    assert report["ok"] is True
