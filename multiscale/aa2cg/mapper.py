"""Native AA-to-CG single-frame structure mapper.

The mapper consumes an AA ``.gro`` frame, topology metadata read from either a
GROMACS text dump or a binary TPR, and an in-memory bead mapping/topology
definition loaded from the formal JSON configuration. It produces a
CGSystem with mass-weighted bead coordinates, mapping diagnostics, optional
visual-validation artifacts (``cg.xmol`` and ``view_mapping.vmd``), and the
backend-independent ``cg_system_mapped.json``. No trajectory processing or
potential-table generation is performed here.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
from typing import Any

from multiscale.aa2cg.model import MappingSpec, TopologySpec
from multiscale.core.cg_system import CGAngle, CGBead, CGBond, CGSystem
from multiscale.aa2cg.gromacs import parse_gro
from multiscale.aa2cg.tpr_dump import TprDumpAtom, TprDumpTopology, TprTopologyInput, load_tpr_topology


@dataclass(frozen=True)
class StructureMapReport:
    outdir: str
    files: dict[str, str]
    counts: dict[str, Any]
    checks: dict[str, bool | None]
    warnings: list[str]
    topology_source: dict[str, Any] | None = None
    entrypoint: str = "run"
    backend: str = "internal"

    @property
    def ok(self) -> bool:
        return all(value is not False for value in self.checks.values())

    def to_dict(self) -> dict[str, Any]:
        return {
            "outdir": self.outdir,
            "files": self.files,
            "counts": self.counts,
            "checks": self.checks,
            "warnings": self.warnings,
            "topology_source": self.topology_source,
            "entrypoint": self.entrypoint,
            "backend": self.backend,
            "ok": self.ok,
        }


def _invert_cgtypes(groups: dict[str, list[str]]) -> dict[str, str]:
    return {name: group for group, names in groups.items() for name in names}


def _template_by_index(atoms: list[TprDumpAtom]) -> dict[int, TprDumpAtom]:
    return {atom.index: atom for atom in atoms}


def _topology_input_path(source: TprTopologyInput) -> str:
    if isinstance(source, TprDumpTopology):
        return source.path
    return str(source)


def _cg_name(mapping_name: str | None, topology_name: str | None) -> str:
    if mapping_name:
        return mapping_name.replace(".CG", "")
    if topology_name:
        return f"{topology_name}_CG"
    return "cg_system"


def _topology_group_counts(topology: TopologySpec | None) -> dict[str, list[int] | int]:
    if topology is None:
        return {"bond_groups": [], "angle_groups": [], "bond_total_per_molecule": 0, "angle_total_per_molecule": 0}
    bond_groups = [len(group.pairs) for group in topology.bond_groups]
    angle_groups = [len(group.triples) for group in topology.angle_groups]
    return {
        "bond_groups": bond_groups,
        "angle_groups": angle_groups,
        "bond_total_per_molecule": sum(bond_groups),
        "angle_total_per_molecule": sum(angle_groups),
    }


def _expand_bonds_angles(
    topology: TopologySpec,
    *,
    bead_name_to_local_index: dict[str, int],
    name: str,
    nmol: int,
    beads_per_molecule: int,
) -> tuple[list[CGBond], list[CGAngle]]:
    bonds: list[CGBond] = []
    angles: list[CGAngle] = []

    def local_index(bead_name: str) -> int:
        try:
            return bead_name_to_local_index[bead_name]
        except KeyError as exc:
            known = ", ".join(sorted(bead_name_to_local_index))
            raise ValueError(f"Configured topology refers to unknown bead {bead_name!r}; known beads: {known}") from exc

    for mol_index in range(1, nmol + 1):
        offset = (mol_index - 1) * beads_per_molecule
        for group in topology.bond_groups:
            bond_type = f"{name}:{group.type_id}"
            for left, right in group.pairs:
                bonds.append(
                    CGBond(
                        i=offset + local_index(left),
                        j=offset + local_index(right),
                        bond_type=bond_type,
                        molecule_index=mol_index,
                    )
                )
        for group in topology.angle_groups:
            angle_type = f"{name}:{group.type_id}"
            for left, center, right in group.triples:
                angles.append(
                    CGAngle(
                        i=offset + local_index(left),
                        j=offset + local_index(center),
                        k=offset + local_index(right),
                        angle_type=angle_type,
                        molecule_index=mol_index,
                    )
                )
    return bonds, angles


def _length_to_angstrom(value: float, unit: str | None) -> float:
    normalized = (unit or "angstrom").lower()
    if normalized in {"angstrom", "angstroms", "a", "å"}:
        return float(value)
    if normalized in {"nm", "nanometer", "nanometers"}:
        return float(value) * 10.0
    raise ValueError(f"Unsupported length unit for XMOL/VMD export: {unit!r}")


def _system_box_angstrom(system: CGSystem) -> tuple[float, float, float] | None:
    if system.box is None:
        return None
    return tuple(_length_to_angstrom(float(v), system.box_unit) for v in system.box)  # type: ignore[return-value]


def _bead_position_angstrom(system: CGSystem, bead: CGBead) -> tuple[float, float, float]:
    if bead.position is None:
        return (0.0, 0.0, 0.0)
    return tuple(_length_to_angstrom(float(v), system.coordinate_unit) for v in bead.position)  # type: ignore[return-value]


def render_xmol(system: CGSystem, *, step: int = 0) -> str:
    """Render an XMOL coordinate frame for visual validation.

    The coordinates are written in Angstrom.  This is intentional because XMOL
    does not carry robust unit metadata and is mainly used here for VMD overlay
    checks against an AA GRO structure.
    """
    box = _system_box_angstrom(system)
    if box:
        comment = f" after {step:12d}. BOX: {box[0]:9.3f} {box[1]:9.3f} {box[2]:9.3f}"
    else:
        comment = f" after {step:12d}. BOX:    0.000    0.000    0.000"
    lines = [f"{len(system.beads):12d}", comment]
    for bead in sorted(system.beads, key=lambda b: b.index):
        x, y, z = _bead_position_angstrom(system, bead)
        name = str(bead.name or bead.bead_type or f"B{bead.index}")[:16]
        lines.append(f"{name:<16}{x:15.5f}{y:15.5f}{z:15.5f}")
    return "\n".join(lines) + "\n"


def _tcl_quote(path: str | Path) -> str:
    text = str(path).replace("\\", "/")
    return "{" + text.replace("}", "\\}") + "}"


def render_vmd_script(
    *,
    aa_gro: str | Path,
    cg_xmol: str | Path,
    cg_box_angstrom: tuple[float, float, float] | None = None,
) -> str:
    lines = [
        "# Generated by multiscale visual validation support.",
        "# Load AA GRO and mapped CG XMOL in VMD to review AA -> CG geometry.",
        "# cg.xmol is an XYZ/XMOL-style visual-validation file.",
        "# VMD loads it with the xyz reader; its BOX line is a comment, not a standard VMD unit-cell record.",
        f"set aa_gro {_tcl_quote(aa_gro)}",
        f"set cg_xmol {_tcl_quote(cg_xmol)}",
        "",
        "mol new $aa_gro type gro waitfor all",
        "set aa_mol [molinfo top get id]",
        "mol rename $aa_mol AA_reference_GRO",
        "mol delrep 0 $aa_mol",
        "mol representation Lines 1.0",
        "mol color Name",
        "mol selection all",
        "mol material Transparent",
        "mol addrep $aa_mol",
        "",
        "mol new $cg_xmol type xyz waitfor all",
        "set cg_mol [molinfo top get id]",
        "mol rename $cg_mol CG_mapped_XMOL",
    ]
    if cg_box_angstrom is not None:
        bx, by, bz = cg_box_angstrom
        lines.extend([
            f"molinfo $cg_mol set a {bx:.6f}",
            f"molinfo $cg_mol set b {by:.6f}",
            f"molinfo $cg_mol set c {bz:.6f}",
            "molinfo $cg_mol set alpha 90",
            "molinfo $cg_mol set beta 90",
            "molinfo $cg_mol set gamma 90",
        ])
    lines.extend([
        "mol delrep 0 $cg_mol",
        "mol representation VDW 1.2 20",
        "mol color Name",
        "mol selection all",
        "mol material Opaque",
        "mol addrep $cg_mol",
        "",
        "# To obtain a GROMACS .gro of the mapped CG structure, run:",
        "#   multiscale export-gromacs <cg_system_mapped.json> --outdir <gromacs_outdir>",
        "# and load <gromacs_outdir>/cg.gro in VMD if needed.",
        "",
        "catch {package require pbctools; pbc box -molid $aa_mol -color white}",
        "display resetview",
        "axes location Off",
        "puts \"Loaded AA GRO: $aa_gro\"",
        "puts \"Loaded mapped CG XMOL: $cg_xmol\"",
        "puts \"Visual check: CG beads should sit on corresponding AA molecular regions; exact atom overlap is not expected.\"",
    ])
    return "\n".join(lines) + "\n"


def map_structure_to_cgsystem(
    *,
    topology_dump: TprTopologyInput,
    coordinates: str | Path,
    mapping: MappingSpec,
    cg_topology: TopologySpec | None = None,
) -> tuple[CGSystem, StructureMapReport]:
    tpr = load_tpr_topology(topology_dump)
    gro = parse_gro(coordinates)
    cfg = mapping
    topology = cg_topology
    mapping_source = cfg.source
    topology_input_source = topology.source if topology else None
    cgtypes = topology.cgtypes if topology else {}
    cgtype_by_name = _invert_cgtypes(cgtypes)

    warnings: list[str] = []
    if not cfg.beads:
        raise ValueError(f"No bead mappings found in {mapping_source or 'AA-to-CG configuration'}")
    max_mapped_atom = max(max(bead.atom_indices_1based or [0]) for bead in cfg.beads)
    atoms_per_molecule = tpr.atoms_per_molecule or max_mapped_atom
    if atoms_per_molecule < max_mapped_atom:
        raise ValueError(
            f"Mapping refers to atom index {max_mapped_atom}, but TPR topology has only {atoms_per_molecule} atoms per molecule"
        )
    if gro.natoms % atoms_per_molecule != 0:
        raise ValueError(
            f"GRO atom count {gro.natoms} is not divisible by atoms_per_molecule={atoms_per_molecule}"
        )
    nmol = gro.natoms // atoms_per_molecule
    if tpr.molecules is not None and tpr.molecules != nmol:
        warnings.append(f"TPR topology molecule count {tpr.molecules} differs from GRO-derived molecule count {nmol}; using GRO-derived count.")
    if tpr.total_atoms is not None and tpr.total_atoms != gro.natoms:
        warnings.append(f"TPR topology total atom count {tpr.total_atoms} differs from GRO atom count {gro.natoms}; using GRO coordinates.")

    atom_template = _template_by_index(tpr.atoms)
    if len(atom_template) < max_mapped_atom:
        raise ValueError(
            f"TPR topology atom template has {len(atom_template)} parsed atoms, but mapping needs atom index {max_mapped_atom}"
        )

    name = _cg_name(cfg.name, tpr.name)
    beads_per_molecule = len(cfg.beads)
    beads: list[CGBead] = []
    mapping_debug: list[dict[str, Any]] = []
    for mol_index in range(1, nmol + 1):
        atom_offset = (mol_index - 1) * atoms_per_molecule
        for mol_bead_index, bead_map in enumerate(cfg.beads, start=1):
            total_mass = 0.0
            total_charge = 0.0
            weighted = [0.0, 0.0, 0.0]
            mapped_names: list[str | None] = []
            global_atom_indices: list[int] = []
            for local_atom_index in bead_map.atom_indices_1based:
                atom_meta = atom_template[local_atom_index]
                global_atom_index = atom_offset + local_atom_index
                gro_atom = gro.atoms[global_atom_index - 1]
                mass = float(atom_meta.mass)
                total_mass += mass
                total_charge += float(atom_meta.charge)
                weighted[0] += mass * gro_atom.x
                weighted[1] += mass * gro_atom.y
                weighted[2] += mass * gro_atom.z
                mapped_names.append(atom_meta.name)
                global_atom_indices.append(global_atom_index)
            if total_mass <= 0:
                raise ValueError(f"Mapped bead {bead_map.name} in molecule {mol_index} has non-positive total mass")
            index = (mol_index - 1) * beads_per_molecule + mol_bead_index
            position = (weighted[0] / total_mass, weighted[1] / total_mass, weighted[2] / total_mass)
            bead = CGBead(
                index=index,
                name=bead_map.name,
                bead_type=cgtype_by_name.get(bead_map.name, bead_map.name),
                molecule_index=mol_index,
                molecule_bead_index=mol_bead_index,
                mass=total_mass,
                charge=total_charge,
                position=position,
            )
            beads.append(bead)
            mapping_debug.append(
                {
                    "cg_index": index,
                    "molecule_index": mol_index,
                    "molecule_bead_index": mol_bead_index,
                    "bead_name": bead_map.name,
                    "bead_type": bead.bead_type,
                    "local_atom_indices": bead_map.atom_indices_1based,
                    "global_atom_indices": global_atom_indices,
                    "aa_atom_names": mapped_names,
                    "mass": total_mass,
                    "charge": total_charge,
                    "x_nm": position[0],
                    "y_nm": position[1],
                    "z_nm": position[2],
                    "x_angstrom": position[0] * 10.0,
                    "y_angstrom": position[1] * 10.0,
                    "z_angstrom": position[2] * 10.0,
                }
            )

    bead_name_to_local_index = {bead.name: idx for idx, bead in enumerate(cfg.beads, start=1)}
    topology_source = "none"
    topology_counts = _topology_group_counts(topology)

    if topology and topology.has_bond_angle_topology:
        bonds, angles = _expand_bonds_angles(
            topology,
            bead_name_to_local_index=bead_name_to_local_index,
            name=name,
            nmol=nmol,
            beads_per_molecule=beads_per_molecule,
        )
        topology_source = "configuration"
    else:
        bonds, angles = [], []
        warnings.append(
            "No CG bond/angle topology was configured; the mapped CGSystem contains coordinates, masses, and charges but no CG bonds or angles."
        )

    system = CGSystem(
        name=name,
        box=gro.box or tpr.box,
        box_unit="nm" if gro.box else ("nm" if tpr.box else None),
        coordinate_unit="nm",
        index_base=1,
        nmol=nmol,
        beads_per_molecule=beads_per_molecule,
        beads=beads,
        bonds=bonds,
        angles=angles,
        metadata={
            "source_files": {
                "topology_dump": _topology_input_path(topology_dump),
                "coordinates": str(coordinates),
                "mapping": mapping_source,
                "cg_topology": topology_input_source,
            },
            "mapping": cfg.to_dict(),
            "mapping_debug": mapping_debug,
            "bead_type_groups": cgtypes,
            "cg_topology": topology.to_dict() if topology else None,
            "cg_topology_source": topology_source,
            "cg_topology_counts": topology_counts,
            "aa_topology": {
                "name": tpr.name,
                "total_atoms": tpr.total_atoms,
                "molecules": tpr.molecules,
                "atoms_per_molecule": atoms_per_molecule,
                "parsed_atom_template_atoms": len(tpr.atoms),
                "reader": tpr.reader,
            },
            "aa_coordinates": {
                "title": gro.title,
                "natoms": gro.natoms,
                "box_nm": list(gro.box) if gro.box else None,
            },
            "visual_validation": {
                "primary_overlay": "Load AA GRO and generated cg.xmol in VMD; CG beads should sit on corresponding AA molecular regions.",
                "xmol_coordinate_unit": "angstrom",
            },
            "notes": [
                "Coordinates are mass-weighted centers of mapped AA atoms from a single GRO frame.",
                "Masses and charges are summed from the TPR atom template.",
                "CG bonds and angles are built from the topology declared in the AA-to-CG configuration.",
                "Generated XMOL is intended for VMD visual overlay validation before selecting GALAMOST/GROMACS/LAMMPS export.",
                "The mapping step writes CGSystem JSON; target-specific files are produced by the selected exporters.",
                "No trajectory is read during single-frame mapping. A trajectory is required only when RDF analysis is enabled.",
                "The supplied coordinate file is used as-is. Multiscale does not unwrap, recenter, or otherwise modify periodic coordinates; the user is responsible for coordinate suitability.",
            ],
        },
    )

    checks = {
        "has_beads": len(beads) > 0,
        "gro_atoms_divisible_by_atoms_per_molecule": gro.natoms % atoms_per_molecule == 0,
        "nmol_times_beads_per_molecule_matches_cg_atoms": nmol * beads_per_molecule == len(beads),
        "mapping_counts_consistent": all(bead.is_count_consistent for bead in cfg.beads),
        "tpr_template_covers_mapping": len(atom_template) >= max_mapped_atom,
        "has_box": system.box is not None,
        "has_cg_bonds_angles": bool(bonds or angles) if topology else None,
        "has_configured_bonds_angles": bool(topology and topology.has_bond_angle_topology) if topology else None,
        "box_inherited_from_aa_gro": tuple(system.box) == tuple(gro.box) if system.box and gro.box else None,
    }
    report = StructureMapReport(
        outdir="",
        files={},
        counts={
            "aa_atoms": gro.natoms,
            "atoms_per_molecule": atoms_per_molecule,
            "nmol": nmol,
            "beads_per_molecule": beads_per_molecule,
            "cg_atoms": len(beads),
            "cg_bonds": len(bonds),
            "cg_angles": len(angles),
            "cg_topology_source": topology_source,
            "bond_groups": len(topology.bond_groups) if topology else 0,
            "angle_groups": len(topology.angle_groups) if topology else 0,
            "mapping_debug_rows": len(mapping_debug),
            "xmol_coordinate_unit": "angstrom",
        },
        checks=checks,
        warnings=warnings,
    )
    return system, report


def write_mapping_debug_csv(system: CGSystem, out: str | Path) -> None:
    rows = system.metadata.get("mapping_debug") or []
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "cg_index",
        "molecule_index",
        "molecule_bead_index",
        "bead_name",
        "bead_type",
        "local_atom_indices",
        "global_atom_indices",
        "aa_atom_names",
        "mass",
        "charge",
        "x_nm",
        "y_nm",
        "z_nm",
        "x_angstrom",
        "y_angstrom",
        "z_angstrom",
    ]
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            serializable = dict(row)
            for key in ("local_atom_indices", "global_atom_indices", "aa_atom_names"):
                serializable[key] = ";".join(str(x) for x in serializable.get(key, []))
            writer.writerow({key: serializable.get(key, "") for key in fieldnames})


def write_mapped_structure(
    *,
    topology_dump: TprTopologyInput,
    coordinates: str | Path,
    mapping: MappingSpec,
    outdir: str | Path,
    cg_topology: TopologySpec | None = None,
    topology_source: dict[str, Any] | None = None,
    entrypoint: str = "run",
) -> StructureMapReport:
    system, report = map_structure_to_cgsystem(
        topology_dump=topology_dump,
        coordinates=coordinates,
        mapping=mapping,
        cg_topology=cg_topology,
    )
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    xmol_path = outdir / "cg.xmol"
    json_path = outdir / "cg_system_mapped.json"
    debug_csv_path = outdir / "mapping_debug.csv"
    vmd_path = outdir / "view_mapping.vmd"
    report_path = outdir / "map_structure_report.json"

    xmol_path.write_text(render_xmol(system), encoding="utf-8")
    data = system.to_dict()
    data["schema"] = "multiscale.cgsystem.v0"
    json_path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    write_mapping_debug_csv(system, debug_csv_path)
    vmd_path.write_text(
        render_vmd_script(
            aa_gro=coordinates,
            cg_xmol=xmol_path,
            cg_box_angstrom=tuple(float(x) * 10.0 for x in system.box),
        ),
        encoding="utf-8",
    )

    final_report = StructureMapReport(
        outdir=str(outdir),
        files={
            "xmol": str(xmol_path),
            "cgsystem": str(json_path),
            "mapping_debug_csv": str(debug_csv_path),
            "vmd_script": str(vmd_path),
            "report": str(report_path),
        },
        counts=report.counts,
        checks=report.checks,
        warnings=report.warnings,
        topology_source=topology_source
        or {
            "mode": "dump" if not isinstance(topology_dump, TprDumpTopology) else "tpr",
            "dump_path": None if isinstance(topology_dump, TprDumpTopology) else str(topology_dump),
            "tpr_path": topology_dump.path if isinstance(topology_dump, TprDumpTopology) else None,
            "reader": topology_dump.reader if isinstance(topology_dump, TprDumpTopology) else "gmx_dump",
            "dump_generated": False,
        },
        entrypoint=entrypoint,
        backend="internal",
    )
    report_path.write_text(json.dumps(final_report.to_dict(), indent=2, ensure_ascii=False), encoding="utf-8")
    return final_report
