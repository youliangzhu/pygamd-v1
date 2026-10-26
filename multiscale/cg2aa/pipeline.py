from __future__ import annotations

import json
from itertools import permutations
from pathlib import Path

import numpy as np

from .cg_xml import read_cg_xml
from .chem import Chem
from .pipeline_config import load_pipeline_config, resolve_path, resolved_pipeline_document
from .pipeline_forcefield import ff
from .geometry import apply_transform, frame_fit, kabsch_fit, wrap_positions
from .molecules import rebuild_hydrogens, set_missing_h_molecule_ids, split_molecule_with_coords
from .refinement import geometry_audit, refine_reactive_geometry
from .reactions import (
    apply_reactions,
    build_instance_graph,
    build_reaction_events,
    build_reaction_templates,
)
from .templates import TemplateSpec, load_template
from .writer import write_forcefield_xml, write_pdb


def _sequential_instances(config: dict, cg, templates: dict[str, TemplateSpec]) -> list[dict]:
    placement = config.get("placement", {})
    sequence = placement.get("sequence")
    if not sequence:
        raise ValueError("placement.sequence is required for sequential placement")
    cursor = int(placement.get("start_index", 0))
    instances = []
    for entry in sequence:
        template_name = entry["template"]
        template = templates[template_name]
        site_order = entry.get("sites", template.site_names)
        count = int(entry["count"])
        for _ in range(count):
            if cursor + len(site_order) > cg.natoms:
                raise ValueError("sequential placement consumed more CG beads than available")
            site_to_bead = {site: cursor + offset for offset, site in enumerate(site_order)}
            instances.append(
                {
                    "id": len(instances),
                    "template": template_name,
                    "type": entry.get("type", template.type),
                    "site_to_bead": site_to_bead,
                }
            )
            cursor += len(site_order)
    if placement.get("require_all_beads", True) and cursor != cg.natoms:
        raise ValueError(f"sequential placement consumed {cursor} CG beads, but XML has {cg.natoms}")
    return instances


def build_instances(config: dict, cg, templates: dict[str, TemplateSpec]) -> list[dict]:
    placement = config.get("placement", {})
    if placement.get("mode", "sequential") != "sequential":
        raise ValueError("Only placement.mode='sequential' is implemented.")
    return _sequential_instances(config, cg, templates)


def _minimum_image(delta: np.ndarray, box: tuple[float, float, float]) -> np.ndarray:
    adjusted = np.asarray(delta, dtype=float).copy()
    for dim, length in enumerate(box):
        if length > 0.0:
            adjusted[dim] -= length * np.round(adjusted[dim] / length)
    return adjusted


def _target_points(cg, site_to_bead: dict[str, int], sites: list[str], scale: float, unwrap: bool) -> np.ndarray:
    points = []
    for site in sites:
        bead = int(site_to_bead[site])
        if bead < 0 or bead >= cg.natoms:
            raise ValueError(f"CG bead index {bead} out of range")
        points.append(cg.positions[bead].copy())
    points = np.asarray(points, dtype=float)
    if unwrap and len(points) > 1:
        box = tuple(float(v) for v in cg.box[:3])
        unwrapped = [points[0]]
        for point in points[1:]:
            unwrapped.append(unwrapped[-1] + _minimum_image(point - unwrapped[-1], box))
        points = np.asarray(unwrapped, dtype=float)
    return points * scale


def _source_points(template: TemplateSpec, sites: list[str]) -> np.ndarray:
    return np.asarray([template.site_centers[site] for site in sites], dtype=float)


def _fit_template(template: TemplateSpec, source: np.ndarray, target: np.ndarray):
    if template.fit_mode == "frame":
        return frame_fit(source, target)
    return kabsch_fit(source, target)


def _fit_error(source: np.ndarray, target: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> float:
    fitted = apply_transform(source, rotation, translation)
    delta = fitted - target
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def _site_permutation_candidates(template: TemplateSpec, sites: list[str]) -> list[list[str]]:
    if not template.auto_permute_equivalent_sites:
        return [list(sites)]

    groups: dict[str, list[int]] = {}
    for idx, site in enumerate(sites):
        groups.setdefault(template.site_labels.get(site, site), []).append(idx)
    candidates = [list(sites)]
    for indices in groups.values():
        if len(indices) < 2:
            continue
        updated = []
        values = [sites[idx] for idx in indices]
        for base in candidates:
            for permuted in permutations(values):
                candidate = list(base)
                for idx, site in zip(indices, permuted):
                    candidate[idx] = site
                if candidate not in updated:
                    updated.append(candidate)
        candidates = updated
    return candidates


def _best_fit_for_instance(template: TemplateSpec, cg, site_to_bead: dict[str, int], cg_position_scale: float, unwrap_sites: bool):
    target_sites = template.anchor_sites
    target = _target_points(cg, site_to_bead, target_sites, cg_position_scale, unwrap_sites)
    best = None
    for source_sites in _site_permutation_candidates(template, template.anchor_sites):
        source = _source_points(template, source_sites)
        rotation, translation = _fit_template(template, source, target)
        error = _fit_error(source, target, rotation, translation)
        if best is None or error < best["error"]:
            best = {
                "source_sites": source_sites,
                "rotation": rotation,
                "translation": translation,
                "error": error,
            }
    return best


def _copy_instance(template: TemplateSpec, coords: np.ndarray, instance: dict, molecule, global_coords: list[np.ndarray]) -> None:
    local_to_global = {}
    atom_names = template.atom_names
    for local_idx, atom in enumerate(template.molecule.GetAtoms()):
        global_idx = molecule.AddAtom(atom)
        global_atom = molecule.GetAtomWithIdx(global_idx)
        global_atom.SetIntProp("molecule_id", int(instance["id"]))
        global_atom.SetProp("TemplateName", template.name)
        global_atom.SetProp("AtomName", atom_names[local_idx])
        try:
            global_atom.GetProp("ResidueName")
        except KeyError:
            global_atom.SetProp("ResidueName", "MOL")
        local_to_global[local_idx] = global_idx
        global_coords.append(np.asarray(coords[local_idx], dtype=float))

    for bond in template.molecule.GetBonds():
        molecule.AddBond(
            local_to_global[bond.GetBeginAtomIdx()],
            local_to_global[bond.GetEndAtomIdx()],
            bond.GetBondType(),
        )
    instance["local_to_global"] = local_to_global


def _build_template_system(config: dict, cg, templates: dict[str, TemplateSpec], instances: list[dict], cg_position_scale: float):
    molecule = Chem.RWMol()
    coords: list[np.ndarray] = []
    bead_lookup = {}
    unwrap_sites = bool(config.get("placement", {}).get("unwrap_sites", True))
    for instance in instances:
        template = templates[instance["template"]]
        missing = [site for site in template.anchor_sites if site not in instance["site_to_bead"]]
        if missing:
            raise ValueError(f"instance {instance['id']} missing anchor sites: {missing}")
        best_fit = _best_fit_for_instance(template, cg, instance["site_to_bead"], cg_position_scale, unwrap_sites)
        rotation = best_fit["rotation"]
        translation = best_fit["translation"]
        instance["fit_source_sites"] = best_fit["source_sites"]
        instance["fit_target_sites"] = list(template.anchor_sites)
        instance["target_to_source_site"] = dict(zip(template.anchor_sites, best_fit["source_sites"]))
        instance["site_chem_atoms"] = {
            site: list(atoms)
            for site, atoms in template.site_chem_atoms.items()
        }
        instance["fit_error"] = best_fit["error"]
        placed_coords = apply_transform(template.coords, rotation, translation)
        _copy_instance(template, placed_coords, instance, molecule, coords)
        for site, bead in instance["site_to_bead"].items():
            source_site = instance["target_to_source_site"].get(site, site)
            bead_lookup[int(bead)] = {
                "instance": instance["id"],
                "site": source_site,
                "target_site": site,
                "source_site": source_site,
                "template": template.name,
                "type": instance["type"],
            }
    return molecule, coords, bead_lookup


def _reactants_meta(templates: dict[str, TemplateSpec]) -> dict[str, dict]:
    meta = {}
    for template in templates.values():
        if template.smiles:
            meta[template.type] = {"smiles": template.smiles}
    return meta


def _prepare_hydrogens(
    config: dict,
    molecule,
    coords: list[np.ndarray],
    reaction_count: int,
    *,
    angstrom_to_coordinate: float,
    box,
):
    mode = config.get("hydrogens", {}).get("mode")
    if mode is None:
        mode = "rebuild" if reaction_count else "template"
    if mode == "rebuild":
        return rebuild_hydrogens(
            molecule,
            coords,
            angstrom_to_coordinate=angstrom_to_coordinate,
            box=box,
        )
    if mode in {"template", "preserve"}:
        return set_missing_h_molecule_ids(molecule), coords
    raise ValueError("hydrogens.mode must be 'template'/'preserve' or 'rebuild'")


def _output_fragments(config: dict, molecule, coords: list[np.ndarray]):
    output_cfg = config.get("output", {})
    if output_cfg.get("split_molecules", True):
        return split_molecule_with_coords(molecule, coords)
    return [(molecule, coords)]


def run_pipeline(config_file: str | Path) -> dict:
    config = load_pipeline_config(config_file)
    cg = read_cg_xml(resolve_path(config, config["xml"]))
    units = config.get("units", {})
    cg_position_scale = float(units.get("cg_position_scale", 10.0))
    xml_position_scale = float(units.get("xml_position_scale", 0.1))
    output_box_scale = float(units.get("box_scale", 1.0))
    angstrom_to_coordinate = float(units.get("angstrom_to_coordinate", 1.0))
    coordinate_to_nm = float(units.get("coordinate_to_nm", 0.1))
    coordinate_to_angstrom = float(units.get("coordinate_to_angstrom", 1.0))
    coordinate_box = tuple(float(value) * cg_position_scale for value in cg.box)
    output_box = tuple(float(value) * output_box_scale for value in cg.box)

    templates = {
        name: load_template(name, raw, config)
        for name, raw in config["templates"].items()
    }
    instances = build_instances(config, cg, templates)
    molecule, coords, bead_lookup = _build_template_system(
        config,
        cg,
        templates,
        instances,
        cg_position_scale,
    )

    instance_graph = build_instance_graph(cg, instances, bead_lookup)
    reaction_templates = build_reaction_templates(config)
    reaction_events = build_reaction_events(config, cg, instance_graph, bead_lookup)
    refinement_report = {
        "applied": False,
        "iterations": 0,
    }
    if reaction_events:
        reactants_meta = _reactants_meta(templates)
        missing = sorted(
            {
                reactant_type
                for meta in reaction_templates.values()
                for reactant_type in meta["reactants_tuple"]
                if reactant_type not in reactants_meta
            }
        )
        if missing:
            raise ValueError(f"Reaction reactants missing template SMILES: {missing}")
        molecule, coords = apply_reactions(
            molecule,
            coords,
            instances,
            reaction_templates,
            reaction_events,
            reactants_meta,
        )
        coords, refinement_report = refine_reactive_geometry(
            molecule,
            coords,
            coordinate_box,
            angstrom_to_coordinate=angstrom_to_coordinate,
            coordinate_to_nm=coordinate_to_nm,
        )

    molecule, coords = _prepare_hydrogens(
        config,
        molecule,
        coords,
        len(reaction_events),
        angstrom_to_coordinate=angstrom_to_coordinate,
        box=coordinate_box,
    )
    geometry_report = geometry_audit(
        molecule,
        coords,
        coordinate_box,
        coordinate_to_nm=coordinate_to_nm,
        angstrom_to_coordinate=angstrom_to_coordinate,
    )
    geometry_report["refinement"] = refinement_report
    fragments = _output_fragments(config, molecule, coords)

    output_cfg = config.get("output", {})
    output_dir = resolve_path(config, output_cfg.get("directory", "outputs"))
    output_dir.mkdir(parents=True, exist_ok=True)
    resolved_config_path = output_dir / "cg2aa_resolved.json"
    with resolved_config_path.open("w", encoding="utf-8") as handle:
        json.dump(resolved_pipeline_document(config), handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    geometry_report_path = output_dir / "cg2aa_geometry_report.json"
    with geometry_report_path.open("w", encoding="utf-8") as handle:
        json.dump(geometry_report, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

    prefix = output_cfg.get("prefix", "aa")
    program = output_cfg.get("program", "galamost")
    version = output_cfg.get("version", "1.3")

    force_field = config.get("force_field", {})
    defaults = force_field.get("defaults", {})
    use_smarts = force_field.get("use_smarts", False)
    neutralize = force_field.get("neutralize", True)

    outputs = []
    inbox_outputs = []
    pdb_outputs = []
    for idx, (fragment, fragment_coords) in enumerate(fragments):
        Chem.SanitizeMol(fragment)
        typed, bonds, angles, dihedrals = ff(
            fragment,
            defaults=defaults,
            use_smarts=use_smarts,
            neutralize=neutralize,
        )
        suffix = "%06d" % idx if output_cfg.get("split_molecules", True) else ""
        stem = f"{prefix}_{suffix}" if suffix else prefix
        xml_path = output_dir / f"{stem}.xml"
        outputs.append(
            write_forcefield_xml(
                xml_path,
                typed,
                fragment_coords,
                output_box,
                bonds,
                angles,
                dihedrals,
                position_scale=xml_position_scale,
                program=program,
                version=version,
            )
        )
        wrapped_coords, images = wrap_positions(fragment_coords, coordinate_box)
        inbox_path = output_dir / f"{stem}_inbox.xml"
        inbox_outputs.append(
            write_forcefield_xml(
                inbox_path,
                typed,
                wrapped_coords,
                output_box,
                bonds,
                angles,
                dihedrals,
                position_scale=xml_position_scale,
                program=program,
                version=version,
                images=images,
            )
        )
        if output_cfg.get("pdb", True):
            pdb_outputs.append(
                write_pdb(
                    output_dir / f"{stem}.pdb",
                    typed,
                    fragment_coords,
                    box=coordinate_box,
                    position_scale=coordinate_to_angstrom,
                )
            )

    summary = {
        "config": str(Path(config_file).resolve()),
        "xml": str(resolve_path(config, config["xml"])),
        "cg_atoms": cg.natoms,
        "instances": len(instances),
        "instance_bonds": instance_graph.number_of_edges(),
        "reaction_events": len(reaction_events),
        "aa_fragments": len(fragments),
        "resolved_config": str(resolved_config_path),
        "geometry_report": str(geometry_report_path),
        "geometry_ok": geometry_report["ok"],
        "aa_outputs": len(outputs),
        "outputs": outputs,
        "inbox_outputs": inbox_outputs,
        "pdb_outputs": pdb_outputs,
    }
    with (output_dir / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    return summary
