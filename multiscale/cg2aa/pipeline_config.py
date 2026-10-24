"""Formal ``multiscale.cg2aa.v1`` configuration loading and compilation."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any

SCHEMA = "multiscale.cg2aa.v1"
RESOLVED_SCHEMA = "multiscale.cg2aa.resolved.v1"

_UNIT_TO_NM = {
    "nm": 1.0,
    "nanometer": 1.0,
    "nanometers": 1.0,
    "angstrom": 0.1,
    "angstroms": 0.1,
    "a": 0.1,
}

_EVENT_KEYS = {
    "source",
    "bond_types",
    "allow_reverse",
    "section",
    "columns",
    "index_space",
    "center_type",
    "neighbor_types",
    "cases",
    "event",
}


def _object(value: Any, label: str) -> dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object")
    return value


def _reject_unknown_keys(data: dict[str, Any], allowed: set[str], label: str) -> None:
    unknown = sorted(set(data) - allowed)
    if unknown:
        fields = ", ".join(unknown)
        raise ValueError(f"{label} contains unknown field(s): {fields}")


def _string(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty string")
    return value.strip()


def _boolean(value: Any, label: str, *, default: bool) -> bool:
    if value is None:
        return default
    if not isinstance(value, bool):
        raise ValueError(f"{label} must be true or false")
    return value


def _unit(value: Any, label: str, default: str) -> str:
    name = default if value is None else _string(value, label).lower()
    if name not in _UNIT_TO_NM:
        raise ValueError(f"{label} must be 'nm' or 'angstrom'")
    return "angstrom" if name in {"angstrom", "angstroms", "a"} else "nm"


def _factor(source: str, target: str) -> float:
    return _UNIT_TO_NM[source] / _UNIT_TO_NM[target]


def _absolute(base: Path, value: str) -> str:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    return str(path.resolve())


def _compile_units(value: Any) -> tuple[dict[str, float], dict[str, str]]:
    data = _object(value, "units")
    _reject_unknown_keys(data, {"cg_xml", "template_pdb", "output_xml"}, "units")
    cg_unit = _unit(data.get("cg_xml"), "units.cg_xml", "nm")
    template_unit = _unit(data.get("template_pdb"), "units.template_pdb", "angstrom")
    output_unit = _unit(data.get("output_xml"), "units.output_xml", "nm")
    numeric = {
        "cg_position_scale": _factor(cg_unit, template_unit),
        "xml_position_scale": _factor(template_unit, output_unit),
        "box_scale": _factor(cg_unit, output_unit),
        "angstrom_to_coordinate": _factor("angstrom", template_unit),
        "coordinate_to_nm": _UNIT_TO_NM[template_unit],
        "coordinate_to_angstrom": _factor(template_unit, "angstrom"),
    }
    declared = {
        "cg_xml": cg_unit,
        "template_pdb": template_unit,
        "output_xml": output_unit,
    }
    return numeric, declared


def _compile_templates(value: Any, *, base: Path) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    templates = _object(value, "templates")
    if not templates:
        raise ValueError("templates must contain at least one template")

    compiled: dict[str, dict[str, Any]] = {}
    resolved: dict[str, dict[str, Any]] = {}
    for name, raw_value in templates.items():
        raw = _object(raw_value, f"templates.{name}")
        _reject_unknown_keys(raw, {"cg_type", "pdb", "smiles", "add_hydrogens", "smiles_to_pdb_atoms", "sites", "fit"}, f"templates.{name}")
        pdb = _string(raw.get("pdb"), f"templates.{name}.pdb")
        item: dict[str, Any] = {
            "type": _string(raw.get("cg_type", name), f"templates.{name}.cg_type"),
            "pdb": pdb,
            "atom_index_base": 1,
        }
        resolved_item: dict[str, Any] = {
            "cg_type": item["type"],
            "pdb": _absolute(base, pdb),
        }
        for key in ("smiles", "add_hydrogens", "smiles_to_pdb_atoms"):
            if key in raw:
                item[key] = deepcopy(raw[key])
                resolved_item[key] = deepcopy(raw[key])

        sites = raw.get("sites")
        if not isinstance(sites, list) or not sites:
            raise ValueError(f"templates.{name}.sites must be a non-empty array")
        compiled_sites: list[dict[str, Any]] = []
        resolved_sites: list[dict[str, Any]] = []
        for index, site_value in enumerate(sites):
            site = _object(site_value, f"templates.{name}.sites[{index}]")
            _reject_unknown_keys(site, {"name", "atoms", "group"}, f"templates.{name}.sites[{index}]")
            site_name = _string(site.get("name"), f"templates.{name}.sites[{index}].name")
            atoms = site.get("atoms")
            if not isinstance(atoms, list) or not atoms:
                raise ValueError(f"templates.{name}.sites[{index}].atoms must be a non-empty integer array")
            if any(isinstance(atom, bool) or not isinstance(atom, int) or atom < 1 for atom in atoms):
                raise ValueError(f"templates.{name}.sites[{index}].atoms must contain positive 1-based integers")
            compiled_site = {"name": site_name, "atoms": list(atoms)}
            resolved_site = {"name": site_name, "atoms": list(atoms)}
            if site.get("group") is not None:
                group = _string(site["group"], f"templates.{name}.sites[{index}].group")
                compiled_site["label"] = group
                resolved_site["group"] = group
            compiled_sites.append(compiled_site)
            resolved_sites.append(resolved_site)
        item["sites"] = compiled_sites
        resolved_item["sites"] = resolved_sites

        fit = _object(raw.get("fit"), f"templates.{name}.fit")
        _reject_unknown_keys(fit, {"mode", "sites", "permute_equivalent"}, f"templates.{name}.fit")
        if fit:
            mode = _string(fit.get("mode", "kabsch"), f"templates.{name}.fit.mode")
            compiled_fit: dict[str, Any] = {"mode": mode}
            resolved_fit: dict[str, Any] = {"mode": mode}
            if fit.get("sites") is not None:
                if not isinstance(fit["sites"], list) or not fit["sites"]:
                    raise ValueError(f"templates.{name}.fit.sites must be a non-empty string array")
                compiled_fit["sites"] = [_string(v, f"templates.{name}.fit.sites") for v in fit["sites"]]
                resolved_fit["sites"] = list(compiled_fit["sites"])
            permute = _boolean(
                fit.get("permute_equivalent"),
                f"templates.{name}.fit.permute_equivalent",
                default=False,
            )
            compiled_fit["auto_permute_equivalent_sites"] = permute
            resolved_fit["permute_equivalent"] = permute
            item["fit"] = compiled_fit
            resolved_item["fit"] = resolved_fit

        compiled[str(name)] = item
        resolved[str(name)] = resolved_item
    return compiled, resolved


def _compile_placement(value: Any, templates: dict[str, dict[str, Any]]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if not isinstance(value, list) or not value:
        raise ValueError("placement must be a non-empty JSON array")
    sequence: list[dict[str, Any]] = []
    resolved: list[dict[str, Any]] = []
    for index, entry_value in enumerate(value):
        entry = _object(entry_value, f"placement[{index}]")
        _reject_unknown_keys(entry, {"template", "count", "sites", "cg_type"}, f"placement[{index}]")
        template_name = _string(entry.get("template"), f"placement[{index}].template")
        if template_name not in templates:
            raise ValueError(f"placement[{index}] references unknown template {template_name!r}")
        count = entry.get("count")
        if isinstance(count, bool) or not isinstance(count, int) or count < 1:
            raise ValueError(f"placement[{index}].count must be a positive integer")
        item = deepcopy(entry)
        item["template"] = template_name
        item["count"] = count
        item["type"] = _string(entry.get("cg_type", templates[template_name]["type"]), f"placement[{index}].cg_type")
        item.pop("cg_type", None)
        if item.get("sites") is not None:
            if not isinstance(item["sites"], list) or not item["sites"]:
                raise ValueError(f"placement[{index}].sites must be a non-empty string array")
            item["sites"] = [_string(site, f"placement[{index}].sites") for site in item["sites"]]
        sequence.append(item)
        resolved.append(deepcopy(entry))
    return {
        "mode": "sequential",
        "unwrap_sites": True,
        "require_all_beads": True,
        "sequence": sequence,
    }, resolved


def _compile_reactions(value: Any) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    if value is None:
        raw_reactions: list[Any] = []
    elif isinstance(value, list):
        raw_reactions = value
    else:
        raise ValueError("reactions must be a JSON array")

    templates: list[dict[str, Any]] = []
    events: list[dict[str, Any]] = []
    resolved: list[dict[str, Any]] = []
    for index, reaction_value in enumerate(raw_reactions):
        reaction = _object(reaction_value, f"reactions[{index}]")
        _reject_unknown_keys(
            reaction,
            {"id", "reactants", "smarts", "probability", "keep_products", "templates", *_EVENT_KEYS},
            f"reactions[{index}]",
        )
        reaction_id = _string(reaction.get("id"), f"reactions[{index}].id")
        if reaction.get("smarts") is None:
            raise ValueError(f"reactions[{index}].smarts is required")
        template = {key: deepcopy(item) for key, item in reaction.items() if key not in _EVENT_KEYS}
        templates.append(template)

        event_block = _object(reaction.get("event"), f"reactions[{index}].event")
        _reject_unknown_keys(event_block, (_EVENT_KEYS - {"event"}) | {"reactants", "reaction"}, f"reactions[{index}].event")
        event = deepcopy(event_block)
        for key in _EVENT_KEYS - {"event"}:
            if key in reaction and key not in event:
                event[key] = deepcopy(reaction[key])
        event.setdefault("source", "bonds")
        event["reaction"] = reaction_id
        if event["source"] == "bonds" and "reactants" not in event and reaction.get("reactants") is not None:
            event["reactants"] = deepcopy(reaction["reactants"])
        events.append(event)

        resolved_item = deepcopy(reaction)
        resolved_item.setdefault("source", event["source"])
        resolved_item.setdefault("probability", 1.0)
        resolved.append(resolved_item)
    return templates, events, resolved


def load_pipeline_config(path: str | Path) -> dict[str, Any]:
    config_path = Path(path).expanduser().resolve()
    if not config_path.is_file():
        raise ValueError(f"CG-to-AA JSON configuration not found: {config_path}")
    try:
        raw = json.loads(config_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in {config_path}: {exc}") from exc
    if not isinstance(raw, dict):
        raise ValueError("CG-to-AA configuration root must be a JSON object")
    if raw.get("schema") != SCHEMA:
        raise ValueError(f"schema must be {SCHEMA!r}, got {raw.get('schema')!r}")
    _reject_unknown_keys(raw, {"schema", "name", "input", "units", "templates", "placement", "reactions", "hydrogens", "force_field", "output"}, "configuration")

    base = config_path.parent
    input_data = _object(raw.get("input"), "input")
    _reject_unknown_keys(input_data, {"cg_xml"}, "input")
    cg_xml = _string(input_data.get("cg_xml"), "input.cg_xml")
    numeric_units, declared_units = _compile_units(raw.get("units"))
    templates, resolved_templates = _compile_templates(raw.get("templates"), base=base)
    placement, resolved_placement = _compile_placement(raw.get("placement"), templates)
    reactions, reaction_events, resolved_reactions = _compile_reactions(raw.get("reactions"))

    hydrogens = raw.get("hydrogens")
    if hydrogens is None:
        hydrogen_mode = "rebuild" if reactions else "template"
    else:
        hydrogen_mode = _string(hydrogens, "hydrogens")
    if hydrogen_mode not in {"template", "preserve", "rebuild"}:
        raise ValueError("hydrogens must be 'template', 'preserve', or 'rebuild'")

    force_field = deepcopy(_object(raw.get("force_field"), "force_field"))
    _reject_unknown_keys(force_field, {"name", "atom_defaults", "neutralize", "use_smarts"}, "force_field")
    if "atom_defaults" in force_field:
        force_field["defaults"] = force_field.pop("atom_defaults")
    force_field.setdefault("name", "opls")
    force_field.setdefault("neutralize", True)

    output_raw = _object(raw.get("output"), "output")
    _reject_unknown_keys(output_raw, {"directory", "prefix", "split", "pdb"}, "output")
    output = {
        "directory": output_raw.get("directory", "outputs/cg2aa"),
        "prefix": output_raw.get("prefix", "aa"),
        "split_molecules": _boolean(output_raw.get("split"), "output.split", default=False),
        "pdb": _boolean(output_raw.get("pdb"), "output.pdb", default=True),
    }

    resolved_configuration: dict[str, Any] = {
        "schema": SCHEMA,
        "name": raw.get("name", "CG-to-AA backmapping"),
        "input": {"cg_xml": _absolute(base, cg_xml)},
        "units": declared_units,
        "templates": resolved_templates,
        "placement": resolved_placement,
        "reactions": resolved_reactions,
        "hydrogens": hydrogen_mode,
        "force_field": deepcopy(raw.get("force_field", {})) or {"name": "opls", "neutralize": True},
        "output": {
            "directory": _absolute(base, _string(output["directory"], "output.directory")),
            "prefix": _string(output["prefix"], "output.prefix"),
            "split": output["split_molecules"],
            "pdb": output["pdb"],
        },
    }

    config: dict[str, Any] = {
        "schema": "multiscale.cg2aa.internal.v1",
        "name": raw.get("name", "CG-to-AA backmapping"),
        "xml": cg_xml,
        "units": numeric_units,
        "templates": templates,
        "placement": placement,
        "reactions": reactions,
        "reaction_events": reaction_events,
        "hydrogens": {"mode": hydrogen_mode},
        "force_field": force_field,
        "output": output,
        "_base_dir": str(base),
        "_config_path": str(config_path),
        "_resolved_configuration": resolved_configuration,
    }
    return config


def resolve_path(config: dict[str, Any], value: str | Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path.resolve()
    return (Path(config["_base_dir"]) / path).resolve()


def resolved_pipeline_document(config: dict[str, Any]) -> dict[str, Any]:
    """Return the complete formal configuration with resolved paths/defaults."""

    return {
        "schema": RESOLVED_SCHEMA,
        "resolved": True,
        "source_config": config.get("_config_path"),
        "configuration": deepcopy(config["_resolved_configuration"]),
    }
