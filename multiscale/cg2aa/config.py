"""CG-to-AA configuration summary and validation helpers."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from multiscale.forcefield.opls import validate_config_forcefield

from .pipeline_config import SCHEMA, load_pipeline_config, resolve_path


@dataclass(frozen=True)
class PathCheck:
    key: str
    path: str
    exists: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class CG2AAConfigSummary:
    path: str
    schema: str
    xml: str
    xml_exists: bool
    template_count: int
    template_paths: list[PathCheck]
    reaction_count: int
    output_directory: str
    force_field_name: str
    force_field_neutralize: bool
    force_field_use_smarts: bool

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["template_paths"] = [item.to_dict() for item in self.template_paths]
        return data


def load_config(path: str | Path) -> dict[str, Any]:
    return load_pipeline_config(path)


def summarize_config(path: str | Path) -> CG2AAConfigSummary:
    config_path = Path(path).expanduser().resolve()
    config = load_config(config_path)
    xml = resolve_path(config, config["xml"])
    templates = config["templates"]
    template_paths = [
        PathCheck(
            key=str(name),
            path=str(resolve_path(config, spec["pdb"])),
            exists=resolve_path(config, spec["pdb"]).is_file(),
        )
        for name, spec in templates.items()
    ]
    output = config["output"]
    force_field = config["force_field"]
    return CG2AAConfigSummary(
        path=str(config_path),
        schema=SCHEMA,
        xml=str(xml),
        xml_exists=xml.is_file(),
        template_count=len(templates),
        template_paths=template_paths,
        reaction_count=len(config.get("reactions", [])),
        output_directory=str(resolve_path(config, output["directory"])),
        force_field_name=str(force_field.get("name", "opls")),
        force_field_neutralize=bool(force_field.get("neutralize", True)),
        force_field_use_smarts=bool(force_field.get("use_smarts", False)),
    )


def validate_config(path: str | Path, *, check_forcefield: bool = True) -> dict[str, Any]:
    config_path = Path(path).expanduser().resolve()
    config = load_config(config_path)
    summary = summarize_config(config_path)
    checks: dict[str, bool] = {
        "config_exists": config_path.is_file(),
        "schema": summary.schema == SCHEMA,
        "has_xml": bool(summary.xml),
        "xml_exists": summary.xml_exists,
        "has_templates": summary.template_count > 0,
        "templates_exist": all(item.exists for item in summary.template_paths),
        "has_output_directory": bool(summary.output_directory),
        "has_placement": bool(config.get("placement", {}).get("sequence")),
    }
    forcefield_data: dict[str, Any] | None = None
    if check_forcefield:
        forcefield_data = validate_config_forcefield(config_path)
        checks["forcefield_valid"] = bool(forcefield_data.get("ok"))
    return {
        "ok": all(checks.values()),
        "config": str(config_path),
        "summary": summary.to_dict(),
        "checks": checks,
        "forcefield": forcefield_data,
    }
