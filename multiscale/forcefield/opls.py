"""OPLS force-field asset discovery, parsing, and config validation."""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from .models import ForceFieldAssets, OPLSAtomType, OPLSSummary

REQUIRED_DEFAULT_TYPES = {"C": "@atom:84", "H": "@atom:85", "N": "@atom:732", "O": "@atom:122"}


def packaged_opls_asset_dir() -> Path:
    """Return the bundled OPLS asset directory inside the Multiscale package."""
    return Path(__file__).resolve().parent / "assets" / "opls"


def find_opls_assets() -> ForceFieldAssets:
    """Discover bundled OPLS assets.

    OPLS assets are self-contained inside :mod:`multiscale.forcefield`.
    """
    asset_dir = packaged_opls_asset_dir()
    template = asset_dir / "STaGE_opls_tomoltemplate_opls.txt"
    bonded = asset_dir / "bonded.itp"
    missing: list[str] = []
    if not asset_dir.is_dir():
        missing.append("asset_dir")
    if not template.is_file():
        missing.append("STaGE_opls_tomoltemplate_opls.txt")
    if not bonded.is_file():
        missing.append("bonded.itp")
    return ForceFieldAssets("opls", "internal-multiscale", asset_dir, template, bonded, tuple(missing))


def parse_opls_template(path: str | Path) -> list[OPLSAtomType]:
    rows: list[OPLSAtomType] = []
    for raw in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("*") or "|" not in line:
            continue
        parts = [p.strip() for p in line.split("|")]
        if len(parts) < 7:
            continue
        element, atom_name, gromacs_type, smarts, moltemplate_type, charge, desc = parts[:7]
        try:
            charge_f = float(charge)
        except ValueError:
            charge_f = 0.0
        rows.append(OPLSAtomType(element, atom_name, gromacs_type, smarts, moltemplate_type, charge_f, desc.strip().strip('"')))
    return rows


def parse_bonded_section_counts(path: str | Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    current: str | None = None
    section_re = re.compile(r"^\s*\[\s*([^\]]+)\s*\]")
    for raw in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        m = section_re.match(line)
        if m:
            current = m.group(1).strip()
            counts.setdefault(current, 0)
            continue
        if current is None or line.startswith(";") or line.startswith("#"):
            continue
        counts[current] = counts.get(current, 0) + 1
    return counts


def summarize_opls() -> OPLSSummary:
    assets = find_opls_assets()
    if not assets.ok or assets.template_path is None or assets.bonded_path is None:
        return OPLSSummary(assets, 0, 0, (), {}, {k: False for k in REQUIRED_DEFAULT_TYPES})
    rows = parse_opls_template(assets.template_path)
    moltemplate_types = {r.moltemplate_type for r in rows}
    elements = tuple(sorted({r.element for r in rows if r.element}))
    default_presence = {el: typ in moltemplate_types for el, typ in REQUIRED_DEFAULT_TYPES.items()}
    return OPLSSummary(assets, len(rows), len(moltemplate_types), elements, parse_bonded_section_counts(assets.bonded_path), default_presence)


def atom_type_lookup() -> dict[str, dict[str, Any]]:
    summary = summarize_opls()
    if not summary.assets.template_path:
        return {}
    return {row.moltemplate_type: row.to_dict() for row in parse_opls_template(summary.assets.template_path)}


def validate_config_forcefield(config_path: str | Path) -> dict[str, Any]:
    config_path = Path(config_path).resolve()
    cfg = json.loads(config_path.read_text(encoding="utf-8"))
    ff_cfg = cfg.get("force_field", {})
    name = ff_cfg.get("name", "opls")
    defaults = ff_cfg.get("atom_defaults", ff_cfg.get("defaults", {})) or {}
    if name != "opls":
        return {"ok": False, "config": str(config_path), "force_field": ff_cfg, "error": f"unsupported force_field.name: {name}"}
    summary = summarize_opls()
    lookup = atom_type_lookup()
    default_checks = {symbol: value in lookup for symbol, value in defaults.items()}
    ok = summary.ok and all(default_checks.values())
    return {
        "ok": ok,
        "config": str(config_path),
        "force_field": {
            "name": name,
            "defaults": defaults,
            "neutralize": bool(ff_cfg.get("neutralize", True)),
            "use_smarts": bool(ff_cfg.get("use_smarts", False)),
        },
        "assets": summary.assets.to_dict(),
        "default_checks": default_checks,
        "opls_summary": summary.to_dict(),
    }
