from __future__ import annotations

import json
import re
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from multiscale.cg2aa.config import summarize_config, validate_config
from multiscale.cg2aa.internal import internal_module_status
from multiscale.cg2aa.pipeline_config import load_pipeline_config, resolved_pipeline_document
from multiscale.cg2aa.runner import run_cg2aa
from multiscale.cli import main


def _write_minimal_project(tmp_path: Path, *, cg_distance: float = 0.154) -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    (tmp_path / "carbon_template.pdb").write_text(
        "REMARK minimal carbon template\n"
        "ATOM      1  C   MOL A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "END\n",
        encoding="utf-8",
    )
    (tmp_path / "minimal_reaction_cg.xml").write_text(
        f"""<?xml version="1.0" encoding="UTF-8"?>
<galamost_xml version="1.3">
<configuration time_step="0" dimensions="3" natoms="2" >
<box lx="5.0" ly="5.0" lz="5.0"/>
<position num="2">
0.0 0.0 0.0
{cg_distance:.6f} 0.0 0.0
</position>
<type num="2">
M
M
</type>
<bond num="1">
M-M 0 1
</bond>
</configuration>
</galamost_xml>
""",
        encoding="utf-8",
    )
    config = {
        "schema": "multiscale.cg2aa.v1",
        "name": "Minimal methane-to-ethane backmapping",
        "input": {"cg_xml": "minimal_reaction_cg.xml"},
        "units": {
            "cg_xml": "nm",
            "template_pdb": "angstrom",
            "output_xml": "nm",
        },
        "templates": {
            "METHANE": {
                "cg_type": "M",
                "pdb": "carbon_template.pdb",
                "smiles": "C",
                "add_hydrogens": False,
                "sites": [{"name": "CENTER", "atoms": [1]}],
            }
        },
        "placement": [{"template": "METHANE", "count": 2}],
        "reactions": [
            {
                "id": "couple_carbons",
                "reactants": ["M", "M"],
                "bond_types": ["M-M"],
                "smarts": "[C:1].[C:2]>>[C:1][C:2]",
            }
        ],
        "hydrogens": "rebuild",
        "force_field": {"name": "opls", "neutralize": True},
        "output": {
            "directory": "output",
            "prefix": "ethane",
            "split": False,
        },
    }
    path = tmp_path / "minimal_reaction.json"
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return path


def _xml_rows(path: Path, tag: str) -> list[str]:
    root = ET.parse(path).getroot()
    element = root.find(f".//{tag}")
    assert element is not None
    return [line.strip() for line in (element.text or "").splitlines() if line.strip()]


def _bond_lengths_by_elements(path: Path) -> dict[tuple[str, str], list[float]]:
    positions = np.asarray([[float(value) for value in row.split()] for row in _xml_rows(path, "position")])
    elements = [label[0] for label in _xml_rows(path, "opls_type")]
    result: dict[tuple[str, str], list[float]] = {}
    for row in _xml_rows(path, "bond"):
        fields = row.split()
        i, j = int(fields[-2]), int(fields[-1])
        key = tuple(sorted((elements[i], elements[j])))
        result.setdefault(key, []).append(float(np.linalg.norm(positions[j] - positions[i])))
    return result


def test_internal_status_lists_packaged_modules():
    data = internal_module_status()
    assert data["ok"] is True
    assert data["release_status"] == "public-self-contained"
    assert data["checks"]["modules_inside_package"] is True


def test_formal_config_validates(tmp_path: Path):
    path = _write_minimal_project(tmp_path)
    summary = summarize_config(path)
    assert summary.schema == "multiscale.cg2aa.v1"
    assert summary.xml_exists is True
    assert all(item.exists for item in summary.template_paths)
    assert validate_config(path)["ok"] is True


def test_resolved_document_contains_only_formal_configuration(tmp_path: Path):
    path = _write_minimal_project(tmp_path)
    config = load_pipeline_config(path)
    resolved = resolved_pipeline_document(config)
    assert resolved["schema"] == "multiscale.cg2aa.resolved.v1"
    assert resolved["configuration"]["schema"] == "multiscale.cg2aa.v1"
    assert Path(resolved["configuration"]["input"]["cg_xml"]).is_absolute()
    assert isinstance(resolved["configuration"]["placement"], list)
    assert "pipeline" not in resolved
    assert "instances" not in resolved["configuration"]


def test_minimal_reaction_run(tmp_path: Path):
    config_path = _write_minimal_project(tmp_path)
    result = run_cg2aa(config_path)
    output = tmp_path / "output"
    assert result["reaction_events"] == 1
    assert result["geometry_ok"] is True
    assert (output / "ethane.pdb").is_file()
    assert "CRYST1" in (output / "ethane.pdb").read_text(encoding="utf-8")
    xml_path = output / "ethane.xml"
    inbox_path = output / "ethane_inbox.xml"
    xml_text = xml_path.read_text(encoding="utf-8")
    assert re.search(r'natoms="8"', xml_text)
    assert '<bond num="7">' in xml_text
    assert inbox_path.is_file()
    assert len(_xml_rows(inbox_path, "image")) == 8

    bond_lengths = _bond_lengths_by_elements(xml_path)
    assert 0.14 <= bond_lengths[("C", "C")][0] <= 0.17
    assert min(bond_lengths[("C", "H")]) >= 0.10
    assert max(bond_lengths[("C", "H")]) <= 0.115

    original = np.asarray([[float(value) for value in row.split()] for row in _xml_rows(xml_path, "position")])
    wrapped = np.asarray([[float(value) for value in row.split()] for row in _xml_rows(inbox_path, "position")])
    images = np.asarray([[int(value) for value in row.split()] for row in _xml_rows(inbox_path, "image")])
    box_element = ET.parse(inbox_path).getroot().find(".//box")
    assert box_element is not None
    box = np.asarray([float(box_element.attrib[key]) for key in ("lx", "ly", "lz")])
    assert np.all(wrapped >= -box / 2.0 - 1.0e-8)
    assert np.all(wrapped < box / 2.0 + 1.0e-8)
    assert np.max(np.abs(original - (wrapped + images * box))) < 2.0e-6

    geometry = json.loads((output / "cg2aa_geometry_report.json").read_text(encoding="utf-8"))
    assert geometry["ok"] is True
    assert geometry["checks"]["reaction_bond_lengths_reasonable"] is True
    assert geometry["checks"]["reaction_center_angles_reasonable"] is True
    assert geometry["checks"]["hydrogen_bond_lengths_reasonable"] is True
    assert geometry["checks"]["no_severe_interinstance_heavy_atom_clashes"] is True
    assert geometry["checks"]["no_severe_nonbonded_all_atom_clashes"] is True
    assert (output / "summary.json").is_file()
    resolved = json.loads((output / "cg2aa_resolved.json").read_text(encoding="utf-8"))
    assert resolved["schema"] == "multiscale.cg2aa.resolved.v1"
    assert resolved["configuration"]["hydrogens"] == "rebuild"


def test_reaction_geometry_relaxes_a_long_initial_cg_bond(tmp_path: Path):
    config_path = _write_minimal_project(tmp_path, cg_distance=0.80)
    result = run_cg2aa(config_path)
    assert result["geometry_ok"] is True

    output = tmp_path / "output"
    bond_lengths = _bond_lengths_by_elements(output / "ethane.xml")
    assert 0.14 <= bond_lengths[("C", "C")][0] <= 0.17

    report = json.loads((output / "cg2aa_geometry_report.json").read_text(encoding="utf-8"))
    assert report["refinement"]["before"]["reaction_bonds"]["max_nm"] > 0.50
    assert report["refinement"]["after"]["reaction_bonds"]["max_nm"] < 0.20


def test_unknown_and_removed_cg2aa_fields_are_rejected(tmp_path: Path):
    config_path = _write_minimal_project(tmp_path)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    config["instances"] = []
    config_path.write_text(json.dumps(config), encoding="utf-8")
    with pytest.raises(ValueError, match="unknown field.*instances"):
        load_pipeline_config(config_path)

    config.pop("instances")
    config["placement"][0]["test_count"] = 2
    config_path.write_text(json.dumps(config), encoding="utf-8")
    with pytest.raises(ValueError, match="unknown field.*test_count"):
        load_pipeline_config(config_path)


def test_nonformal_cg2aa_document_is_rejected(tmp_path: Path):
    path = tmp_path / "old.json"
    path.write_text(
        json.dumps(
            {
                "xml": "system.xml",
                "templates": {},
                "placement": {"sequence": []},
                "output": {"directory": "out"},
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="schema"):
        validate_config(path, check_forcefield=False)


def test_cli_summary_validate_and_run(tmp_path: Path, capsys):
    config_path = _write_minimal_project(tmp_path)
    assert main(["cg2aa", "summary", str(config_path), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["ok"] is True
    assert main(["cg2aa", "validate", str(config_path), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["ok"] is True
    assert main(["cg2aa", "run", str(config_path)]) == 0
    assert json.loads(capsys.readouterr().out)["ok"] is True


def test_removed_cg2aa_command_is_not_an_api(capsys):
    assert main(["cg2aa", "prepare-" + "config"]) == 1
    data = json.loads(capsys.readouterr().out)
    assert data["ok"] is False
