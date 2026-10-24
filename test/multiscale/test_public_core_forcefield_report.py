from __future__ import annotations

import json
from pathlib import Path

from multiscale.cli import main
from multiscale.aa2cg.exporters.galamost import export_galamost_topology
from multiscale.core import CGBead, CGBond, CGSystem, compare_galamost_counts, read_galamost_xml, write_galamost_xml
from multiscale.core.report import REPORT_SCHEMA, make_report, validate_report
from multiscale.forcefield import assign_forcefield, internal_module_status, list_forcefields
from multiscale.forcefield.opls import parse_bonded_section_counts, parse_opls_template, summarize_opls


def _system() -> CGSystem:
    return CGSystem(
        name="tiny",
        box=(10.0, 10.0, 10.0),
        box_unit="angstrom",
        coordinate_unit="angstrom",
        beads=[
            CGBead(index=1, name="A", bead_type="A", mass=12.0, charge=0.0, position=(1.0, 2.0, 3.0)),
            CGBead(index=2, name="B", bead_type="B", mass=14.0, charge=0.0, position=(4.0, 5.0, 6.0)),
        ],
        bonds=[CGBond(i=1, j=2, bond_type="A-B:1")],
    )


def test_galamost_roundtrip_and_cli(tmp_path: Path, capsys):
    path = tmp_path / "topology.xml"
    write_galamost_xml(_system(), path)
    xml = read_galamost_xml(path)
    assert xml.counts()["natoms"] == 2
    assert xml.bonds == [("A-B:1", 0, 1)]
    assert compare_galamost_counts(path, path)["ok"] is True
    assert main(["core", "inspect-galamost", str(path), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["galamost"]["counts"]["natoms"] == 2


def test_report_helpers_and_cli(tmp_path: Path, capsys):
    report = make_report(module="core", command="inspect", payload={"ok": True, "counts": {"natoms": 2}})
    assert report["schema"] == REPORT_SCHEMA
    assert validate_report(report)["ok"] is True

    payload = tmp_path / "payload.json"
    out = tmp_path / "report.json"
    payload.write_text(json.dumps({"ok": True, "counts": {"n": 1}}), encoding="utf-8")
    assert main(["report", "wrap", str(payload), "--module", "unit", "--command", "test", "--out", str(out), "--json"]) == 0
    capsys.readouterr()
    assert out.exists()
    assert main(["report", "validate", str(out), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["ok"] is True


def test_forcefield_assets_and_cli(tmp_path: Path, capsys):
    template = tmp_path / "opls.txt"
    template.write_text('* comment\nC | CT | opls_139 | [C] | @atom:84 | 0.0 | "Alkane >C<"\n', encoding="utf-8")
    bonded = tmp_path / "bonded.itp"
    bonded.write_text('[ bondtypes ]\nC H 1 0.1 100\n[ angletypes ]\nC C H 1 109 50\n', encoding="utf-8")
    assert parse_opls_template(template)[0].moltemplate_type == "@atom:84"
    assert parse_bonded_section_counts(bonded) == {"bondtypes": 1, "angletypes": 1}

    listed = list_forcefields()
    assert listed[0]["name"] == "opls"
    assert listed[0]["source"] == "internal-multiscale"
    summary = summarize_opls().to_dict()
    assert summary["ok"] is True
    assert summary["assets"]["source"] == "internal-multiscale"
    assert callable(assign_forcefield)
    assert internal_module_status()["ok"] is True
    assert main(["forcefield", "list"]) == 0
    text_output = capsys.readouterr().out
    assert text_output.startswith("ok: True\n")
    assert "forcefields:" in text_output

    assert main(["forcefield", "list", "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["forcefields"][0]["source"] == "internal-multiscale"


def test_core_reader_understands_galamost_images(tmp_path: Path):
    report = export_galamost_topology(_system(), tmp_path)
    assert report.topology_inbox["ok"] is True
    inbox = read_galamost_xml(tmp_path / "topology_inbox.xml")
    assert inbox.counts()["images"] == inbox.natoms
    assert len(inbox.images) == inbox.natoms
    assert inbox.consistency_checks()["images_match_natoms_or_absent"] is True
