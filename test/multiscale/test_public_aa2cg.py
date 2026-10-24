from __future__ import annotations

import hashlib
import json
from pathlib import Path

from multiscale.aa2cg.exporters.galamost import export_galamost_auxiliary, export_galamost_topology
from multiscale.aa2cg.exporters.gromacs import export_gromacs_scaffold
from multiscale.aa2cg.exporters.lammps import export_lammps_scaffold
from multiscale.aa2cg.internal import aa2cg_internal_status
from multiscale.aa2cg.json_config import load_aa2cg_config, validate_aa2cg_config
from multiscale.aa2cg.mapper import map_structure_to_cgsystem, write_mapped_structure
from multiscale.cli import main
from multiscale.core.cg_system import CGAngle, CGBead, CGBond, CGSystem, write_cgsystem_json


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_project(tmp_path: Path, *, mapping_artifacts: bool = False, rdf: object = False) -> dict[str, Path]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    dump = tmp_path / "topology.dump"
    dump.write_text(
        """md.tpr:
   natoms = 8
topology:
   name=\"MOL\"
   #atoms                         = 8
   #molblock                      = 1
   molblock (0):
      moltype              = 0 \"MOL\"
      #molecules                     = 2
   moltype (0):
      name=\"MOL\"
      atoms:
         atom (4):
            atom[     0]={type=  0, typeB=  0, ptype=    Atom, m= 1.20000e+01, q= 1.00000e-01, mB= 1.20000e+01, qB= 1.00000e-01, resind=    0, atomnumber= -1}
            atom[     1]={type=  1, typeB=  1, ptype=    Atom, m= 1.00000e+00, q=-1.00000e-01, mB= 1.00000e+00, qB=-1.00000e-01, resind=    0, atomnumber= -1}
            atom[     2]={type=  0, typeB=  0, ptype=    Atom, m= 2.00000e+00, q= 2.00000e-01, mB= 2.00000e+00, qB= 2.00000e-01, resind=    0, atomnumber= -1}
            atom[     3]={type=  1, typeB=  1, ptype=    Atom, m= 2.00000e+00, q=-2.00000e-01, mB= 2.00000e+00, qB=-2.00000e-01, resind=    0, atomnumber= -1}
         atom (4):
            atom[0]={name=\"C1\"}
            atom[1]={name=\"H1\"}
            atom[2]={name=\"C2\"}
            atom[3]={name=\"H2\"}
box (3x3):
   box[    0]={ 4.00000e+00,  0.00000e+00,  0.00000e+00}
   box[    1]={ 0.00000e+00,  5.00000e+00,  0.00000e+00}
   box[    2]={ 0.00000e+00,  0.00000e+00,  6.00000e+00}
""",
        encoding="utf-8",
    )
    gro = tmp_path / "input.gro"
    gro.write_text(
        """AA coordinate frame
    8
    1MOL     C1    1   0.000   0.000   0.000
    1MOL     H1    2   1.300   0.000   0.000
    1MOL     C2    3   0.000   2.000   0.000
    1MOL     H2    4   0.000   4.000   0.000
    2MOL     C1    5   2.000   0.000   0.000
    2MOL     H1    6   3.300   0.000   0.000
    2MOL     C2    7   2.000   2.000   0.000
    2MOL     H2    8   2.000   4.000   0.000
   4.00000   5.00000   6.00000
""",
        encoding="utf-8",
    )
    tpr = tmp_path / "topology.tpr"
    tpr.write_bytes(b"MULTISCALE_TPR_INPUT\n")
    config = {
        "schema": "multiscale.aa2cg.v1",
        "name": "three-bead mapping",
        "input": {"tpr": tpr.name, "gro": gro.name},
        "model": {
            "name": "MOL_CG",
            "parent": "MOL",
            "type_order": ["A", "B", "C"],
            "beads": [
                {"name": "B1", "type": "A", "atoms": [1, 2]},
                {"name": "B2", "type": "B", "atoms": [3]},
                {"name": "B3", "type": "C", "atoms": [4]},
            ],
            "bonds": {"1": [["B1", "B2"], ["B2", "B3"]]},
            "angles": {"2": [["B1", "B2", "B3"]]},
        },
        "output": {
            "directory": "outputs/aa2cg",
            "formats": ["galamost", "gromacs", "lammps"],
            "mapping_artifacts": mapping_artifacts,
        },
        "rdf": rdf,
        "runtime": {"gmx_bin": str(tmp_path / "gmx")},
    }
    config_path = tmp_path / "aa2cg.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    gmx = tmp_path / "gmx"
    gmx.write_text("#!/usr/bin/env bash\ncat \"$MULTISCALE_TPR_DUMP\"\n", encoding="utf-8")
    gmx.chmod(0o755)
    return {"dump": dump, "gro": gro, "tpr": tpr, "config": config_path, "gmx": gmx}


def _export_system() -> CGSystem:
    return CGSystem(
        name="polymer",
        box=(5.0, 6.0, 7.0),
        box_unit="nm",
        coordinate_unit="nm",
        index_base=1,
        nmol=1,
        beads_per_molecule=3,
        beads=[
            CGBead(index=1, name="B1", bead_type="BB", molecule_index=1, molecule_bead_index=1, mass=12.0, charge=0.1, position=(1.0, 1.0, 1.0)),
            CGBead(index=2, name="B2", bead_type="BBE", molecule_index=1, molecule_bead_index=2, mass=13.0, charge=-0.1, position=(2.0, 2.0, 2.0)),
            CGBead(index=3, name="B3", bead_type="CO", molecule_index=1, molecule_bead_index=3, mass=14.0, charge=0.0, position=(3.0, 3.0, 3.0)),
        ],
        bonds=[CGBond(i=1, j=2, bond_type="polymer:1", molecule_index=1)],
        angles=[CGAngle(i=1, j=2, k=3, angle_type="polymer:5", molecule_index=1)],
        metadata={"bead_type_groups": {"BBE": ["B2"], "BB": ["B1"], "CO": ["B3"]}},
    )


def test_internal_status_lists_packaged_modules():
    data = aa2cg_internal_status()
    assert data["ok"] is True
    assert data["release_status"] == "public-self-contained"
    assert data["checks"]["modules_inside_package"] is True


def test_formal_config_maps_structure(tmp_path: Path):
    paths = _write_project(tmp_path)
    config = load_aa2cg_config(paths["config"])
    system, report = map_structure_to_cgsystem(
        topology_dump=paths["dump"],
        coordinates=config.coordinates_gro,
        mapping=config.mapping,
        cg_topology=config.topology,
    )
    assert report.ok
    assert system.summary()["natoms"] == 6
    assert system.summary()["nbonds"] == 4
    assert system.summary()["nangles"] == 2
    assert abs(system.beads[0].position[0] - 0.1) < 1e-12

    written = write_mapped_structure(
        topology_dump=paths["dump"],
        coordinates=config.coordinates_gro,
        mapping=config.mapping,
        cg_topology=config.topology,
        outdir=tmp_path / "details",
    )
    assert written.ok
    assert (tmp_path / "details" / "cg_system_mapped.json").is_file()


def test_run_uses_formal_config_and_hides_mapping_details_by_default(tmp_path: Path, capsys, monkeypatch):
    paths = _write_project(tmp_path)
    monkeypatch.setenv("MULTISCALE_TPR_DUMP", str(paths["dump"]))
    assert main(["aa2cg", "run", str(paths["config"]), "--gmx-bin", str(paths["gmx"]), "--json"]) == 0
    result = json.loads(capsys.readouterr().out)
    root = tmp_path / "outputs" / "aa2cg"
    assert result["ok"] is True
    assert result["mapping"]["artifacts"]["kept"] is False
    assert not (root / "mapped").exists()
    assert (root / "galamost" / "topology.xml").is_file()
    assert (root / "galamost" / "topology_inbox.xml").is_file()
    assert result["exports"]["galamost"]["topology_inbox"]["ok"] is True
    assert result["files"]["galamost_topology_inbox"].endswith("topology_inbox.xml")
    assert (root / "gromacs" / "cg.gro").is_file()
    assert (root / "lammps" / "system.data").is_file()
    resolved = json.loads((root / "aa2cg_resolved.json").read_text(encoding="utf-8"))
    assert resolved["schema"] == "multiscale.aa2cg.resolved.v1"
    assert resolved["configuration"]["schema"] == "multiscale.aa2cg.v1"
    assert resolved["configuration"]["output"]["mapping_artifacts"] is False


def test_mapping_details_do_not_change_final_exports(tmp_path: Path, capsys, monkeypatch):
    left = tmp_path / "left"
    right = tmp_path / "right"
    left.mkdir()
    right.mkdir()
    a = _write_project(left, mapping_artifacts=False)
    b = _write_project(right, mapping_artifacts=True)
    for paths in (a, b):
        monkeypatch.setenv("MULTISCALE_TPR_DUMP", str(paths["dump"]))
        assert main(["aa2cg", "run", str(paths["config"]), "--gmx-bin", str(paths["gmx"]), "--json"]) == 0
        capsys.readouterr()
    for rel in [
        "galamost/topology.xml",
        "galamost/topology_inbox.xml",
        "galamost/exclusionsSR.inc.py",
        "galamost/exclusionsEL.inc.py",
        "galamost/tables.inc.py",
        "gromacs/cg.gro",
        "gromacs/cg.itp",
        "gromacs/topol.top",
        "lammps/system.data",
        "lammps/in.cg",
    ]:
        assert _sha256(left / "outputs" / "aa2cg" / rel) == _sha256(right / "outputs" / "aa2cg" / rel)
    assert not (left / "outputs" / "aa2cg" / "mapped").exists()
    assert (right / "outputs" / "aa2cg" / "mapped" / "cg_system_mapped.json").is_file()


def test_rdf_requires_a_trajectory_only_when_enabled(tmp_path: Path):
    disabled = _write_project(tmp_path / "disabled")
    assert validate_aa2cg_config(disabled["config"])["ok"] is True

    enabled_root = tmp_path / "enabled"
    enabled_root.mkdir()
    enabled = _write_project(enabled_root, rdf={})
    data = validate_aa2cg_config(enabled["config"])
    assert data["ok"] is False
    assert "trajectory" in data["error"]


def test_nonformal_aa2cg_document_is_rejected(tmp_path: Path):
    path = tmp_path / "old.json"
    path.write_text(json.dumps({"input" + "s": {}, "mole" + "cule": {}, "topo" + "logy": {}}), encoding="utf-8")
    data = validate_aa2cg_config(path, require_files=False)
    assert data["ok"] is False
    assert "schema" in data["error"]


def test_unknown_aa2cg_field_is_rejected(tmp_path: Path):
    paths = _write_project(tmp_path)
    config = json.loads(paths["config"].read_text(encoding="utf-8"))
    config["output"]["test_only"] = True
    paths["config"].write_text(json.dumps(config), encoding="utf-8")
    data = validate_aa2cg_config(paths["config"], require_files=False)
    assert data["ok"] is False
    assert "unknown field" in data["error"]
    assert "test_only" in data["error"]


def test_removed_aa2cg_command_is_not_an_api(capsys):
    assert main(["aa2cg", "map"]) == 2
    assert "failed" in capsys.readouterr().err


def test_exporters_accept_retained_cgsystem(tmp_path: Path, capsys):
    source = tmp_path / "cgsystem.json"
    write_cgsystem_json(_export_system(), source)
    for command in ["export-galamost", "export-gromacs", "export-lammps"]:
        outdir = tmp_path / command
        assert main(["aa2cg", command, str(source), "--outdir", str(outdir), "--json"]) == 0
        assert json.loads(capsys.readouterr().out)["backend"] == "internal"
    assert (tmp_path / "export-galamost" / "topology.xml").is_file()
    assert (tmp_path / "export-galamost" / "topology_inbox.xml").is_file()
    assert (tmp_path / "export-gromacs" / "cg.gro").is_file()
    assert (tmp_path / "export-lammps" / "system.data").is_file()


def test_exporter_functions_write_files(tmp_path: Path):
    system = _export_system()
    assert export_galamost_topology(system, tmp_path / "gal").ok
    assert export_galamost_auxiliary(system, tmp_path / "gal").ok
    assert export_gromacs_scaffold(system, tmp_path / "gmx").ok
    assert export_lammps_scaffold(system, tmp_path / "lmp").ok


def test_galamost_inbox_wraps_coordinates_and_preserves_topology(tmp_path: Path):
    system = CGSystem(
        name="unwrapped",
        box=(10.0, 12.0, 14.0),
        box_unit="angstrom",
        coordinate_unit="angstrom",
        index_base=1,
        beads=[
            CGBead(index=1, name="A", bead_type="A", mass=12.0, charge=0.1, position=(16.0, -2.0, 5.0)),
            CGBead(index=2, name="B", bead_type="B", mass=14.0, charge=-0.1, position=(4.0, 19.0, 33.0)),
            CGBead(index=3, name="C", bead_type="C", mass=16.0, charge=0.0, position=(5.0, 6.0, 7.0)),
        ],
        bonds=[CGBond(i=1, j=2, bond_type="A-B:1")],
        angles=[CGAngle(i=1, j=2, k=3, angle_type="A-B-C:1")],
    )
    report = export_galamost_topology(system, tmp_path / "galamost")
    assert report.ok is True
    assert report.topology_ok is True
    inbox = report.topology_inbox
    assert inbox["ok"] is True
    assert inbox["counts"]["images"] == 3
    assert all(inbox["checks"].values())

    original = (tmp_path / "galamost" / "topology.xml").read_text(encoding="utf-8")
    wrapped = (tmp_path / "galamost" / "topology_inbox.xml").read_text(encoding="utf-8")
    assert "<image" not in original
    assert '<image num="3">' in wrapped
