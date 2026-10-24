from __future__ import annotations

import json
from pathlib import Path

import pytest

from multiscale import doctor as doctor_mod
from multiscale.audit.component import build_component_matrix, scan_external_dependency_references, self_contained_status
from multiscale.clean import clean_cache, collect_cleanable_paths
from multiscale.cli import main


def test_version_json(capsys):
    assert main(["version", "--json"]) == 0
    assert json.loads(capsys.readouterr().out) == {"ok": True, "version": "0.1.0", "release": "public"}


def test_help_lists_formal_workflows_only(capsys):
    assert main(["--help"]) == 0
    out = capsys.readouterr().out
    assert "multiscale aa2cg run <aa2cg.json>" in out
    assert "multiscale cg2aa run <cg2aa.json>" in out
    assert "--work" + "dir" not in out
    assert "prepare-" + "config" not in out


def test_audit_status_and_matrix():
    refs = scan_external_dependency_references()
    assert refs["dependency_reference_count"] == 0
    status = self_contained_status()
    assert status["self_contained_complete"] is True
    assert status["release_status"] == "public-self-contained"
    matrix = build_component_matrix()
    assert matrix["schema"] == "multiscale.component_matrix.v1"
    aa = next(item for item in matrix["components"] if item["component"] == "aa2cg")
    cg = next(item for item in matrix["components"] if item["component"] == "cg2aa")
    assert "multiscale aa2cg map" not in aa["commands"]
    assert "multiscale cg2aa prepare-" + "config" not in cg["commands"]


def test_audit_cli_json(capsys):
    assert main(["audit", "status", "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["self_contained_complete"] is True


def test_clean_cache_removes_runtime_noise(tmp_path: Path):
    (tmp_path / "pkg" / "__pycache__").mkdir(parents=True)
    (tmp_path / "pkg" / "__pycache__" / "mod.cpython-311.pyc").write_bytes(b"x")
    (tmp_path / ".pytest_cache").mkdir()
    (tmp_path / "multiscale_md.egg-info").mkdir()
    (tmp_path / "run.log").write_text("log\n", encoding="utf-8")
    (tmp_path / "multiscale_validation_run_001").mkdir()
    (tmp_path / "build").mkdir()

    dry_names = {path.name for path in collect_cleanable_paths(tmp_path)}
    assert {"__pycache__", ".pytest_cache", "multiscale_md.egg-info", "run.log", "multiscale_validation_run_001"} <= dry_names
    assert "build" not in dry_names

    result = clean_cache(tmp_path)
    assert result.ok is True
    assert len(result.removed) == 5
    assert not (tmp_path / "pkg" / "__pycache__").exists()
    assert not (tmp_path / "multiscale_validation_run_001").exists()
    assert (tmp_path / "build").exists()


def test_unknown_command_is_rejected():
    with pytest.raises(SystemExit):
        main(["forcefield", "external-module", "opls", "--json"])


def test_doctor_supports_python311_and_optional_gromacs(monkeypatch):
    monkeypatch.setattr(doctor_mod, "_gmx_info", lambda gmx_bin: {"ok": False, "required": False, "error": "not found"})
    data = doctor_mod.collect_doctor()
    assert data["python"]["requires_python"] == ">=3.11"
    assert data["checks"]["python_supported"] is True
    assert data["optional_tools"]["gromacs_available"] is False
    assert set(data["required_packages"]) == {"numpy", "networkx"}
    assert data["ok"] is True
