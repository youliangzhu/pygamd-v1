from __future__ import annotations

import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_release_metadata():
    metadata = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    project = metadata["project"]
    assert project["name"] == "multiscale-md"
    assert project["version"] == "0.1.0"
    assert project["requires-python"] == ">=3.11"
    assert project["license"] == "GPL-3.0-or-later"
    assert set(project["license-files"]) >= {"LICENSE", "NOTICE.md"}
    assert set(project["dependencies"]) == {"numpy>=1.26", "networkx>=3.0"}
    assert project["optional-dependencies"]["rdf"] == ["MDAnalysis>=2.8"]
    assert metadata["build-system"]["requires"][0] == "setuptools>=77"
    assert metadata["tool"]["setuptools"]["packages"]["find"]["where"] == ["."]


def _public_text() -> str:
    paths = [
        ROOT / "README.md",
        ROOT / "NOTICE.md",
        ROOT / "MANIFEST.in",
        ROOT / ".gitignore",
        ROOT / "pyproject.toml",
        *sorted((ROOT / "docs").glob("*.md")),
        *sorted((ROOT / "scripts").glob("*.sh")),
        *sorted((ROOT / "multiscale").rglob("*.py")),
    ]
    return "\n".join(path.read_text(encoding="utf-8", errors="ignore") for path in paths)


def test_public_material_describes_only_formal_json_workflows():
    text = _public_text()
    forbidden = [
        "AA2CG_" + "WORK" + "DIR",
        "prepare-" + "config",
        "prepare-" + "run",
        "magic_" + "cg" + "traj_parser",
        "cg" + "traj_",
        "source_" + "style",
        "MIGRA" + "TED_" + "MODULES",
        "multiscale_" + "smo" + "ke",
        "smo" + "ke-test",
        "compati" + "bility parser",
        "verbose configuration",
        "PLACE" + "HOLDER",
        "St" + "ep-",
        "Multiscale_CG",
        "Multiscale_Linux",
        "/home/gsh",
        "src/multiscale",
    ]
    for marker in forbidden:
        assert marker not in text
    assert ("m" + "dx") not in text.lower()
    assert "multiscale.aa2cg.v1" in text
    assert "multiscale.cg2aa.v1" in text
    assert "MagiC-3" in text
    assert "in-house" in text
    assert "TEST" + "!!!!" not in text
    assert "covered by " + "tests for a later backend swap" not in text
    assert "GROMACS OPLS-AA bonded parameters" in text
    assert "asteeves/rdlt" in text
    assert "authorized secondary developer" in text


def test_repository_layout_contains_only_current_entrypoints():
    assert (ROOT / "scripts" / "multiscale_validate.sh").is_file()
    assert not (ROOT / "scripts" / ("multiscale_" + "smo" + "ke.sh")).exists()
    for path in [
        ROOT / "multiscale" / "aa2cg" / "config.py",
        ROOT / "multiscale" / "aa2cg" / "discovery.py",
        ROOT / "multiscale" / "aa2cg" / ("magic_" + "cg" + "traj_parser.py"),
        ROOT / "multiscale" / "aa2cg" / "schema.py",
    ]:
        assert not path.exists()
    egg_info = {path.name for path in ROOT.glob("*.egg-info")}
    assert egg_info <= {"multiscale_md.egg-info"}
    assert not (ROOT / "examples").exists()
    assert not (ROOT / "multiscale" / "cg2aa" / "examples").exists()
