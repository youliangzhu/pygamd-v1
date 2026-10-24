from __future__ import annotations

import json
from pathlib import Path

from multiscale.aa2cg.model import TopologySpec
from multiscale.aa2cg.rdf import (
    CGFrame,
    RDFPair,
    RDFSpec,
    compute_rdf_histograms,
    make_rdf_spec_from_topology,
    resolve_rdf_mode_from_spec,
    write_rdf_outputs,
)
from multiscale.cli import main


def _spec() -> RDFSpec:
    return RDFSpec(
        cgtypes={"A": ["B1"], "B": ["B2"]},
        pairs=[RDFPair("A", "A"), RDFPair("A", "B")],
        rmax_input=10.0,
        dr_input=5.0,
        input_distance_unit="angstrom",
        output_distance_unit="nm",
        rmax_nm=1.0,
        dr_nm=0.5,
        source="sample",
        topology_exclusions=frozenset({("B1", "B2")}),
        has_bond_angle_topology=True,
        has_topology_exclusion_controls=True,
    )


def test_resolve_and_compute_rdf_histograms():
    spec = _spec()
    mode, source, warnings = resolve_rdf_mode_from_spec(spec)
    assert mode == "topology"
    assert source == "configuration"
    assert warnings == []
    data = compute_rdf_histograms(
        [CGFrame(positions_nm=[(0.0, 0.0, 0.0), (0.1, 0.0, 0.0), (0.9, 0.0, 0.0)], box_nm=(3.0, 3.0, 3.0))],
        bead_types=["A", "B", "A"],
        molecule_ids=[1, 1, 2],
        bead_names=["B1", "B2", "B1"],
        spec=spec,
        mode="topology",
        rdf_engine="python",
    )
    assert data["n_frames"] == 1
    assert data["topology_exclusion_count"] == 1
    assert data["histograms"]["A-A"] == [0, 1]


def test_write_rdf_outputs_are_public(tmp_path: Path):
    rdf_data = {
        "radii_nm": [0.25, 0.75],
        "histograms": {"A-A": [0, 1], "A-B": [1, 0]},
        "g_r": {"A-A": [0.0, 1.0], "A-B": [2.0, 0.0]},
        "n_frames": 1,
        "box_volume_average": 27.0,
        "bead_counts_per_type": {"A": 2, "B": 1},
        "rdf_engine": "python",
        "topology_exclusion_count": 1,
        "topology_exclusion_pairs": [["B1", "B2"]],
    }
    result = write_rdf_outputs(
        outdir=tmp_path / "rdf",
        spec=_spec(),
        rdf_data=rdf_data,
        mode="topology",
        source_trajectory="traj.xtc",
        inputs={"rdf_mode_requested": "auto", "rdf_mode_effective": "topology", "rdf_mode_source": "configuration", "warnings": []},
    )
    assert result.report["entrypoint"] == "rdf"
    header = (tmp_path / "rdf" / "AA_rdf_A-A.dat").read_text(encoding="utf-8")
    assert "# Multiscale RDF" in header
    assert ("multiscale" + "_cg") not in header
    assert json.loads((tmp_path / "rdf" / "rdf_report.json").read_text(encoding="utf-8"))["entrypoint"] == "rdf"


def test_cg2aa_help_describes_formal_schema(capsys):
    try:
        main(["cg2aa", "--help"])
    except SystemExit as exc:
        assert exc.code == 0
    out = capsys.readouterr().out
    assert "multiscale.cg2aa.v1" in out
    assert "multiscale.cg2aa.v1" in out


def _topology() -> TopologySpec:
    return TopologySpec(
        cgtypes={"A": ["B1"], "B": ["B2"], "C": ["B3"]},
        bond_groups=[],
        angle_groups=[],
        source="test",
    )


def test_make_rdf_spec_normalizes_and_deduplicates_explicit_pairs():
    spec = make_rdf_spec_from_topology(
        _topology(),
        pairs=[["B", "A"], ["A", "B"], ["C", "C"], ["C", "A"]],
    )

    assert spec.pairs == [RDFPair("A", "B"), RDFPair("C", "C"), RDFPair("A", "C")]


def test_make_rdf_spec_rejects_unknown_explicit_pair_type():
    try:
        make_rdf_spec_from_topology(_topology(), pairs=[["A", "UNKNOWN"]])
    except ValueError as exc:
        assert str(exc) == "RDF pair refers to unknown bead type(s): UNKNOWN"
    else:
        raise AssertionError("unknown RDF pair type was accepted")
