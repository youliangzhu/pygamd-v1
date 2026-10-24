"""Shared core namespace for Multiscale.

Shared data/IO utilities for aa2cg/cg2aa
algorithm delegation unchanged.
"""

from .runtime import RuntimeLayout, discover_runtime_layout
from .report import REPORT_SCHEMA, make_report, validate_report
from .cg_system import CGAngle, CGBead, CGBond, CGSystem, load_cgsystem_json, write_cgsystem_json
from .galamost_xml import (
    GalamostXml,
    compare_galamost_counts,
    inspect_galamost_xml,
    read_galamost_xml,
    render_galamost_xml,
    write_galamost_xml,
)

__all__ = [
    "RuntimeLayout",
    "discover_runtime_layout",
    "CGAngle",
    "CGBead",
    "CGBond",
    "CGSystem",
    "load_cgsystem_json",
    "write_cgsystem_json",
    "GalamostXml",
    "read_galamost_xml",
    "inspect_galamost_xml",
    "compare_galamost_counts",
    "render_galamost_xml",
    "write_galamost_xml",
    "REPORT_SCHEMA",
    "make_report",
    "validate_report",
]
