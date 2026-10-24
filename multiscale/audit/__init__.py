"""Self-contained source/runtime audit helpers for Multiscale."""

from .component import (
    COMPONENTS,
    build_component_matrix,
    classify_source_file,
    matrix_to_markdown,
    scan_external_dependency_references,
    self_contained_status,
    source_file_inventory,
    source_project_summary,
    write_json,
    write_markdown_matrix,
)

__all__ = [
    "COMPONENTS",
    "build_component_matrix",
    "classify_source_file",
    "matrix_to_markdown",
    "scan_external_dependency_references",
    "self_contained_status",
    "source_file_inventory",
    "source_project_summary",
    "write_json",
    "write_markdown_matrix",
]
