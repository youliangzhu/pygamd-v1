from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np


@dataclass
class CgXml:
    timestep: int
    dimensions: int
    natoms: int
    box: tuple[float, float, float, float, float, float]
    positions: np.ndarray
    types: list[str]
    bonds: list[tuple[str, int, int]]
    sections: dict[str, list[list[str]]]

    def rows(self, section: str) -> list[list[str]]:
        return self.sections.get(section, [])


def _rows(node) -> list[list[str]]:
    if node is None or node.text is None:
        return []
    rows = []
    for line in node.text.splitlines():
        fields = line.split()
        if fields:
            rows.append(fields)
    return rows


def read_cg_xml(path: str | Path) -> CgXml:
    root = ET.parse(path).getroot()
    conf = root.find("configuration")
    if conf is None:
        raise ValueError(f"No <configuration> found in {path}")

    natoms = int(conf.attrib["natoms"])
    timestep = int(conf.attrib.get("time_step", 0))
    dimensions = int(conf.attrib.get("dimensions", 3))

    box_node = conf.find("box")
    if box_node is None:
        raise ValueError(f"No <box> found in {path}")
    box = (
        float(box_node.attrib.get("lx", 0.0)),
        float(box_node.attrib.get("ly", 0.0)),
        float(box_node.attrib.get("lz", 0.0)),
        float(box_node.attrib.get("xy", 0.0)),
        float(box_node.attrib.get("xz", 0.0)),
        float(box_node.attrib.get("yz", 0.0)),
    )

    sections = {child.tag: _rows(child) for child in conf}
    positions = np.asarray([[float(v) for v in row[:3]] for row in sections.get("position", [])], dtype=float)
    if positions.shape != (natoms, 3):
        raise ValueError(f"position count mismatch: expected {natoms}, got {len(positions)}")

    types = [str(row[0]) for row in sections.get("type", [])]
    if len(types) != natoms:
        raise ValueError(f"type count mismatch: expected {natoms}, got {len(types)}")

    bonds = []
    for row in sections.get("bond", []):
        if len(row) >= 3:
            bonds.append((str(row[0]), int(row[1]), int(row[2])))

    return CgXml(
        timestep=timestep,
        dimensions=dimensions,
        natoms=natoms,
        box=box,
        positions=positions,
        types=types,
        bonds=bonds,
        sections=sections,
    )
