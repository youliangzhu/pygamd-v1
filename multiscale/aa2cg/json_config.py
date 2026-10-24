"""Formal ``multiscale.aa2cg.v1`` JSON configuration loading and validation."""

from __future__ import annotations

from dataclasses import dataclass
import importlib.util
import json
from pathlib import Path
import shutil
from typing import Any

from .model import AngleGroup, BeadMapping, BondGroup, MappingSpec, TopologySpec

SCHEMA = "multiscale.aa2cg.v1"
RESOLVED_SCHEMA = "multiscale.aa2cg.resolved.v1"


@dataclass(frozen=True)
class ExportTarget:
    enabled: bool
    directory: Path
    epsilon: float = 1.0


@dataclass(frozen=True)
class RDFRunConfig:
    enabled: bool
    trajectory: Path | None
    directory: Path
    start: int = 0
    stop: int | None = None
    stride: int = 1
    max_frames: int | None = None
    progress_every: int = 1
    progress_pairs: bool = False
    trajectory_engine: str = "mdanalysis"
    histogram_engine: str = "numpy"
    mode: str = "auto"
    pairs: str | list[list[str]] = "all"
    rmax_nm: float = 2.0
    bin_width_nm: float = 0.01
    exclude_bonded_pairs: bool = True
    exclude_angle_endpoints: bool = True


@dataclass(frozen=True)
class AA2CGConfig:
    path: Path
    name: str
    topology_tpr: Path
    coordinates_gro: Path
    mapping: MappingSpec
    topology: TopologySpec
    output_root: Path
    mapped_directory: Path
    keep_mapping_artifacts: bool
    galamost: ExportTarget
    gromacs: ExportTarget
    lammps: ExportTarget
    rdf: RDFRunConfig
    gmx_bin: str
    keep_topology_dump: bool

    def to_summary(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "config": str(self.path),
            "name": self.name,
            "inputs": {
                "topology_tpr": str(self.topology_tpr),
                "coordinates_gro": str(self.coordinates_gro),
                "trajectory": str(self.rdf.trajectory) if self.rdf.trajectory else None,
            },
            "mapping": {
                "molecule": self.mapping.name,
                "parent_type": self.mapping.parent_type,
                "beads_per_molecule": len(self.mapping.beads),
                "bead_types": list(self.topology.cgtypes),
                "bond_groups": len(self.topology.bond_groups),
                "angle_groups": len(self.topology.angle_groups),
            },
            "output": {
                "directory": str(self.output_root),
                "mapping_artifacts": {
                    "kept": self.keep_mapping_artifacts,
                    "directory": str(self.mapped_directory) if self.keep_mapping_artifacts else None,
                },
                "galamost": {"enabled": self.galamost.enabled, "directory": str(self.galamost.directory)},
                "gromacs": {"enabled": self.gromacs.enabled, "directory": str(self.gromacs.directory)},
                "lammps": {"enabled": self.lammps.enabled, "directory": str(self.lammps.directory)},
            },
            "rdf": {
                "enabled": self.rdf.enabled,
                "directory": str(self.rdf.directory),
                "trajectory_required": self.rdf.enabled,
                "trajectory": str(self.rdf.trajectory) if self.rdf.trajectory else None,
                "mode": self.rdf.mode,
                "rmax_nm": self.rdf.rmax_nm,
                "bin_width_nm": self.rdf.bin_width_nm,
            },
            "runtime": {
                "gmx_bin": self.gmx_bin,
                "gmx_available": _command_available(self.gmx_bin),
                "mdanalysis_available": _mdanalysis_available(),
                "topology_reader": (
                    "gmx_dump"
                    if _command_available(self.gmx_bin)
                    else ("mdanalysis" if _mdanalysis_available() else None)
                ),
                "keep_topology_dump": self.keep_topology_dump,
            },
            "input_contract": {
                "coordinates_used_as_is": True,
                "periodic_coordinate_preprocessing_performed": False,
                "user_responsible_for_coordinate_suitability": True,
            },
        }

    def to_resolved_dict(self, *, gmx_bin: str | None = None) -> dict[str, Any]:
        bead_type_by_name = {
            bead_name: bead_type
            for bead_type, bead_names in self.topology.cgtypes.items()
            for bead_name in bead_names
        }
        formats = [
            name
            for name, target in (
                ("galamost", self.galamost),
                ("gromacs", self.gromacs),
                ("lammps", self.lammps),
            )
            if target.enabled
        ]
        rdf: bool | dict[str, Any]
        if not self.rdf.enabled:
            rdf = False
        else:
            rdf = {
                "trajectory": str(self.rdf.trajectory),
                "directory": str(self.rdf.directory),
                "frames": {
                    "start": self.rdf.start,
                    "stop": self.rdf.stop,
                    "stride": self.rdf.stride,
                    "max": self.rdf.max_frames,
                },
                "range": {
                    "max_nm": self.rdf.rmax_nm,
                    "bin_nm": self.rdf.bin_width_nm,
                },
                "pairs": self.rdf.pairs,
                "mode": self.rdf.mode,
                "exclusions": {
                    "bonded_pairs": self.rdf.exclude_bonded_pairs,
                    "angle_endpoints": self.rdf.exclude_angle_endpoints,
                },
                "trajectory_engine": self.rdf.trajectory_engine,
                "histogram_engine": self.rdf.histogram_engine,
                "progress_every": self.rdf.progress_every,
                "progress_pairs": self.rdf.progress_pairs,
            }
        return {
            "schema": RESOLVED_SCHEMA,
            "resolved": True,
            "source_config": str(self.path),
            "configuration": {
                "schema": SCHEMA,
                "name": self.name,
                "input": {
                    "tpr": str(self.topology_tpr),
                    "gro": str(self.coordinates_gro),
                },
                "model": {
                    "name": self.mapping.name,
                    "parent": self.mapping.parent_type,
                    "type_order": list(self.topology.cgtypes),
                    "beads": [
                        {
                            "name": bead.name,
                            "type": bead_type_by_name[bead.name],
                            "atoms": list(bead.atom_indices_1based),
                        }
                        for bead in self.mapping.beads
                    ],
                    "bonds": {
                        group.type_id: [list(pair) for pair in group.pairs]
                        for group in self.topology.bond_groups
                    },
                    "angles": {
                        group.type_id: [list(triple) for triple in group.triples]
                        for group in self.topology.angle_groups
                    },
                },
                "output": {
                    "directory": str(self.output_root),
                    "mapping_artifacts": (
                        {"directory": str(self.mapped_directory)}
                        if self.keep_mapping_artifacts
                        else False
                    ),
                    "formats": formats,
                    "galamost": {
                        "directory": str(self.galamost.directory),
                        "epsilon": self.galamost.epsilon,
                    },
                    "gromacs": {"directory": str(self.gromacs.directory)},
                    "lammps": {"directory": str(self.lammps.directory)},
                },
                "rdf": rdf,
                "runtime": {
                    "gmx_bin": gmx_bin or self.gmx_bin,
                    "keep_topology_dump": self.keep_topology_dump,
                },
            },
            "input_contract": {
                "coordinates_used_as_is": True,
                "periodic_coordinate_preprocessing_performed": False,
                "user_responsible_for_coordinate_suitability": True,
            },
        }


def _object(value: Any, label: str) -> dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object")
    return value


def _reject_unknown_keys(data: dict[str, Any], allowed: set[str], label: str) -> None:
    unknown = sorted(set(data) - allowed)
    if unknown:
        fields = ", ".join(unknown)
        raise ValueError(f"{label} contains unknown field(s): {fields}")


def _string(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty string")
    return value.strip()


def _boolean(value: Any, label: str, *, default: bool) -> bool:
    if value is None:
        return default
    if not isinstance(value, bool):
        raise ValueError(f"{label} must be true or false")
    return value


def _integer(value: Any, label: str, *, default: int | None = None, minimum: int | None = None) -> int | None:
    if value is None:
        result = default
    elif isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{label} must be an integer or null")
    else:
        result = value
    if result is not None and minimum is not None and result < minimum:
        raise ValueError(f"{label} must be >= {minimum}")
    return result


def _number(value: Any, label: str, *, default: float, positive: bool = False) -> float:
    if value is None:
        result = float(default)
    elif isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be numeric")
    else:
        result = float(value)
    if positive and result <= 0.0:
        raise ValueError(f"{label} must be > 0")
    return result


def _resolve(base: Path, value: Any, label: str) -> Path:
    path = Path(_string(value, label)).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _output_path(output_root: Path, value: Any, default: str, label: str) -> Path:
    path = Path(default if value is None else _string(value, label)).expanduser()
    if not path.is_absolute():
        path = output_root / path
    return path.resolve()


def _command_available(command: str) -> bool:
    candidate = Path(command).expanduser()
    if candidate.is_absolute() or candidate.parent != Path("."):
        return candidate.is_file()
    return shutil.which(command) is not None


def _mdanalysis_available() -> bool:
    try:
        return importlib.util.find_spec("MDAnalysis") is not None
    except (ImportError, ValueError):
        return False


def _formats(value: Any) -> list[str]:
    if value is None:
        items: list[Any] = ["galamost"]
    elif isinstance(value, str):
        items = [value]
    elif isinstance(value, list):
        items = value
    else:
        raise ValueError("output.formats must be a string or string array")
    result: list[str] = []
    for index, item in enumerate(items):
        name = _string(item, f"output.formats[{index}]").lower()
        if name not in {"galamost", "gromacs", "lammps"}:
            raise ValueError("output.formats entries must be one of: galamost, gromacs, lammps")
        if name not in result:
            result.append(name)
    return result


def _parse_mapping(model: dict[str, Any], *, source: Path) -> tuple[MappingSpec, dict[str, list[str]]]:
    name = _string(model.get("name"), "model.name")
    parent = _string(model.get("parent", name), "model.parent")
    raw_beads = model.get("beads")
    if not isinstance(raw_beads, list) or not raw_beads:
        raise ValueError("model.beads must be a non-empty JSON array")

    beads: list[BeadMapping] = []
    type_groups: dict[str, list[str]] = {}
    bead_names: set[str] = set()
    used_atoms: set[int] = set()
    for index, raw_value in enumerate(raw_beads):
        bead = _object(raw_value, f"model.beads[{index}]")
        _reject_unknown_keys(bead, {"name", "type", "atoms"}, f"model.beads[{index}]")
        bead_name = _string(bead.get("name"), f"model.beads[{index}].name")
        bead_type = _string(bead.get("type"), f"model.beads[{index}].type")
        if bead_name in bead_names:
            raise ValueError(f"duplicate bead name: {bead_name}")
        bead_names.add(bead_name)
        raw_atoms = bead.get("atoms")
        if not isinstance(raw_atoms, list) or not raw_atoms:
            raise ValueError(f"model.beads[{index}].atoms must be a non-empty integer array")
        atoms: list[int] = []
        for atom_index, value in enumerate(raw_atoms):
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"model.beads[{index}].atoms[{atom_index}] must be a positive 1-based integer")
            atoms.append(value)
        if len(atoms) != len(set(atoms)):
            raise ValueError(f"model.beads[{index}].atoms contains duplicate indices")
        overlap = sorted(set(atoms) & used_atoms)
        if overlap:
            raise ValueError(f"AA atom indices are assigned to multiple beads: {overlap}")
        used_atoms.update(atoms)
        beads.append(BeadMapping(name=bead_name, atom_indices_1based=atoms))
        type_groups.setdefault(bead_type, []).append(bead_name)

    raw_order = model.get("type_order")
    if raw_order is not None:
        if not isinstance(raw_order, list) or not raw_order:
            raise ValueError("model.type_order must be a non-empty string array")
        order = [_string(item, f"model.type_order[{index}]") for index, item in enumerate(raw_order)]
        if len(order) != len(set(order)):
            raise ValueError("model.type_order contains duplicate types")
        if set(order) != set(type_groups):
            missing = sorted(set(type_groups) - set(order))
            extra = sorted(set(order) - set(type_groups))
            raise ValueError(f"model.type_order must list every bead type exactly once; missing={missing}, extra={extra}")
        type_groups = {bead_type: type_groups[bead_type] for bead_type in order}

    return MappingSpec(name=name, parent_type=parent, beads=beads, source=str(source)), type_groups


def _parse_topology(model: dict[str, Any], *, mapping: MappingSpec, type_groups: dict[str, list[str]], source: Path) -> TopologySpec:
    known = {bead.name for bead in mapping.beads}
    raw_bonds = model.get("bonds", {})
    raw_angles = model.get("angles", {})
    if not isinstance(raw_bonds, dict):
        raise ValueError("model.bonds must be a JSON object keyed by bond type")
    if not isinstance(raw_angles, dict):
        raise ValueError("model.angles must be a JSON object keyed by angle type")

    bond_groups: list[BondGroup] = []
    for type_id, rows in raw_bonds.items():
        if not isinstance(rows, list):
            raise ValueError(f"model.bonds.{type_id} must be an array")
        pairs: list[tuple[str, str]] = []
        for index, row in enumerate(rows):
            if not isinstance(row, list) or len(row) != 2:
                raise ValueError(f"model.bonds.{type_id}[{index}] must contain two bead names")
            left = _string(row[0], f"model.bonds.{type_id}[{index}][0]")
            right = _string(row[1], f"model.bonds.{type_id}[{index}][1]")
            unknown = [item for item in (left, right) if item not in known]
            if unknown:
                raise ValueError(f"model bond refers to unknown bead(s): {', '.join(unknown)}")
            pairs.append((left, right))
        bond_groups.append(BondGroup(mapping.name, str(type_id), pairs))

    angle_groups: list[AngleGroup] = []
    for type_id, rows in raw_angles.items():
        if not isinstance(rows, list):
            raise ValueError(f"model.angles.{type_id} must be an array")
        triples: list[tuple[str, str, str]] = []
        for index, row in enumerate(rows):
            if not isinstance(row, list) or len(row) != 3:
                raise ValueError(f"model.angles.{type_id}[{index}] must contain three bead names")
            names = tuple(_string(item, f"model.angles.{type_id}[{index}]") for item in row)
            unknown = [item for item in names if item not in known]
            if unknown:
                raise ValueError(f"model angle refers to unknown bead(s): {', '.join(unknown)}")
            triples.append((names[0], names[1], names[2]))
        angle_groups.append(AngleGroup(mapping.name, str(type_id), triples))

    return TopologySpec(type_groups, bond_groups, angle_groups, source=str(source))


def _parse_export_target(
    options: Any,
    *,
    enabled: bool,
    output_root: Path,
    default_directory: str,
    label: str,
    supports_epsilon: bool = False,
) -> ExportTarget:
    data = _object(options, label)
    allowed = {"directory", "epsilon"} if supports_epsilon else {"directory"}
    _reject_unknown_keys(data, allowed, label)
    directory = _output_path(output_root, data.get("directory"), default_directory, f"{label}.directory")
    epsilon = _number(data.get("epsilon"), f"{label}.epsilon", default=1.0, positive=True) if supports_epsilon else 1.0
    return ExportTarget(enabled=enabled, directory=directory, epsilon=epsilon)


def _parse_mapping_artifacts(output: dict[str, Any], output_root: Path) -> tuple[bool, Path]:
    value = output.get("mapping_artifacts", False)
    if isinstance(value, bool):
        return value, (output_root / "mapped").resolve()
    data = _object(value, "output.mapping_artifacts")
    _reject_unknown_keys(data, {"enabled", "directory"}, "output.mapping_artifacts")
    enabled = _boolean(data.get("enabled"), "output.mapping_artifacts.enabled", default=True)
    directory = _output_path(output_root, data.get("directory"), "mapped", "output.mapping_artifacts.directory")
    return enabled, directory


def _parse_rdf(value: Any, *, base: Path, output_root: Path, type_groups: dict[str, list[str]]) -> RDFRunConfig:
    if value in (None, False):
        return RDFRunConfig(enabled=False, trajectory=None, directory=(output_root / "rdf").resolve())
    if value is True:
        raise ValueError("rdf must be false or a JSON object containing trajectory")
    data = _object(value, "rdf")
    _reject_unknown_keys(
        data,
        {"trajectory", "directory", "frames", "range", "pairs", "mode", "exclusions", "trajectory_engine", "histogram_engine", "progress_every", "progress_pairs"},
        "rdf",
    )
    trajectory = _resolve(base, data.get("trajectory"), "rdf.trajectory")
    directory = _output_path(output_root, data.get("directory"), "rdf", "rdf.directory")

    frames = _object(data.get("frames"), "rdf.frames")
    _reject_unknown_keys(frames, {"start", "stop", "stride", "max"}, "rdf.frames")
    start = _integer(frames.get("start"), "rdf.frames.start", default=0, minimum=0)
    stop = _integer(frames.get("stop"), "rdf.frames.stop", default=None, minimum=0)
    stride = _integer(frames.get("stride"), "rdf.frames.stride", default=1, minimum=1)
    max_frames = _integer(frames.get("max"), "rdf.frames.max", default=None, minimum=1)
    assert start is not None and stride is not None
    if stop is not None and stop < start:
        raise ValueError("rdf.frames.stop must be greater than or equal to rdf.frames.start")

    range_data = _object(data.get("range"), "rdf.range")
    _reject_unknown_keys(range_data, {"max_nm", "bin_nm"}, "rdf.range")
    rmax_nm = _number(range_data.get("max_nm"), "rdf.range.max_nm", default=2.0, positive=True)
    bin_width_nm = _number(range_data.get("bin_nm"), "rdf.range.bin_nm", default=0.01, positive=True)
    if bin_width_nm >= rmax_nm:
        raise ValueError("rdf.range.bin_nm must be smaller than rdf.range.max_nm")

    pairs = data.get("pairs", "all")
    if pairs != "all":
        if not isinstance(pairs, list):
            raise ValueError("rdf.pairs must be 'all' or an array of [type1, type2] pairs")
        normalized: list[list[str]] = []
        known_types = set(type_groups)
        for index, row in enumerate(pairs):
            if not isinstance(row, list) or len(row) != 2:
                raise ValueError(f"rdf.pairs[{index}] must contain two bead type names")
            left = _string(row[0], f"rdf.pairs[{index}][0]")
            right = _string(row[1], f"rdf.pairs[{index}][1]")
            unknown = [item for item in (left, right) if item not in known_types]
            if unknown:
                raise ValueError(f"rdf pair refers to unknown bead type(s): {', '.join(unknown)}")
            normalized.append([left, right])
        pairs = normalized

    mode = _string(data.get("mode", "auto"), "rdf.mode").lower().replace("-", "_")
    if mode not in {"auto", "topology", "intermolecular", "all"}:
        raise ValueError("rdf.mode must be one of: auto, topology, intermolecular, all")
    trajectory_engine = _string(data.get("trajectory_engine", "mdanalysis"), "rdf.trajectory_engine").lower()
    if trajectory_engine != "mdanalysis":
        raise ValueError("rdf.trajectory_engine currently supports only 'mdanalysis'")
    histogram_engine = _string(data.get("histogram_engine", "numpy"), "rdf.histogram_engine").lower()
    histogram_engine = histogram_engine.replace("-", "").replace("_", "")
    if histogram_engine not in {"numpy", "python", "kdtree"}:
        raise ValueError("rdf.histogram_engine must be 'numpy', 'python', or 'kdtree'")

    exclusions = _object(data.get("exclusions"), "rdf.exclusions")
    _reject_unknown_keys(exclusions, {"bonded_pairs", "angle_endpoints"}, "rdf.exclusions")
    progress_every = _integer(data.get("progress_every"), "rdf.progress_every", default=1, minimum=0)
    assert progress_every is not None
    return RDFRunConfig(
        enabled=True,
        trajectory=trajectory,
        directory=directory,
        start=start,
        stop=stop,
        stride=stride,
        max_frames=max_frames,
        progress_every=progress_every,
        progress_pairs=_boolean(data.get("progress_pairs"), "rdf.progress_pairs", default=False),
        trajectory_engine=trajectory_engine,
        histogram_engine=histogram_engine,
        mode=mode,
        pairs=pairs,
        rmax_nm=rmax_nm,
        bin_width_nm=bin_width_nm,
        exclude_bonded_pairs=_boolean(exclusions.get("bonded_pairs"), "rdf.exclusions.bonded_pairs", default=True),
        exclude_angle_endpoints=_boolean(exclusions.get("angle_endpoints"), "rdf.exclusions.angle_endpoints", default=True),
    )


def load_aa2cg_config(path: str | Path, *, require_files: bool = True) -> AA2CGConfig:
    config_path = Path(path).expanduser().resolve()
    if not config_path.is_file():
        raise ValueError(f"AA-to-CG JSON configuration not found: {config_path}")
    try:
        raw = json.loads(config_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in {config_path}: {exc}") from exc
    if not isinstance(raw, dict):
        raise ValueError("AA-to-CG configuration root must be a JSON object")
    if raw.get("schema") != SCHEMA:
        raise ValueError(f"schema must be {SCHEMA!r}, got {raw.get('schema')!r}")
    _reject_unknown_keys(raw, {"schema", "name", "input", "model", "output", "rdf", "runtime"}, "configuration")

    base = config_path.parent
    input_data = _object(raw.get("input"), "input")
    _reject_unknown_keys(input_data, {"tpr", "gro"}, "input")
    topology_tpr = _resolve(base, input_data.get("tpr"), "input.tpr")
    coordinates_gro = _resolve(base, input_data.get("gro"), "input.gro")

    model = _object(raw.get("model"), "model")
    _reject_unknown_keys(model, {"name", "parent", "type_order", "beads", "bonds", "angles"}, "model")
    mapping, type_groups = _parse_mapping(model, source=config_path)
    topology = _parse_topology(model, mapping=mapping, type_groups=type_groups, source=config_path)
    name = _string(raw.get("name", mapping.name), "name")

    output = _object(raw.get("output"), "output")
    _reject_unknown_keys(output, {"directory", "formats", "mapping_artifacts", "galamost", "gromacs", "lammps"}, "output")
    output_root = _resolve(base, output.get("directory", "outputs/aa2cg"), "output.directory")
    formats = _formats(output.get("formats"))
    keep_mapping_artifacts, mapped_directory = _parse_mapping_artifacts(output, output_root)
    galamost = _parse_export_target(
        output.get("galamost"),
        enabled="galamost" in formats,
        output_root=output_root,
        default_directory="galamost",
        label="output.galamost",
        supports_epsilon=True,
    )
    gromacs = _parse_export_target(
        output.get("gromacs"),
        enabled="gromacs" in formats,
        output_root=output_root,
        default_directory="gromacs",
        label="output.gromacs",
    )
    lammps = _parse_export_target(
        output.get("lammps"),
        enabled="lammps" in formats,
        output_root=output_root,
        default_directory="lammps",
        label="output.lammps",
    )

    rdf = _parse_rdf(raw.get("rdf", False), base=base, output_root=output_root, type_groups=type_groups)
    runtime = _object(raw.get("runtime"), "runtime")
    _reject_unknown_keys(runtime, {"gmx_bin", "keep_topology_dump"}, "runtime")
    gmx_bin = _string(runtime.get("gmx_bin", "gmx"), "runtime.gmx_bin")
    keep_topology_dump = _boolean(runtime.get("keep_topology_dump"), "runtime.keep_topology_dump", default=False)

    if require_files:
        missing = [item for item in (topology_tpr, coordinates_gro) if not item.is_file()]
        if rdf.enabled and rdf.trajectory is not None and not rdf.trajectory.is_file():
            missing.append(rdf.trajectory)
        if missing:
            joined = "\n".join(f"  - {item}" for item in missing)
            raise ValueError(f"AA-to-CG input file(s) not found:\n{joined}")

    return AA2CGConfig(
        path=config_path,
        name=name,
        topology_tpr=topology_tpr,
        coordinates_gro=coordinates_gro,
        mapping=mapping,
        topology=topology,
        output_root=output_root,
        mapped_directory=mapped_directory,
        keep_mapping_artifacts=keep_mapping_artifacts,
        galamost=galamost,
        gromacs=gromacs,
        lammps=lammps,
        rdf=rdf,
        gmx_bin=gmx_bin,
        keep_topology_dump=keep_topology_dump,
    )


def validate_aa2cg_config(path: str | Path, *, require_files: bool = True) -> dict[str, Any]:
    try:
        config = load_aa2cg_config(path, require_files=require_files)
    except Exception as exc:
        return {
            "ok": False,
            "config": str(Path(path).expanduser().resolve()),
            "error_type": type(exc).__name__,
            "error": str(exc),
        }
    gmx_available = _command_available(config.gmx_bin) if require_files else True
    mdanalysis_available = _mdanalysis_available() if require_files else True
    checks = {
        "schema": True,
        "topology_tpr_exists": config.topology_tpr.is_file() if require_files else True,
        "coordinates_gro_exists": config.coordinates_gro.is_file() if require_files else True,
        "has_beads": bool(config.mapping.beads),
        "has_bead_types": bool(config.topology.cgtypes),
        "topology_references_known_beads": True,
        "trajectory_present_when_rdf_enabled": (
            (not config.rdf.enabled or bool(config.rdf.trajectory and config.rdf.trajectory.is_file()))
            if require_files
            else True
        ),
        # These two are diagnostic capabilities.  Either topology reader is
        # sufficient; gmx itself is no longer mandatory.
        "gmx_available": gmx_available,
        "mdanalysis_available": mdanalysis_available,
        "topology_reader_available": gmx_available or mdanalysis_available,
        "mdanalysis_available_when_rdf_enabled": (not config.rdf.enabled) or mdanalysis_available,
    }
    required_checks = {
        key: value
        for key, value in checks.items()
        if key not in {"gmx_available", "mdanalysis_available"}
    }
    warnings: list[str] = []
    if not gmx_available and mdanalysis_available:
        warnings.append(
            f"GROMACS executable is not available: {config.gmx_bin}; binary TPR will be read with MDAnalysis"
        )
    elif not gmx_available:
        warnings.append(
            f"Neither GROMACS executable {config.gmx_bin!r} nor MDAnalysis is available to read the TPR"
        )
    return {
        "ok": all(required_checks.values()),
        "summary": config.to_summary(),
        "checks": checks,
        "warnings": warnings,
    }
