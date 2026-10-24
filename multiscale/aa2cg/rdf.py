"""Internal aa2cg RDF analysis implementation."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import math
import time
from collections.abc import Callable, Iterable, Sequence

from multiscale.aa2cg.model import MappingSpec, TopologySpec
from multiscale.aa2cg.gromacs import parse_gro
from multiscale.aa2cg.tpr_dump import TprDumpAtom, TprTopologyInput, load_tpr_topology



@dataclass(frozen=True)
class RDFPair:
    left: str
    right: str

    @property
    def name(self) -> str:
        return f"{self.left}-{self.right}"


def _unique_pairs(
    pairs: Iterable[RDFPair],
    cgtype_order: Sequence[str],
) -> list[RDFPair]:
    """Normalize symmetric RDF pairs and remove duplicates."""
    type_index = {name: index for index, name in enumerate(cgtype_order)}
    unique: list[RDFPair] = []
    seen: set[tuple[str, str]] = set()

    for pair in pairs:
        unknown = [name for name in (pair.left, pair.right) if name not in type_index]
        if unknown:
            names = ", ".join(dict.fromkeys(unknown))
            raise ValueError(f"RDF pair refers to unknown bead type(s): {names}")

        if type_index[pair.left] <= type_index[pair.right]:
            normalized = RDFPair(pair.left, pair.right)
        else:
            normalized = RDFPair(pair.right, pair.left)

        key = (normalized.left, normalized.right)
        if key not in seen:
            seen.add(key)
            unique.append(normalized)

    return unique


@dataclass(frozen=True)
class RDFSpec:
    cgtypes: dict[str, list[str]]
    pairs: list[RDFPair]
    rmax_input: float
    dr_input: float
    input_distance_unit: str
    output_distance_unit: str
    rmax_nm: float
    dr_nm: float
    source: str
    n_pair_bonds_exclude: int = 1
    n_angle_bonds_exclude: int = 1
    topology_exclusions: frozenset[tuple[str, str]] = frozenset()
    topology_exclusion_source: str = "configuration"
    has_bond_angle_topology: bool = False
    has_topology_exclusion_controls: bool = False

    def infer_default_mode(self) -> str:
        """Infer RDF pair filtering from the configured CG topology."""
        if self.has_bond_angle_topology:
            if self.n_pair_bonds_exclude > 0 or self.n_angle_bonds_exclude > 0:
                return "topology"
            return "all"
        return "intermolecular"




def make_rdf_spec_from_topology(
    topology: TopologySpec,
    *,
    pairs: str | list[list[str]] = "all",
    rmax_nm: float = 2.0,
    dr_nm: float = 0.01,
    exclude_bonded_pairs: bool = True,
    exclude_angle_endpoints: bool = True,
    source: str = "aa2cg.json",
) -> RDFSpec:
    """Build an RDF specification from the unified AA-to-CG JSON topology."""
    if not topology.cgtypes:
        raise ValueError("RDF requires at least one configured bead type")
    if rmax_nm <= 0.0 or dr_nm <= 0.0 or dr_nm >= rmax_nm:
        raise ValueError("RDF requires 0 < bin_width_nm < rmax_nm")
    cgtype_order = list(topology.cgtypes)
    if pairs == "all":
        pair_list = [
            RDFPair(cgtype_order[i], cgtype_order[j])
            for i in range(len(cgtype_order))
            for j in range(i, len(cgtype_order))
        ]
    else:
        pair_list = _unique_pairs((RDFPair(left, right) for left, right in pairs), cgtype_order)
    n_pair = 1 if exclude_bonded_pairs else 0
    n_angle = 1 if exclude_angle_endpoints else 0
    exclusions = _build_topology_exclusions(
        topology,
        n_pair_bonds_exclude=n_pair,
        n_angle_bonds_exclude=n_angle,
    )
    return RDFSpec(
        cgtypes=topology.cgtypes,
        pairs=pair_list,
        rmax_input=float(rmax_nm),
        dr_input=float(dr_nm),
        input_distance_unit="nm",
        output_distance_unit="nm",
        rmax_nm=float(rmax_nm),
        dr_nm=float(dr_nm),
        source=source,
        n_pair_bonds_exclude=n_pair,
        n_angle_bonds_exclude=n_angle,
        topology_exclusions=exclusions,
        topology_exclusion_source="configuration",
        has_bond_angle_topology=topology.has_bond_angle_topology,
        has_topology_exclusion_controls=True,
    )


@dataclass(frozen=True)
class CGFrame:
    positions_nm: list[tuple[float, float, float]]
    box_nm: tuple[float, float, float]


@dataclass(frozen=True)
class CGMappingTemplate:
    bead_names: list[str]
    bead_types: list[str]
    bead_molecule_ids: list[int]
    atom_indices_by_bead: list[list[int]]
    atom_masses_by_bead: list[list[float]]
    atoms_per_molecule: int
    nmol: int


@dataclass(frozen=True)
class RDFResult:
    outdir: Path
    files: dict[str, str]
    report: dict


def _pair_key(left: str, right: str) -> tuple[str, str]:
    return tuple(sorted((left, right)))


def _build_topology_exclusions(
    topology,
    *,
    n_pair_bonds_exclude: int,
    n_angle_bonds_exclude: int,
) -> frozenset[tuple[str, str]]:
    """Build local bead-pair exclusions from configured bonds and angles."""
    exclusions: set[tuple[str, str]] = set()
    if n_pair_bonds_exclude > 0:
        for group in topology.bond_groups:
            for left, right in group.pairs:
                if left != right:
                    exclusions.add(_pair_key(left, right))
    if n_angle_bonds_exclude > 0:
        for group in topology.angle_groups:
            for left, _center, right in group.triples:
                if left != right:
                    exclusions.add(_pair_key(left, right))
    return frozenset(exclusions)


def resolve_rdf_mode_from_spec(spec: RDFSpec) -> tuple[str, str, list[str]]:
    """Return ``(effective, source, warnings)`` for RDF pair filtering."""
    if spec.has_topology_exclusion_controls:
        if spec.topology_exclusions:
            return "topology", "configuration", []
        return spec.infer_default_mode(), "configuration", []
    return (
        "intermolecular",
        "fallback",
        ["no topology exclusion information found; using fallback intermolecular rule"],
    )


def _template_by_index(atoms: list[TprDumpAtom]) -> dict[int, TprDumpAtom]:
    return {atom.index: atom for atom in atoms}


def build_cg_mapping_template(
    *,
    topology_dump: TprTopologyInput,
    coordinates: str | Path,
    mapping: MappingSpec,
    cg_topology: TopologySpec | None = None,
    rdf_spec: RDFSpec | None = None,
) -> CGMappingTemplate:
    tpr = load_tpr_topology(topology_dump)
    gro = parse_gro(coordinates)
    cfg = mapping
    if rdf_spec is not None:
        spec = rdf_spec
    elif cg_topology is not None:
        spec = make_rdf_spec_from_topology(cg_topology, source=cg_topology.source or "configuration")
    else:
        raise ValueError("CG topology or RDF specification is required to build an RDF mapping template")
    bead_type_by_name = {bead: cg_type for cg_type, beads in spec.cgtypes.items() for bead in beads}
    max_mapped_atom = max(max(bead.atom_indices_1based or [0]) for bead in cfg.beads)
    atoms_per_molecule = tpr.atoms_per_molecule or max_mapped_atom
    if gro.natoms % atoms_per_molecule != 0:
        raise ValueError(f"GRO atom count {gro.natoms} is not divisible by atoms_per_molecule={atoms_per_molecule}")
    nmol = gro.natoms // atoms_per_molecule
    atom_template = _template_by_index(tpr.atoms)

    bead_names: list[str] = []
    bead_types: list[str] = []
    bead_molecule_ids: list[int] = []
    atom_indices_by_bead: list[list[int]] = []
    atom_masses_by_bead: list[list[float]] = []
    for mol_index in range(1, nmol + 1):
        offset = (mol_index - 1) * atoms_per_molecule
        for bead in cfg.beads:
            global_indices = [offset + idx for idx in bead.atom_indices_1based]
            masses = [float(atom_template[idx].mass) for idx in bead.atom_indices_1based]
            bead_names.append(bead.name)
            bead_types.append(bead_type_by_name.get(bead.name, bead.name))
            bead_molecule_ids.append(mol_index)
            atom_indices_by_bead.append(global_indices)
            atom_masses_by_bead.append(masses)

    return CGMappingTemplate(
        bead_names=bead_names,
        bead_types=bead_types,
        bead_molecule_ids=bead_molecule_ids,
        atom_indices_by_bead=atom_indices_by_bead,
        atom_masses_by_bead=atom_masses_by_bead,
        atoms_per_molecule=atoms_per_molecule,
        nmol=nmol,
    )


def map_aa_frame_to_cg(
    aa_positions_nm: Sequence[Sequence[float]],
    box_nm: Sequence[float],
    template: CGMappingTemplate,
) -> CGFrame:
    cg_positions: list[tuple[float, float, float]] = []
    for atom_indices, masses in zip(template.atom_indices_by_bead, template.atom_masses_by_bead):
        total_mass = sum(masses)
        if total_mass <= 0:
            raise ValueError("Cannot map bead with non-positive total mass")
        weighted = [0.0, 0.0, 0.0]
        for atom_index, mass in zip(atom_indices, masses):
            pos = aa_positions_nm[atom_index - 1]
            weighted[0] += mass * float(pos[0])
            weighted[1] += mass * float(pos[1])
            weighted[2] += mass * float(pos[2])
        cg_positions.append((weighted[0] / total_mass, weighted[1] / total_mass, weighted[2] / total_mass))
    return CGFrame(positions_nm=cg_positions, box_nm=(float(box_nm[0]), float(box_nm[1]), float(box_nm[2])))


def _minimum_image_delta(value: float, box_length: float) -> float:
    if box_length <= 0:
        raise ValueError("RDF requires positive orthorhombic box lengths")
    return value - box_length * round(value / box_length)


def _distance_pbc(a: Sequence[float], b: Sequence[float], box: Sequence[float]) -> float:
    dx = _minimum_image_delta(float(a[0]) - float(b[0]), float(box[0]))
    dy = _minimum_image_delta(float(a[1]) - float(b[1]), float(box[1]))
    dz = _minimum_image_delta(float(a[2]) - float(b[2]), float(box[2]))
    return math.sqrt(dx * dx + dy * dy + dz * dz)



def _require_numpy():
    try:
        import numpy as np  # type: ignore[import-not-found]
    except ImportError as exc:
        raise RuntimeError(
            "NumPy is required for RDF histogram accumulation. "
            "Install Multiscale with its runtime dependencies using `pip install -e .`, "
            "or set rdf.histogram_engine to 'python' for the slower pure-Python implementation."
        ) from exc
    return np


def _require_kdtree():
    np = _require_numpy()
    try:
        from scipy.spatial import cKDTree  # type: ignore[import-not-found]
    except ImportError as exc:
        raise RuntimeError(
            "SciPy is required for rdf.histogram_engine='kdtree'. "
            "Install SciPy or select the 'numpy' or 'python' RDF backend."
        ) from exc
    return np, cKDTree


def _format_duration(seconds: float) -> str:
    if seconds < 0 or not math.isfinite(seconds):
        return "unknown"
    if seconds < 60:
        return f"{seconds:.0f}s"
    minutes, sec = divmod(int(seconds), 60)
    if minutes < 60:
        return f"{minutes}m{sec:02d}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h{minutes:02d}m{sec:02d}s"


def _validate_box_np(box: Sequence[float], np):
    box_arr = np.asarray(box, dtype=np.float64)
    if box_arr.shape[0] < 3:
        raise ValueError("RDF requires three orthorhombic box lengths")
    box_arr = box_arr[:3]
    if np.any(box_arr <= 0.0):
        raise ValueError("RDF requires positive orthorhombic box lengths")
    return box_arr


def _histogram_distances_np(delta, box_arr, *, rmax_nm: float, dr_nm: float, n_bins: int, np):
    if delta.size == 0:
        return np.zeros(n_bins, dtype=np.int64)
    delta = np.asarray(delta, dtype=np.float64)
    delta -= box_arr * np.rint(delta / box_arr)
    dist2 = np.einsum("ij,ij->i", delta, delta)
    in_range = dist2 < (rmax_nm * rmax_nm)
    if not np.any(in_range):
        return np.zeros(n_bins, dtype=np.int64)
    distances = np.sqrt(dist2[in_range])
    bin_ids = np.floor(distances / dr_nm).astype(np.int64)
    valid = (bin_ids >= 0) & (bin_ids < n_bins)
    return np.bincount(bin_ids[valid], minlength=n_bins)[:n_bins]



def _build_numpy_pair_plans(
    *,
    spec: RDFSpec,
    indices_by_type: dict[str, list[int]],
    molecule_ids: Sequence[int],
    bead_names: Sequence[str] | None,
    mode: str,
    np,
) -> dict[str, dict]:
    molecule_ids_np = np.asarray(molecule_ids, dtype=np.int64)
    bead_names_np = np.asarray(bead_names, dtype=object) if bead_names is not None else None
    plans: dict[str, dict] = {}
    for pair in spec.pairs:
        left_indices = np.asarray(indices_by_type.get(pair.left, []), dtype=np.int64)
        right_indices = np.asarray(indices_by_type.get(pair.right, []), dtype=np.int64)
        if pair.left == pair.right:
            candidate_pairs = left_indices.size * (left_indices.size - 1) // 2
            if left_indices.size < 2:
                pair_i = np.asarray([], dtype=np.int64)
                pair_j = np.asarray([], dtype=np.int64)
            else:
                local_i, local_j = np.triu_indices(left_indices.size, k=1)
                pair_i = left_indices[local_i]
                pair_j = left_indices[local_j]
        else:
            candidate_pairs = left_indices.size * right_indices.size
            if left_indices.size == 0 or right_indices.size == 0:
                pair_i = np.asarray([], dtype=np.int64)
                pair_j = np.asarray([], dtype=np.int64)
            else:
                pair_i = np.repeat(left_indices, right_indices.size)
                pair_j = np.tile(right_indices, left_indices.size)

        if pair_i.size:
            same_molecule = molecule_ids_np[pair_i] == molecule_ids_np[pair_j]
            if mode == "intermolecular":
                keep = ~same_molecule
                pair_i = pair_i[keep]
                pair_j = pair_j[keep]
            elif mode == "topology" and bead_names_np is not None and spec.topology_exclusions:
                excluded_same_molecule = np.zeros(pair_i.shape, dtype=bool)
                same_positions = np.nonzero(same_molecule)[0]
                for pos in same_positions:
                    local_key = _pair_key(str(bead_names_np[pair_i[pos]]), str(bead_names_np[pair_j[pos]]))
                    if local_key in spec.topology_exclusions:
                        excluded_same_molecule[pos] = True
                keep = ~excluded_same_molecule
                pair_i = pair_i[keep]
                pair_j = pair_j[keep]

        plans[pair.name] = {
            "pair_i": pair_i,
            "pair_j": pair_j,
            "candidate_pairs": int(candidate_pairs),
            "evaluated_pairs": int(pair_i.size),
        }
    return plans


def _histogram_pair_plan_numpy(
    positions,
    box: Sequence[float],
    pair_i,
    pair_j,
    *,
    rmax_nm: float,
    dr_nm: float,
    n_bins: int,
    np,
):
    if pair_i.size == 0:
        return np.zeros(n_bins, dtype=np.int64)
    box_arr = _validate_box_np(box, np)
    delta = positions[pair_i] - positions[pair_j]
    return _histogram_distances_np(
        delta,
        box_arr,
        rmax_nm=rmax_nm,
        dr_nm=dr_nm,
        n_bins=n_bins,
        np=np,
    )

def _add_numpy_counts_to_list(counts: list[int], increments) -> None:
    for idx, value in enumerate(increments):
        if value:
            counts[idx] += int(value)


def _build_kdtree_plan(
    *,
    spec: RDFSpec,
    bead_types: Sequence[str],
    molecule_ids: Sequence[int],
    bead_names: Sequence[str] | None,
    np,
) -> dict:
    """Prepare immutable type and topology lookups for cutoff neighbor queries."""
    type_index = {name: index for index, name in enumerate(spec.cgtypes)}
    n_types = len(type_index)
    pair_name_by_code = {
        type_index[pair.left] * n_types + type_index[pair.right]: pair.name
        for pair in spec.pairs
    }
    requested_pair_codes = np.asarray(list(pair_name_by_code), dtype=np.int64)
    plan = {
        "bead_type_codes": np.asarray([type_index[name] for name in bead_types], dtype=np.int64),
        "molecule_ids": np.asarray(molecule_ids, dtype=np.int64),
        "n_types": n_types,
        "pair_name_by_code": pair_name_by_code,
        "requested_pair_codes": requested_pair_codes,
        "bead_name_codes": None,
        "n_bead_names": 0,
        "topology_exclusion_codes": np.asarray([], dtype=np.int64),
    }
    if bead_names is None or not spec.topology_exclusions:
        return plan

    name_index = {name: index for index, name in enumerate(dict.fromkeys(bead_names))}
    n_bead_names = len(name_index)
    exclusion_codes = []
    for left, right in spec.topology_exclusions:
        left_code = name_index.get(left)
        right_code = name_index.get(right)
        if left_code is not None and right_code is not None:
            low, high = sorted((left_code, right_code))
            exclusion_codes.append(low * n_bead_names + high)
    plan["bead_name_codes"] = np.asarray([name_index[name] for name in bead_names], dtype=np.int64)
    plan["n_bead_names"] = n_bead_names
    plan["topology_exclusion_codes"] = np.asarray(exclusion_codes, dtype=np.int64)
    return plan


def _histogram_kdtree_frame(
    positions,
    box: Sequence[float],
    *,
    plan: dict,
    spec: RDFSpec,
    mode: str,
    hist: dict[str, list[int]],
    n_bins: int,
    np,
    cKDTree,
):
    """Accumulate one RDF frame from unique periodic neighbors inside rmax."""
    box_arr = _validate_box_np(box, np)
    wrapped = np.mod(positions, box_arr)
    pair_indices = cKDTree(wrapped, boxsize=box_arr).query_pairs(spec.rmax_nm, output_type="ndarray")
    if pair_indices.size == 0:
        empty = np.asarray([], dtype=np.int64)
        return empty, empty

    pair_i = pair_indices[:, 0]
    pair_j = pair_indices[:, 1]
    bead_type_codes = plan["bead_type_codes"]
    low_types = np.minimum(bead_type_codes[pair_i], bead_type_codes[pair_j])
    high_types = np.maximum(bead_type_codes[pair_i], bead_type_codes[pair_j])
    type_pair_codes = low_types * plan["n_types"] + high_types
    requested = np.isin(type_pair_codes, plan["requested_pair_codes"])
    neighbor_pair_codes = type_pair_codes[requested]
    if not np.any(requested):
        return neighbor_pair_codes, np.asarray([], dtype=np.int64)

    keep = requested.copy()
    if mode == "intermolecular":
        molecule_ids = plan["molecule_ids"]
        keep &= molecule_ids[pair_i] != molecule_ids[pair_j]
    elif mode == "topology" and plan["topology_exclusion_codes"].size:
        molecule_ids = plan["molecule_ids"]
        bead_name_codes = plan["bead_name_codes"]
        low_names = np.minimum(bead_name_codes[pair_i], bead_name_codes[pair_j])
        high_names = np.maximum(bead_name_codes[pair_i], bead_name_codes[pair_j])
        local_pair_codes = low_names * plan["n_bead_names"] + high_names
        excluded = (molecule_ids[pair_i] == molecule_ids[pair_j]) & np.isin(
            local_pair_codes, plan["topology_exclusion_codes"]
        )
        keep &= ~excluded

    if not np.any(keep):
        return neighbor_pair_codes, np.asarray([], dtype=np.int64)

    pair_i = pair_i[keep]
    pair_j = pair_j[keep]
    kept_type_pair_codes = type_pair_codes[keep]
    delta = wrapped[pair_i] - wrapped[pair_j]
    delta -= box_arr * np.rint(delta / box_arr)
    distances = np.sqrt(np.einsum("ij,ij->i", delta, delta))
    in_range = distances < spec.rmax_nm
    if not np.any(in_range):
        return neighbor_pair_codes, np.asarray([], dtype=np.int64)

    evaluated_pair_codes = kept_type_pair_codes[in_range]
    bin_ids = np.floor(distances[in_range] / spec.dr_nm).astype(np.int64)
    valid_bins = (bin_ids >= 0) & (bin_ids < n_bins)
    if not np.any(valid_bins):
        return neighbor_pair_codes, np.asarray([], dtype=np.int64)

    evaluated_pair_codes = evaluated_pair_codes[valid_bins]
    bin_ids = bin_ids[valid_bins]
    for pair_code, pair_name in plan["pair_name_by_code"].items():
        selected = evaluated_pair_codes == pair_code
        if np.any(selected):
            _add_numpy_counts_to_list(
                hist[pair_name],
                np.bincount(bin_ids[selected], minlength=n_bins)[:n_bins],
            )
    return neighbor_pair_codes, evaluated_pair_codes


def compute_rdf_histograms(
    frames: Iterable[CGFrame],
    *,
    bead_types: Sequence[str],
    molecule_ids: Sequence[int],
    spec: RDFSpec,
    bead_names: Sequence[str] | None = None,
    mode: str = "intermolecular",
    rdf_engine: str = "numpy",
    progress_every: int = 0,
    total_frames: int | None = None,
    progress_pairs: bool = False,
    log: Callable[[str], None] | None = None,
) -> dict:
    rdf_engine = rdf_engine.lower().replace("-", "").replace("_", "")
    if mode not in {"intermolecular", "all", "topology"}:
        raise ValueError("mode must be 'intermolecular', 'all', or 'topology'")
    if mode == "topology" and bead_names is None:
        raise ValueError("mode='topology' requires bead_names so local topology exclusions can be applied")
    if rdf_engine not in {"numpy", "python", "kdtree"}:
        raise ValueError("rdf_engine must be 'numpy', 'python', or 'kdtree'")
    if progress_every < 0:
        raise ValueError("progress_every must be >= 0")
    n_bins = int(math.ceil(spec.rmax_nm / spec.dr_nm))
    hist = {pair.name: [0 for _ in range(n_bins)] for pair in spec.pairs}
    frame_count = 0
    volume_sum = 0.0
    counts_by_type = {cg_type: bead_types.count(cg_type) for cg_type in spec.cgtypes}

    indices_by_type = {
        cg_type: [idx for idx, bead_type in enumerate(bead_types) if bead_type == cg_type]
        for cg_type in spec.cgtypes
    }

    np = None
    np_pair_plans = None
    kdtree = None
    kdtree_plan = None
    if rdf_engine == "numpy":
        np = _require_numpy()
        np_pair_plans = _build_numpy_pair_plans(
            spec=spec,
            indices_by_type=indices_by_type,
            molecule_ids=molecule_ids,
            bead_names=bead_names,
            mode=mode,
            np=np,
        )
    elif rdf_engine == "kdtree":
        np, kdtree = _require_kdtree()
        kdtree_plan = _build_kdtree_plan(
            spec=spec,
            bead_types=bead_types,
            molecule_ids=molecule_ids,
            bead_names=bead_names,
            np=np,
        )

    started_all = time.perf_counter()
    for frame in frames:
        frame_count += 1
        frame_started = time.perf_counter()
        show_progress = bool(log and progress_every and (frame_count == 1 or frame_count % progress_every == 0))
        total_label = str(total_frames) if total_frames is not None else "?"
        box = frame.box_nm
        volume_sum += box[0] * box[1] * box[2]
        if show_progress:
            log(
                f"[rdf] frame {frame_count}/{total_label}: "
                f"box=({box[0]:.5f}, {box[1]:.5f}, {box[2]:.5f}) nm; "
                f"starting {len(spec.pairs)} RDF pair histogram(s) with {rdf_engine} backend"
            )

        positions_np = None
        if rdf_engine in {"numpy", "kdtree"}:
            positions_np = np.asarray(frame.positions_nm, dtype=np.float64)
            if positions_np.ndim != 2 or positions_np.shape[1] != 3:
                raise ValueError(f"CG frame positions must have shape (n_beads, 3), got {positions_np.shape}")

        if rdf_engine == "kdtree":
            neighbor_pair_codes, evaluated_pair_codes = _histogram_kdtree_frame(
                positions_np,
                box,
                plan=kdtree_plan,
                spec=spec,
                mode=mode,
                hist=hist,
                n_bins=n_bins,
                np=np,
                cKDTree=kdtree,
            )
            if show_progress and progress_pairs:
                for pair_index, pair in enumerate(spec.pairs, start=1):
                    pair_code = next(
                        code for code, pair_name in kdtree_plan["pair_name_by_code"].items() if pair_name == pair.name
                    )
                    neighbor_count = int(np.count_nonzero(neighbor_pair_codes == pair_code))
                    evaluated_count = int(np.count_nonzero(evaluated_pair_codes == pair_code))
                    log(
                        f"[rdf] frame {frame_count}/{total_label}: "
                        f"pair {pair_index}/{len(spec.pairs)} {pair.name} "
                        f"neighbor_candidates={neighbor_count}, evaluated={evaluated_count}"
                    )
            if show_progress:
                elapsed = time.perf_counter() - frame_started
                message = f"[rdf] frame {frame_count}/{total_label}: done in {elapsed:.2f} s"
                if total_frames is not None and frame_count > 0:
                    avg_elapsed = (time.perf_counter() - started_all) / frame_count
                    remaining = max(total_frames - frame_count, 0)
                    message += f"; avg={avg_elapsed:.2f} s/frame; eta={_format_duration(avg_elapsed * remaining)}"
                log(message)
            continue

        for pair_index, pair in enumerate(spec.pairs, start=1):
            left_indices = indices_by_type.get(pair.left, [])
            right_indices = indices_by_type.get(pair.right, [])
            counts = hist[pair.name]
            if show_progress and progress_pairs:
                if rdf_engine == "numpy":
                    plan = np_pair_plans[pair.name]
                    candidate_pairs = plan["candidate_pairs"]
                    evaluated_pairs = plan["evaluated_pairs"]
                    log(
                        f"[rdf] frame {frame_count}/{total_label}: "
                        f"pair {pair_index}/{len(spec.pairs)} {pair.name} "
                        f"candidates={candidate_pairs}, evaluated={evaluated_pairs}"
                    )
                else:
                    if pair.left == pair.right:
                        candidate_pairs = len(left_indices) * (len(left_indices) - 1) // 2
                    else:
                        candidate_pairs = len(left_indices) * len(right_indices)
                    log(
                        f"[rdf] frame {frame_count}/{total_label}: "
                        f"pair {pair_index}/{len(spec.pairs)} {pair.name} "
                        f"candidates={candidate_pairs}"
                    )

            if rdf_engine == "numpy":
                plan = np_pair_plans[pair.name]
                increments = _histogram_pair_plan_numpy(
                    positions_np,
                    box,
                    plan["pair_i"],
                    plan["pair_j"],
                    rmax_nm=spec.rmax_nm,
                    dr_nm=spec.dr_nm,
                    n_bins=n_bins,
                    np=np,
                )
                _add_numpy_counts_to_list(counts, increments)
                continue

            if pair.left == pair.right:
                iterable = ((left_indices[i], left_indices[j]) for i in range(len(left_indices)) for j in range(i + 1, len(left_indices)))
            else:
                iterable = ((i, j) for i in left_indices for j in right_indices)
            for i, j in iterable:
                if i == j:
                    continue
                if molecule_ids[i] == molecule_ids[j]:
                    if mode == "intermolecular":
                        continue
                    if mode == "topology" and bead_names is not None:
                        if _pair_key(bead_names[i], bead_names[j]) in spec.topology_exclusions:
                            continue
                distance = _distance_pbc(frame.positions_nm[i], frame.positions_nm[j], box)
                if 0.0 <= distance < spec.rmax_nm:
                    counts[int(distance / spec.dr_nm)] += 1

        if show_progress:
            elapsed = time.perf_counter() - frame_started
            message = f"[rdf] frame {frame_count}/{total_label}: done in {elapsed:.2f} s"
            if total_frames is not None and frame_count > 0:
                avg_elapsed = (time.perf_counter() - started_all) / frame_count
                remaining = max(total_frames - frame_count, 0)
                message += f"; avg={avg_elapsed:.2f} s/frame; eta={_format_duration(avg_elapsed * remaining)}"
            log(message)

    if frame_count == 0:
        raise ValueError("No frames were processed for RDF")
    avg_volume = volume_sum / frame_count
    radii = [(i + 0.5) * spec.dr_nm for i in range(n_bins)]
    shell_volumes = [
        (4.0 / 3.0) * math.pi * (((i + 1) * spec.dr_nm) ** 3 - (i * spec.dr_nm) ** 3)
        for i in range(n_bins)
    ]
    gr: dict[str, list[float]] = {}
    for pair in spec.pairs:
        n_ref = counts_by_type.get(pair.left, 0)
        n_target = counts_by_type.get(pair.right, 0)
        rho_target = n_target / avg_volume if avg_volume > 0 else 0.0
        values: list[float] = []
        for raw, shell_volume in zip(hist[pair.name], shell_volumes):
            effective = raw * 2 if pair.left == pair.right else raw
            denom = frame_count * n_ref * rho_target * shell_volume
            values.append(effective / denom if denom > 0 else 0.0)
        gr[pair.name] = values

    return {
        "radii_nm": radii,
        "histograms": hist,
        "g_r": gr,
        "n_frames": frame_count,
        "box_volume_average": avg_volume,
        "bead_counts_per_type": counts_by_type,
        "progress_every": progress_every,
        "rdf_engine": rdf_engine,
        "topology_exclusion_count": len(spec.topology_exclusions),
        "topology_exclusion_pairs": [list(pair) for pair in sorted(spec.topology_exclusions)],
    }

def write_rdf_outputs(
    *,
    outdir: str | Path,
    spec: RDFSpec,
    rdf_data: dict,
    mode: str,
    source_trajectory: str | Path,
    inputs: dict,
) -> RDFResult:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    files: dict[str, str] = {}
    for pair in spec.pairs:
        path = out / f"AA_rdf_{pair.name}.dat"
        lines = [
            "# Multiscale RDF",
            f"# pair: {pair.name}",
            "# r_unit: nm",
            "# g_r_unit: dimensionless",
            f"# mode: {mode}",
            f"# rdf_mode_requested: {inputs.get('rdf_mode_requested', 'auto')}",
            f"# rdf_mode_effective: {inputs.get('rdf_mode_effective', mode)}",
            f"# rdf_mode_source: {inputs.get('rdf_mode_source', 'unknown')}",
            f"# topology_exclusion_count: {rdf_data.get('topology_exclusion_count', 0)}",
            f"# source_trajectory: {source_trajectory}",
            "# generated_by: multiscale aa2cg run",
        ]
        for radius, value in zip(rdf_data["radii_nm"], rdf_data["g_r"][pair.name]):
            if not math.isfinite(value):
                raise ValueError(f"Non-finite RDF value for {pair.name} at r={radius}")
            lines.append(f"{radius:.6f}  {value:.6f}")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        files[pair.name] = str(path)

    report = {
        "entrypoint": "rdf",
        "inputs": inputs,
        "outputs": files,
        "rmax_input": spec.rmax_input,
        "dr_input": spec.dr_input,
        "input_distance_unit": spec.input_distance_unit,
        "output_distance_unit": spec.output_distance_unit,
        "rmax_nm": spec.rmax_nm,
        "dr_nm": spec.dr_nm,
        "normalization": "3d_number_density",
        "mode": mode,
        "rdf_mode_requested": inputs.get("rdf_mode_requested", "auto"),
        "rdf_mode_effective": inputs.get("rdf_mode_effective", mode),
        "rdf_mode_source": inputs.get("rdf_mode_source", "unknown"),
        "rdf_engine": rdf_data.get("rdf_engine", inputs.get("rdf_engine")),
        "n_frames": rdf_data["n_frames"],
        "box_volume_average": rdf_data["box_volume_average"],
        "pair_counts": {name: sum(values) for name, values in rdf_data["histograms"].items()},
        "bead_counts_per_type": rdf_data["bead_counts_per_type"],
        "topology_exclusions": {
            "n_pair_bonds_exclude": spec.n_pair_bonds_exclude,
            "n_angle_bonds_exclude": spec.n_angle_bonds_exclude,
            "count": rdf_data.get("topology_exclusion_count", len(spec.topology_exclusions)),
            "pairs": rdf_data.get("topology_exclusion_pairs", [list(pair) for pair in sorted(spec.topology_exclusions)]),
        },
        "warnings": inputs.get("warnings", []),
        "histogram_counts_summary": {
            name: {"sum": sum(values), "max": max(values) if values else 0}
            for name, values in rdf_data["histograms"].items()
        },
        "notes": [
            "Native on-the-fly RDF; no CG.xtc is generated.",
            "RDF topology and exclusion rules are derived from the configured CG model.",
            "RDF results are reproducible when trajectory preprocessing, mapping, and topology settings are unchanged.",
            "Current implementation writes nonbonded CG bead-type RDF only and does not perform IBI, IMC, NG, or potential-table generation.",
        ],
    }
    report_path = out / "rdf_report.json"
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    files["report"] = str(report_path)
    return RDFResult(outdir=out, files=files, report=report)
