from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .chem import Chem
from .pipeline_config import resolve_path
from .pdb_template import PdbTemplate, read_pdb


def _atom_indices(values, atom_index_base: int, n_atoms: int) -> list[int]:
    indices = [int(v) - atom_index_base for v in values]
    for idx in indices:
        if idx < 0 or idx >= n_atoms:
            raise ValueError(f"atom index {idx + atom_index_base} is outside template atom range")
    return indices


def _pdb_molecule(pdb: PdbTemplate):
    molecule = Chem.RWMol()
    for atom in pdb.atoms:
        new_atom = Chem.Atom(atom.element)
        new_atom.SetProp("AtomName", atom.name)
        new_atom.SetProp("ResidueName", atom.residue)
        molecule.AddAtom(new_atom)
    for i, j in pdb.bonds:
        molecule.AddBond(i, j, Chem.BondType.SINGLE)
    return molecule


def _adjacency_from_bonds(n_atoms: int, bonds: list[tuple[int, int]]) -> list[set[int]]:
    adjacency = [set() for _ in range(n_atoms)]
    for i, j in bonds:
        adjacency[i].add(j)
        adjacency[j].add(i)
    return adjacency


def _molecule_adjacency(molecule) -> list[set[int]]:
    adjacency = [set() for _ in range(molecule.GetNumAtoms())]
    for bond in molecule.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        adjacency[i].add(j)
        adjacency[j].add(i)
    return adjacency


def _auto_atom_map(molecule, pdb: PdbTemplate) -> dict[int, int]:
    n_chem = molecule.GetNumAtoms()
    n_pdb = len(pdb.atoms)
    if n_chem > n_pdb:
        raise ValueError(f"SMILES molecule has {n_chem} atoms but PDB has only {n_pdb}")

    chem_adj = _molecule_adjacency(molecule)
    pdb_adj = _adjacency_from_bonds(n_pdb, pdb.bonds)
    order = sorted(range(n_chem), key=lambda idx: len(chem_adj[idx]), reverse=True)
    candidates = {}
    for cidx in range(n_chem):
        atom = molecule.GetAtomWithIdx(cidx)
        element = atom.GetSymbol()
        degree = len(chem_adj[cidx])
        candidates[cidx] = [
            pidx for pidx, pdb_atom in enumerate(pdb.atoms)
            if pdb_atom.element == element and len(pdb_adj[pidx]) >= degree
        ]
        if not candidates[cidx]:
            raise ValueError(f"No PDB atom candidate for SMILES atom {cidx} ({element})")

    mapping: dict[int, int] = {}
    used_pdb = set()

    def backtrack(pos: int) -> bool:
        if pos == len(order):
            return True
        cidx = order[pos]
        for pidx in candidates[cidx]:
            if pidx in used_pdb:
                continue
            ok = True
            for mapped_cidx, mapped_pidx in mapping.items():
                chem_has_edge = mapped_cidx in chem_adj[cidx]
                pdb_has_edge = mapped_pidx in pdb_adj[pidx]
                if chem_has_edge != pdb_has_edge:
                    ok = False
                    break
            if not ok:
                continue
            mapping[cidx] = pidx
            used_pdb.add(pidx)
            if backtrack(pos + 1):
                return True
            used_pdb.remove(pidx)
            del mapping[cidx]
        return False

    if not backtrack(0):
        raise ValueError("Could not map SMILES atoms onto PDB atoms; provide smiles_to_pdb_atoms explicitly")
    return dict(sorted(mapping.items()))


def _explicit_atom_map(raw_map, n_chem: int, n_pdb: int, atom_index_base: int) -> dict[int, int]:
    if isinstance(raw_map, list):
        if len(raw_map) != n_chem:
            raise ValueError("smiles_to_pdb_atoms list length must match SMILES atom count")
        pairs = enumerate(raw_map)
    elif isinstance(raw_map, dict):
        pairs = ((int(k), v) for k, v in raw_map.items())
    else:
        raise ValueError("smiles_to_pdb_atoms must be a list or dict")

    mapping = {}
    for chem_idx, pdb_value in pairs:
        pdb_idx = int(pdb_value) - atom_index_base
        if chem_idx < 0 or chem_idx >= n_chem:
            raise ValueError(f"SMILES atom index {chem_idx} out of range")
        if pdb_idx < 0 or pdb_idx >= n_pdb:
            raise ValueError(f"PDB atom index {pdb_value} out of range")
        mapping[chem_idx] = pdb_idx
    if len(mapping) != n_chem:
        raise ValueError("smiles_to_pdb_atoms must map every SMILES atom")
    return mapping


@dataclass
class TemplateSpec:
    name: str
    type: str
    pdb: PdbTemplate
    molecule: object
    chem_to_pdb: dict[int, int]
    atom_index_base: int
    site_names: list[str]
    site_labels: dict[str, str]
    site_atoms: dict[str, list[int]]
    site_chem_atoms: dict[str, list[int]]
    site_centers: dict[str, np.ndarray]
    anchor_sites: list[str]
    auto_permute_equivalent_sites: bool
    fit_mode: str
    smiles: str | None

    @property
    def coords(self) -> np.ndarray:
        return np.asarray([self.pdb.atoms[self.chem_to_pdb[idx]].coord for idx in range(self.molecule.GetNumAtoms())], dtype=float)

    @property
    def atom_names(self) -> list[str]:
        return [self.pdb.atoms[self.chem_to_pdb[idx]].name for idx in range(self.molecule.GetNumAtoms())]


def load_template(name: str, raw: dict, config: dict) -> TemplateSpec:
    path = resolve_path(config, raw["pdb"])
    pdb = read_pdb(path)
    atom_index_base = int(raw.get("atom_index_base", 1))
    smiles = raw.get("smiles")
    add_hydrogens = bool(raw.get("add_hydrogens", True))
    if smiles:
        molecule = Chem.MolFromSmiles(smiles)
        if add_hydrogens:
            molecule = Chem.AddHs(molecule)
        if raw.get("smiles_to_pdb_atoms") is not None:
            chem_to_pdb = _explicit_atom_map(raw["smiles_to_pdb_atoms"], molecule.GetNumAtoms(), len(pdb.atoms), atom_index_base)
        else:
            chem_to_pdb = _auto_atom_map(molecule, pdb)
    else:
        molecule = _pdb_molecule(pdb)
        chem_to_pdb = {idx: idx for idx in range(molecule.GetNumAtoms())}

    sites = raw.get("sites") or raw.get("beads")
    if isinstance(sites, dict):
        sites = [{"name": key, "atoms": value} for key, value in sites.items()]
    if not sites:
        raise ValueError(f"template {name} must define sites/beads")

    site_names = [site["name"] for site in sites]
    site_labels = {site["name"]: site.get("label", site["name"]) for site in sites}
    site_atoms = {
        site["name"]: _atom_indices(site["atoms"], atom_index_base, len(pdb.atoms))
        for site in sites
    }
    pdb_to_chem = {pdb_idx: chem_idx for chem_idx, pdb_idx in chem_to_pdb.items()}
    site_chem_atoms = {
        site: [pdb_to_chem[idx] for idx in indices if idx in pdb_to_chem]
        for site, indices in site_atoms.items()
    }
    site_centers = {
        site: pdb.coords[indices].mean(axis=0)
        for site, indices in site_atoms.items()
    }

    fit_cfg = raw.get("fit", {})
    fit_mode = fit_cfg.get("mode", raw.get("fit_mode", "kabsch"))
    if fit_mode not in {"kabsch", "frame"}:
        raise ValueError(f"template {name} has unsupported fit mode {fit_mode!r}")
    auto_permute_equivalent_sites = bool(fit_cfg.get("auto_permute_equivalent_sites", raw.get("auto_permute_equivalent_sites", False)))
    anchor_sites = fit_cfg.get("sites", raw.get("anchor_sites", raw.get("anchor_beads", site_names)))
    if fit_mode == "frame" and len(anchor_sites) != 3:
        raise ValueError(f"template {name} fit.mode='frame' requires exactly three anchor sites")
    for site in anchor_sites:
        if site not in site_centers:
            raise ValueError(f"template {name} anchor site {site!r} is not defined")

    return TemplateSpec(
        name=name,
        type=raw.get("type", name),
        pdb=pdb,
        molecule=molecule,
        chem_to_pdb=chem_to_pdb,
        atom_index_base=atom_index_base,
        site_names=site_names,
        site_labels=site_labels,
        site_atoms=site_atoms,
        site_chem_atoms=site_chem_atoms,
        site_centers=site_centers,
        anchor_sites=list(anchor_sites),
        auto_permute_equivalent_sites=auto_permute_equivalent_sites,
        fit_mode=fit_mode,
        smiles=smiles,
    )
