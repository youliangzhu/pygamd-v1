'''
CuPy-native particle device-data layer (MVP-1 / CuInfo).

`CuInfo` is the single owner of the dynamics main-line particle data on the GPU.
It turns the host NumPy arrays produced by the snapshot reader into CuPy device
arrays (`cupy.asarray`), enforces the dtype contract, and hands out cached 1-D
CuPy views for RawKernel launches.  This replaces the previously-scattered
`cupy.asarray(info.d_xxx)` wraps inside RawPair / RawNVT: those backends now ask
CuInfo for their device arrays.

Dtype contract (normalized on upload):
    position / velocity / force / virial_potential / box / mass : float32
    type / image / tag / rtag                                   : int32

Layout (unchanged from the engine's existing convention):
    d_pos               (npa, 4)  x, y, z, type_id(as float)
    d_vel               (npa, 4)  vx, vy, vz, mass
    d_velo              (npa, 4)  DPDGWVV predictor velocity, same layout
    d_force             (npa, 3)
    d_virial_potential  (npa, 2)
    d_image             (npa, 3)
    d_box               (3,) or (6,)
    d_type              (npa,)    int32   (derived from pos[:, 3])
    d_mass              (npa,)    float32 (derived from vel[:, 3])
    d_tag, d_rtag       (npa,)    int32

Import safety: this module imports only NumPy at load time.  CuPy is imported
lazily (inside __init__, via the central bridge), so `import pygamd` and even
`import pygamd.backend.cupy` never require a CUDA driver.
'''

import numpy as np

from pygamd.backend.cupy.context import import_cupy
from pygamd.backend.cupy.runtime import CORE_FIELDS

_F32 = np.float32
_I32 = np.int32

# dtype contract for optional per-particle fields (real -> float32, index -> int32)
_FIELD_DTYPE = {
    "charge": _F32, "diameter": _F32, "rotangle": _F32, "rotation": _F32,
    "orientation": _F32, "quaternion": _F32, "inert": _F32, "multipole": _F32,
    "body": _I32, "molecule": _I32, "init": _I32, "cris": _I32,
}


class CuTopology:
    '''CuPy-native topology block for one kind (bond / angle / dihedral / vsite).

    Stores connectivity in BOTH tag-based and index-based form, plus a per-entry
    type id with a local type-name -> id map:
        host_tags / d_tags        (count, members) int32  -- tags as read (XML/MST)
        host_indices / d_indices  (count, members) int32  -- current particle index
        host_type / d_type        (count,)         int32  -- bonded type id
        typemap                   list[str]                -- type id -> name
        entries                   raw [name, tag0, ...] host list (dump roundtrip)

    index = rtag[tag].  At read time rtag is identity so index == tag; after a
    future Sort reorders particles, call refresh_indices(new_rtag) to rebuild the
    index-based arrays (tags never change).  RawBond/RawAngle/RawDihedral will
    consume the *index-based* arrays.
    '''

    _MEMBERS = {"bond": 2, "angle": 3, "dihedral": 4, "vsite": 4}

    def __init__(self, cp, kind, entries, npa, rtag_host):
        self._cp = cp
        self.kind = kind
        self.entries = list(entries)
        self.count = len(self.entries)
        self.members = self._MEMBERS.get(
            kind, (len(self.entries[0]) - 1) if self.entries else 0)
        self.typemap = []
        host_type = np.zeros(self.count, dtype=_I32)
        host_tags = np.zeros((self.count, self.members), dtype=_I32)
        for i, e in enumerate(self.entries):
            host_type[i] = self._type_id(e[0])
            for j in range(self.members):
                t = int(e[1 + j])
                if t < 0 or t >= npa:
                    raise RuntimeError(
                        "Error! %s entry %d references particle tag %d outside "
                        "[0, %d)." % (kind, i, t, npa))
                host_tags[i, j] = t
        self.host_type = host_type
        self.host_tags = host_tags
        self.host_indices = self._tags_to_indices(host_tags, rtag_host)
        self.d_type = cp.asarray(host_type)
        self.d_tags = cp.asarray(host_tags)
        self.d_indices = cp.asarray(self.host_indices)
        self._host_indices_valid = True

    def _type_id(self, name):
        if name not in self.typemap:
            self.typemap.append(name)
        return self.typemap.index(name)

    @staticmethod
    def _tags_to_indices(host_tags, rtag_host):
        rt = np.ascontiguousarray(rtag_host, dtype=_I32)
        return rt[host_tags].astype(_I32, copy=False)

    def refresh_indices(self, rtag):
        '''Recompute index-based connectivity from the (unchanged) tags using a
        new rtag.  `rtag` may be a host NumPy array or a CuPy device array.

        Device input is the sort hot path: update d_indices IN PLACE from
        d_rtag[d_tags], so bonded-force graph signatures keep the same device
        array identity across sort operations.  host_indices are refreshed lazily
        only when to_numpy("indices") is requested.
        '''
        cp = self._cp
        if self.count == 0:
            return
        try:
            xp = cp.get_array_module(rtag)
        except Exception:
            xp = None
        if xp is cp:
            idx = rtag[self.d_tags].astype(_I32, copy=False)
            if self.d_indices.shape == idx.shape:
                self.d_indices[...] = idx
            else:
                # Shape changes should not happen for normal sort, but keep the
                # method robust for reconstructed topology blocks.
                self.d_indices = cp.asarray(idx, dtype=_I32)
            self._host_indices_valid = False
            return
        self.host_indices = self._tags_to_indices(self.host_tags, rtag)
        if self.d_indices.shape == self.host_indices.shape:
            self.d_indices[...] = cp.asarray(self.host_indices, dtype=_I32)
        else:
            self.d_indices = cp.asarray(self.host_indices, dtype=_I32)
        self._host_indices_valid = True

    def to_numpy(self, which="indices"):
        if which == "tags":
            return self.host_tags.copy()
        if which == "type":
            return self.host_type.copy()
        if not getattr(self, "_host_indices_valid", True):
            self.host_indices = self._cp.asnumpy(self.d_indices).astype(_I32, copy=False)
            self._host_indices_valid = True
        return self.host_indices.copy()


class CuInfo:
    '''CuPy-native device-data container for the core particle arrays.

    Construct from host NumPy arrays; CuInfo uploads them to CuPy with the dtype
    contract above and owns the resulting device buffers.  `dev(name)` returns a
    cached, id-guarded 1-D CuPy view suitable for a RawKernel launch.
    '''

    # marker so consumers can detect the CuPy-native device-data layer
    is_cuinfo = True

    def __init__(self, npa, pos, vel, force, virial_potential, image, box,
                 tag=None, rtag=None):
        cp = import_cupy("CuInfo CuPy device-data layer")
        self._cp = cp
        self.npa = int(npa)

        # --- canonical CuPy device arrays (host NumPy -> cupy, dtype-normalized) -
        self.d_pos = self._upload(cp, pos, _F32)                 # (npa, 4)
        self.d_vel = self._upload(cp, vel, _F32)                 # (npa, 4)
        # DPDGWVV predictor velocity.  It starts as a copy of vel and is updated
        # by the RawDPDGWVV first step.  Keep it as a core stable CuPy array so
        # RawDPD can graph-capture a fixed pointer when velocity_source="velo".
        self.d_velo = self.d_vel.copy()                           # (npa, 4)
        self.d_force = self._upload(cp, force, _F32)             # (npa, 3)
        self.d_virial_potential = self._upload(cp, virial_potential, _F32)  # (npa, 2)
        self.d_image = self._upload(cp, image, _I32)             # (npa, 3)
        self.d_box = self._upload(cp, box, _F32)                 # (3,) / (6,)
        self.d_tag = self._upload(
            cp, tag if tag is not None else np.arange(self.npa), _I32)
        self.d_rtag = self._upload(
            cp, rtag if rtag is not None else np.arange(self.npa), _I32)

        # derived per-particle scalars (CuPy-native; not consumed by the current
        # RawKernels, which read type from pos[:,3] and mass from vel[:,3], but
        # required by the CuInfo field contract and useful for future kernels).
        host_pos = np.ascontiguousarray(pos)
        host_vel = np.ascontiguousarray(vel)
        self.d_type = self._upload(cp, host_pos[:, 3], _I32)     # (npa,)
        self.d_mass = self._upload(cp, host_vel[:, 3], _F32)     # (npa,)

        # cache of raveled (1-D) views, id-guarded against array replacement
        self._views = {}

        # optional per-particle fields (name -> cupy array), beyond the core set
        self._optional_fields = {}
        # topology blocks (kind -> CuTopology), e.g. "bond"/"angle"/"dihedral"
        self.topology = {}

    # --- helpers --------------------------------------------------------------
    @staticmethod
    def _upload(cp, host, dtype):
        '''host array-like -> contiguous CuPy device array of `dtype`.'''
        return cp.asarray(np.ascontiguousarray(host, dtype=dtype))

    @property
    def cp(self):
        '''The cupy module (already imported).'''
        return self._cp

    def dev(self, name, ravel=True):
        '''Return a CuPy view of core array `name` (e.g. "pos", "vel", "force",
        "virial_potential", "image", "box").  Default `ravel=True` returns the
        cached 1-D view used by the RawKernel launches; the underlying buffer is
        owned by CuInfo so the view identity is stable (safe for graph capture).
        '''
        arr = getattr(self, "d_" + name)
        key = (name, ravel)
        cached = self._views.get(key)
        if cached is not None and cached[0] is arr:
            return cached[1]
        view = arr.ravel() if ravel else arr
        self._views[key] = (arr, view)
        return view

    # core fields always present after construction
    _CORE_FIELDS = CORE_FIELDS

    def to_numpy(self, name):
        '''Copy field `name` (core, optional, or topology accessor like "bond")
        back to a host NumPy array (cupy.asnumpy).'''
        arr = getattr(self, "d_" + name, None)
        if arr is None:
            raise KeyError("CuInfo has no field %r" % (name,))
        return self._cp.asnumpy(arr)

    # --- optional per-particle fields -----------------------------------------
    def add_field(self, name, host, dtype=None):
        '''Upload an optional per-particle field to a CuPy device array and
        register it as info.d_<name>.  dtype defaults from the field contract
        (real -> float32, index/id -> int32).'''
        if dtype is None:
            dtype = _FIELD_DTYPE.get(name, _F32)
        arr = self._upload(self._cp, host, dtype)
        setattr(self, "d_" + name, arr)
        self._optional_fields[name] = arr
        return arr

    def has_field(self, name):
        '''True if a (core or optional) per-particle field is present.'''
        return getattr(self, "d_" + name, None) is not None

    def fields(self):
        '''List of available field names (core + optional).'''
        return list(self._CORE_FIELDS) + sorted(self._optional_fields.keys())

    # --- topology -------------------------------------------------------------
    def add_topology(self, kind, entries, rtag_host=None):
        '''Build and register a CuTopology block (bond/angle/dihedral/vsite).
        `entries` are the raw [type_name, tag0, ...] rows from the snapshot.'''
        if rtag_host is None:
            # Topology construction is outside the sort hot path.  Keep the
            # constructor host-based so validation errors and dump roundtrip data
            # remain simple and explicit.
            rtag_host = self._cp.asnumpy(self.d_rtag)
        topo = CuTopology(self._cp, kind, entries, self.npa, rtag_host)
        self.topology[kind] = topo
        return topo

    def has_topology(self, kind):
        return kind in self.topology

    # --- particle reordering (Sort) -------------------------------------------
    def reorder(self, order):
        '''Reorder ALL per-particle device arrays by `order` (new index k takes
        old particle order[k]), IN PLACE so array identity / cached views / the
        snapshot d_* aliases stay valid.  Then rebuild d_rtag (rtag[tag]=index)
        and refresh topology index-based connectivity.  `order` is a length-npa
        permutation (host or device int array).'''
        cp = self._cp
        order = cp.asarray(order).astype(_I32)
        # per-particle arrays to gather (NOT box; NOT rtag -- rtag is rebuilt below)
        names = [n for n in self._CORE_FIELDS if n not in ("box", "rtag")]
        names += list(self._optional_fields.keys())
        for name in names:
            arr = getattr(self, "d_" + name, None)
            if arr is None:
                continue
            arr[:] = arr[order]            # in-place gather (RHS fully evaluated first)
        # rebuild rtag on device: rtag[tag_at_index_k] = k
        idx = cp.arange(self.npa, dtype=_I32)
        self.d_rtag[self.d_tag] = idx
        # Keep cached raveled views: all arrays were updated in place, so view
        # identities remain valid and graph signatures stay stable across sort.
        self.refresh_topology_indices(self.d_rtag)

    def refresh_topology_indices(self, rtag=None):
        '''Rebuild index-based connectivity for ALL topology blocks (call after a
        Sort reorders particles, once d_rtag has been updated).  `rtag` may be
        omitted, host NumPy, or CuPy device; omitted now uses the device d_rtag
        to keep the sort path free of host roundtrips.'''
        if rtag is None:
            rtag = self.d_rtag
        for topo in self.topology.values():
            topo.refresh_indices(rtag)

    # convenience accessors (None / 0 when absent)
    def has_bond(self):
        return "bond" in self.topology

    def has_angle(self):
        return "angle" in self.topology

    def has_dihedral(self):
        return "dihedral" in self.topology

    @property
    def bond_count(self):
        return self.topology["bond"].count if "bond" in self.topology else 0

    @property
    def angle_count(self):
        return self.topology["angle"].count if "angle" in self.topology else 0

    @property
    def dihedral_count(self):
        return self.topology["dihedral"].count if "dihedral" in self.topology else 0

    @property
    def d_bond(self):
        return self.topology["bond"].d_indices if "bond" in self.topology else None

    @property
    def d_angle(self):
        return self.topology["angle"].d_indices if "angle" in self.topology else None

    @property
    def d_dihedral(self):
        return self.topology["dihedral"].d_indices if "dihedral" in self.topology else None
