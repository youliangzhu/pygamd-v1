'''
Bonded-topology consumption interface (read-only) for the CuPy main line.

This is the stable data interface that the FUTURE RawBond / RawAngle /
RawDihedral force kernels will consume.  It does NOT implement any force; it
only exposes the CuInfo topology CuPy arrays under fixed field names so the
bonded-force step can be implemented later without reworking the data layer.

    get_bond_topology(info)     -> BondedTopologyView (members=2)
    get_angle_topology(info)    -> BondedTopologyView (members=3)
    get_dihedral_topology(info) -> BondedTopologyView (members=4)

Each view exposes:
    kind            "bond" / "angle" / "dihedral"
    members         particles per entry (2 / 3 / 4)
    count           number of entries (0 if the topology is absent)
    has_topology    bool
    d_tags          CuPy int32 (count, members)  -- tag-based (sort-invariant)
    d_indices       CuPy int32 (count, members)  -- current particle index (kernel input)
    d_typeid        CuPy int32 (count,)          -- per-entry bonded type id
    typemap         list[str]                    -- type id -> name
    refresh_indices()  -- rebuild d_indices from tags via the current rtag

Conventions for the future kernels:
  * XML/MST OUTPUT uses tag-based connectivity (`d_tags` / topology.entries):
    stable across a Sort.
  * KERNELS consume index-based connectivity (`d_indices`): valid for the
    current particle ordering.  After a Sort reorders particles you MUST call
    info.cu.refresh_topology_indices() (or view.refresh_indices()) before the
    next bonded-force evaluation.

Pure Python at import time (no cupy / no numba); the CuPy arrays it returns are
the ones CuInfo already allocated.
'''


class BondedTopologyView:
    '''Stable, read-only view over one CuInfo topology block (or an empty view
    when the topology is absent).'''

    def __init__(self, kind, members, cu, topo):
        self.kind = kind
        self.members = members
        self._cu = cu
        self._topo = topo
        self.has_topology = topo is not None

    # Live properties so a post-sort info.cu.refresh_topology_indices() update
    # is immediately visible here.  The sort hot path now updates d_indices in
    # place to preserve graph signatures.
    @property
    def count(self):
        return int(self._topo.count) if self._topo is not None else 0

    @property
    def d_tags(self):
        return self._topo.d_tags if self._topo is not None else None

    @property
    def d_indices(self):
        return self._topo.d_indices if self._topo is not None else None

    @property
    def d_typeid(self):
        return self._topo.d_type if self._topo is not None else None

    @property
    def typemap(self):
        return self._topo.typemap if self._topo is not None else []

    def refresh_indices(self, rtag=None):
        '''Rebuild index-based connectivity from the (unchanged) tags using the
        current rtag.  Call after a Sort reorders particles.  No-op when the
        topology is absent.  If omitted, use the device d_rtag to avoid a host
        copy in the sort path.'''
        if self._topo is None:
            return
        if rtag is None and self._cu is not None:
            rtag = self._cu.d_rtag
        self._topo.refresh_indices(rtag)   # d_indices is a live property

    def __repr__(self):
        return ("BondedTopologyView(kind=%s, members=%d, count=%d, has=%s)"
                % (self.kind, self.members, self.count, self.has_topology))


def _view(info, kind, members):
    cu = getattr(info, "cu", None)
    topo = cu.topology.get(kind) if (cu is not None) else None
    return BondedTopologyView(kind, members, cu, topo)


def get_bond_topology(info):
    '''Bond topology view (2 particles/entry); count==0 if no bonds.'''
    return _view(info, "bond", 2)


def get_angle_topology(info):
    '''Angle topology view (3 particles/entry); count==0 if no angles.'''
    return _view(info, "angle", 3)


def get_dihedral_topology(info):
    '''Dihedral topology view (4 particles/entry); count==0 if no dihedrals.'''
    return _view(info, "dihedral", 4)
