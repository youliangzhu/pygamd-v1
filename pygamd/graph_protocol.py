'''
PYGAMD eager-graph protocol helpers.

This module is intentionally small and CUDA-driver-free.  It defines the
placeholder protocol used by the eager_graph primary execution path to discover
force backends that can participate in CUDA Graph capture.  The goal of this protocol is to
keep the custom-force extension seam open without introducing custom-force
benchmarks or a full public custom-force contract yet.

Force wrappers are duck-typed.  A wrapper may expose either
``wrapper.data.backend`` or ``wrapper.data`` as the graph backend.  Built-in and
future custom backends are treated the same way: no concrete force class is
hard-coded in the graph planner.

Minimal graph-force protocol, placeholder v0:

    graph_ready() -> bool
    graph_unready_reason() -> str | None       # optional but recommended
    graph_signature() -> hashable/tuple        # stable capture signature
    graph_prepare_replay(timestep)             # optional dynamic device inputs
    graph_force_add()                          # additive force contribution
    graph_force_write()                        # exclusive force-buffer writer

``graph_force_add`` is the preferred extension point because most force objects
contribute additively after the graph planner clears the global force/virial
buffers.  ``graph_force_write`` exists for legacy single-force fast paths that
overwrite the global force buffer directly.
'''

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Tuple

FORCE_GRAPH_PROTOCOL_VERSION = "p4-placeholder-v0"

GRAPH_READY = "graph_ready"
GRAPH_UNREADY_REASON = "graph_unready_reason"
GRAPH_SIGNATURE = "graph_signature"
GRAPH_PREPARE_REPLAY = "graph_prepare_replay"
GRAPH_FORCE_ADD = "graph_force_add"
GRAPH_FORCE_WRITE = "graph_force_write"


@dataclass(frozen=True)
class ForceGraphCapabilities:
    backend_name: str
    backend_kind: str
    has_ready: bool
    has_unready_reason: bool
    has_signature: bool
    has_prepare_replay: bool
    can_add: bool
    can_write: bool

    @property
    def graph_capable(self) -> bool:
        return self.has_ready and self.has_signature and (self.can_add or self.can_write)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "backend_name": self.backend_name,
            "backend_kind": self.backend_kind,
            "has_ready": self.has_ready,
            "has_unready_reason": self.has_unready_reason,
            "has_signature": self.has_signature,
            "has_prepare_replay": self.has_prepare_replay,
            "can_add": self.can_add,
            "can_write": self.can_write,
            "graph_capable": self.graph_capable,
        }


def _callable_attr(obj: Any, name: str) -> bool:
    return callable(getattr(obj, name, None))


def force_graph_backend_from_wrapper(force_obj: Any) -> Any:
    '''Return the object that should satisfy the graph-force protocol.

    Public force wrappers in pygamd usually hold implementation state in
    ``.data`` and the CuPy raw implementation under ``.data.backend``.  Future
    custom forces can expose either shape; this resolver deliberately does not
    check concrete classes.
    '''
    data = getattr(force_obj, "data", force_obj)
    backend = getattr(data, "backend", None)
    return backend if backend is not None else data


def force_graph_backend_name(backend: Any) -> str:
    return str(getattr(backend, "graph_name", backend.__class__.__name__))


def force_graph_backend_kind(backend: Any) -> str:
    mod = getattr(backend.__class__, "__module__", "") or ""
    cls = getattr(backend.__class__, "__name__", backend.__class__.__name__)
    if mod:
        return "%s.%s" % (mod, cls)
    return str(cls)


def force_graph_capabilities(backend: Any) -> ForceGraphCapabilities:
    return ForceGraphCapabilities(
        backend_name=force_graph_backend_name(backend),
        backend_kind=force_graph_backend_kind(backend),
        has_ready=_callable_attr(backend, GRAPH_READY),
        has_unready_reason=_callable_attr(backend, GRAPH_UNREADY_REASON),
        has_signature=_callable_attr(backend, GRAPH_SIGNATURE),
        has_prepare_replay=_callable_attr(backend, GRAPH_PREPARE_REPLAY),
        can_add=_callable_attr(backend, GRAPH_FORCE_ADD),
        can_write=_callable_attr(backend, GRAPH_FORCE_WRITE),
    )


def force_graph_unready_reason(backend: Any) -> str:
    caps = force_graph_capabilities(backend)
    if not caps.has_ready:
        return "force_graph_missing_graph_ready"
    if not caps.has_signature:
        return "force_graph_missing_graph_signature"
    if not (caps.can_add or caps.can_write):
        return "force_graph_missing_graph_force_method"
    fn = getattr(backend, GRAPH_UNREADY_REASON, None)
    if callable(fn):
        try:
            reason = fn()
            if reason:
                return str(reason)
        except Exception:
            return "force_graph_unready_reason_failed"
    return "force_graph_not_ready"


def force_graph_ready(backend: Any) -> bool:
    caps = force_graph_capabilities(backend)
    if not caps.graph_capable:
        return False
    fn = getattr(backend, GRAPH_READY, None)
    try:
        return bool(fn()) if callable(fn) else False
    except Exception:
        return False


def force_graph_signature(backend: Any) -> Tuple[Any, ...]:
    fn = getattr(backend, GRAPH_SIGNATURE, None)
    if callable(fn):
        try:
            sig = fn()
        except Exception:
            return ("force_graph_signature_failed", force_graph_backend_kind(backend), id(backend))
        if isinstance(sig, tuple):
            return sig
        return (sig,)
    return ("force_graph_missing_graph_signature", force_graph_backend_kind(backend), id(backend))


def force_graph_prepare_replay(backend: Any, timestep: int) -> None:
    fn = getattr(backend, GRAPH_PREPARE_REPLAY, None)
    if callable(fn):
        fn(timestep)


def force_graph_can_add(backend: Any) -> bool:
    return force_graph_capabilities(backend).can_add


def force_graph_can_write(backend: Any) -> bool:
    return force_graph_capabilities(backend).can_write


def describe_force_graph_backend(backend: Any) -> Dict[str, Any]:
    '''Return a JSON-serializable protocol description for diagnostics/tests.'''
    return force_graph_capabilities(backend).as_dict()
