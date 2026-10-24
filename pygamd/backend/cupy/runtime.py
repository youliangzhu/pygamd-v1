'''
Small CuPy-runtime utilities shared by RawKernel backends.

This module is intentionally CUDA-driver-free at import time.  It contains only
plain-Python constants and helpers used to keep the CuPy-native branch coherent:

* stable names for core CuInfo fields;
* DPD velocity-source normalization for ordinary DPD vs DPDGWVV;
* one-dimensional RawKernel launch geometry;
* a tiny id-guarded device-view cache for graph-stable wrappers.

It does not import CuPy directly.  Callers pass CuPy arrays/modules that were
created through the existing lazy backend.cupy.context bridge.
'''

from __future__ import annotations

import math
from typing import Any, Dict, Tuple

# Core per-particle / global fields owned by CuInfo.  The values are public to
# internal backends only; user-facing code should still access fields through
# CuInfo.dev()/CuInfo.to_numpy().
CORE_FIELDS = (
    "pos", "vel", "velo", "force", "virial_potential", "image",
    "box", "type", "mass", "tag", "rtag",
)

# DPD velocity sources.  Ordinary DPD reads d_vel.  DPDGWVV reads d_velo, the
# graph-stable predictor velocity written by integration.dpdgwvv.
VELOCITY_SOURCE_VEL = "vel"
VELOCITY_SOURCE_VELO = "velo"
VALID_DPD_VELOCITY_SOURCES = (VELOCITY_SOURCE_VEL, VELOCITY_SOURCE_VELO)
DPDGWVV_DEFAULT_LAMBDA = 0.65


def normalize_dpd_velocity_source(source: Any) -> str:
    '''Return a canonical RawDPD velocity-source name or raise ValueError.'''
    source = str(source).strip().lower()
    if source not in VALID_DPD_VELOCITY_SOURCES:
        raise ValueError(
            "DPD velocity_source must be one of %s, got %r"
            % (VALID_DPD_VELOCITY_SOURCES, source)
        )
    return source


def require_cuinfo(info: Any, owner: str = "CuPy backend") -> Any:
    '''Return info.cu or raise a concise error for non-CuInfo snapshots.'''
    cu = getattr(info, "cu", None)
    if cu is None:
        raise RuntimeError("%s requires a CuInfo snapshot" % owner)
    return cu


def launch_config_1d(n_items: int, block_size: int) -> Tuple[Tuple[int], Tuple[int]]:
    '''Return ``(grid, block)`` for a one-dimensional RawKernel launch.'''
    n_items = int(n_items)
    block_size = int(block_size)
    if block_size < 1:
        raise ValueError("block_size must be >= 1")
    nblocks = max(1, int(math.ceil(max(n_items, 0) / float(block_size))))
    return (nblocks,), (block_size,)


class DeviceViewCache:
    '''Id-guarded cache for raveled device views used by RawKernel arguments.'''

    def __init__(self, as_cupy_array_func):
        self._as_cupy_array = as_cupy_array_func
        self._cache: Dict[str, Tuple[Any, Any]] = {}

    def get(self, key: str, arr: Any, ravel: bool = True) -> Any:
        cached = self._cache.get(key)
        if cached is not None and cached[0] is arr:
            return cached[1]
        view = self._as_cupy_array(arr, ravel=ravel)
        self._cache[key] = (arr, view)
        return view

    def clear(self) -> None:
        self._cache.clear()
